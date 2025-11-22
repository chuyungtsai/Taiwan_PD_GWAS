## ================================================================
## 08_beta_comparison.R
##
## Compare Taiwan PD GWAS effect sizes with:
##   - Nalls et al. 2019 (European PD GWAS)
##   - Foo et al. 2020 (Asian PD GWAS)
##
## Inputs:
##   - TWB GWAS: 20251004_TWB_imputed_logistic_PLINK2_frequency.csv
##   - Nalls:    20250906_Nalls_PD_PGS_hg38.xlsx (sheet 'Nalls_PD_hg_38')
##               near_gens_PD_GWAS.xlsx
##               Nalls_PD_Loci_Frequency.xlsx
##   - Foo:      20251004_Foo_PD_GWAS_location.xlsx
##               Foo_GWAS_frequency.xlsx
##
## Outputs:
##   - 20251004_Nalls_PD_lookup.csv
##   - 20251005_Nalls_sig_loci_for_PRS.csv
##   - 20251005_Foo_lookup_table.csv
##   - 20251005_Foo_output.csv
##   - 20251004_Nalls_beta_beta_plot.jpg
##   - 20251004_Foo_beta_beta_plot.jpg
##   - 20251004_Nalls_Foo_beta_plot_sharedLegend.jpg
## ================================================================

setwd("/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PLINK_analysis/20241225_TWB_imputed_analysis")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

## ---------------- USER CONFIGURATION ----------------------------

## TWB GWAS summary statistics (PLINK2 logistic + frequency)
sumstats_file <- "20251004_TWB_imputed_logistic_PLINK2_frequency.csv"

## Nalls 2019
nalls_file        <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/20250906_Nalls_PD_PGS_hg38.xlsx"
nalls_sheet       <- "Nalls_PD_hg_38"
nalls_gene_file   <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/near_gens_PD_GWAS.xlsx"
nalls_freq_file   <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/Nalls_PD_Loci_Frequency.xlsx"

## Foo 2020
foo_file          <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/20251004_Foo_PD_GWAS_location.xlsx"
foo_freq_file     <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/Foo_GWAS_frequency.xlsx"

## ================================================================
## 1. LOAD TAIWAN GWAS SUMMARY STATISTICS
## ================================================================

result_imputed_2 <- fread(sumstats_file)

## Expect columns: CHR, BP, ID, REF, ALT, A1_FREQ, OR, P, ...
## Harmonize and create numeric CHROM/POS for joining
result_imputed_2b <- result_imputed_2 %>%
  mutate(
    CHROM = as.integer(CHR),
    POS   = as.integer(BP),
    REF   = toupper(REF),
    ALT   = toupper(ALT)
  )

## Keep key fields plus everything else (for flexibility)
## (No select() here so we preserve original columns.)

## ================================================================
## 2. NALLS 2019 – LOOKUP & BETA–BETA COMPARISON
## ================================================================

## 2.1 Load Nalls loci + annotations
Nalls_PD           <- read_excel(nalls_file, sheet = nalls_sheet)
Nalls_PD_gene_names <- read_excel(nalls_gene_file)
Nalls_PD_gene_freq  <- read_excel(nalls_freq_file)

## Ensure integer positions and uppercase allele
Nalls_PD2 <- Nalls_PD %>%
  mutate(
    chr_name        = as.integer(chr_name),
    chr_position_38 = as.integer(chr_position_38),
    effect_allele   = toupper(effect_allele)
  )

## 2.2 Join with Taiwan GWAS by chr/pos
Nalls_PD_lookup <- Nalls_PD2 %>%
  left_join(
    result_imputed_2b,
    by = c("chr_name" = "CHROM", "chr_position_38" = "POS")
  ) %>%
  mutate(
    logOR = log(OR)
  ) %>%
  left_join(Nalls_PD_gene_names, by = "rsID") %>%
  left_join(Nalls_PD_gene_freq,  by = "rsID") %>%
  mutate(
    source              = "Nalls 2019",
    Nalls_effect_weight = effect_weight,
    ## orient Taiwan logOR to Nalls effect allele
    logOR_strand = if_else(effect_allele == ALT, logOR, -logOR)
  )

## Export full lookup table
Nalls_output <- Nalls_PD_lookup %>%
  select(
    source,
    rsID,
    Nearest_Gene,
    chr_name,
    chr_position_38,
    A1,
    REF,
    A1_FREQ,
    OR,
    L95,
    U95,
    logOR,
    P,
    `EAS Frequency`,
    `EUR Frequency`,
    Nalls_effect_weight
  )

fwrite(Nalls_output,
       "20251004_Nalls_PD_lookup.csv",
       quote = FALSE, row.names = FALSE)

## Restrict to variants actually present in TWB GWAS
Nalls_PD_lookup_present <- Nalls_PD_lookup %>%
  filter(!is.na(ID), is.finite(logOR_strand), is.finite(effect_weight))

## Nalls significance categories (based on TWB P)
Nalls_PD_lookup_present <- Nalls_PD_lookup_present %>%
  mutate(
    sig_status = case_when(
      P < 1e-5 ~ "Significant (p < 1e-5)",
      P < 0.05 ~ "Nominally Significant",
      TRUE     ~ "Not Significant"
    )
  )

## Loci with P < 0.05 in Taiwan (for PRS table)
Nalls_sig_loci <- Nalls_PD_lookup_present %>%
  filter(P < 0.05)

fwrite(Nalls_sig_loci,
       "20251005_Nalls_sig_loci_for_PRS.csv",
       quote = FALSE, row.names = FALSE)

## 2.3 Correlation Nalls vs TWB (beta-beta)
cor_test_Nalls <- cor.test(
  Nalls_PD_lookup_present$effect_weight,
  Nalls_PD_lookup_present$logOR_strand,
  method = "pearson"
)

R_val_Nalls <- round(cor_test_Nalls$estimate, 2)
p_val_Nalls <- signif(cor_test_Nalls$p.value, 2)

## 2.4 Nalls beta–beta plot
Nalls_plot <- ggplot(Nalls_PD_lookup_present,
                     aes(x = logOR_strand,
                         y = effect_weight,
                         color = sig_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text_repel(
    data = subset(Nalls_PD_lookup_present, P < 1e-5),
    aes(label = rsID),
    size = 3,
    max.overlaps = 20
  ) +
  annotate(
    "text",
    x = min(Nalls_PD_lookup_present$logOR_strand, na.rm = TRUE),
    y = max(Nalls_PD_lookup_present$effect_weight, na.rm = TRUE),
    label = paste0("R = ", R_val_Nalls, ", p = ", p_val_Nalls),
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  labs(
    x     = "Taiwan log(OR)",
    y     = "Nalls effect size (β)",
    color = "P-value (Taiwan)",
    title = "A. Taiwan vs. Europe (Nalls 2019)"
  ) +
  theme_minimal(base_size = 14)

ggsave("20251004_Nalls_beta_beta_plot.jpg",
       Nalls_plot, width = 6, height = 6, dpi = 400)

## ================================================================
## 3. FOO 2020 – LOOKUP & BETA–BETA COMPARISON
## ================================================================

Foo_PD          <- read_excel(foo_file)
Foo_GWAS_freq   <- read_excel(foo_freq_file)

Foo_PD2 <- Foo_PD %>%
  mutate(
    chr_name        = as.integer(Chromosome),
    chr_position_38 = as.integer(Position),
    ALT_Foo         = toupper(ALT_Foo),
    REF_Foo         = toupper(REF_Foo)
  )

Foo_PD_lookup <- Foo_PD2 %>%
  left_join(
    result_imputed_2b,
    by = c("chr_name" = "CHROM", "chr_position_38" = "POS")
  ) %>%
  mutate(
    logOR = log(OR),
    ALT   = toupper(ALT),
    REF   = toupper(REF)
  ) %>%
  ## orient TWB effect in same direction as Foo (ALT_Foo)
  mutate(
    logOR_strand = if_else(ALT_Foo == ALT, logOR, -logOR),
    OR_strand    = if_else(ALT_Foo == ALT, OR, 1 / OR)
  ) %>%
  left_join(Foo_GWAS_freq, by = "rsID") %>%
  mutate(
    source     = "Foo 2020",
    sig_status = case_when(
      P < 1e-5 ~ "Significant (p < 1e-5)",
      P < 0.05 ~ "Nominally Significant",
      TRUE     ~ "Not Significant"
    )
  )

## Full lookup & compact output table
Foo_output_full <- Foo_PD_lookup

Foo_output <- Foo_PD_lookup %>%
  select(
    source,
    rsID,
    Gene,
    chr_name,
    chr_position_38,
    ALT_Foo,
    REF_Foo,
    A1_FREQ,
    OR_strand,
    logOR_strand,
    P,
    `EAS Frequency`,
    `EUR Frequency`,
    Beta_Foo
  )

fwrite(Foo_output_full,
       "20251005_Foo_lookup_table.csv",
       quote = FALSE, row.names = FALSE)
fwrite(Foo_output,
       "20251005_Foo_output.csv",
       quote = FALSE, row.names = FALSE)

## 3.1 Correlation Foo vs TWB (beta-beta)
Foo_PD_lookup_present <- Foo_PD_lookup %>%
  filter(is.finite(logOR_strand), is.finite(Beta_Foo))

cor_test_Foo <- cor.test(
  Foo_PD_lookup_present$Beta_Foo,
  Foo_PD_lookup_present$logOR_strand,
  method = "pearson"
)

R_val_Foo <- round(cor_test_Foo$estimate, 2)
p_val_Foo <- signif(cor_test_Foo$p.value, 2)

## 3.2 Foo beta–beta plot
Foo_plot <- ggplot(Foo_PD_lookup_present,
                   aes(x = logOR_strand,
                       y = Beta_Foo,
                       color = sig_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text_repel(
    data = subset(Foo_PD_lookup_present, P < 1e-5),
    aes(label = rsID),
    size = 3,
    max.overlaps = 20
  ) +
  annotate(
    "text",
    x = min(Foo_PD_lookup_present$logOR_strand, na.rm = TRUE),
    y = max(Foo_PD_lookup_present$Beta_Foo,    na.rm = TRUE),
    label = paste0("R = ", R_val_Foo, ", p = ", p_val_Foo),
    hjust = 0,
    vjust = 1,
    size  = 3
  ) +
  labs(
    x     = "Taiwan log(OR)",
    y     = "Foo effect size (β)",
    color = "P-value (Taiwan)",
    title = "B. Taiwan vs. Asia (Foo 2020)"
  ) +
  theme_minimal(base_size = 14)

ggsave("20251004_Foo_beta_beta_plot.jpg",
       Foo_plot, width = 6, height = 6, dpi = 400)

## ================================================================
## 4. COMBINED 2-PANEL PLOT WITH SHARED LEGEND
## ================================================================

combined <- Nalls_plot + Foo_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("20251004_Nalls_Foo_beta_plot_sharedLegend.jpg",
       combined, width = 10, height = 5.5, dpi = 400)

message("08_beta_comparison finished:
  - Nalls lookup:                  20251004_Nalls_PD_lookup.csv
  - Nalls PRS loci (P<0.05):       20251005_Nalls_sig_loci_for_PRS.csv
  - Foo lookup + output tables:    20251005_Foo_lookup_table.csv, 20251005_Foo_output.csv
  - Beta–beta plots:               20251004_Nalls_beta_beta_plot.jpg,
                                    20251004_Foo_beta_beta_plot.jpg,
                                    20251004_Nalls_Foo_beta_plot_sharedLegend.jpg")
