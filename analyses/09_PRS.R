## ================================================================
## 09_PRS.R
##
## Polygenic Risk Score (PRS) workflow for Taiwan PD GWAS:
##   A. Derive Nalls 2019 overlap (76–77 SNPs) with TWB
##      → score file 20251006_Nalls_77.txt
##      → PRS_Nalls77 (per individual)
##   B. Build hybrid Nalls + Foo panel
##      → score file 20251006_Nalls_Foo_score.txt
##      → PRS_Nalls_Foo
##   C. Add 4 Taiwan-specific variants
##      → score file 2025_TW_4_variants_score.txt
##      → PRS_TW_only
##   D. Combine PRS and evaluate:
##      - AUC: Nalls vs Nalls+Foo vs Nalls+Foo+TW
##      - ROC figure
##      - Decile OR plot (default: PRS_Nalls_Foo_TW)
## ================================================================

## --------- 0. Setup ------------------------------------------------
setwd("/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PLINK_analysis/20241225_TWB_imputed_analysis")

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(glue)
  library(pROC)
  library(broom)
  library(forcats)
  library(tibble)
  library(tidyr)
  library(purrr)
})

## Helper: write PLINK2 --score file (ID, effect_allele, effect_weight; no header)
write_for_plink <- function(df, file) {
  df_out <- df %>%
    select(1:3) %>%
    mutate(across(everything(), as.character))
  write.table(
    df_out,
    file = file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

## Helper: strip leading "#" in column names (e.g. "#FID" → "FID")
normalize_cols <- function(df) {
  names(df) <- gsub("^#", "", names(df))
  df
}

## Paths & PLINK settings
sumstats_file   <- "20251004_TWB_imputed_logistic_PLINK2_frequency.csv"
plink_exec      <- "plink2"
bfile_prefix    <- "chr"
bfile_suffix    <- "_TWB_imputed"
out_dir         <- "PRS_chr_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

## Nalls 2019 files
nalls_file      <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/20250906_Nalls_PD_PGS_hg38.xlsx"
nalls_sheet     <- "Nalls_PD_hg_38"
nalls_gene_file <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/near_gens_PD_GWAS.xlsx"
nalls_freq_file <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/Nalls_PD_Loci_Frequency.xlsx"

## Foo 2020 files
foo_file        <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/20251004_Foo_PD_GWAS_location.xlsx"
foo_freq_file   <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PD_PGS/Foo_GWAS_frequency.xlsx"

## ================================================================
## 1. Load TWB summary statistics
## ================================================================
cat("\n[STEP 1] Loading Taiwan GWAS summary statistics...\n")

result_imputed_2 <- fread(sumstats_file)

## Expect PLINK2-style columns: CHR, BP, ID, REF, ALT, OR, P, ...
result_imputed_2b <- result_imputed_2 %>%
  mutate(
    CHROM = as.integer(CHR),
    POS   = as.integer(BP),
    REF   = toupper(REF),
    ALT   = toupper(ALT)
  )

## ================================================================
## 2. Nalls 2019 overlap and Nalls-only PRS (PRS_Nalls77)
## ================================================================
cat("\n[STEP 2] Building Nalls 2019 overlap and Nalls-only PRS...\n")

## 2.1 Load Nalls loci + annotation/frequency
Nalls_PD            <- read_excel(nalls_file, sheet = nalls_sheet)
Nalls_PD_gene_names <- read_excel(nalls_gene_file)
Nalls_PD_gene_freq  <- read_excel(nalls_freq_file)

Nalls_PD2 <- Nalls_PD %>%
  mutate(
    chr_name        = as.integer(chr_name),
    chr_position_38 = as.integer(chr_position_38),
    effect_allele   = toupper(effect_allele)
  )

## 2.2 Join to TWB by hg38 coordinates
Nalls_PD_lookup <- Nalls_PD2 %>%
  left_join(
    result_imputed_2b,
    by = c("chr_name" = "CHROM", "chr_position_38" = "POS")
  ) %>%
  mutate(
    logOR = log(OR)
  ) %>%
  left_join(Nalls_PD_gene_names, by = "rsID") %>%
  left_join(Nalls_PD_gene_freq,  by = "rsID")

## Keep variants present in TWB GWAS
Nalls_PD_lookup_present <- Nalls_PD_lookup %>%
  filter(!is.na(ID))

## Orient Taiwan logOR to Nalls effect_allele (for QC/plots; not required for PRS)
Nalls_PD_lookup_present <- Nalls_PD_lookup_present %>%
  mutate(
    logOR_strand = if_else(effect_allele == ALT, logOR, -logOR)
  )

## Export lookup table and simple PRS table
fwrite(
  Nalls_PD_lookup_present,
  "20251006_Nalls_PD_lookup_present.csv",
  quote = FALSE,
  row.names = FALSE
)

Nalls_PD_PRS_table <- Nalls_PD_lookup_present %>%
  select(rsID, chr_name, chr_position_38, effect_allele, effect_weight)

fwrite(
  Nalls_PD_PRS_table,
  "20251010_Nalls_PD_76.csv",
  quote = FALSE,
  row.names = FALSE
)

## 2.3 Build Nalls-77 score file for PLINK2 --score
## NOTE: original comment: "need to remove one row chr3:28664199:T:C"
nalls_score_df <- Nalls_PD_lookup_present %>%
  select(ID, effect_allele, effect_weight) %>%
  mutate(
    ID            = as.character(ID),
    effect_allele = as.character(effect_allele),
    effect_weight = as.numeric(effect_weight)
  ) %>%
  filter(ID != "chr3:28664199:T:C")  ## remove problematic SNP explicitly

write_for_plink(nalls_score_df, "20251006_Nalls_77.txt")

## 2.4 Run PLINK2 scoring for Nalls-77 across chr1–22
score_file_Nalls <- "20251006_Nalls_77.txt"

for (chr in 1:22) {
  bfile_path <- sprintf("%s%d%s", bfile_prefix, chr, bfile_suffix)
  out_path   <- file.path(out_dir, sprintf("chr%d_Nalls_77", chr))
  
  cmd <- sprintf(
    '%s --bfile %s --score %s cols=+scoresums --out %s',
    plink_exec, bfile_path, score_file_Nalls, out_path
  )
  cat("\nRunning:", cmd, "\n")
  system(cmd)
}

## 2.5 Merge .sscore files into PRS_Nalls77
sscore_files_nalls <- list.files(
  out_dir,
  pattern = "_Nalls_77.sscore$",
  full.names = TRUE
)

PRS_all <- map_dfr(sscore_files_nalls, fread) %>%
  normalize_cols() %>%
  group_by(FID, IID, PHENO1) %>%
  summarise(
    PRS_Nalls77 = sum(SCORE1_SUM, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  PRS_all,
  file.path(out_dir, "20251006_PRS_Nalls77_combined.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n✅ Finished computing PRS_Nalls77.\n")

## ================================================================
## 3. Hybrid Nalls + Foo PRS (PRS_Nalls_Foo)
## ================================================================
cat("\n[STEP 3] Building Nalls+Foo hybrid PRS...\n")

## 3.1 Foo GWAS summary (location + frequency)
Foo_PD        <- read_excel(foo_file)
Foo_GWAS_freq <- read_excel(foo_freq_file)

Foo_PD2 <- Foo_PD %>%
  mutate(
    chr_name        = as.integer(Chromosome),
    chr_position_38 = as.integer(Position),
    ALT_Foo         = toupper(ALT_Foo),
    REF_Foo         = toupper(REF_Foo)
  )

## Join Foo to TWB to inherit ID and allele (for QC / linking, not for weight derivation)
Foo_PD_lookup <- Foo_PD2 %>%
  left_join(
    result_imputed_2b,
    by = c("chr_name" = "CHROM", "chr_position_38" = "POS")
  ) %>%
  left_join(Foo_GWAS_freq, by = "rsID") %>%
  mutate(
    ALT   = toupper(ALT),
    REF   = toupper(REF),
    logOR = log(OR)
  )

## Output Foo PRS table (for record)
Foo_PRS_table <- Foo_PD_lookup %>%
  select(rsID, chr_name, chr_position_38, ID, ALT_Foo, Beta_Foo)

fwrite(
  Foo_PRS_table,
  "20251010_Foo_PRS.csv",
  quote = FALSE,
  row.names = FALSE
)

## 3.2 Build combined Nalls + Foo score file
Foo_effect <- Foo_PD_lookup %>%
  select(ID, ALT_Foo, Beta_Foo) %>%
  filter(!is.na(ID)) %>%
  mutate(
    ID       = as.character(ID),
    ALT_Foo  = as.character(ALT_Foo),
    Beta_Foo = as.numeric(Beta_Foo)
  )

colnames(Foo_effect) <- c("ID", "effect_allele", "effect_weight")

Foo_ID <- Foo_effect$ID

Nalls_wo_Foo <- Nalls_PD_lookup_present %>%
  filter(!is.na(ID)) %>%
  filter(!(ID %in% Foo_ID)) %>%
  filter(ID != "chr3:28664199:T:C") %>%  ## same problematic variant
  select(ID, effect_allele, effect_weight) %>%
  mutate(
    ID            = as.character(ID),
    effect_allele = as.character(effect_allele),
    effect_weight = as.numeric(effect_weight)
  )

Nalls_Foo_score <- bind_rows(Nalls_wo_Foo, Foo_effect) %>%
  distinct(ID, effect_allele, .keep_all = TRUE)

write_for_plink(Nalls_Foo_score, "20251006_Nalls_Foo_score.txt")

## 3.3 Run PLINK2 scoring for Nalls_Foo across chr1–22
score_file_Nalls_Foo <- "20251006_Nalls_Foo_score.txt"

for (chr in 1:22) {
  bfile_path <- sprintf("%s%d%s", bfile_prefix, chr, bfile_suffix)
  out_path   <- file.path(out_dir, sprintf("chr%d_Nalls_Foo", chr))
  
  cmd <- sprintf(
    '%s --bfile %s --score %s cols=+scoresums --out %s',
    plink_exec, bfile_path, score_file_Nalls_Foo, out_path
  )
  cat("\nRunning:", cmd, "\n")
  system(cmd)
}

sscore_files_foo <- list.files(
  out_dir,
  pattern = "_Nalls_Foo.sscore$",
  full.names = TRUE
)

PRS_Nalls_Foo <- map_dfr(sscore_files_foo, fread) %>%
  normalize_cols() %>%
  group_by(FID, IID, PHENO1) %>%
  summarise(
    PRS_Nalls_Foo = sum(SCORE1_SUM, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  PRS_Nalls_Foo,
  file.path(out_dir, "20251006_PRS_Nalls_Foo_combined.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n✅ Finished computing PRS_Nalls_Foo.\n")

## ================================================================
## 4. AUC: Nalls-77 vs Nalls+Foo (with DeLong test)
## ================================================================
cat("\n[STEP 4] AUC comparison: Nalls-77 vs Nalls+Foo...\n")

df_comp_12 <- PRS_all %>%
  inner_join(PRS_Nalls_Foo, by = c("FID", "IID", "PHENO1"))

df_comp_12 <- df_comp_12 %>%
  mutate(
    CaseCtrl = case_when(
      PHENO1 == 2 ~ 1,
      PHENO1 == 1 ~ 0,
      TRUE        ~ NA_real_
    )
  ) %>%
  filter(!is.na(CaseCtrl))

df_comp_12 <- df_comp_12 %>%
  mutate(
    PRS_Nalls77_z = as.numeric(scale(PRS_Nalls77)),
    PRS_Foo_z     = as.numeric(scale(PRS_Nalls_Foo))
  )

roc_nalls77 <- roc(df_comp_12$CaseCtrl, df_comp_12$PRS_Nalls77_z, quiet = TRUE, direction = "auto")
roc_foo     <- roc(df_comp_12$CaseCtrl, df_comp_12$PRS_Foo_z,     quiet = TRUE, direction = "auto")

auc_tbl_12 <- tibble(
  Model = c("Nalls-77 PRS", "PRS_Nalls_Foo"),
  AUC   = c(as.numeric(auc(roc_nalls77)), as.numeric(auc(roc_foo))),
  CI_lo = c(ci.auc(roc_nalls77)[1],       ci.auc(roc_foo)[1]),
  CI_hi = c(ci.auc(roc_nalls77)[3],       ci.auc(roc_foo)[3])
)

print(auc_tbl_12)
fwrite(auc_tbl_12, "AUC_compare_Nalls77_vs_Foo.tsv", sep = "\t", quote = FALSE)

cat("\nDeLong test (Nalls-77 vs Nalls+Foo):\n")
print(roc.test(roc_nalls77, roc_foo, method = "delong"))

## Optional base R ROC plot
jpeg("20251006_ROC_Nalls77_vs_NallsFoo.jpg", width = 900, height = 900, res = 150)
plot(roc_nalls77, col = "#1f78b4", lwd = 2,
     main = "ROC: Nalls-77 vs Nalls+Foo")
plot(roc_foo, add = TRUE, col = "#e31a1c", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey70")
legend("bottomright",
       legend = c(
         sprintf("Nalls-77 (AUC=%.3f)", auc(roc_nalls77)),
         sprintf("Nalls+Foo (AUC=%.3f)", auc(roc_foo))
       ),
       col = c("#1f78b4", "#e31a1c"), lwd = 2, bty = "n")
dev.off()

## ================================================================
## 5. Add 4 Taiwan-specific variants → PRS_TW_only
## ================================================================
cat("\n[STEP 5] Adding 4 Taiwan-specific variants...\n")

## LRRK2 G2385R, R1628P, SCARB2, VPS13C
TW <- data.frame(
  ID            = c("chr12:40363526:G:A",
                    "chr12:40320043:G:C",
                    "chr4:76217199:C:G",
                    "chr15:61994009:A:G"),
  effect_allele = c("A", "C", "G", "G"),
  effect_weight = c(0.47, 0.398, 0.183, -0.189)
)

TW_PRS_table <- data.frame(
  rsID          = c("rs34778348", "rs33949390", "rs76591264", "rs12900645"),
  chromosome    = c(12, 12, 4, 15),
  position      = c(40363526, 40320043, 76217199, 61994009),
  effect_allele = c("A", "C", "G", "G"),
  effect_weight = c(0.47, 0.398, 0.183, -0.189)
)

fwrite(
  TW_PRS_table,
  "20251010_TW_PRS_table.csv",
  quote = FALSE,
  row.names = FALSE
)

write_for_plink(TW, "2025_TW_4_variants_score.txt")

## Run PLINK2 only on chr4, chr12, chr15
score_file_TW <- "2025_TW_4_variants_score.txt"
chr_TW        <- c(4, 12, 15)

for (chr in chr_TW) {
  bfile_path <- sprintf("%s%d%s", bfile_prefix, chr, bfile_suffix)
  out_path   <- file.path(out_dir, sprintf("chr%d_TW_4", chr))
  
  cmd <- sprintf(
    '%s --bfile %s --score %s cols=+scoresums --out %s',
    plink_exec, bfile_path, score_file_TW, out_path
  )
  cat("\nRunning:", cmd, "\n")
  system(cmd)
}

sscore_files_tw <- list.files(
  out_dir,
  pattern = "_TW_4.sscore$",
  full.names = TRUE
)

PRS_TW_only <- map_dfr(sscore_files_tw, fread) %>%
  normalize_cols() %>%
  group_by(FID, IID, PHENO1) %>%
  summarise(
    PRS_TW_only = sum(SCORE1_SUM, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n✅ Finished computing PRS_TW_only.\n")

## ================================================================
## 6. Combine PRS and 3-model ROC (Nalls / Nalls+Foo / Nalls+Foo+TW)
## ================================================================
cat("\n[STEP 6] Combining PRS and generating 3-model ROC...\n")

PRS_total <- PRS_all %>%
  inner_join(PRS_Nalls_Foo, by = c("FID", "IID", "PHENO1")) %>%
  inner_join(PRS_TW_only,   by = c("FID", "IID", "PHENO1"))

PRS_total_short <- PRS_total %>%
  mutate(PRS_Nalls_Foo_TW = PRS_Nalls_Foo + PRS_TW_only)

fwrite(
  PRS_total_short,
  "20251006_PRS_total.csv",
  quote = FALSE,
  row.names = FALSE
)

## ROC for 3 models
dat <- PRS_total_short %>%
  mutate(CaseCtrl = ifelse(PHENO1 == 2, 1, 0)) %>%
  drop_na(PRS_Nalls77, PRS_Nalls_Foo, PRS_Nalls_Foo_TW)

roc_nalls   <- roc(dat$CaseCtrl, dat$PRS_Nalls77,      quiet = TRUE, direction = "auto")
roc_foo2    <- roc(dat$CaseCtrl, dat$PRS_Nalls_Foo,    quiet = TRUE, direction = "auto")
roc_foo_tw  <- roc(dat$CaseCtrl, dat$PRS_Nalls_Foo_TW, quiet = TRUE, direction = "auto")

auc_tbl_3 <- tibble(
  Model = c("Nalls", "Nalls+Foo", "Nalls+Foo+TW"),
  AUC   = c(auc(roc_nalls), auc(roc_foo2), auc(roc_foo_tw)),
  CI_lo = c(ci.auc(roc_nalls)[1], ci.auc(roc_foo2)[1], ci.auc(roc_foo_tw)[1]),
  CI_hi = c(ci.auc(roc_nalls)[3], ci.auc(roc_foo2)[3], ci.auc(roc_foo_tw)[3])
)

print(auc_tbl_3)
fwrite(auc_tbl_3, "AUC_PRSTotal_models.tsv", sep = "\t", quote = FALSE)

## Build ROC dataframe for ggplot
roc_df <- bind_rows(
  tibble(Model = "Nalls",
         FPR   = 1 - roc_nalls$specificities,
         TPR   = roc_nalls$sensitivities),
  tibble(Model = "Nalls+Foo",
         FPR   = 1 - roc_foo2$specificities,
         TPR   = roc_foo2$sensitivities),
  tibble(Model = "Nalls+Foo+TW",
         FPR   = 1 - roc_foo_tw$specificities,
         TPR   = roc_foo_tw$sensitivities)
)

cols <- c(
  "Nalls"        = "#f1a340",
  "Nalls+Foo"    = "#0571b0",
  "Nalls+Foo+TW" = "#2ca25f"
)

annot_tbl <- tibble(
  label = c(
    glue("Nalls: AUC={sprintf('%.3f', auc_tbl_3$AUC[1])} (95%CI {sprintf('%.3f', auc_tbl_3$CI_lo[1])}-{sprintf('%.3f', auc_tbl_3$CI_hi[1])})"),
    glue("Nalls+Foo: AUC={sprintf('%.3f', auc_tbl_3$AUC[2])} (95%CI {sprintf('%.3f', auc_tbl_3$CI_lo[2])}-{sprintf('%.3f', auc_tbl_3$CI_hi[2])})"),
    glue("Nalls+Foo+TW: AUC={sprintf('%.3f', auc_tbl_3$AUC[3])} (95%CI {sprintf('%.3f', auc_tbl_3$CI_lo[3])}-{sprintf('%.3f', auc_tbl_3$CI_hi[3])})")
  ),
  ypos    = c(0.27, 0.20, 0.13),
  dot_col = c(cols["Nalls"], cols["Nalls+Foo"], cols["Nalls+Foo+TW"])
)

box_xmin <- 0.52; box_xmax <- 0.995
box_ymin <- 0.05; box_ymax <- 0.32
x_text   <- box_xmin + 0.04
x_dots   <- box_xmin + 0.02

prs_AUC <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  annotate("rect", xmin = box_xmin, xmax = box_xmax, ymin = box_ymin, ymax = box_ymax,
           fill = "white", color = "grey80", linewidth = 0.4) +
  geom_point(
    data = annot_tbl,
    aes(x = x_dots, y = ypos),
    inherit.aes = FALSE,
    shape = 15, size = 3, color = annot_tbl$dot_col
  ) +
  geom_text(
    data = annot_tbl,
    aes(x = x_text, y = ypos, label = label),
    inherit.aes = FALSE,
    hjust = 0, size = 3.2, color = "black"
  ) +
  scale_color_manual(values = cols) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "ROC Curves for Polygenic Risk Score Models",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave(
  "20251006_ROC_PRSTotal_noDelong.jpg",
  plot = prs_AUC,
  width = 8,
  height = 6,
  dpi = 400
)

cat("\n✅ ROC figure (3 models) saved as 20251006_ROC_PRSTotal_noDelong.jpg\n")

## ================================================================
## 7. Decile OR plot for PRS_Nalls_Foo_TW
## ================================================================
cat("\n[STEP 7] Decile OR plot for PRS_Nalls_Foo_TW...\n")

PRS_COL <- "PRS_Nalls_Foo_TW"

df_dec <- PRS_total_short %>%
  select(FID, IID, PHENO1, all_of(PRS_COL)) %>%
  rename(PRS = !!PRS_COL) %>%
  filter(!is.na(PRS), !is.na(PHENO1)) %>%
  mutate(
    PRS_decile = ntile(PRS, 10),
    PRS_decile = factor(PRS_decile,
                        levels = 1:10,
                        labels = paste0("D", 1:10)),
    CaseCtrl = ifelse(PHENO1 == 2, 1, 0)
  ) %>%
  filter(!is.na(CaseCtrl))

fit_dec <- glm(CaseCtrl ~ PRS_decile, data = df_dec, family = binomial())

or_tbl <- tidy(fit_dec, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Decile = gsub("PRS_decile", "", term),
    Decile = factor(Decile, levels = paste0("D", 2:10)),
    OR     = estimate,
    CI_low = conf.low,
    CI_high = conf.high
  )

ref_row <- tibble(Decile = "D1", OR = 1, CI_low = 1, CI_high = 1)

or_tbl <- bind_rows(ref_row, or_tbl) %>%
  mutate(Decile = factor(Decile, levels = paste0("D", 1:10)))

p_decile <- ggplot(or_tbl, aes(x = Decile, y = OR)) +
  geom_point(size = 3, color = "#2c7bb6") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, color = "#2c7bb6") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_log10(
    limits = c(0.5, 7),
    breaks = c(0.5, 1, 2, 3, 4, 5),
    labels = c("0.5", "1", "2", "3", "4", "5")
  ) +
  labs(
    title    = sprintf("Odds of PD across PRS deciles (%s)", PRS_COL),
    subtitle = "Reference = D1 (lowest decile); Odds Ratio with 95% CI",
    x        = "PRS Decile",
    y        = "Odds Ratio (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x   = element_text(size = 12, color = "black"),
    axis.text.y   = element_text(size = 12, color = "black")
  )

out_decile <- sprintf("20251006_PRS_decile_stratification_%s.jpg", PRS_COL)
ggsave(out_decile, plot = p_decile, width = 7.5, height = 5.5, dpi = 400)
fwrite(or_tbl, "20251006_OR_decile.csv", quote = FALSE, row.names = FALSE)

cat("\n✅ Decile OR plot saved as", out_decile, "\n")
cat("✅ Decile OR table saved as 20251006_OR_decile.csv\n")

cat("\nAll PRS steps completed.\n")
