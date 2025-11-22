## ================================================================
## SNCA haplotype analysis and visualization (chr4 region)
## - PLINK2 / PLINK1.9 pre-processing (region subset, LD, VCF)
## - geneHapR haplotypes and haplotype table plot
## - LD heatmap for 5 key SNCA variants
## - haplo.cc (haplo.stats) case–control association
## - publication-style combined haplotype + forest plot
## ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(geneHapR)
  library(haplo.stats)
  library(forcats)
  library(patchwork)
})

## ---------------- USER CONFIGURATION ----------------------------

# Working directory for this analysis
setwd("/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PLINK_analysis/20241225_TWB_imputed_analysis")

# Input genotype prefix (TWB imputed, chr4)
bfile_chr4 <- "chr4_TWB_imputed"

# Covariate file (PCs etc.)
covar_file <- "~/Desktop/Bioinfo/NTU_GWAS/PLINK_analysis/20241003_GWAS/20241119_covariate_PC15.txt"

# PLINK / tabix executables (assumes they are in $PATH)
plink2_bin  <- "plink2"
plink1_bin  <- "plink19"
tabix_bin   <- "tabix"

# SNCA region (hg38)
chr        <- 4
start_bp   <- 89500000L
end_bp     <- 89900000L

# 5 key SNCA variants (PLINK IDs must match your .bim/.vcf)
SNCA_variants <- c(
  "chr4:89550094:T:A",  # rs2870004, SNCA 3'
  "chr4:89704960:G:A",  # rs356182, Nalls 2019, SNCA 3'
  "chr4:89744890:C:T",  # rs356203, European signal
  "chr4:89761323:T:C",  # rs6826785, Foo 2020
  "chr4:89765968:G:A"   # rs3857061 (Taiwan SNCA intron 4 top hit)
)

SNCA_labels <- c(
  "rs2870004 (3')",
  "rs356182 (Nalls 2019)",
  "rs356203 (EUR)",
  "rs6826785 (Foo 2020)",
  "chr4:89765968 (TW top hit)"
)

# Output file prefixes
region_prefix   <- "SNCA_region"
ld_prefix       <- "SNCA_5_variants_LD"
hap_vcf_prefix  <- "SNCA_region_5_variants"
block_prefix    <- "20251006_SNCA_block"
ped_prefix      <- "20251006_SNCA"
hap_table_plot  <- "SNCA_haplotype_table_plot.jpg"
ld_heatmap_file <- "SNCA_LD_heatmap_lower_triangle.jpg"
hap_unadj_out   <- "SNCA_haplotype_results_unadjusted.txt"
hap_adj_out     <- "SNCA_haplotype_results_adjusted.txt"
hap_or_table    <- "SNCA_haplotype_OR_table.csv"
hap_or_common   <- "Common_SNCA_haplotype_OR_table.csv"
combo_plot_file <- "SNCA_haplotype_combined_final.jpg"

## ================================================================
## 1. PLINK PRE-PROCESSING (REGION, HAPLOTYPE VCF, LD)
## ================================================================

# Write SNCA variant list used by PLINK --extract
write.table(
  SNCA_variants,
  file      = "SNCA_variants.txt",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 1.1 Subset chr4 SNCA region
system(glue::glue(
  "{plink2_bin} --bfile {bfile_chr4} ",
  "--chr {chr} --from-bp {start_bp} --to-bp {end_bp} ",
  "--geno 0.02 --make-bed --out {region_prefix}"
))

# 1.2 Extract 5 SNCA variants and export VCF (for geneHapR)
system(glue::glue(
  "{plink1_bin} --bfile {region_prefix} ",
  "--extract SNCA_variants.txt ",
  "--recode vcf bgz --out {hap_vcf_prefix}"
))

system(glue::glue(
  "{tabix_bin} -p vcf {hap_vcf_prefix}.vcf.gz"
))

# 1.3 Compute LD among the 5 SNCA variants
system(glue::glue(
  "{plink1_bin} --bfile {region_prefix} ",
  "--extract SNCA_variants.txt ",
  "--r2 --ld-window 99999 --ld-window-kb 250 --ld-window-r2 0 ",
  "--out {ld_prefix}"
))

## ================================================================
## 2. HAPLOTYPES VIA geneHapR (HOMOZYGOUS ONLY) + TABLE PLOT
## ================================================================

hapdat <- import_vcf(paste0(hap_vcf_prefix, ".vcf.gz"))

# Include heterozygotes for frequency; use hapResult_homo as in your script
hapResult_hetero <- vcf2hap(hapdat, hetero_remove = FALSE)
hapResult_homo   <- vcf2hap(hapdat, hetero_remove = TRUE)

hapSummary <- hap_summary(hapResult_homo)
write.hap(hapResult_homo, file = "20251005_SNCA_homo.hapResult")

p_hap <- plotHapTable(hapSummary)
ggsave(
  filename = hap_table_plot,
  plot     = p_hap,
  width    = 8,
  height   = 6,
  dpi      = 600
)

## ================================================================
## 3. LD HEATMAP FOR 5 SNCA VARIANTS
## ================================================================

ld <- read.table(paste0(ld_prefix, ".ld"), header = TRUE)

ids <- SNCA_variants   # same order for matrix & labels
snp_labels <- SNCA_labels

# Build symmetric R^2 matrix
mat <- matrix(
  1,
  nrow = length(ids),
  ncol = length(ids),
  dimnames = list(ids, ids)
)

ld_sub <- ld %>%
  filter(SNP_A %in% ids & SNP_B %in% ids)

for (i in seq_len(nrow(ld_sub))) {
  mat[ld_sub$SNP_A[i], ld_sub$SNP_B[i]] <- ld_sub$R2[i]
  mat[ld_sub$SNP_B[i], ld_sub$SNP_A[i]] <- ld_sub$R2[i]
}

# Lower triangle only, for plotting
df_ld <- as.data.frame(mat) %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "r2") %>%
  mutate(
    Var1 = factor(Var1, levels = ids),
    Var2 = factor(Var2, levels = ids),
    i    = as.integer(Var1),
    j    = as.integer(Var2)
  ) %>%
  filter(i >= j)

p_ld_lower <- ggplot(df_ld, aes(x = Var1, y = Var2, fill = r2)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 4, color = "black") +
  scale_x_discrete(labels = snp_labels, expand = c(0, 0)) +
  scale_y_discrete(labels = snp_labels, expand = c(0, 0)) +
  scale_fill_gradient(
    name   = expression(r^2),
    limits = c(0, 1),
    low    = "gray95",
    high   = "firebrick3"
  ) +
  coord_fixed() +
  labs(
    title    = "LD among SNCA variants",
    subtitle = "Five key SNCA variants from TW and literature sources",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y   = element_text(size = 12, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, lineheight = 1.2),
    legend.position = "right",
    panel.grid    = element_blank()
  )

ggsave(
  ld_heatmap_file,
  p_ld_lower,
  width  = 6.5,
  height = 6,
  dpi    = 600
)

## ================================================================
## 4. haplo.cc ANALYSIS (haplo.stats)
## ================================================================

# 4.1 Define blocks (optional) and export genotype
system(glue::glue(
  "{plink1_bin} --bfile {region_prefix} --blocks --out {block_prefix}"
))

system(glue::glue(
  "{plink1_bin} --bfile {region_prefix} --recode --out {ped_prefix}"
))

# 4.2 Load PED/MAP into R
PED <- fread(paste0(ped_prefix, ".ped"), data.table = FALSE)
MAP <- fread(paste0(ped_prefix, ".map"), data.table = FALSE)
colnames(MAP) <- c("CHR", "SNP", "CM", "POS")

# Duplicate SNPs for two allele columns and order by POS
MAP_1 <- MAP
MAP_2 <- MAP
MAP_1$SNP <- paste0(MAP_1$SNP, "_1")
MAP_2$SNP <- paste0(MAP_2$SNP, "_2")

MAP_2alleles <- rbind(MAP_1, MAP_2)
MAP_pos      <- MAP_2alleles[order(MAP_2alleles$POS), ]
SNCA_alleles <- MAP_pos$SNP

# Assign column names to PED
colnames(PED) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO", SNCA_alleles)
write.table(
  PED,
  file      = "Geno_matrix_SNCA.tab",
  quote     = FALSE,
  row.names = FALSE,
  sep       = "\t"
)

# 4.3 Covariates (AGE, SEX, PC1–PC5)
SampleInfo_Adjustment <- read.table(covar_file, header = TRUE)
colnames(SampleInfo_Adjustment) <- c(
  "Case", "IID", "SEX", "AGE",
  paste0("PC", 1:15)
)
adj <- data.frame(SampleInfo_Adjustment[, c("SEX", "AGE", paste0("PC", 1:5))])

# 4.4 Extract haplotype SNPs (5 SNCA markers)
SNCA_allele_list <- c(
  "chr4:89550094",
  "chr4:89704960",
  "chr4:89744890",
  "chr4:89761323",
  "chr4:89765968"
)

Geno_matrix_SNCA <- read.delim("Geno_matrix_SNCA.tab", check.names = FALSE)

H1_SNCA <- Geno_matrix_SNCA[, c(
  "FID", "IID", "PAT", "MAT", "SEX", "PHENO",
  grep(paste(SNCA_allele_list, collapse = "|"),
       colnames(Geno_matrix_SNCA),
       value = TRUE)
)]

# 4.5 Genotype matrix (two columns per SNP)
geno <- data.frame(H1_SNCA[, 7:ncol(H1_SNCA)], check.names = FALSE)

# 4.6 Build haplo.stats geno object
myGeno <- setupGeno(geno)
summaryGeno(myGeno)

# 4.7 Locus labels (for haplo.stats output)
label <- c("rs2870004", "rs356182", "rs356203", "rs6826785", "rs3857061")

# 4.8 Binary phenotype (0 = control, 1 = case)
H1_SNCA$PHENO_01 <- H1_SNCA$PHENO - 1
y.bin <- 1L * (H1_SNCA$PHENO_01 == 1)

# 4.9 Non-adjusted haplotype association
H1_SNCA_res <- haplo.cc(
  y           = y.bin,
  geno        = geno,
  locus.label = label,
  control     = haplo.glm.control(haplo.freq.min = 0.01)
)

H1_SNCA_unadj_df      <- H1_SNCA_res$cc.df
H1_SNCA_unadj_df_sort <- H1_SNCA_unadj_df[order(H1_SNCA_unadj_df$`p-val`), ]

## ================================================================
## 5. ADJUSTED haplo.cc (AGE, SEX, PC1–PC5)
## ================================================================

# Align IDs across genotype, phenotype and covariates
geno_IDs <- H1_SNCA$IID
adj_IDs  <- SampleInfo_Adjustment$IID

rownames(geno) <- geno_IDs
rownames(adj)  <- adj_IDs

common_ids <- intersect(rownames(geno), adj_IDs)

geno_aligned <- geno[common_ids, , drop = FALSE]
adj_aligned  <- adj[common_ids, , drop = FALSE]

y_aligned <- y.bin[match(common_ids, H1_SNCA$IID)]

stopifnot(identical(rownames(geno_aligned), rownames(adj_aligned)))
stopifnot(length(y_aligned) == nrow(geno_aligned))

H1_SNCA_adj <- haplo.cc(
  y           = y_aligned,
  geno        = geno_aligned,
  locus.label = label,
  x.adj       = adj_aligned,
  control     = haplo.glm.control(haplo.freq.min = 0.01)
)

H1_SNCA_adj_df      <- H1_SNCA_adj$cc.df
H1_SNCA_adj_df_sort <- H1_SNCA_adj_df[order(H1_SNCA_adj_df$`p-val`), ]

## ================================================================
## 6. EXPORT HAPLOTYPE RESULT TABLES
## ================================================================

write.table(
  H1_SNCA_unadj_df_sort,
  hap_unadj_out,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  H1_SNCA_adj_df_sort,
  hap_adj_out,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

hap_results <- H1_SNCA_adj_df_sort %>%
  select(
    rs2870004, rs356182, rs356203, rs6826785, rs3857061,
    glm.eff, control.hf, case.hf, OR, OR.lower, OR.upper, `p-val`
  ) %>%
  rename(
    P_value      = `p-val`,
    Control_Freq = control.hf,
    Case_Freq    = case.hf
  ) %>%
  mutate(
    Significant = ifelse(P_value < 0.05, "Yes", "No")
  )

write.csv(hap_results, hap_or_table, row.names = FALSE)

# Keep only common haplotypes (Control_Freq > 0.01)
SNCA_hap_results_common <- hap_results %>%
  filter(Control_Freq > 0.01)

write.csv(SNCA_hap_results_common, hap_or_common, row.names = FALSE)

## ================================================================
## 7. PUBLICATION-STYLE FIGURE:
##    LEFT: allele grid + case/control frequencies
##    RIGHT: forest plot for OR (AGCCA reference; AGTTG/AATTG highlighted)
## ================================================================

df <- SNCA_hap_results_common %>%
  mutate(
    Haplotype = paste0(rs2870004, rs356182, rs356203, rs6826785, rs3857061),
    Haplotype = toupper(gsub("\\s+", "", Haplotype)),
    Ctrl_pct  = 100 * Control_Freq,
    Case_pct  = 100 * Case_Freq
  )

# Order haplotypes: reference (AGCCA) first
ord_after  <- df %>% arrange(desc(Control_Freq)) %>% pull(Haplotype) %>% unique()
ord_levels <- c("AGCCA", setdiff(ord_after, "AGCCA"))
df$Haplotype <- factor(df$Haplotype, levels = ord_levels)

# ---- Allele grid (left panel) ----
pos_order   <- c("89550094", "89704960", "89744890", "89761323", "89765968")
refalt_labs <- c("A/T", "G/A", "C/T", "C/T", "A/G")
col_headers <- paste0(pos_order, "\n", refalt_labs)

allele_long <- df %>%
  select(Haplotype) %>%
  distinct() %>%
  mutate(
    a1 = substr(Haplotype, 1, 1),
    a2 = substr(Haplotype, 2, 2),
    a3 = substr(Haplotype, 3, 3),
    a4 = substr(Haplotype, 4, 4),
    a5 = substr(Haplotype, 5, 5)
  ) %>%
  pivot_longer(a1:a5, names_to = "col", values_to = "Allele") %>%
  mutate(
    POS = factor(
      c("89550094", "89704960", "89744890", "89761323", "89765968")[
        as.integer(substr(col, 2, 2))
      ],
      levels = pos_order
    )
  ) %>%
  select(Haplotype, POS, Allele) %>%
  mutate(Highlight = ifelse(Haplotype == "AGCCA", "Highlight", "Normal"))

freq_df <- df %>%
  select(Haplotype, Ctrl_pct, Case_pct) %>%
  distinct()

base_palette <- c(
  A = "#E15759",
  C = "#59A14F",
  G = "#4E79A7",
  T = "#B07AA1"
)

p_grid <- ggplot(allele_long, aes(x = POS, y = Haplotype, fill = Allele)) +
  geom_tile(
    aes(alpha = ifelse(Highlight == "Highlight", 1, 0.8),
        colour = ifelse(Highlight == "Highlight", "#000000", NA)),
    linewidth = 0.8,
    width     = 0.98,
    height    = 0.98
  ) +
  geom_text(aes(label = Allele), size = 4.2, color = "black") +
  scale_fill_manual(values = base_palette, drop = FALSE) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_x_discrete(
    labels = setNames(col_headers, pos_order),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Haplotype") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 10, color = "black"),
    axis.text.y     = element_text(size = 11, color = "black"),
    panel.grid      = element_blank(),
    plot.margin     = margin(t = 5, r = 0, b = 5, l = 5)
  )

# Case + control frequencies (middle narrow column)
p_freq <- ggplot(freq_df, aes(x = 1, y = Haplotype)) +
  geom_tile(fill = "grey90", width = 0.98, height = 0.98) +
  geom_text(
    aes(label = sprintf("Ctrl %.2f%%\nCase %.2f%%", Ctrl_pct, Case_pct)),
    size      = 3.7,
    lineheight = 1.0
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = 1, labels = "freq") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    panel.grid    = element_blank(),
    axis.text.x   = element_text(color = "grey20"),
    plot.margin   = margin(t = 5, r = 5, b = 5, l = 0)
  )

left_panel <- p_grid + p_freq + plot_layout(widths = c(5, 1.2))

# ---- Forest plot (right panel) ----
forest_df <- df %>%
  mutate(
    label_txt = ifelse(
      Haplotype == "AGCCA",
      sprintf("OR %.2f (%.2f–%.2f)", OR, OR.lower, OR.upper),
      sprintf("OR %.2f (%.2f–%.2f), p=%.3g",
              OR, OR.lower, OR.upper, P_value)
    ),
    color_group = case_when(
      Haplotype %in% c("AGTTG", "AATTG") ~ "highlight",
      TRUE                               ~ "normal"
    ),
    Haplotype = factor(Haplotype, levels = levels(df$Haplotype))
  )

x_min <- min(forest_df$OR.lower, na.rm = TRUE)
x_max <- min(3, max(forest_df$OR.upper, na.rm = TRUE) * 1.3)

col_map <- c(
  highlight = "firebrick3",
  normal    = "black"
)

p_forest <- ggplot(forest_df, aes(x = OR, y = Haplotype)) +
  geom_point(aes(color = color_group), size = 3) +
  geom_errorbarh(
    aes(xmin = OR.lower, xmax = OR.upper, color = color_group),
    height   = 0.25,
    linewidth = 0.6
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_text(
    aes(label = label_txt),
    vjust = -1.1,
    size  = 3.3,
    color = "black"
  ) +
  scale_color_manual(values = col_map) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = c(0.5, 1, 1.5, 2, 3),
    labels = c("0.5", "1", "1.5", "2", "3")
  ) +
  labs(
    x     = "Odds Ratio (log scale, 95% CI)",
    y     = NULL,
    title = "SNCA haplotype association"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    axis.text.y        = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(t = 15, r = 30, b = 10, l = 10)
  )

combo <- left_panel | p_forest

ggsave(
  combo_plot_file,
  combo,
  width  = 13,
  height = 6,
  dpi    = 600
)

message("SNCA haplotype pipeline completed:
  - Haplotype table plot:              ", hap_table_plot, "
  - LD heatmap:                        ", ld_heatmap_file, "
  - Unadjusted haplo.cc results:       ", hap_unadj_out, "
  - Adjusted haplo.cc results:         ", hap_adj_out, "
  - Adjusted OR table (all):           ", hap_or_table, "
  - Adjusted OR table (common only):   ", hap_or_common, "
  - Combined figure (grid + forest):   ", combo_plot_file)
