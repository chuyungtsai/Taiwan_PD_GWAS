# 1. QC and PCA
## ================================================================
## Taiwan PD GWAS – SNP & sample QC, PCA, and covariate file
## Starting point: 20241117_merged.{bed,bim,fam} (Genomestudio → PLINK)
## ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(gtsummary)
  library(writexl)
})

## ---------------- USER CONFIGURATION ----------------------------

plink_bin <- "plink19"

# Input raw merged genotype
raw_prefix <- "20241117_merged"

# SNP list of known underperforming GP2 markers
underperf_snp_file <- "underperforming_GP2_SNPs_name"

# QC prefixes
qc0_prefix <- "20241117_merged_2"   # after removing underperforming SNPs
qc1_prefix <- "20241117_merged_QC_1" # after geno/mind/maf/HWE
qc2_prefix <- "20241117_merged_QC_2" # after differential missingness filters
qc3_prefix <- "20241117_merged_QC_3" # after removing IBD / curated list

# Manually curated list of samples to remove (duplicates/relatives etc.)
# Fill this file with FID IID pairs before running, or leave empty.
ibd_remove_file <- "20241003_remove.txt"

# PCA/IBD/pruning output prefixes
prune_prefix <- "20241117_prune"
ibd_prefix   <- "20241117_IBD"
pca_prefix   <- "20241117_pca"

# Phenotype / covariate info
# total_profile must exist in the session:
#   data.frame with at least: IID, sex, age, PD (0/1)
# e.g.:
# total_profile <- read.table("total_profile.txt", header = TRUE)
# Here we assume it's loaded before sourcing this script.

covar_outfile <- "20241119_covariate_PC15.txt"  # used later by PLINK2 --covar
table1_xlsx   <- "20241117_Table1_demographics.xlsx"

## ================================================================
## 1. BASIC SNP & SAMPLE QC
## ================================================================

# 1.1 Remove poor performing SNPs
system(glue::glue(
  "{plink_bin} --bfile {raw_prefix} ",
  "--exclude {underperf_snp_file} ",
  "--make-bed --out {qc0_prefix}"
))

# Optionally remove original PLINK files after QC stage 0
system(glue::glue("rm {raw_prefix}.bed {raw_prefix}.bim {raw_prefix}.fam"), ignore.stdout = TRUE, ignore.stderr = TRUE)

# 1.2 Call rate, MAF, HWE
system(glue::glue(
  "{plink_bin} --bfile {qc0_prefix} ",
  "--geno 0.05 --mind 0.05 --maf 0.01 --hwe 1e-4 ",
  "--make-bed --out {qc1_prefix}"
))

## ================================================================
## 2. DIFFERENTIAL MISSINGNESS (CASE/CONTROL + HAPLOTYPE)
## ================================================================

# 2.1 test-missing (by case/control)
system(glue::glue(
  "{plink_bin} --bfile {qc1_prefix} ",
  "--test-missing --out 20241117_test_missing"
))

test_missing <- read.table("20241117_test_missing.missing", header = TRUE)
missing_SNP1 <- test_missing %>%
  filter(P < 1e-4) %>%
  distinct(SNP)

# 2.2 test-mishap (haplotype-based)
system(glue::glue(
  "{plink_bin} --bfile {qc1_prefix} ",
  "--test-mishap --out 20241117_test_mishap"
))

test_mishap <- read.table("20241117_test_mishap.missing.hap", header = TRUE)
missing_SNP2 <- test_mishap %>%
  filter(P < 1e-4) %>%
  distinct(SNP)

# 2.3 Combine SNPs to remove and write for PLINK --exclude
missing_SNP_all <- bind_rows(missing_SNP1, missing_SNP2) %>%
  distinct(SNP)

write.table(
  missing_SNP_all,
  file      = "20241117_missing_SNP.txt",
  quote     = FALSE,
  sep       = " ",
  row.names = FALSE,
  col.names = FALSE
)

# 2.4 Exclude differential missing SNPs
system(glue::glue(
  "{plink_bin} --bfile {qc1_prefix} ",
  "--exclude 20241117_missing_SNP.txt ",
  "--make-bed --out {qc2_prefix}"
))

# Optionally remove intermediate QC_1 files
system(glue::glue("rm {qc1_prefix}.bed {qc1_prefix}.bim {qc1_prefix}.fam"),
       ignore.stdout = TRUE, ignore.stderr = TRUE)

## ================================================================
## 3. SAMPLE-LEVEL IBD & FINAL QC DATASET
## ================================================================

# Ensure a remove file exists (may be empty)
if (!file.exists(ibd_remove_file)) {
  writeLines(character(0), ibd_remove_file)
}

# 3.1 Remove pre-identified duplicates/relatives (if any)
system(glue::glue(
  "{plink_bin} --bfile {qc2_prefix} ",
  "--remove {ibd_remove_file} ",
  "--make-bed --out {qc3_prefix}"
))

# 3.2 LD-pruning for IBD/PCA
system(glue::glue(
  "{plink_bin} --bfile {qc3_prefix} ",
  "--indep-pairwise 50 5 0.2 ",
  "--out {prune_prefix}"
))

# 3.3 IBD on pruned SNPs
system(glue::glue(
  "{plink_bin} --bfile {qc3_prefix} ",
  "--extract {prune_prefix}.prune.in ",
  "--genome --out {ibd_prefix}"
))

IBD <- read.table(paste0(ibd_prefix, ".genome"), header = TRUE)

# Users can inspect IBD$PI_HAT offline and update ibd_remove_file if needed,
# then re-run from section 3.1.

## ================================================================
## 4. PCA ON FINAL QC DATASET
## ================================================================

# 4.1 PCA on pruned SNPs
system(glue::glue(
  "{plink_bin} --bfile {qc3_prefix} ",
  "--extract {prune_prefix}.prune.in ",
  "--pca 20 --out {pca_prefix}"
))

# 4.2 Load FAM + PCA eigenvalues/eigenvectors
fam <- read.table(paste0(qc3_prefix, ".fam"), header = FALSE)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "Sex_plink", "PD")

eigenval <- read.table(paste0(pca_prefix, ".eigenval"), header = FALSE)
colnames(eigenval) <- "eigval"

eigenvec <- read.table(paste0(pca_prefix, ".eigenvec"), header = FALSE)
colnames(eigenvec) <- c(
  "FID", "IID",
  paste0("PC", 1:20)
)

# 4.3 Scree plot (variance explained)
var_explained <- eigenval$eigval^2 / sum(eigenval$eigval^2)
scree_df <- data.frame(
  PC   = seq_along(var_explained),
  Var  = var_explained
)

p_scree <- ggplot(scree_df[1:20, ], aes(x = PC, y = Var)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Principal Component",
    y = "Variance explained",
    title = "PCA scree plot (first 20 PCs)"
  ) +
  theme_bw(base_size = 12)

ggsave("20241117_PCA_scree_PC1_PC20.png", p_scree, width = 5, height = 4, dpi = 300)

# 4.4 Merge PCA with FAM
pca <- fam %>%
  left_join(eigenvec, by = c("FID", "IID"))

# Optional: hospital code from FID (if FID encodes hospital)
# Adjust as appropriate for your data.
pca <- pca %>%
  mutate(hospital = as.factor(FID))

# 4.5 PC1 vs PC2 plot colored by PD status
p_pca_pd <- ggplot(pca, aes(x = PC1, y = PC2, color = factor(PD))) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("0" = "black", "1" = "red"),
                     name  = "PD status",
                     labels = c("0" = "Control", "1" = "PD")) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA: PC1 vs PC2 by PD status"
  ) +
  theme_bw(base_size = 12)

ggsave("20241117_PCA_PC1_PC2_PD.png", p_pca_pd, width = 5, height = 4, dpi = 300)

# 4.6 Controls only, colored by hospital (optional QC figure)
controls <- pca %>% filter(PD == 1)
p_pca_ctrl_hosp <- ggplot(controls, aes(x = PC1, y = PC2, color = hospital)) +
  geom_point(alpha = 0.8, size = 1.2) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA: controls only (PC1 vs PC2, colored by hospital)"
  ) +
  theme_bw(base_size = 12)

ggsave("20241117_PCA_PC1_PC2_controls_hospital.png",
       p_pca_ctrl_hosp, width = 5, height = 4, dpi = 300)

## ================================================================
## 5. BUILD COVARIATE FILE (PC1–PC15, AGE, SEX)
## ================================================================

# total_profile must already be loaded in the environment.
# It should contain at least:
#   IID: sample ID (matches PLINK IID)
#   sex: 0/1 or 1/2 coding (will be used as covariate)
#   age: numeric age
#   PD : case/control indicator (0/1)

if (!exists("total_profile")) {
  stop("total_profile is not loaded. Please load it before running this script.")
}

pca_sex_age <- pca %>%
  left_join(total_profile, by = "IID")

# Ensure age is numeric
pca_sex_age$age <- suppressWarnings(as.numeric(pca_sex_age$age))

# Covariate file: FID, IID, sex, age, PC1–PC15
PC15 <- pca_sex_age %>%
  select(FID, IID, sex, age, PC1:PC15)

# PLINK2 --covar expects header unless you add 'noheader';
# here we keep column names.
write.table(
  PC15,
  file      = covar_outfile,
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE,
  col.names = TRUE
)

## ================================================================
## 6. TABLE 1 – DEMOGRAPHICS (AGE, SEX BY PD STATUS)
## ================================================================

# Filter: non-missing age, valid sex coding
table1_data <- pca_sex_age %>%
  filter(!is.na(age)) %>%
  filter(!is.na(sex))   # adjust or refine (e.g. sex != 2) based on coding

table1 <- table1_data %>%
  tbl_summary(
    include = c(age, sex),
    by      = PD,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 2
  ) %>%
  add_p() %>%
  add_overall() %>%
  bold_labels()

# Export Table 1 to Excel as plain table
table1_tibble <- gtsummary::as_tibble(table1)
write_xlsx(list("Table1_demographics" = table1_tibble), table1_xlsx)

message("QC/PCA/covariate pipeline completed:
  - Final QC dataset:               ", qc3_prefix, ".bed/.bim/.fam
  - Pruned SNP list:                ", prune_prefix, ".prune.in
  - IBD results:                    ", ibd_prefix, ".genome
  - PCA eigenvalues/eigenvectors:   ", pca_prefix, ".eigenval/.eigenvec
  - Covariate file for PLINK2:      ", covar_outfile, "
  - Demographic Table 1 (Excel):    ", table1_xlsx, "
  - PCA figures:                    20241117_PCA_*.png")
