## ================================================================
## 08_LRRK2_analysis.R
##
## Taiwan PD GWAS – LRRK2 Asian variants (G2385R, R1628P)
## 1) Extract LRRK2 genotypes from PLINK (0/1/2)
## 2) Case–control association for LRRK2 dosage (0 vs 1 vs 2)
## 3) Age at onset (AAO) comparison across LRRK2 genotypes
##
## Required inputs:
##   - PLINK bed/bim/fam:   LRRK2_region_5_variants
##   - Covariate file:      20241119_covariate_PC15.txt
##   - AAO file (Excel):    20250907_PD_AAO_total.xlsx
##
## Outputs (in working directory):
##   - 20251013_LRRK2_variants.csv        (IID, PHENOTYPE, G2385R, R1628P, LRRK2)
##   - 20251013_LRRK2_count.csv           (2x3 contingency table)
##   - 20251013_LRRK2_OR.csv              (ORs: 1 vs 0, 2 vs 0, 2 vs 1)
##   - 20251013_LRRK2_AAO_summary.csv     (n, mean, SD AAO by LRRK2)
##   - 20251013_LRRK2_AAO_tests.txt       (ANOVA / Welch / Kruskal–Wallis / pairwise tests)
## ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(broom)
  library(readxl)
})

plink_bin <- "plink19"

## ================================================================
## 1. EXTRACT LRRK2 GENOTYPES (G2385R, R1628P) FROM PLINK
## ================================================================

# PLINK recode (additive dosage, .raw format)
system(paste(
  plink_bin,
  "--bfile LRRK2_region_5_variants",
  "--recode",
  "--out LRRK2_5_variants"
))

# Read PLINK .raw (FID IID PAT MAT SEX PHENOTYPE SNP1_A SNP2_A ...)
LRRK2_RAW <- fread("LRRK2_5_variants.raw")

# Select G2385R and R1628P dosages (A1 allele count)
# G2385R: chr12:40363526:G:A  → column name "chr12:40363526:G:A_A"
# R1628P: chr12:40320043:G:C  → column name "chr12:40320043:G:C_C"
LRRK2_2var <- LRRK2_RAW %>%
  select(IID, PHENOTYPE,
         `chr12:40363526:G:A_A`,
         `chr12:40320043:G:C_C`) %>%
  rename(
    G2385R = `chr12:40363526:G:A_A`,
    R1628P = `chr12:40320043:G:C_C`
  ) %>%
  mutate(
    # LRRK2 dosage = sum of the two risk-allele dosages (0/1/2 each)
    LRRK2 = G2385R + R1628P
  )

# Quick check: counts by total LRRK2 dosage
message("LRRK2 dosage distribution (0/1/2):")
print(table(LRRK2_2var$LRRK2))

# Save variant table for downstream use (and for 1000G scripts)
fwrite(LRRK2_2var, "20251013_LRRK2_variants.csv", quote = FALSE)

# 2x3 contingency table: PHENOTYPE vs LRRK2 dosage
xtab <- with(LRRK2_2var, table(PHENOTYPE, LRRK2))
write.csv(xtab, "20251013_LRRK2_count.csv", quote = FALSE, row.names = FALSE)

## ================================================================
## 2. CASE–CONTROL ASSOCIATION: LRRK2 0/1/2 (ADJUSTED)
## ================================================================

# 2.1 Load covariates (sex, age, PC1–PC15)
# This should come from 01_QC_PCA_covariates.R
covar_file <- "20241119_covariate_PC15.txt"
PC15 <- fread(covar_file, header = TRUE)

# Expect columns: FID, IID, sex, age, PC1..PC15
# Merge LRRK2 genotypes with covariates
LRRK2_covar <- LRRK2_2var %>%
  left_join(PC15, by = "IID")

# PHENOTYPE in PLINK is typically 1 = control, 2 = case
LRRK2_covar <- LRRK2_covar %>%
  mutate(PD = PHENOTYPE - 1L)

df <- LRRK2_covar

# 2.2 Define outcome Y (binary 0/1)
outcome_candidates <- c("PD")
outcome_col <- outcome_candidates[outcome_candidates %in% names(df)][1]
if (is.na(outcome_col)) {
  stop("No binary outcome column found. Expected 'PD'.")
}
y_raw <- df[[outcome_col]]
if (is.numeric(y_raw)) {
  Y <- ifelse(y_raw > 0, 1L, 0L)
} else {
  y_chr <- tolower(as.character(y_raw))
  Y <- ifelse(y_chr %in% c("case","1","pd","patient","yes","true"), 1L, 0L)
}
df$Y <- Y

# 2.3 Predictor: LRRK2 0/1/2 as factor (reference = 0)
if (!"LRRK2" %in% names(df)) {
  stop("Missing column 'LRRK2' (expected values 0/1/2).")
}
df$LRRK2_f <- factor(df$LRRK2, levels = c(0,1,2), labels = c("0","1","2"))

# 2.4 Covariates: sex, age, PCs
sex_candidates <- c("SEX","Sex","sex","GENDER","Gender","gender")
age_candidates <- c("AGE","Age","age","AGE_AT_ENROLLMENT","age_at_enrollment",
                    "AGE_AT_ONSET","AAO","age_at_onset")

sex_col <- sex_candidates[sex_candidates %in% names(df)][1]
age_col <- age_candidates[age_candidates %in% names(df)][1]
pcs_present <- paste0("PC", 1:15)
pcs_present <- pcs_present[pcs_present %in% names(df)]

rhs_terms <- c(
  "LRRK2_f",
  if (!is.na(sex_col)) sex_col,
  if (!is.na(age_col)) age_col,
  pcs_present
)

form <- as.formula(paste("Y ~", paste(rhs_terms, collapse = " + ")))
message("Logistic regression formula:")
print(form)

# 2.5 Fit logistic regression
fit <- glm(form, family = binomial(), data = df)

# 2.6 Wald OR/CI helper for arbitrary contrasts
wald_or <- function(L, fit, label) {
  b <- coef(fit); V <- vcov(fit)
  est <- sum(L * b)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))
  data.frame(
    contrast = label,
    OR       = exp(est),
    CI_low   = exp(est - 1.96 * se),
    CI_high  = exp(est + 1.96 * se),
    p_value  = p,
    row.names = NULL
  )
}

coef_names <- names(coef(fit))
b1 <- grep("^LRRK2_f1$", coef_names, value = TRUE)  # LRRK2=1 vs 0
b2 <- grep("^LRRK2_f2$", coef_names, value = TRUE)  # LRRK2=2 vs 0
if (length(b1) != 1 || length(b2) != 1) {
  stop("Could not find LRRK2 coefficients (levels 1/2 vs 0). Check that LRRK2 has levels 0,1,2 and all levels are present.")
}

L_10 <- setNames(rep(0, length(coef(fit))), coef_names); L_10[b1] <- 1                 # 1 vs 0
L_20 <- setNames(rep(0, length(coef(fit))), coef_names); L_20[b2] <- 1                 # 2 vs 0
L_21 <- setNames(rep(0, length(coef(fit))), coef_names); L_21[b2] <- 1; L_21[b1] <- -1 # 2 vs 1

res_OR <- dplyr::bind_rows(
  wald_or(L_10, fit, "LRRK2: 1 vs 0"),
  wald_or(L_20, fit, "LRRK2: 2 vs 0"),
  wald_or(L_21, fit, "LRRK2: 2 vs 1")
)

print(res_OR, row.names = FALSE)
write.csv(res_OR, "20251013_LRRK2_OR.csv", quote = FALSE, row.names = FALSE)

## ================================================================
## 3. AGE AT ONSET (AAO) ANALYSIS BY LRRK2 GENOTYPE
## ================================================================

# 3.1 Load AAO data
AAO_file <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/patient_list/Total_PD_list/20250907_PD_AAO_total.xlsx"
AAO_df <- read_excel(AAO_file)

# Expect column "IID" and "Onset_age" at least
LRRK2_AAO <- AAO_df %>%
  left_join(LRRK2_2var, by = "IID")

# Optional subsets (compound hetero or homozygotes)
LRRK2_com_hetero   <- LRRK2_AAO %>% filter(G2385R == 1, R1628P == 1)
LRRK2_G2385R_homo  <- LRRK2_AAO %>% filter(G2385R == 2)
LRRK2_R1628P_homo  <- LRRK2_AAO %>% filter(R1628P == 2)

# 3.2 Prepare dataset for global LRRK2 0/1/2 comparison
dat2 <- LRRK2_AAO %>%
  mutate(
    LRRK2 = factor(LRRK2, levels = c(0, 1, 2), labels = c("0", "1", "2"))
  ) %>%
  filter(!is.na(Onset_age), !is.na(LRRK2))

# 3.3 Mean ± SD AAO by LRRK2 group
summary_tbl <- dat2 %>%
  group_by(LRRK2) %>%
  summarise(
    n        = n(),
    mean_AAO = mean(Onset_age),
    sd_AAO   = sd(Onset_age),
    .groups  = "drop"
  )

summary_fmt <- summary_tbl %>%
  mutate(mean_sd = sprintf("%.2f ± %.2f", mean_AAO, sd_AAO)) %>%
  select(LRRK2, n, mean_sd)

message("\nAAO summary by LRRK2 (mean ± SD):")
print(summary_fmt)

write.csv(summary_tbl, "20251013_LRRK2_AAO_summary.csv", quote = FALSE, row.names = FALSE)

# 3.4 Global group comparison
sink("20251013_LRRK2_AAO_tests.txt")

cat("=== Age at onset (AAO) – LRRK2 0/1/2 ===\n\n")

cat("Summary (n, mean, SD):\n")
print(summary_tbl)
cat("\n")

# One-way ANOVA
fit_aov <- aov(Onset_age ~ LRRK2, data = dat2)
cat("\n[One-way ANOVA]\n")
print(broom::tidy(fit_aov))

# Welch ANOVA
welch <- oneway.test(Onset_age ~ LRRK2, data = dat2)
cat("\n[Welch ANOVA]\n")
print(welch)

# Kruskal–Wallis
kw <- kruskal.test(Onset_age ~ LRRK2, data = dat2)
cat("\n[Kruskal–Wallis test]\n")
print(kw)

# 3.5 Pairwise comparisons
# Tukey HSD
tukey <- TukeyHSD(fit_aov, "LRRK2")
cat("\n[Tukey HSD (parametric, equal variances)]\n")
print(tukey)

# Welch pairwise t-tests (Holm-adjusted)
pair_t <- pairwise.t.test(dat2$Onset_age, dat2$LRRK2,
                          p.adjust.method = "holm", pool.sd = FALSE)
cat("\n[Pairwise Welch t-tests (Holm-adjusted)]\n")
print(pair_t)

# Pairwise Wilcoxon (Holm-adjusted)
pair_w <- pairwise.wilcox.test(dat2$Onset_age, dat2$LRRK2,
                               p.adjust.method = "holm",
                               exact = FALSE)
cat("\n[Pairwise Wilcoxon tests (Holm-adjusted)]\n")
print(pair_w)

sink()

message("LRRK2 analysis completed:
  - LRRK2 variant table:          20251013_LRRK2_variants.csv
  - LRRK2 dosage 2x3 table:       20251013_LRRK2_count.csv
  - LRRK2 ORs (0/1/2):            20251013_LRRK2_OR.csv
  - AAO summary by LRRK2:         20251013_LRRK2_AAO_summary.csv
  - AAO global & pairwise tests:  20251013_LRRK2_AAO_tests.txt")
