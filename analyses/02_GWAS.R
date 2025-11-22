#!/usr/bin/env Rscript

## ================================================================
## GWAS on imputed data (chr1–22) + QQ plot & genomic control
## - Runs PLINK2 per chromosome using a common covariate file
## - Merges .glm.logistic.hybrid outputs
## - Computes λGC and draws QQ plot
## ================================================================

suppressPackageStartupMessages({
  library(glue)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

## ---------------- USER CONFIGURATION ----------------------------

# PLINK2 binary (assumed in PATH; change if needed)
plink2_bin <- "plink2"

# Covariate file (PC1–PC15 etc.)
covar_file <- "20241119_covariate_PC15.txt"

# Base prefix for PLINK2 output (per chromosome)
out_prefix <- "20241228_imputed_chr"

# Genotype prefix pattern (per chromosome)
# Expecting files like: chr1_TWB_imputed.{bed,bim,fam}
bfile_pattern <- "chr{chr}B_imputed"

# Output QQ plot file
qq_outfile <- "QQplot_lambdaGC.png"

# Chromosomes to run
chroms <- 1:22

## ================================================================
## 1. RUN PLINK2 GWAS PER CHROMOSOME
## ================================================================

for (chr in chroms) {
  fam_path <- glue("{bfile_pattern}.fam", chr = chr)
  if (!file.exists(fam_path)) {
    message("Skip chr", chr, ": FAM not found (", fam_path, ")")
    next
  }
  
  message("Running PLINK2 on chr", chr, "...")
  
  cmd <- glue(
    "{plink2_bin} ",
    "--bfile {bfile_pattern} ",
    "--glm hide-covar cols=+a1count,+a1freq ",
    "--covar {covar_file} ",
    "--covar-variance-standardize ",
    "--ci 0.95 ",
    "--out {out_prefix}{chr}",
    chr = chr
  )
  
  system(cmd)
}

## ================================================================
## 2. READ & MERGE PLINK2 OUTPUTS
## ================================================================

gwas_list <- list()

for (chr in chroms) {
  f_glm <- glue("{out_prefix}{chr}.PHENO1.glm.logistic.hybrid")
  if (!file.exists(f_glm)) {
    message("Result not found for chr", chr, ": ", f_glm)
    next
  }
  
  dt <- fread(f_glm)
  
  # Keep additive model only; drop covariate rows etc.
  if ("TEST" %in% names(dt)) {
    dt <- dt[TEST == "ADD"]
  }
  
  # Basic cleaning: drop rows with error codes if ERRCODE is present
  if ("ERRCODE" %in% names(dt)) {
    dt <- dt[is.na(ERRCODE) | ERRCODE == "NO_ERROR"]
  }
  
  # Standardize chromosome / position column names if using PLINK2 defaults
  if ("#CHROM" %in% names(dt)) setnames(dt, "#CHROM", "CHR")
  if ("POS"    %in% names(dt)) setnames(dt, "POS",    "BP")
  
  # Add chr if missing (some PLINK versions have CHR, some #CHROM)
  if (!"CHR" %in% names(dt)) {
    dt[, CHR := chr]
  }
  
  gwas_list[[as.character(chr)]] <- dt
}

if (length(gwas_list) == 0L) {
  stop("No GWAS result files were found. Check paths and PLINK output.")
}

gwas_all <- rbindlist(gwas_list, use.names = TRUE, fill = TRUE)

# Optional: reorder columns (if they exist) for clarity
col_order <- c(
  "CHR", "BP", "ID", "REF", "ALT", "A1",
  "A1_CT", "A1_FREQ", "OBS_CT",
  "OR", "SE", "L95", "U95", "Z_STAT", "STAT", "P",
  "TEST", "FIRTH", "ERRCODE"
)
col_order <- intersect(col_order, names(gwas_all))
setcolorder(gwas_all, c(col_order, setdiff(names(gwas_all), col_order)))

# Save merged GWAS table (tab-delimited)
fwrite(gwas_all, "20241228_GWAS_imputed_allchr.tsv", sep = "\t")

## ================================================================
## 3. GENOMIC CONTROL (λGC) & QQ PLOT
## ================================================================

# Extract valid p-values
pvals <- gwas_all$P
pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]

if (length(pvals) < 1000) {
  warning("Fewer than 1000 P-values after filtering; QQ plot may be unstable.")
}

# λGC: median(χ²)/median(χ² under null)
chi2       <- qchisq(1 - pvals, df = 1)
lambda_gc  <- median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1)

message(sprintf("Genomic control λGC = %.3f", lambda_gc))

# QQ data
p_sorted <- sort(pvals)
exp_vals <- -log10((seq_along(p_sorted)) / (length(p_sorted) + 1))
obs_vals <- -log10(p_sorted)

qq_df <- data.frame(
  exp = exp_vals,
  obs = obs_vals
)

qq_plot <- ggplot(qq_df, aes(x = exp, y = obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray40") +
  geom_point(size = 0.8, color = "navy", alpha = 0.6) +
  labs(
    title = sprintf("QQ Plot of GWAS P-values (λGC = %.3f)", lambda_gc),
    x     = "Expected -log10(P)",
    y     = "Observed -log10(P)"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = qq_outfile,
  plot     = qq_plot,
  width    = 6,
  height   = 6,
  dpi      = 300
)

message("Done:
  - Merged GWAS: 20241228_GWAS_imputed_allchr.tsv
  - QQ plot:     ", qq_outfile,
        "\n  - λGC:       ", sprintf('%.3f', lambda_gc))
