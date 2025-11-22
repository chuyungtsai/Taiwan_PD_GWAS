## 04a_SNCA_conditional.R
## Conditional analyses at the SNCA locus and LocusZoom-style plots
## - Model 1: conditional on Taiwan lead variant (rs5860181, chr4:89759478:A:ATGCATATT)
## - Model 2: conditional on rs356182 (Nalls 2019, chr4:89704960:G:A)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(glue)
})

## ------------------------------------------------------------------
## 0. Paths and configuration
## ------------------------------------------------------------------
proj_dir   <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PLINK_analysis/20241225_TWB_imputed_analysis"
lz_dir     <- "/Users/chuyungtsai/Desktop/Bioinfo/NTU_GWAS/PLINK_analysis/20251005_locusZoom"
gtf_path   <- "/Users/chuyungtsai/Desktop/Bioinfo/hg38/Homo_sapiens.GRCh38.109.gtf.gz"

plink2_bin <- "plink2"
bfile_chr4 <- file.path(proj_dir, "chr4_TWB_imputed")
covar_file <- file.path(proj_dir, "20241119_covariate_PC15.txt")
keep_ids   <- file.path(proj_dir, "20251010_control_list.txt")   # FID IID of controls (from PCA script)

dir.create(lz_dir, showWarnings = FALSE)
setwd(proj_dir)

## ------------------------------------------------------------------
## 1. Helpers
## ------------------------------------------------------------------

## PLINK 2.0 hybrid logistic header
plink_header_2 <- c(
  "CHR", "BP", "ID", "REF", "ALT", "A1",
  "A1_CT", "A1_FREQ", "FIRTH", "TEST",
  "OBS_CT", "OR", "[LOG(OR)_]SE", "L95", "U95",
  "STAT", "P", "ERRCODE"
)

read_plink_glm <- function(path) {
  df <- fread(path, data.table = FALSE)
  if (ncol(df) != length(plink_header_2)) {
    stop("Unexpected number of columns in ", basename(path),
         " (", ncol(df), " vs ", length(plink_header_2), ").")
  }
  colnames(df) <- plink_header_2
  df
}

## Load GTF once; locuszoom_preloaded_reLD will do regional filtering
gtf <- fread(
  gtf_path,
  sep = "\t",
  header = FALSE,
  quote = "",
  na.strings = c("", "NA"),
  data.table = FALSE
)
colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attr"
)

## ------------------------------------------------------------------
## 2. Conditional on Taiwan lead rs5860181 (chr4:89759478:A:ATGCATATT)
##    → second hit rs1442145 (chr4:89815772:G:A)
## ------------------------------------------------------------------

lead_tw      <- "chr4:89759478:A:ATGCATATT"   # Taiwan SNCA lead
glm_prefix_1 <- file.path(proj_dir, "20251005_TWB_SNCA")
out_prefix_1 <- file.path(lz_dir,  "20251011_chr4_SNCA_conditional_TWlead")

if (!file.exists(paste0(glm_prefix_1, ".PHENO1.glm.logistic.hybrid"))) {
  cmd1 <- glue(
    "{plink2_bin} --bfile {bfile_chr4} ",
    "--glm 'hide-covar' cols=+a1count,+a1freq ",
    "--covar {covar_file} --covar-variance-standardize ",
    "--condition {lead_tw} --ci 0.95 ",
    "--out {glm_prefix_1}"
  )
  cat("\n[INFO] Running conditional GWAS on SNCA (condition on Taiwan lead):\n", cmd1, "\n")
  system(cmd1)
}

SNCA_cond1 <- read_plink_glm(paste0(glm_prefix_1, ".PHENO1.glm.logistic.hybrid"))

## LocusZoom-style plot around the conditional secondary signal rs1442145
locuszoom_preloaded_reLD(
  lead_var       = "chr4:89815772:G:A",   # rs1442145
  chr            = 4,
  lead_bp        = 89815772,
  out_prefix     = out_prefix_1,
  ss_df          = SNCA_cond1,
  gtf_df         = gtf,
  source_bfile   = bfile_chr4,
  keep_ids_path  = keep_ids,
  window_kb      = 300,
  forced_labels  = c(
    "chr4:89815772:G:A" = "rs1442145"
  )
)

## ------------------------------------------------------------------
## 3. Conditional on rs356182 (chr4:89704960:G:A, Nalls 2019)
##    → second hit rs3775437 (chr4:89784792:G:A)
## ------------------------------------------------------------------

lead_rs356182 <- "chr4:89704960:G:A"
glm_prefix_2  <- file.path(proj_dir, "20251104_TWB_SNCA_rs356182")
out_prefix_2  <- file.path(lz_dir,  "20251104_chr4_SNCA_conditional_rs356182")

if (!file.exists(paste0(glm_prefix_2, ".PHENO1.glm.logistic.hybrid"))) {
  cmd2 <- glue(
    "{plink2_bin} --bfile {bfile_chr4} ",
    "--glm 'hide-covar' cols=+a1count,+a1freq ",
    "--covar {covar_file} --covar-variance-standardize ",
    "--condition {lead_rs356182} --ci 0.95 ",
    "--out {glm_prefix_2}"
  )
  cat("\n[INFO] Running conditional GWAS on SNCA (condition on rs356182):\n", cmd2, "\n")
  system(cmd2)
}

SNCA_cond2 <- read_plink_glm(paste0(glm_prefix_2, ".PHENO1.glm.logistic.hybrid"))

## LocusZoom-style plot around rs3775437 (secondary signal) plus TW lead label
locuszoom_preloaded_reLD(
  lead_var       = "chr4:89784792:G:A",   # rs3775437
  chr            = 4,
  lead_bp        = 89784792,
  out_prefix     = out_prefix_2,
  ss_df          = SNCA_cond2,
  gtf_df         = gtf,
  source_bfile   = bfile_chr4,
  keep_ids_path  = keep_ids,
  window_kb      = 300,
  forced_labels  = c(
    "chr4:89784792:G:A"         = "rs3775437",
    "chr4:89759478:A:ATGCATATT" = "rs5860181"   # Taiwan lead
  )
)

cat("\n[DONE] Conditional SNCA analyses and plots written under:\n  ", lz_dir, "\n")
