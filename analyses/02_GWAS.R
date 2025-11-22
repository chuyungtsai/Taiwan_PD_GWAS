# 2 GWAS
# after imputation
library(glue)

# run the GWAS in separate chromosomes
for (i in seq(1,22)){
  print(i)
  c0 <- glue('chr{i}_TWB_imputed.fam')
  if (file.exists(c0)){
    print('exists!')
    c1 <- glue("plink2 --bfile chr{i}_imputed --glm 'hide-covar' cols=+a1count,+a1freq --covar 20241119_covariate_PC15.txt --covar-variance-standardize --ci 0.95 --out 20241228_imputed_chr{i}")
    system(c1)
  }
}

# read the results
result_chr1 <- read.table('imputed_chr1.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr2 <- read.table('imputed_chr2.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr3 <- read.table('imputed_chr3.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr4 <- read.table('imputed_chr4.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr5 <- read.table('imputed_chr5.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr6 <- read.table('imputed_chr6.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr7 <- read.table('imputed_chr7.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr8 <- read.table('imputed_chr8.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr9 <- read.table('imputed_chr9.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr10 <- read.table('imputed_chr10.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr11 <- read.table('imputed_chr11.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr12 <- read.table('imputed_chr12.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr13 <- read.table('imputed_chr13.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr14 <- read.table('imputed_chr14.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr15 <- read.table('imputed_chr15.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr16 <- read.table('imputed_chr16.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr17 <- read.table('imputed_chr17.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr18 <- read.table('imputed_chr18.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr19 <- read.table('imputed_chr19.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr20 <- read.table('imputed_chr20.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr21 <- read.table('imputed_chr21.PHENO1.glm.logistic.hybrid', head=FALSE)
result_chr22 <- read.table('imputed_chr22.PHENO1.glm.logistic.hybrid', head=FALSE)

# merge the results
result_imputed_2 <- rbind(result_chr1, result_chr2, result_chr3, result_chr4, result_chr5, result_chr6, result_chr7, result_chr8, result_chr9, result_chr10, result_chr11, result_chr12, result_chr13, result_chr14, result_chr15, result_chr16, result_chr17, result_chr18, result_chr19, result_chr20, result_chr21, result_chr22)

plink_header_2 <- c('CHR', 'BP', 'ID', 'REF', 'ALT', 'A1', 'A1_CT', 'A1_FREQ', 'FIRTH', 'TEST', 'OBS_CT', 'OR', '[LOG(OR)_]SE', 'L95', 'U95', 'STAT', 'P', 'ERRCODE')
colnames(result_imputed_2) <- plink_header_2

### qqplot and genomic control
library(ggplot2)

# Suppose you already have your GWAS p-values
pvals <- result_imputed_2$P
pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]

# Compute λGC
chi2 <- qchisq(1 - pvals, df = 1)
lambda_gc <- median(chi2, na.rm = TRUE) / 0.4559364
lambda_gc 

# Create QQ plot
qq_df <- data.frame(
  exp = -log10((1:length(pvals)) / (length(pvals) + 1)),
  obs = -log10(sort(pvals))
)

qq_plot <- ggplot(qq_df, aes(x = exp, y = obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray40") +
  geom_point(size = 0.8, color = "navy", alpha = 0.6) +
  labs(
    title = sprintf("QQ Plot of GWAS P-values (λGC = %.3f)", lambda_gc),
    x = "Expected -log10(P)",
    y = "Observed -log10(P)"
  ) +
  theme_bw(base_size = 12)

# Output to JPG
jpeg("QQplot_lambdaGC.jpg", width = 1200, height = 1200, res = 150)
print(qq_plot)
dev.off()

