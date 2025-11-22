# Taiwan PD GWAS Script
# 1. Quality control from plink file (already converted from idat using GenomeStudio)

library(dplyr)
library(readxl)
library(ggplot2)
library(data.table)
library(gtsummary)
library(writexl) # write_xlsx
library(gtsummary)

# basic QC
# remove poor performing SNP first!!
# file: underperforming_GP2_SNPs_name
system('plink19 --bfile 20241117_merged --exclude underperforming_GP2_SNPs_name --make-bed --out 20241117_merged_2')
system('rm 20241117_merged.*')

system('plink19 --bfile 20241117_merged_2 --geno 0.05 --mind 0.05 --maf 0.01 --hwe 1e-4 --make-bed --out 20241117_merged_QC_1')

#Variant Pruning. Missingness by: 
# 1. case/control (1e-4) 
# 2. haplotype (1e-4)
# test-missing
# test-mishap

system('plink19 --bfile 20241117_merged_QC_1 --test-missing --out 20241117_test_missing')

test_missing <- read.table('20241117_test_missing.missing', header = TRUE)
missing_SNP <- test_missing %>% filter(P<1e-4) %>% select(SNP) %>% unique()
write_for_plink(missing_SNP, '20241117_missing_SNP.txt')

system('plink19 --bfile 20241117_merged_QC_1 --test-mishap --out 20241117_test_missing') 
test_mishap <- read.table('20241117_test_missing.missing.hap', header = TRUE)
missing_SNP_2 <- test_mishap %>%filter(P<1e-4) %>% select(SNP) %>% unique()

write.table(missing_SNP_2,  '20241117_missing_SNP.txt', append=TRUE, quote=FALSE, sep=' ', col.names = FALSE, row.names=FALSE)

# exclude test-missing and test-mishap
system('plink19 --bfile 20241117_merged_QC_1 --exclude 20241117_missing_SNP.txt --make-bed --out 20241117_merged_QC_2')

system('rm 202410117_merged_QC_1.*')

# IBD
txt <- c() # sample ID to be removed.
writeLines(txt, "20241003_remove.txt")

system('plink19 --bfile 20241117_merged_QC_2 --remove 20241003_remove.txt --make-bed --out 20241117_merged_QC_3')

# then re-run pruning and IBD
system('plink19 --bfile 20241117_merged_QC_4 --indep-pairwise 50 5 0.2 --out 20241117_prune')

# IBD
system('plink19 --bfile 20241117_merged_QC_4 --extract 20241117_prune.prune.in --genome --out 20241117_IBD')

# read IBD result
IBD <- read.table('20241117_IBD.genome', header = TRUE)
unique(IBD$RT)

# selecte related samples
IBD %>% select(PI_HAT)%>% hist()
IBD %>% filter(0.6> PI_HAT & PI_HAT >0.2) %>% select(FID1, IID1, FID2, IID2, PI_HAT)

# identical
IBD %>% filter(PI_HAT >= 0.6) %>% select(FID1, IID1, FID2, IID2, PI_HAT)

# PCA
system('plink19 --bfile 20241117_merged_QC_4 --indep-pairwise 50 5 0.2 --out 20241117_prune')

system('plink19 --bfile 20241117_merged_QC_4 --extract 20241117_prune.prune.in --pca --out 20241117_pca')

final_pt <- read.table('20241117_merged_QC_4.fam')
View(final_pt)

# Scree plot
pca1 <- read.table('20241117_pca.eigenval')
plot(pca1$V1^2/sum(pca1$V1^2), ylab= 'Variance Explained', xlab='PC', main='First 20 PC')

pca2 <- read.table('20241117_pca.eigenvec')
colnames(pca2)= c('FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20')

colnames(final_pt)= c('FID', 'IID', 'F', 'M', 'Sex', 'PD')
View(final_pt)

pca <- merge(final_pt, pca2, by='IID')
pca$hospital <- pca$FID.x+1

# take a look 
pairs(pca[, 8:12])
plot(pca$PC1, pca$PC2, col=pca$PD, xlab='PC1', ylab='PC2')
plot(pca$PC1, pca$PC2, col=pca$hospital, xlab='PC1', ylab='PC2')

plot(pca$PC1, pca$PC2, col=pca$PD, xlab='PC1', ylab='PC2')
legend("topright", legend = c('CON', 'PD'), fill =c('black', 'red'))
legend("topright", legend = c('CON', 'PD'), col = unique(pca$PD),lty=1)

plot(pca$PC1, pca$PC2, col=pca$hospital, xlab='PC1', ylab='PC2')

# only look at controls
controls <- pca%>% filter(PD==1)
plot(controls$PC1, controls$PC2, col=controls$hospital, xlab='PC1', ylab='PC2')

# covariate files
# previous files
profile <- read.table('sex_age.txt')
colnames(profile) <- c('FID', 'IID', 'sex_plink', 'age', 'sex')

# merge profile
pca_sex_age <- merge(pca, total_profile, by='IID')
View(pca_sex_age)

View(final_pt$IID)
View(total_profile$IID)

PC15 <- pca_sex_age %>% select(FID.x, IID, sex, age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15)

write_for_plink(PC15, 'covariate_PC15.txt')

# table 1
pca_sex_age$age <- as.numeric(pca_sex_age$age)

table1 <- pca_sex_age %>% filter(!is.na(age)) %>% filter(sex !=2) %>% tbl_summary(include=c(age, sex), by = PD, statistic = list(
  all_continuous() ~ "{mean} ({sd})",
  all_categorical() ~ "{n} / {N} ({p}%)"
), digits = all_continuous() ~ 2,) %>% 
  add_p() %>% 
  add_overall() %>% 
  bold_labels()

