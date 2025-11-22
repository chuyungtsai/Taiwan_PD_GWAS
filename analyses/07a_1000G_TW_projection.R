
## ================================================================
## 07a_1000G_TW_projection.R
## 1000 Genomes + Taiwan PCA projection 
## Produces:
##   - 20251019_TW_1000G_PCA.jpeg
## ================================================================

setwd("./1000G_vcf")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

## ---------------- USER CONFIGURATION ----------------------------

# 1000G pedigree / meta file (tab-delimited)
ped_file <- "20130606_g1k.ped"

# PCA eigenvectors for merged Taiwan + 1000G (PLINK output)
eigenvec_file <- file.path(
  "1000G_pruned", "20251019_pruned",
  "20251019_1000G_TW_PCA.eigenvec"
)

# Output figure
pca_global_png <- "20251019_TW_1000G_PCA.jpeg"

## ================================================================
## 1. LOAD PCA (TW + 1000G) AND 1000G POPULATION LABELS
## ================================================================

# PLINK eigenvec: FID IID PC1 PC2 ...
eigenvec_raw <- read.table(eigenvec_file, header = FALSE, sep = " ")
colnames(eigenvec_raw) <- c("FID", "IID", paste0("PC", seq_len(ncol(eigenvec_raw) - 2)))

TW_1000G_PCA <- eigenvec_raw
rownames(TW_1000G_PCA) <- TW_1000G_PCA$IID

# 1000G meta / population
PED <- read.table(ped_file, header = TRUE, sep = "\t")

# Keep only 1000G individuals present in PCA
PED <- PED[PED$Individual.ID %in% rownames(TW_1000G_PCA), ]

# Order to match eigenvec IIDs
PED <- PED[match(rownames(TW_1000G_PCA), PED$Individual.ID), ]
stopifnot(all(PED$Individual.ID == rownames(TW_1000G_PCA)))

PED <- PED %>%
  mutate(IID = Individual.ID)

Population_1000G <- PED %>%
  select(IID, Population)

## ================================================================
## 2. DEFINE TAIWAN vs 1000G AND BUILD POP TABLE
## ================================================================

TW_1000G_PCA$IID <- rownames(TW_1000G_PCA)

# Taiwan samples = IIDs present in PCA but not in 1000G PED
TW_IID <- setdiff(TW_1000G_PCA$IID, Population_1000G$IID)
TW_df  <- data.frame(IID = TW_IID, Population = "TW", stringsAsFactors = FALSE)

Total_population <- dplyr::bind_rows(TW_df, Population_1000G)

pca_tab <- Total_population %>%
  left_join(TW_1000G_PCA, by = "IID")

## ================================================================
## 3. SUPERPOPULATION MAPPING (TW / AFR / AMR / EAS / EUR / SAS)
## ================================================================

pcaped <- pca_tab %>%
  mutate(
    Superpopulation = dplyr::case_when(
      Population %in% c("JPT","CHB","CDX","CHS","KHV","CHD") ~ "EAS",
      Population %in% c("BEB","GIH","ITU","PJL","STU")       ~ "SAS",
      Population %in% c("CLM","MXL","PEL","PUR")             ~ "AMR",
      Population %in% c("CEU","FIN","GBR","IBS","TSI")       ~ "EUR",
      Population %in% c("ACB","ASW","ESN","GWD","LWK","MSL","YRI") ~ "AFR",
      TRUE ~ NA_character_
    ),
    POP = if_else(is.na(Superpopulation), Population, Superpopulation),
    POP = factor(POP, levels = c("TW","AFR","AMR","EAS","EUR","SAS"))
  )

pal_super <- c(
  TW  = "grey60",
  AFR = "#1f77b4",
  AMR = "#ff7f0e",
  EAS = "#2ca02c",
  EUR = "#d62728",
  SAS = "#9467bd"
)

## ================================================================
## 4. GLOBAL PCA PLOT (PC1–PC2) – 1000G + TAIWAN
## ================================================================

p_global <- ggplot(pcaped, aes(PC1, PC2, color = POP)) +
  geom_point(size = 1.6, alpha = 0.9) +
  scale_color_manual(values = pal_super, drop = FALSE, name = "POP") +
  labs(
    x = "PC1",
    y = "PC2",
    title = "1000 Genome PCA projection"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.title.position = "plot"
  )

ggsave(pca_global_png, p_global, width = 6.5, height = 5, dpi = 300)

message("06a finished:
  - Global 1000G + TW PCA (PC1–PC2): ", pca_global_png)
