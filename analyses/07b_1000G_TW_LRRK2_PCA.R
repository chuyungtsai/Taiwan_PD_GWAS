## ================================================================
## 07b_1000G_TW_LRRK2_PCA.R
## East Asian + Taiwan PCA  with LRRK2 carriers highlighted
## Produces:
##   - 20251019_LRRK2_EAS_TW.jpeg
##
## Requires:
##   - 20130606_g1k.ped
##   - 20251018_TW_1000G_Merge_QC.{bed,bim,fam}
##   - 1000G_pruned/20251019_pruned/20251019_1000G_TW_PCA.eigenvec
##   - LRRK2 carrier list for TW:
##       20251013_LRRK2_variants.csv (columns IID, G2385R, R1628P)
##   - (Optional) G2385R_1000G, R1628P_1000G data.frames with IID column
## ================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

## ---------------- USER CONFIGURATION ----------------------------

ped_file <- "20130606_g1k.ped"
fam_file <- "20251018_TW_1000G_Merge_QC.fam"

eigenvec_file <- file.path(
  "1000G_pruned", "20251019_pruned",
  "20251019_1000G_TW_PCA.eigenvec"
)

lrrk2_carrier_file <- "~/20251013_LRRK2_variants.csv"

keep_file <- file.path("1000G_pruned", "20251019_pruned", "20251018_TW_EAS_keep.txt")
tw_eas_prefix <- "20251019_TW_EAS_PCA"

pca_lrrk2_png <- "20251019_LRRK2_EAS_TW_PC2_PC3.jpeg"

plink_bin <- "plink19"

## ================================================================
## 1. RECREATE TW + 1000G META (SAME LOGIC AS 06a)
## ================================================================

# 1.1 PCA (TW + 1000G)
eigenvec_raw <- read.table(eigenvec_file, header = FALSE, sep = " ")
colnames(eigenvec_raw) <- c("FID", "IID", paste0("PC", seq_len(ncol(eigenvec_raw) - 2)))

TW_1000G_PCA <- eigenvec_raw
rownames(TW_1000G_PCA) <- TW_1000G_PCA$IID
TW_1000G_PCA$IID <- rownames(TW_1000G_PCA)

# 1.2 1000G population information
PED <- read.table(ped_file, header = TRUE, sep = "\t")
PED <- PED[PED$Individual.ID %in% rownames(TW_1000G_PCA), ]
PED <- PED[match(rownames(TW_1000G_PCA), PED$Individual.ID), ]
stopifnot(all(PED$Individual.ID == rownames(TW_1000G_PCA)))
PED <- PED %>% mutate(IID = Individual.ID)

Population_1000G <- PED %>%
  select(IID, Population)

# 1.3 Taiwan vs 1000G
TW_IID <- setdiff(TW_1000G_PCA$IID, Population_1000G$IID)
TW_df  <- data.frame(IID = TW_IID, Population = "TW", stringsAsFactors = FALSE)

Total_population <- bind_rows(TW_df, Population_1000G)

pca_tab <- Total_population %>%
  left_join(TW_1000G_PCA, by = "IID")

pcaped <- pca_tab %>%
  mutate(
    Superpopulation = case_when(
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

## ================================================================
## 2. BUILD TW + EAS SUBSET AND RUN PCA (PLINK)
## ================================================================

TW_EAS_meta <- pcaped %>%
  filter(POP %in% c("TW", "EAS")) %>%
  select(IID, POP, Population)

# Merged TW + 1000G FAM (after QC/merge)
TW_1000G_fam <- read.table(fam_file, header = FALSE)
colnames(TW_1000G_fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

TW_EAS_keep <- TW_1000G_fam %>%
  filter(IID %in% TW_EAS_meta$IID) %>%
  select(FID, IID)

# Write keep file (FID IID) for PLINK
write.table(
  TW_EAS_keep,
  file      = keep_file,
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# Run PCA on TW + EAS only
pruned_dir <- file.path(getwd(), "1000G_pruned", "20251019_pruned")
oldwd <- getwd()
setwd(pruned_dir)

system(paste(
  plink_bin,
  "--bfile 20251019_1000G_TW_merge",
  "--keep 20251018_TW_EAS_keep.txt",
  "--make-bed --out 20251019_TW_EAS"
))

system(paste(
  plink_bin,
  "--bfile 20251019_TW_EAS",
  "--pca 20 --out", tw_eas_prefix
))

setwd(oldwd)

# Read TW + EAS PCA
TW_EAS_eigenvec <- read.table(
  file.path(pruned_dir, paste0(tw_eas_prefix, ".eigenvec")),
  header = FALSE, sep = " "
)
colnames(TW_EAS_eigenvec) <- c("FID", "IID", paste0("PC", seq_len(ncol(TW_EAS_eigenvec) - 2)))
TW_EAS_eigenvec$IID <- as.character(TW_EAS_eigenvec$IID)

TW_EAS_PCA <- TW_EAS_meta %>%
  left_join(TW_EAS_eigenvec, by = "IID")

## ================================================================
## 3. LRRK2 CARRIER LABELING (TW + 1000G EAS)
## ================================================================

# 3.1 Taiwan carriers from TWB LRRK2 file
LRRK2_carrier <- read.csv(lrrk2_carrier_file, header = TRUE)

G2385R_TW <- LRRK2_carrier %>%
  filter(G2385R >= 1) %>%
  select(IID)

R1628P_TW <- LRRK2_carrier %>%
  filter(R1628P >= 1) %>%
  select(IID)

# 3.2 1000G carriers — if not already loaded, create empty placeholders
if (!exists("G2385R_1000G")) {
  G2385R_1000G <- data.frame(IID = character(0), stringsAsFactors = FALSE)
}

if (!exists("R1628P_1000G")) {
  R1628P_1000G <- data.frame(IID = character(0), stringsAsFactors = FALSE)
}

G2385R_carrier <- bind_rows(G2385R_TW, G2385R_1000G) %>%
  mutate(IID = str_trim(as.character(IID))) %>%
  distinct()

R1628P_carrier <- bind_rows(R1628P_TW, R1628P_1000G) %>%
  mutate(IID = str_trim(as.character(IID))) %>%
  distinct()

G2385R_vec <- G2385R_carrier$IID
R1628P_vec <- R1628P_carrier$IID

TW_EAS_PCA <- TW_EAS_PCA %>%
  mutate(
    G2385R_carrier = IID %in% G2385R_vec,
    R1628P_carrier = IID %in% R1628P_vec
  )

## ================================================================
## 4. PCA PLOTS (PC2–PC3) WITH LRRK2 CARRIERS HIGHLIGHTED
## ================================================================

pal_eas <- c(
  "TW" = "grey80",   # Taiwan cohort
  "JPT" = "#56B4E9",
  "CHB" = "#009E73",
  "CDX" = "#E69F00",
  "CHS" = "#0072B2",
  "KHV" = "#CC79A7"
)

p_base_pc23 <- ggplot(TW_EAS_PCA, aes(PC2, PC3, color = Population)) +
  geom_point(size = 1.6, alpha = 0.9) +
  scale_color_manual(values = pal_eas, drop = FALSE, name = "Population") +
  labs(
    x = "PC2",
    y = "PC3"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_G2385R <- p_base_pc23 +
  geom_point(
    data        = dplyr::filter(TW_EAS_PCA, G2385R_carrier),
    aes(PC2, PC3),
    inherit.aes = FALSE,
    shape       = 21,
    fill        = NA,
    color       = "#d62728",
    size        = 3.2,
    stroke      = 1.1
  ) +
  labs(title = "G2385R") +
  theme(
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.title.position = "plot"
  )

p_R1628P <- p_base_pc23 +
  geom_point(
    data        = dplyr::filter(TW_EAS_PCA, R1628P_carrier),
    aes(PC2, PC3),
    inherit.aes = FALSE,
    shape       = 21,
    fill        = NA,
    color       = "#d62728",
    size        = 3.2,
    stroke      = 1.1
  ) +
  labs(title = "R1628P") +
  theme(
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.title.position = "plot"
  )

xlim_pc23 <- range(TW_EAS_PCA$PC2, na.rm = TRUE)
ylim_pc23 <- range(TW_EAS_PCA$PC3, na.rm = TRUE)

p1 <- p_R1628P +
  coord_cartesian(xlim = xlim_pc23, ylim = ylim_pc23) +
  scale_color_manual(
    values = pal_eas,
    name   = "Population (1000G EAS; TW = Taiwan)"
  )

p2 <- p_G2385R +
  coord_cartesian(xlim = xlim_pc23, ylim = ylim_pc23) +
  scale_color_manual(
    values = pal_eas,
    name   = "Population (1000G EAS; TW = Taiwan)"
  )

combined_lrrk2 <- (p1 + p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "LRRK2 p.R1628P and p.G2385R on East Asian PCA (PC2–PC3)",
    subtitle = "1000 Genomes East Asian populations and Taiwan cohort; carriers highlighted with red rings",
    caption  = "LD-pruned autosomal SNPs; unrelated individuals; GRCh38"
  ) &
  theme(
    legend.position   = "right",
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle     = element_text(hjust = 0.5),
    plot.caption      = element_text(size = 9, colour = "grey40")
  )

ggsave(pca_lrrk2_png, combined_lrrk2, width = 12, height = 6, dpi = 300)

message("06b finished:
  - LRRK2 carriers (EAS + TW): ", pca_lrrk2_png)
