# Genome-Wide Association and Population-Tailored Polygenic Risk for Parkinson’s Disease in Taiwan
GP2 ❤️ Open Science 😍  
Last Updated: November 2025

## Summary
This repository accompanies the manuscript **“Genome-Wide Association and Population-Tailored Polygenic Risk for Parkinson’s Disease in Taiwan.”** It contains analysis scripts and summary resources for the first and largest genome-wide investigation of Parkinson’s disease in Taiwanese populations. Genotype imputation was performed using the [Taiwan Biobank Imputation Server](https://twbimpute.biobank.org.tw/index.html), leveraging a population-specific whole-genome reference panel to optimize variant coverage and accuracy.


## Data Statement
All data used in this repository were generated using GP2-funded resources from the East Asian Parkinson’s Disease Genomics Consortium (EAPDGC), as part of the Global Parkinson’s Genetics Program (GP2). The underlying genotype and phenotype data are not distributed through this repository and can be accessed upon reasonable request to the EAPDGC/GP2 data access committees, in accordance with their data-sharing policies.

### Helpful links
- [EAPDGC](https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(21)00373-2/fulltext)
- [GP2 website](https://gp2.org/)
  - [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
  - [Other GP2 publication](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)
- [Taiwan Biobank Imputation Server](https://twbimpute.biobank.org.tw/index.html)
  - [Taiwan Biobank: A rich biomedical research database of the Taiwanese population](https://www.sciencedirect.com/science/article/pii/S2666979X2200146X)

## Citation
(pending pubication)

## Tables
(pending pubication)

## Figures
(pending pubication)

## Repository Orientation
The analysis/ directory includes all analyses discussed in the manuscript.

<pre> THIS_REPO/ 
  ├── analyses/ 
  |     ├── 01_QC_PCA.R
  |     ├── 02_GWAS.R
  |     ├── 03_manhattan.R
  |     ├── 04_regional_plot.R
  |     ├── 05_SNCA_conditional.R
  |     ├── 06_haplotype.R
  |     ├── 07a_1000G_TW_projection.R
  |     ├── 07b_1000G_TW_LRRK2_PCA.R
  |     ├── 08_LRRK2.R
  |     ├── 09_beta_comparison.R
  |     └── 10_PRS.R
  ├── figures/
  ├── tables/
  ├── LICENSE
  └── README.md 
</pre>

## Analysis Scripts
Primary language: R

| Directory   | Script                        | Description                                                                 |
|------------|-------------------------------|-----------------------------------------------------------------------------|
| `analyses/`| `01_QC_PCA.R`                 | Genotyping QC and principal component analysis. Input genotypes are PLINK files converted from IDAT via Illumina GenomeStudio. |
| `analyses/`| `02_GWAS.R`                   | Genome-wide association analysis using imputed data and PLINK 2.0.         |
| `analyses/`| `03_manhattan.R`              | Generation of Manhattan for GWAS results.                      |
| `analyses/`| `04_regional_plot.R`          | LocusZoom-style regional association plots for key loci.                    |
| `analyses/`| `05_SNCA_conditional.R`       | Conditional association analyses at the SNCA locus and regional plots.     |
| `analyses/`| `06_haplotype.R`              | Haplotype estimation and visualization for SNCA and LRRK2.                  |
| `analyses/`| `07a_1000G_TW_projection.R`   | PCA-based ancestry projection using 1000 Genomes and the Taiwan cohort.     |
| `analyses/`| `07b_1000G_TW_LRRK2_PCA.R`    | East Asian + Taiwan PCA with LRRK2 variant carriers highlighted.            |
| `analyses/`| `08_LRRK2.R`                  | Case–control and age-at-onset analyses stratified by LRRK2 Asian variants.  |
| `analyses/`| `09_beta_comparison.R`        | Effect size comparison with prior PD GWAS (e.g., Nalls 2019, Foo 2020) and beta–beta plots. |
| `analyses/`| `10_PRS.R`                    | Construction and evaluation of polygenic risk scores (AUC, ROC, deciles).   |


## Software / Packages

| Software                               | Version(s) | Resource URL                                                                 | RRID        | Notes |
|----------------------------------------|------------|------------------------------------------------------------------------------|-------------|-------|
| `GenomeStudio`                         | 2.0        | https://www.illumina.com/products/by-type/informatics-products/microarray-software/genomestudio.html | SCR_010973 | Used to convert IDAT intensity files to PLINK-formatted genotype data. |
| `PLINK`                                | 2.0        | https://www.cog-genomics.org/plink/2.0/                                     | SCR_001757 | Used for GWAS, scoring (PRS), and basic genotype processing. |
| `PLINK`                                | 1.9        | https://www.cog-genomics.org/plink/                                         | SCR_001757 | Used for LD calculations, pruning, IBD, and legacy QC steps. |
| `R Project for Statistical Computing`  | 4.3        | https://www.r-project.org/                                                  | SCR_001905 | Used for downstream analyses and visualization. Key packages include: `dplyr`, `tidyr`, `data.table`, `ggplot2`, `haplo.stats`, `geneHapR`, `pROC`, and `patchwork`. |
