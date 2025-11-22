# Genome-Wide Association and Population-Tailored Polygenic Risk for Parkinson’s Disease in Taiwan
GP2 ❤️ Open Science 😍
Last Updated: November 2025

## Summary
This is the online repository for the manuscript titled "Genome-Wide Association and Population-Tailored Polygenic Risk for Parkinson’s Disease in Taiwan". This study represents the first and largest genome-wide assessment of Parkinson’s disease in the Taiwanese populations.

## Data Statement
All data used in this repository were generated using GP2-funded resources from the East Asian Parkinson’s Disease Genomics Consortium (EAPDGC), as part of the Global Parkinson’s Genetics Program (GP2). The underlying genotype and phenotype data are not distributed through this repository and can be accessed upon reasonable request to the EAPDGC/GP2 data access committees, in accordance with their data-sharing policies.

### Helpful links
- [EAPDGC](https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(21)00373-2/fulltext)
- [GP2 website](https://gp2.org/)
  - [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)

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
  |     ├── 00_notebook.ipynb
  |     └── 01_notebook.ipynb
  ├── figures/
  |     └── 00_figure.png
  ├── tables/
  |     └── 00_table.csv
  ├── LICENSE
  └── README.md 
</pre>

## Analysis Scripts
Languages: R

| Directory |        Notebooks     |     Description     | 
|-----------|----------------------|---------------------|
|`analyses/`| `01_QC.R`  | Quality Control|
|`analyses/`| `02_GWAS.R`  | GWAS and Logistic Regression and related figures|
|`analyses/`| `03_regional_plot.R`  | regional_plots |
|`analyses/`| `04_haplotype.R`  | SNCA and LRRK2 haplotype staistics and figures|
|`analyses/`| `05_LRRK2.R`  | analysis based on LRRK2 Asian variant status |
|`analyses/`| `06_beta_comparison.R`  | beta-beta plot |
|`analyses/`| `07_PRS.R`  | polygenic risk score |

## Software/Packages

| Software	| Version(s)	| Resource URL | 
|-----------|----------------------|---------------------|
|`Genome Studio`|`2.0`|[https://www.illumina.com/products/by-type/informatics-products/microarray-software/genomestudio.html](https://www.illumina.com/products/by-type/informatics-products/microarray-software/genomestudio.html)|
|`PILNK`| `2.0`  | [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/) |
|`PILNK`| `1.9`  | [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/)|
|`R Project for Statistical Computing`| `4.3`  | [https://www.r-project.org/](https://www.r-project.org/) |
|`R haplo.stats`|`1.9.7`| [https://cran.r-project.org/web/packages/haplo.stats/index.html](https://cran.r-project.org/web/packages/haplo.stats/index.html) |
