# Manuscript inputs and code for manuscript figures

Warmer temperature accelerates senescence by modifying the aging-dependent changes in the mosquito transcriptome, altering immunity, metabolism, and DNA repair
Lindsay E. Martin, Jordyn S. Barr, Jean-Philippe Cartailler, Shristi Shrestha, Tania Y. Estévez-Lao, Julián F. Hillyer
bioRxiv 2025.10.27.684792; doi: https://doi.org/10.1101/2025.10.27.684792 

## Requirements

R packages:
- DESeq2
- ggplot2
- matrixStats
- readr, dplyr, tidyr
- ggpubr, cowplot, forcats, stringr

## Data

Input files are in `data/`:
- `077-For_PCA.Rdata` (multipart, use 7zip to unpack) - DESeq2 object for PCA analysis
- `103-interaction_analysis_HKEC.RData` / `.csv` - Interaction analysis for HKEC treatment
- `107-interaction_analysis_Naive.Rdata` / `.csv` - Interaction analysis for naive samples
- `VectorDb_*.csv` - GSEA results for various comparisons
- `annotation.csv`, `final-eigengenes.txt`, `final-membership.txt`, `gene_module_membership.tsv.zip`, `meta-network-module-metamodule-label-remap.txt`, `normalized_counts.csv.zst` - WGCNA

## Scripts

All scripts are in `src/` and generate figures in corresponding `results/` subdirectories.

### Figure generation

- `001-Figure1.R` - PCA plots colored by treatment, age, and temperature
- `002-Figure2.R` - WGCNA plots
- `003-Figure3.R` - GSEA enrichment plots for HKEC vs Naive comparison
- `004-Figure4.R` - Temperature effects on gene expression by functional category
- `005-Figure5.R` - Age effects on gene expression by functional category
- `006-Figure6.R` - Early vs late aging effects on gene expression
- `007-Figure7.R` - Differential aging responses between treatments

### Helper functions

Functions in `src/functions/`:
- `calcPCs.R` - Compute PCA from DESeq2 variance-stabilized data
- `plot_per_term.R` - Generate enrichment plots for individual GO terms
