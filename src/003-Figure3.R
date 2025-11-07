##
##      __________  _____
##     / ____/ __ \/ ___/
##    / /   / / / /\__ \ 
##   / /___/ /_/ /___/ / 
##   \____/_____//____/  
##
##  Creative Data Solutions
##  Vanderbilt University
##  https://cds.vanderbilt.edu
##
## Date Created: 2025-11-03
## Shristi Shrestha, Jean-Philippe Cartailler
##
## ---------------------------
## Figure3: 
## Immune induction upregulates, downregulates or has a bidirectional effect
## on the expression of genes in distinct functional categories
##
## ---------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Configuration
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(forcats)
library(cowplot)
source("src/functions/plot_per_term.R")

# Input data : GSEA results from differentially expressed genes HKEC__vs__Naive
HKEC_vs_Naive_GSEA <- read_csv("data/VectorDb_HKEC__vs__Naive.csv")

# output path
output_dir <- "results/003-Figure3/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set x-axis limits
x_min <- floor(min(HKEC_vs_Naive_GSEA$NES))
x_max <- ceiling(max(HKEC_vs_Naive_GSEA$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically
output_dir <- "results/003-Figure3/"
plot_list <- lapply(unique(HKEC_vs_Naive_GSEA$New_representative_term), function(term) plot_per_term(term, HKEC_vs_Naive_GSEA, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(HKEC_vs_Naive_GSEA$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}

