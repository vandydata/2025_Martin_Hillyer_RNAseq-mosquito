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
## Figure5: 
## Aging from 1 to 15 days upregulates, 
## downregulates or has a bidirectional effect on the expression of 
## genes in distinct functional categories
##
## ---------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Configuration
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(forcats)
library(cowplot)
library(patchwork)
source("src/functions/plot_per_term.R")

# Input data : GSEA results 
GSEA_15d__vs__1d_Naive <- read_csv("data/VectorDb_15d__vs__1d_Naive.csv")
GSEA_15d__vs__1d_HKEC <- read_csv("data/VectorDb_15d__vs__1d_HKEC.csv")


# output path
output_dir <- "results/005-Figure5/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1.GSEA_15d__vs__1d_Naive ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set x-axis limits
x_min <- floor(min(GSEA_15d__vs__1d_Naive$NES))
x_max <- ceiling(max(GSEA_15d__vs__1d_Naive$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically

plot_list <- lapply(unique(GSEA_15d__vs__1d_Naive$New_representative_term), function(term) plot_per_term(term, GSEA_15d__vs__1d_Naive, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
output_dir <- "results/005-Figure5/Naive/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(GSEA_15d__vs__1d_Naive$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 2.GSEA_15d__vs__1d_HKEC ----

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set x-axis limits
x_min <- floor(min(GSEA_15d__vs__1d_HKEC$NES))
x_max <- ceiling(max(GSEA_15d__vs__1d_HKEC$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically

plot_list <- lapply(unique(GSEA_15d__vs__1d_HKEC$New_representative_term), function(term) plot_per_term(term, GSEA_15d__vs__1d_HKEC, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
output_dir <- "results/005-Figure5/HKEC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(GSEA_15d__vs__1d_HKEC$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}

