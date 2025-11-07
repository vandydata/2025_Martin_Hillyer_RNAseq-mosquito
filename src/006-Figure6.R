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
## Figure6: 
## Early aging (1 to 5 days) and late aging (10 to 15 days) have 
## different effects in the regulation of the expression of genes in 
## distinct functional categories
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
source("src/functions/plot_per_term.R")

# Input data : GSEA results 
GSEA_5d__vs__1d_Naive <- read_csv("data/VectorDb_5d__vs__1d_Naive.csv")
GSEA_15d__vs__10d_Naive <- read_csv("data/VectorDb_15d__vs__10d_Naive.csv")

GSEA_15d__vs__10d_HKEC <- read_csv("data/VectorDb_15d__vs__10d_HKEC.csv")
GSEA_5d__vs__1d_HKEC <- read_csv("data/VectorDb_5d__vs__1d_HKEC.csv")

# output path
output_dir <- "results/006-Figure6/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1.GSEA_5d__vs__1d_Naive ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set x-axis limits
x_min <- floor(min(GSEA_5d__vs__1d_Naive$NES))
x_max <- ceiling(max(GSEA_5d__vs__1d_Naive$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically

plot_list <- lapply(unique(GSEA_5d__vs__1d_Naive$New_representative_term), function(term) plot_per_term(term, GSEA_5d__vs__1d_Naive, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
output_dir <- "results/006-Figure6/5d__vs__1d_Naive/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(GSEA_5d__vs__1d_Naive$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2.GSEA_15d__vs__10d_Naive ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set x-axis limits
x_min <- floor(min(GSEA_15d__vs__10d_Naive$NES))
x_max <- ceiling(max(GSEA_15d__vs__10d_Naive$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically

plot_list <- lapply(unique(GSEA_15d__vs__10d_Naive$New_representative_term), function(term) plot_per_term(term, GSEA_15d__vs__10d_Naive, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
output_dir <- "results/006-Figure6/15d__vs__10d_Naive/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(GSEA_15d__vs__10d_Naive$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3.GSEA_5d__vs__1d_HKEC ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set x-axis limits
x_min <- floor(min(GSEA_5d__vs__1d_HKEC$NES))
x_max <- ceiling(max(GSEA_5d__vs__1d_HKEC$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically

plot_list <- lapply(unique(GSEA_5d__vs__1d_HKEC$New_representative_term), function(term) plot_per_term(term, GSEA_5d__vs__1d_HKEC, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
output_dir <- "results/006-Figure6/5d__vs__1d_HKEC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(GSEA_5d__vs__1d_HKEC$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4.GSEA_15d__vs__10d_HKEC ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Set x-axis limits
x_min <- floor(min(GSEA_15d__vs__10d_HKEC$NES))
x_max <- ceiling(max(GSEA_15d__vs__10d_HKEC$NES))
x_lim <- max(abs(c(x_min, x_max)))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate and save plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Generate a list of plots, one for each term, then align vertically

plot_list <- lapply(unique(GSEA_15d__vs__10d_HKEC$New_representative_term), function(term) plot_per_term(term, GSEA_15d__vs__10d_HKEC, x_lim))
plots <- align_patches(plot_list, direction = "v")

# Save individually as png + pdf
output_dir <- "results/006-Figure6/15d__vs__10d_HKEC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


for (i in 1:length(plot_list)) {
  filename <- paste0(output_dir, unique(GSEA_15d__vs__10d_HKEC$New_representative_term)[i], ".png")
  ggsave(filename, plots[[i]], width=5, height=6, units='in', scale=1.5, dpi=300)
}

