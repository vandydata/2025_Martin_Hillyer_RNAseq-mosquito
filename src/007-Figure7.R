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
## Figure7: 
## The interaction between warmer temperature and aging shapes biological 
## processes within the mosquito
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
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplotify)

# output path
output_dir <- "results/007-Figure7/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Heatmap (Naive) ----
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Input data (interaction analysis: Differentially expressed genes)
interaction_analysis_Naive <- read_csv("data/107-interaction_analysis_Naive.csv")

#rld is rlogtransformed counts
load("data/107-interaction_analysis_Naive.RData")


genes_of_interest <- unique(interaction_analysis_Naive$AGAP_ID)

# matrix for heatmap. rld is rlogtransformed counts
mat <- assay(rld)

# subset matrix for selected genes
mat2 <- mat[match(genes_of_interest, rownames(mat)), ]

# Set up annotations
top_anno_order_df <- as.data.frame(colData(rld)[, c("condition", "Temp", "Age", "Treatment")])
top_anno_order_df$Temp <- factor(top_anno_order_df$Temp, levels = c("27C", "30C", "32C"))
top_anno_order_df$Age <- factor(top_anno_order_df$Age, levels = c("1d", "5d", "10d", "15d"))
top_anno_order_df$Treatment <- factor(top_anno_order_df$Treatment, levels = c("Naive"))


# rearrange order of the counts matrix
top_anno_order_df <- top_anno_order_df %>%
  arrange(Treatment, Temp, Age)
mat3 <- mat2[, rownames(top_anno_order_df)]

# Heatmap top annotatio 
ha <- HeatmapAnnotation(
  Immune_infection = top_anno_order_df$Treatment,
  Temp = top_anno_order_df$Temp,
  Age = top_anno_order_df$Age,
  col = list(
    Immune_infection = c("Naive" = "#568982", "HKEC" = "#B6633F"),
    Temp = c("27C" = "#4E79A7", "30C" = "#59A14F", "32C" = "#F28E2B"),
    Age = c("10d" = "#66C2A5", "15d" = "#8DA0CB", "1d" = "#E78AC3", "5d" = "#E15759")
  )
)

all(colnames(mat3) == rownames(top_anno_order_df))

colors_heatmap <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(3)

# scale
mat3 <- t(scale(t(mat3)))

set.seed(100)
h1 <- Heatmap(mat3,
              name = "z-score",
              row_km = 10,
              border = TRUE,
              cluster_rows = T,
              show_row_names = F,
              row_title_rot = 0,
              row_title_gp = gpar(fontsize = 18, face = "bold"),
              show_row_dend = FALSE,
              column_title = paste0(nrow(mat3), " genes, padj<0.05, log2FC>1"),
              column_title_gp = grid::gpar(fontsize = 16),
              cluster_columns = FALSE,
              clustering_distance_column = function(x, y) 1 - cor(x, y), 
              clustering_method_column = "complete", 
              col = circlize::colorRamp2(c(-2, 0, 2), colors_heatmap), 
              column_names_rot = 90,
              top_annotation = ha,
              heatmap_legend_param = list(legend_direction = "vertical", 
                                          legend_position = "right",
                                          title_gp = gpar(fontsize = 12, face = "bold")
              )
)

p <- as.ggplot(grid.grabExpr(draw(h1))) 
filename <- paste0(output_dir, "Heatmapinteraction_analysis_Naive.png")
ggsave(filename, plot = p, width = 14, height = 12, units = "in", scale = 1, dpi = 300)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Heatmap (HKEC) ----
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Input data (interaction analysis: Differentially expressed genes)
interaction_analysis_HKEC <- read_csv("data/103-interaction_analysis_HKEC.csv")

#rld is rlogtransformed counts
load("data/103-interaction_analysis_HKEC.RData")

genes_of_interest <- unique(interaction_analysis_HKEC$AGAP_ID)

# matrix for heatmap. rld is rlogtransformed counts
mat <- assay(rld)

# subset matrix for selected genes
mat2 <- mat[match(genes_of_interest, rownames(mat)), ]

# Set up annotations
top_anno_order_df <- as.data.frame(colData(rld)[, c("condition", "Temp", "Age", "Treatment")])
top_anno_order_df$Temp <- factor(top_anno_order_df$Temp, levels = c("27C", "30C", "32C"))
top_anno_order_df$Age <- factor(top_anno_order_df$Age, levels = c("1d", "5d", "10d", "15d"))
top_anno_order_df$Treatment <- factor(top_anno_order_df$Treatment, levels = c("HKEC"))

# rearrange order of the counts matrix
top_anno_order_df <- top_anno_order_df %>%
  arrange(Treatment, Temp, Age)
mat3 <- mat2[, rownames(top_anno_order_df)]

# Heatmap top annotatio 
ha <- HeatmapAnnotation(
  Immune_infection = top_anno_order_df$Treatment,
  Temp = top_anno_order_df$Temp,
  Age = top_anno_order_df$Age,
  col = list(
    Immune_infection = c("Naive" = "#568982", "HKEC" = "#B6633F"),
    Temp = c("27C" = "#4E79A7", "30C" = "#59A14F", "32C" = "#F28E2B"),
    Age = c("10d" = "#66C2A5", "15d" = "#8DA0CB", "1d" = "#E78AC3", "5d" = "#E15759")
  )
)

all(colnames(mat3) == rownames(top_anno_order_df))

colors_heatmap <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(3)

# scale
mat3 <- t(scale(t(mat3)))

set.seed(100)
h1 <- Heatmap(mat3,
              name = "z-score",
              row_km = 10,
              border = TRUE,
              cluster_rows = T,
              show_row_names = F,
              row_title_rot = 0,
              row_title_gp = gpar(fontsize = 18, face = "bold"),
              show_row_dend = FALSE,
              column_title = paste0(nrow(mat3), " genes, padj<0.05, log2FC>1"),
              column_title_gp = grid::gpar(fontsize = 16),
              column_names_gp = grid::gpar(fontsize = 11.5),
              cluster_columns = FALSE,
              clustering_distance_column = function(x, y) 1 - cor(x, y), 
              clustering_method_column = "complete", 
              col = circlize::colorRamp2(c(-2, 0, 2), colors_heatmap), 
              column_names_rot = 90,
              top_annotation = ha,
              heatmap_legend_param = list(legend_direction = "vertical", 
                                          legend_position = "right",
                                          title_gp = gpar(fontsize = 12, face = "bold")    # Increase title font size
              )
)

p <- as.ggplot(grid.grabExpr(draw(h1))) 
filename <- paste0(output_dir, "Heatmapinteraction_analysis_HKEC.png")
ggsave(filename, plot = p, width = 14, height = 12, units = "in", scale = 1, dpi = 300)

