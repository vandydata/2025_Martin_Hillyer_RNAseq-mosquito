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
## Figure1: 
## Immune treatment, temperature, and adult age shape the mosquito’s transcriptome
## Principal component analysis across all samples, color-coded by immune treatment
## temperature or age
##
## Script Description:
##   Generates PCA plots of variance-stabilized RNA-seq data
##   using DESeq2 outputs. Plots colored by Treatment, Age, 
##   and Temperature.
## ---------------------------



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Configuration
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#load libraries
library(DESeq2)
library(ggplot2)
library(matrixStats)
source("functions/calcPCs.R")

# Path to your .RData file
rdata_file <- "data/077-For_PCA.Rdata"

# output path
output_dir <- "results/001-Figure1"
dir.create(output_dir, showWarnings = FALSE)

# Load DESeq2 object and metadata
# dds <- DESeq2 object
# sampleTable_ex <- a dataframe of sample metadata
load(rdata_file, verbose = FALSE)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Figure 1: Principal Component Analysis (PCA)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Variance-stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Compute PCA data
pcaData <- calcPCs(vsd, intgroup = c("condition", "Age", "Temp", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# PCA by Treatment
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
pcaData$Treatment <- factor(pcaData$Treatment, levels = c("Naive", "HKEC"))
treatment_colors <- c("Naive" = "#568982", "HKEC" = "#B6633F")

p <- ggplot(pcaData, aes(PC1, PC2, color = Treatment)) +
  geom_point(size = 3, stroke = 1.5) +
  scale_color_manual(values = treatment_colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color = "Immune Treatment") +
  theme_bw(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "001-PCA_by_treatment.png"), p, width = 8, height = 6, units = "in", dpi = 600)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# PCA by Age
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
pcaData$Age <- factor(pcaData$Age, levels = c("1d", "5d", "10d", "15d"))
age_colors <- c("1d" = "#E78AC3", "5d" = "#E15759", "10d" = "#66C2A5", "15d" = "#8DA0CB")

p <- ggplot(pcaData, aes(PC1, PC2, color = Age)) +
  geom_point(size = 3, stroke = 1.5) +
  scale_color_manual(values = age_colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color = "Age") +
  theme_bw(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "001-PCA_by_age.png"), p, width = 8, height = 6, units = "in", dpi = 600)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# PCA by Temperature
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
pcaData$Temp <- factor(pcaData$Temp, levels = c("27C", "30C", "32C"))
temp_colors <- c("27C" = "#4E79A7", "30C" = "#59A14F", "32C" = "#F28E2B")

p <- ggplot(pcaData, aes(PC1, PC2, color = Temp)) +
  geom_point(size = 3, stroke = 1.5) +
  scale_color_manual(values = temp_colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color = "Temperature") +
  theme_bw(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "001-PCA_by_temp.png"), p, width = 8, height = 6, units = "in", dpi = 600)
