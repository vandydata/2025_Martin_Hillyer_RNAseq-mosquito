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
## Figure2: 
## Fig 2. Immune treatment, temperature, and age have unique effects on different
## functional parts of the mosquito’s transcriptome.
##
## ---------------------------


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel A
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# A network was created from data/final-eigengenes.txt
# and data/final-membership.txt. The method for this is 
# closed source. Please contact author for more information.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel B
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Load packages
suppressPackageStartupMessages({
  require(WGCNA)
  require(readr)
  library(ComplexHeatmap)
  library(ggplot2)
  library(circlize)
  library(ggplotify)
  library(dplyr)
  library(mltools)
  library(data.table)
})

options(stringsAsFactors = FALSE)

# Static configuration variables
run_name <- "iwgcna_power_20"
output_dir <- "results/002-Figure2"
data_dir <- "data"
normalized_file <- glue("{data_dir}/normalized_counts.csv")
metadata_file <- glue("{data_dir}/annotations.csv")
remapping_file <- glue("{data_dir}/meta-network-module-metamodule-label-remap.txt")

# Load data
metadata <- read.csv(metadata_file)
data_norm <- read_csv(normalized_file)
remapping <- read.csv(remapping_file, sep="\t")

# Prepare expression data
colnames(data_norm)[1] <- "gene"
datExpr0 <- as.data.frame(t(data_norm[,-1]))
names(datExpr0) <- data_norm$gene
rownames(datExpr0) <- names(data_norm)[-1]

# Quality control
gsg <- goodSamplesGenes(datExpr0, verbose=3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Prepare trait data
df <- data.frame(exp_sample_id=rownames(datExpr))
metadata.merge <- merge(df, metadata, by.x="exp_sample_id", by.y="cds_id", all.x=TRUE)

# Select and encode traits
allSamples <- rownames(datExpr)
traitRows <- match(allSamples, metadata.merge$exp_sample_id)
datTraits <- metadata.merge[traitRows, c("temperature", "age", "treatment")]

# One-hot encoding
datTraits[] <- lapply(datTraits, as.factor)
datTraits <- one_hot(as.data.table(datTraits))
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- allSamples

# Load module data
finalMembership <- read_delim(paste0(data_dir, "/final-membership.txt"), "\t")
moduleColors <- finalMembership$Module

# Calculate module-trait correlations
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
names(MEs) <- substring(names(MEs), 3)

moduleTraitCor <- cor(MEs, datTraits, use="p", method='pearson')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Remove unclassified modules
moduleTraitCor <- moduleTraitCor[!rownames(moduleTraitCor) %in% "UNCLASSIFIED",]
moduleTraitPvalue <- moduleTraitPvalue[!rownames(moduleTraitPvalue) %in% "UNCLASSIFIED",]

# Create text matrix for heatmap
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
rownames(textMatrix) <- rownames(moduleTraitCor)
colnames(textMatrix) <- colnames(moduleTraitCor)

# Generate heatmaps
c <- colorRampPalette(c("#377eb8", "white", "#e41a1c"))(100)

# Basic heatmap
heat <- Heatmap(
  moduleTraitCor,
  col=c,
  cluster_columns=TRUE,
  cluster_rows=TRUE,
  name="Correlation",
  column_names_rot=45,
  rect_gp=gpar(col="white", lwd=1),
  cell_fun=function(j, i, x, y, width, height, fill) {
    color <- if(abs(moduleTraitCor[i, j]) < 0.5) "black" else "white"
    grid.text(textMatrix[i, j], x, y, gp=gpar(fontsize=8, col=color))
  }
)

p <- as.ggplot(grid.grabExpr(draw(heat)))
ggsave(paste0(output_dir, '/trait-module-heatmap.png'), p, height=6, width=3, scale=2)

# Pretty labeled version with remapping
n <- data.frame(Module=rownames(moduleTraitCor))
n <- left_join(n, remapping, by="Module")
n$ModulePrettyLabel <- paste0(n$MetamoduleLabel, " - ", n$ModuleLabel)
rownames(moduleTraitCor) <- n$ModulePrettyLabel

heat <- Heatmap(
  moduleTraitCor,
  col=c,
  cluster_columns=TRUE,
  cluster_rows=TRUE,
  name="Correlation",
  column_names_rot=45,
  rect_gp=gpar(col="white", lwd=1),
  cell_fun=function(j, i, x, y, width, height, fill) {
    color <- if(abs(moduleTraitCor[i, j]) < 0.5) "black" else "white"
    grid.text(textMatrix[i, j], x, y, gp=gpar(fontsize=8, col=color))
  }
)

p <- as.ggplot(grid.grabExpr(draw(heat)))
ggsave(paste0(output_dir, '/trait-module-heatmap-labeled.png'), p, height=6, width=3, scale=2)




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Panel C, D, E
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#!/usr/bin/env Rscript
# Minified GO Enrichment Analysis for WGCNA Modules

# Load packages
suppressPackageStartupMessages({
  require(org.Mm.eg.db)
  require(GO.db)
  require(clusterProfiler)
  require(ggplot2)
  library(dplyr)
  library(cowplot)
})

# Static configuration variables
run_name <- "iwgcna_power_20"
input_file <- "data/gene_module_membership.tsv"
gaf_file <- "data/VectorBase-66_AgambiaePEST_Curated_GO.gaf" # download from veupathdb

# Load data
data <- read.table(input_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
data$ModuleLabel <- factor(data$ModuleLabel)
data$MetamoduleLabel <- factor(data$MetamoduleLabel)

# Load GAF file for GO annotations
gaf_data <- read.delim(gaf_file, header=FALSE, comment.char="!")
colnames(gaf_data) <- c("db", "db_object_id", "db_object_symbol", "qualifier", "go_id", 
                        "db_reference", "evidence_code", "with_from", "aspect", 
                        "db_object_name", "db_object_synonym", "db_object_type", 
                        "taxon", "date", "assigned_by", "annotation_extension", 
                        "gene_product_form_id")

# Subset by ontology (BP = Biological Process)
gaf_BP <- subset(gaf_data, aspect == "P")

# GO enrichment function
perform_enrichment <- function(gaf_subset, module, gene_list, ontology_name) {
  go_terms <- setNames(gaf_subset$go_id, gaf_subset$db_object_symbol)
  gene2go <- split(go_terms, names(go_terms))
  
  # Perform enrichment
  ego <- enricher(
    gene = gene_list,
    TERM2GENE = data.frame(
      term = unlist(gene2go),
      gene = rep(names(gene2go), sapply(gene2go, length))
    ),
    universe = unique(gaf_subset$db_object_symbol),
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if (nrow(ego) == 0) {
    cat("WARNING: No enrichment found for module", module, "\n")
    return(NULL)
  }
  
  # Add GO term descriptions
  go_desc <- AnnotationDbi::select(GO.db, 
                                   keys = unique(ego@result$ID), 
                                   columns = c("GOID", "TERM"), 
                                   keytype = "GOID")
  ego@result <- ego@result %>% 
    left_join(go_desc, by = c("ID" = "GOID")) %>%
    mutate(Description = TERM,
           Description_plot = paste0(TERM, " (", ID, ")"))
  
  # Save results
  dir.create(file.path(output_dir, "04-go_enrichment", ontology_name), 
             showWarnings = FALSE, recursive = TRUE)
  
  result <- ego@result[order(ego@result$p.adjust), ]
  topResult <- if(nrow(result) > 10) result[1:10, ] else result
  
  
  # Create plot
  p <- dotplot(ego, showCategory = 20, 
               title = paste0("GO (", ontology_name, ") Enrichment for Module ", module),
               font.size = 8) +
    theme(
      plot.title = element_text(hjust = 0, face = 'bold', size = 14),
      panel.background = element_rect(fill = "#f0f0f0", color = "grey80"),
      panel.grid.major.y = element_line(color = "grey30", linetype = "dotted"),
      panel.grid.major.x = element_line(color = "white", linetype = "solid"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "04-go_enrichment", ontology_name, 
                   paste0(module, "_enrichment_plot.png")), 
         p, width = 10, height = 8)
  
  return(ggplotGrob(p))
}

# Analyze specific modules
modules_to_analyze <- c(16, 19, 22)
plots <- list()

for(module in modules_to_analyze) {
  gene_list <- data[data$ModuleLabel == module, "Gene"]
  if(length(gene_list) > 0) {
    plots[[as.character(module)]] <- perform_enrichment(gaf_BP, module, gene_list, "BP")
  }
}

# Remove NULL plots
plots <- plots[!sapply(plots, is.null)]

# Create combined plot if we have results
if(length(plots) > 0) {
  combined_plot <- plot_grid(
    plotlist = plots,
    ncol = 1,
    align = 'v',
    axis = 'l'
  )
  
  # Save combined plot
  dir.create(file.path(output_dir, "04-go_enrichment", "BP"), 
             showWarnings = FALSE, recursive = TRUE)
  
  ggsave(file.path(output_dir, "04-go_enrichment", "BP", "combined_enrichment_plot.png"), 
         combined_plot, width = 12, height = 15, dpi = 300)

  cat("Combined plot saved successfully\n")
}

