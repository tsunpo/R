#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 25/07/24; 14/03/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Graphics.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg38.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/scRNA")
wd.rna.raw <- file.path(wd.rna, "10x")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

#samples0 <- readTable(file.path(wd.rna.raw, "scRNA_GRCh38-2020.list"), header=F, rownames=3, sep="\t")
#samples1 <- readTable(file.path(wd.rna.raw, "scRNA_homemade_ref.list"), header=F, rownames=3, sep="\t")
#samples1 <- samples1[rownames(samples0),]

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(pheatmap)

# -----------------------------------------------------------------------------
# To identify six SPG states
# -----------------------------------------------------------------------------
nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_annotated.RData")))

# -----------------------------------------------------------------------------
# DESeq2
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Get the list of unique cell types
cell_types <- levels(Idents(so.integrated))
# Initialize a list to store results for each cell type
results_list <- list()

for (cell_type in cell_types) {
	  # Step 0: Subset the Seurat object to include only the current cell type
	  so_subset <- subset(so.integrated, idents = cell_type)
	
	  # Step 1: Extract count data and metadata
	  counts <- GetAssayData(so_subset, layer = "counts")
	  # Filter out genes with no counts in any samples
	  non_zero_genes <- rowSums(counts > 0) > 0
	  counts_filtered <- counts[non_zero_genes, ]
	  
	  metadata <- so_subset@meta.data
	  # Ensure the metadata columns of interest are present
	  # Assuming 'age' and 'batch' are stored in metadata
	  metadata$age <- as.numeric(metadata$age)  # Convert age to numeric if it isn't already
	  #metadata$batch <- as.factor(metadata$batch)  # Ensure batch is a factor
	  # Center and scale the age variable
	  #metadata$age_scaled <- scale(metadata$age)
	  
	  # Step 2: Create a DESeqDataSet
	  dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
					  																										colData = metadata,
			  																												design = ~ age)
	  # Use the poscounts method for size factor estimation
	  dds <- estimateSizeFactors(dds, type = "poscounts")
	  
  	# Step 4: Run DESeq2 with the LRT to test for the effect of age
  	dds <- DESeq(dds)
	  # Step 5: Extract the results for the age effect
	  results_age <- results(dds, name = "age")
	  # View the top results
	  head(results_age[order(results_age$pvalue), ])
	  # Store the results in the list with the cell type as the name
	  results_list[[cell_type]] <- results_age
	  
	  # Assume results_age is the DESeq2 results object for the effect of age
	  # Convert to a data frame for ggplot2
	  result_df <- as.data.frame(results_age)
	  result_df$gene <- rownames(result_df)
	  result_df$color <- ifelse(result_df$log2FoldChange < -0.05, "blue",
	  																										ifelse(result_df$log2FoldChange > 0.05, "red", "grey"))
	  
	  # Create a volcano plot
	  volcano_plot <- ggplot(result_df, aes(x = log2FoldChange, y = -log10(padj))) +
	  	  geom_point(aes(color = color)) +
	  	  scale_color_manual(values = c(blue, "grey", red)) +
	  	  theme_minimal() +
	  	  labs(title = paste0(cell_type),
	  			  			x = "Log2 Fold Change",
	  				  		y = "-log10 Adjusted P-value") +
	  	  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
	    	geom_vline(xintercept = 0.05, linetype = "dashed", color = red) +
	    	geom_vline(xintercept = -0.05, linetype = "dashed", color = blue) +
	    	geom_text_repel(data = subset(result_df, abs(log2FoldChange) > 0.05),
	  														  			aes(label = gene), size = 3, max.overlaps = 30)
	  
	  # Print the volcano plot
	  pdf(file = file.path(wd.de.plots, paste0("volcano_DeSeq_~age_", cell_type, ".pdf")), width = 5, height = 5)
	  print(volcano_plot)
	  dev.off()
}
save(results_list, file=file.path(wd.de.data, paste0("DESeq_~age.RData")))
