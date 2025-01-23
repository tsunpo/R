#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
cell_type <- args[1]

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
wd.de.data  <- file.path(wd.de, "data_ALL3_500_QC")
wd.de.plots <- file.path(wd.de, "plots_ALL3_500_QC")

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
# DESeq2
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
library(DESeq2)

nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Get the list of unique cell types
#cell_types <- levels(Idents(so.integrated))
# Initialize a list to store results for each cell type
#results_list <- list()

print(cell_type)

#for (cell_type in cell_types) {
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
	  #head(results_age[order(results_age$pvalue), ])
	  # Store the results in the list with the cell type as the name
	  #results_list[[cell_type]] <- results_age
#}
save(dds, file=file.path(wd.de.data, paste0("DESeq_~age_", cell_type, ".RData")))
