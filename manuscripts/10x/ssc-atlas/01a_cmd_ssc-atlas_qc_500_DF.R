# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 09/05/25
# =============================================================================

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/10x")
wd.rna.raw <- file.path(wd.rna, "atlas")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-atlas"))
wd.de.data  <- file.path(wd.de, "data_500")
wd.de.plots <- file.path(wd.de, "plots_500")

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

load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
library(DoubletFinder)

classification_cols <- list()
for (i in seq_along(so.list)) {
	  sample <- so.list[[i]]
	  raw_counts <- GetAssayData(sample, layer = "counts", assay = "RNA")
	  sample[["RNA"]] <- CreateAssayObject(counts = raw_counts)
	  DefaultAssay(sample) <- "RNA"
	  
	  # Preprocessing
	  sample <- NormalizeData(sample, assay = "RNA")
	  sample <- FindVariableFeatures(sample)
	  sample <- ScaleData(sample)
	  sample <- RunPCA(sample)
	
	  # Perform clustering to generate meaningful identities for modelHomotypic
	  sample <- FindNeighbors(sample, dims = 1:10)
	  sample <- FindClusters(sample, resolution = 0.5)  # Adjust resolution as needed
	  
	  # Find optimal pK
	  sweep_res <- paramSweep_v3(sample, PCs = 1:10, sct = FALSE)
	  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
	  # Save the plot
	  png(file.path(wd.de.data, paste0("pK_optimization_", samples.filtered$Sample[i], ".png")), width = 800, height = 600)
	  optimal_pK <- find.pK(sweep_stats)$pK
	  dev.off()
	  optimal_pK <- as.numeric(as.character(optimal_pK))
	  
	  # Estimate homotypic doublet proportion
	  homotypic_prop <- modelHomotypic(Idents(sample))
	  nExp <- round(0.075 * ncol(sample))  # Adjust expected doublet rate
	  nExp_adj <- round(nExp * (1 - homotypic_prop))
	
	  # Run DoubletFinder
	  sample <- doubletFinder_v3(
		    sample,
		    PCs = 1:10,
		    pN = 0.25,
		    pK = optimal_pK,
		    nExp = nExp_adj,
		    reuse.pANN = FALSE,
		    sct = FALSE
	  )
	  
	  # Extract pK corresponding to the maximum BCmetric
	  pK_stats <- find.pK(sweep_stats)
	  optimal_pK <- pK_stats$pK[which.max(pK_stats$BCmetric)]
	  classification_col <- paste0("DF.classifications_0.25_", optimal_pK, "_", nExp_adj)
	  Idents(sample) <- sample@meta.data[[classification_col]]
	  
	  # Store sample ID and classification column
	  classification_cols[[samples.filtered$Sample[i]]] <- classification_col
	  
	  so.list[[i]] <- sample
}

save(samples.filtered, ids, so.list, classification_cols, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.RData")))

# -----------------------------------------------------------------------------
# Last Modified: 14/02/25
# -----------------------------------------------------------------------------
# Check where "PD53621b_2N" appears in sample_ids
#sample_index <- which(sample_ids == "PD53621b_2N")
# Find the corresponding classification column using the same index
#classification_col <- classification_cols[sample_index]
#pANN_col <- gsub("DF.classifications_", "pANN_", classification_col)

# Subset the sample
#sample_cells <- so.integrated@meta.data[so.integrated@meta.data$orig.ident == "PD53621b_2N", ]
# Get classification values for these 10 cells
#original_scores <- sample_cells[1:10, pANN_col]

# Initialize a dataframe to store mapping information
mapping_table <- data.frame(Sample_ID = character(), pANN_Column = character(), stringsAsFactors = FALSE)
sample_ids <- unique(so.integrated@meta.data$sample.id)

# Iterate through each sample
for (i in seq_along(sample_ids)) {
   sample_id <- sample_ids[i]
 
   # Get the corresponding classification column
   classification_col <- classification_cols[i]
   pANN_col <- gsub("DF.classifications_", "pANN_", classification_col)
 
   # Append to mapping table
   mapping_table <- rbind(mapping_table, data.frame(Sample_ID = sample_id, pANN_Column = pANN_col))
}
# Print the mapping table
print(mapping_table)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Step 1: Assign the correct pANN values to each cell in `so.integrated`
so.integrated@meta.data$pANN_score <- NA  # Initialize a column for pANN scores

# Iterate through each sample to assign the correct pANN column values
for (i in seq_along(mapping_table$Sample_ID)) {
   sample_id <- mapping_table$Sample_ID[i]
   pANN_col <- mapping_table$pANN_Column[i]
 
   # Subset cells from the corresponding sample and assign the correct pANN values
   so.integrated@meta.data$pANN_score[so.integrated@meta.data$sample.id == sample_id] <- 
      so.integrated@meta.data[so.integrated@meta.data$sample.id == sample_id, pANN_col]
}

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Identify columns to remove
columns_to_remove <- grep("^pANN_0.25|^DF.classifications_0.25|^integrated_snn_res", colnames(so.integrated@meta.data), value = TRUE)

# Remove the identified columns
so.integrated@meta.data <- so.integrated@meta.data[, !colnames(so.integrated@meta.data) %in% columns_to_remove]
print(colnames(so.integrated@meta.data))

