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
wd.de.data  <- file.path(wd.de, paste0("data_ALL4_500_QC"))
wd.de.plots <- file.path(wd.de, paste0("plots_ALL4_500_QC"))

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
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
#res <- 0.6
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_res=0.6_", nfeatures, "_20_100.RData")))

## Find original raw count data
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_500",           "ssc_filtered_normalised.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500",      "ssc_filtered_normalised.1.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500",      "ssc_filtered_normalised.2.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_June_500", "ssc_filtered_normalised.3.RData"))

# Remove items 2 and 8 from the list (where there are only 130 and 92 cells in the samples)
so.list <- so.list[-c(2, 8)]
so.list <- c(so.list, so.list.1, so.list.2, so.list.3)

# Extract raw counts from so.list and merge them
# Step 1: Find the union of all gene names
all_genes <- Reduce(union, lapply(so.list, function(x) rownames(GetAssayData(x, layer = "counts"))))

# Step 2: Standardize the row names of all datasets
raw_counts_list <- lapply(so.list, function(x) {
	  counts <- GetAssayData(x, layer = "counts")
	  missing_genes <- setdiff(all_genes, rownames(counts))
	  if (length(missing_genes) > 0) {
		    zero_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(length(missing_genes), ncol(counts)))
		    rownames(zero_matrix) <- missing_genes
		    counts <- rbind(counts, zero_matrix)
	  }
	  counts[all_genes, ]  # Reorder rows
})

# Step 3: Combine all standardized matrices
combined_raw_counts <- do.call(cbind, raw_counts_list)

# Step 4: Assign raw counts to so.integrated
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- subset(so.integrated, cells = colnames(combined_raw_counts))
so.integrated[["RNA"]] <- CreateAssayObject(counts = combined_raw_counts)

# Step 5: Verify
slotNames(GetAssay(so.integrated, assay = "RNA"))

# From so.integrated
integrated_samples <- colnames(GetAssayData(so.integrated, layer = "counts"))
# From so.list
list_samples <- lapply(so.list, function(x) colnames(GetAssayData(x, layer = "counts")))
combined_list_samples <- unlist(list_samples)  # Flatten the list of samples from so.list

# Check if the sample orders are identical
identical(integrated_samples, combined_list_samples)

# -----------------------------------------------------------------------------
# 
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
library(DoubletFinder)

# Use clustering information to estimate the proportion of homotypic doublets.
homotypic_prop <- modelHomotypic(Idents(so.integrated))

# Step 1: Find the optimal pK value (can skip if known)
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- UpdateSeuratObject(so.integrated)

sweep_res <- paramSweep_v3(so.integrated, PCs = 1:10, sct = FALSE)  # Adjust PCs
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
optimal_pK <- find.pK(sweep_stats)





# Identify the clusters you want to keep
clusters_to_remove <- c(1)
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)

ssc0 <- as.data.frame(table(so.integrated@meta.data$sample.id))
ssc <- ssc0
write.table(ssc0, file=file.path(wd.de.data, "ssc_samples_-1_>0.txt"), row.names=T, col.names=T, quote=F, sep='\t')
save(ssc0, ssc, so.integrated, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_res=0.6_5000_-1_>0.RData"))

# Access the metadata
metadata <- so.integrated@meta.data
# Count the number of cells per sample
cell_counts <- metadata %>%
	  group_by(sample.id) %>%
	  summarise(cell_count = n()) %>%
	  ungroup()

# Filter samples with more than 100 cells
samples_to_keep <- cell_counts %>%
	  filter(cell_count > 100) %>%
	  pull(sample.id)

# Subset the Seurat object to keep only the samples with more than 100 cells
so.integrated <- subset(so.integrated, cells = rownames(metadata[metadata$sample.id %in% samples_to_keep, ]))

# Check the result
ssc <- as.data.frame(table(so.integrated@meta.data$sample.id))
write.table(ssc, file=file.path(wd.de.data, "ssc_samples_-1_>100.txt"), row.names=T, col.names=T, quote=F, sep='\t')
save(ssc0, ssc, so.integrated, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_res=0.6_5000_-1_>100.RData"))

