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

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
library(DoubletFinder)

samples <- SplitObject(so.integrated, split.by = "sample")

# Pre-processing steps for DoubletFinder
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- NormalizeData(so.integrated)
so.integrated <- FindVariableFeatures(so.integrated)
so.integrated <- ScaleData(so.integrated)
so.integrated <- RunPCA(so.integrated)

# Find optimal pK value (sweep parameter)
sweep_res <- paramSweep_v3(so.integrated, PCs = 1:10, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
optimal_pK <- find.pK(sweep_stats)$pK

# Estimate homotypic doublets (if clustering information is available)
homotypic_prop <- modelHomotypic(Idents(so.integrated))

# Define expected doublet rate and number of artificial doublets
nExp <- round(0.075 * ncol(so.integrated))  # Example: 7.5% doublet rate
nExp_adj <- round(nExp * (1 - homotypic_prop))

# Run DoubletFinder
so.integrated <- doubletFinder_v3(
	  so.integrated,
	  PCs = 1:10,
	  pN = 0.25,
	  pK = optimal_pK,
	  nExp = nExp_adj,
	  reuse.pANN = FALSE,
	  sct = FALSE
)

# Label doublets
Idents(so.integrated) <- "DF.classifications_0.25_0.01_50"






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

