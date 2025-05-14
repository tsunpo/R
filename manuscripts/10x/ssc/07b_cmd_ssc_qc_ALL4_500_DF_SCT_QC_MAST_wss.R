#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
nfeatures <- as.numeric(args[1])
kms <- as.numeric(args[2])

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
handbooks  <- c("Commons.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

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
wd.de.data  <- file.path(wd.de, "data_ALL4_500_QC")
wd.de.plots <- file.path(wd.de, "plots_ALL4_500_QC")

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
# K-means Clustering of Genes
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_markers_MAST_SSC.RData")))

filtered_markers <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
de_gene_list <- unique(filtered_markers$gene)
length(de_gene_list)
# [1] 2706

# Extract expression data from Seurat object
DefaultAssay(so.integrated) <- "RNA"
#so.integrated <- JoinLayers(so.integrated, assays = "RNA")
expression_matrix <- GetAssayData(so.integrated[["RNA"]], layer = "data")
expression_matrix <- as.matrix(expression_matrix)

expression_matrix_de <- expression_matrix[de_gene_list, ]

# Determine optimal number of clusters (Elbow Method)
wss <- sapply(1:15, function(k){
	  kmeans(expression_matrix_de, centers = k, nstart = 10)$tot.withinss
})
save(so.integrated, markers, expression_matrix, expression_matrix_de, wss, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_markers_MAST_SSC_wss.RData")))

# -----------------------------------------------------------------------------
# K-means Clustering of Genes
# -----------------------------------------------------------------------------
# Set the number of clusters (k)
#kms <- 5  # You can adjust this based on your needs
# Perform K-means clustering on genes (rows)
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(expression_matrix_de, centers = kms)
# Assign cluster labels to the genes
gene_clusters <- kmeans_result$cluster

save(so.integrated, markers, expression_matrix, expression_matrix_de, wss, kms, kmeans_result, gene_clusters, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_markers_MAST_SSC_wss_k=", kms, ".RData")))
