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
wd <- "/lustre/scratch127/casm/team294rr/ty2"    ## ty2@farm
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/10x")
wd.rna.raw <- file.path(wd.rna, "atlas")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-atlas"))
wd.de.data  <- file.path(wd.de, "data_500")
wd.de.plots <- file.path(wd.de, "plots_500")
wd.de.plots <- file.path(wd.de.plots, "00_QC")
dir.create(wd.de.plots, showWarnings = FALSE)

samples <- read.table(file.path(wd.rna, "atlas.list"), header=T, sep="\t")
rownames(samples) <- samples$Sample
samples <- samples[grepl("^GSE", samples$Batch), ]

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 00_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

for (s in 1:nrow(samples)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
	  data <- Read10X_h5(file.path(wd.rna.raw, samples$Sample[s], "filtered_feature_bc_matrix.h5"))
	  so <- CreateSeuratObject(counts=data, project=samples$Batch[s], min.cells=3, min.features=200)
	  
	  # QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  
	  file.name <- file.path(wd.de.plots, paste0(samples$Sample[s], "_VlnPlot"))
	  pdf(paste0(file.name, ".pdf"), width=10)
	  VlnPlot(so, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
	  dev.off()
	  
	  plot1 <- FeatureScatter(so, feature1="nCount_RNA", feature2="percent.mt")
	  plot2 <- FeatureScatter(so, feature1="nCount_RNA", feature2="nFeature_RNA")
	  file.name <- file.path(wd.de.plots,  paste0(samples$Sample[s], "_VlnPlot_plot1+plot2"))
	  pdf(paste0(file.name, ".pdf"), width=10)
	  plot1 + plot2
	  dev.off()
}

# -----------------------------------------------------------------------------
# QC and selecting cells for further analysis
# 00_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#setup-the-seurat-object
# -----------------------------------------------------------------------------
colnames <- c("PD_ID", "genes", "cells")
filtered <- toTable(0, length(colnames), nrow(samples), colnames)
filtered$PD_ID <- rownames(samples)
rownames(filtered) <- rownames(samples)

for (s in 1:nrow(samples)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
   data <- Read10X_h5(file.path(wd.rna.raw, samples$Sample[s], "filtered_feature_bc_matrix.h5"))
   so <- CreateSeuratObject(counts=data, project=samples$Batch[s], min.cells=3, min.features=200)
	
  	# QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  so <- subset(so, subset=nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA < 50000 & percent.mt < 5)

	  filtered[s, 2] <- nrow(so)
	  filtered[s, 3] <- ncol(so)
}
writeTable(filtered, file.path(wd.de.data, "ssc-atlas_filtered.txt"), colnames=T, rownames=F, sep="\t")
save(samples, filtered, file=file.path(wd.de.data, "ssc-atlas_filtered.RData"))

# -----------------------------------------------------------------------------
# Standard Seurat pre-processing workflow (SCT)
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered.RData"))
samples.filtered <- samples[subset(filtered, cells > 0)$PD_ID,]

so.list <- c()
ids = c()
genes <- c()
colnames <- c("PD_ID", "genes", "cells")
normalised <- toTable(0, length(colnames), nrow(samples.filtered), colnames)
normalised$PD_ID <- rownames(samples.filtered)
rownames(normalised) <- rownames(samples.filtered)

for (s in 1:nrow(samples.filtered)) {
	  # Initialize the Seurat object with the raw (non-normalized data)
	  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
   data <- Read10X_h5(file.path(wd.rna.raw, samples$Sample[s], "filtered_feature_bc_matrix.h5"))
   so <- CreateSeuratObject(counts=data, project=samples.filtered$Batch[s], min.cells=3, min.features=200)
	
	  # QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA < 50000 & percent.mt < 5)

	  # Apply sctransform normalization
	  # https://satijalab.org/seurat/articles/sctransform_vignette.html
	  #so <- SCTransform(so, vars.to.regress="percent.mt", verbose=F)
	  # Normalizing the data
	  # https://satijalab.org/seurat/articles/pbmc3k_tutorial#normalizing-the-data
	  so <- NormalizeData(so)
	  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
	  
	  normalised[s, 2] <- nrow(so)
	  normalised[s, 3] <- ncol(so)
	
	  so.list <- c(so.list, so)
	  ids = c(ids, samples.filtered$Sample[s])
	
  	if (length(genes) != 0) {
	    	genes <- intersect(genes, rownames(so))
	  } else {
		    genes <- rownames(so)
	  }
}
writeTable(normalised, file.path(wd.de.data, "ssc_filtered_normalised.txt"), colnames=T, rownames=F, sep="\t")
save(filtered, normalised, samples, samples.filtered, so.list, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))

# Merge Based on Normalized Data
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
so.merged <- merge(x=so.list[[1]], y=so.list[-1], add.cell.ids=ids, project="SSC", merge.data=T)

ids <- c()
for (s in 1:nrow(samples.filtered)) {
	  ids <- c(ids, rep(samples.filtered$Sample[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

ages <- c()
for (s in 1:nrow(samples.filtered)) {
	  ages <- c(ages, rep(samples.filtered$Age[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

#n2s <- c()
#for (s in 1:nrow(samples.filtered)) {
#	  n2s <- c(n2s, rep(samples.filtered$V8[s], ncol(so.list[[s]]@assays$RNA$counts)))
#}

batches <- c()
for (s in 1:nrow(samples.filtered)) {
	  batches <- c(batches, rep(samples.filtered$Batch[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

so.merged@meta.data$sample.id <- ids
so.merged@meta.data$age <- ages
#so.merged@meta.data$age <- factor(so.merged@meta.data$age, levels = c("25","37","48","57","60","71"))
#so.merged@meta.data$n2 <- n2s
so.merged@meta.data$batch <- batches
head(so.merged@meta.data)

so.merged@meta.data$orig.ident <- so.merged@meta.data$sample.id
save(filtered, normalised, samples, samples.filtered, so.merged, ids, ages, batches, file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))
