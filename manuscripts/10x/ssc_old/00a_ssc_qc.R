# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 14/03/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"    ## ty2@farm
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/scRNA")
wd.rna.raw <- file.path(wd.rna, "10x")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-de"))
wd.de.data  <- file.path(wd.de, "data_500")
wd.de.plots <- file.path(wd.de, "plots_500")

samples0 <- readTable(file.path(wd.rna, "scRNA_GRCh38-2020.list"), header=F, rownames=3, sep="\t")

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#library(sctransform)

for (s in 1:nrow(samples0)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
	  data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0$V1[s], "filtered_feature_bc_matrix"))
	  so <- CreateSeuratObject(counts=data, project=samples0$V3[s], min.cells=3, min.features=200)
	  
	  # QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  
	  file.name <- file.path(wd.de.plots, "01_QC", paste0(samples0$V3[s], "_VlnPlot"))
	  pdf(paste0(file.name, ".pdf"), width=10)
	  VlnPlot(so, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
	  dev.off()
	  
	  plot1 <- FeatureScatter(so, feature1="nCount_RNA", feature2="percent.mt")
	  plot2 <- FeatureScatter(so, feature1="nCount_RNA", feature2="nFeature_RNA")
	  file.name <- file.path(wd.de.plots, "01_QC", paste0(samples0$V3[s], "_VlnPlot_plot1+plot2"))
	  pdf(paste0(file.name, ".pdf"), width=10)
	  plot1 + plot2
	  dev.off()
}

# -----------------------------------------------------------------------------
# QC and selecting cells for further analysis
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#setup-the-seurat-object
# -----------------------------------------------------------------------------
colnames <- c("PD_ID", "genes", "cells")
filtered <- toTable(0, length(colnames), nrow(samples0), colnames)
filtered$PD_ID <- rownames(samples0)
rownames(filtered) <- rownames(samples0)

for (s in 1:nrow(samples0)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
  	data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0$V1[s], "filtered_feature_bc_matrix"))
  	so <- CreateSeuratObject(counts=data, project=samples0$V3[s], min.cells=3, min.features=200)
	
  	# QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  so <- subset(so, subset=nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA < 50000 & percent.mt < 5)

	  filtered[s, 2] <- nrow(so)
	  filtered[s, 3] <- ncol(so)
}
writeTable(filtered, file.path(wd.de.data, "ssc_filtered.txt"), colnames=T, rownames=F, sep="\t")
save(samples0, filtered, file=file.path(wd.de.data, "ssc_filtered.RData"))

# -----------------------------------------------------------------------------
# Standard Seurat pre-processing workflow (SCT)
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered.RData"))
samples0.filtered <- samples0[subset(filtered, cells > 0)$PD_ID,]
samples0.filtered$V8 <- mapply(x = 1:nrow(samples0.filtered), function(x) unlist(strsplit(samples0.filtered$V3[x], "_"))[2])
samples0.filtered$V9 <- 1
samples0.filtered$V9[grep("M", samples0.filtered$V8)] <- 2

so.list <- c()
ids = c()
genes <- c()
colnames <- c("PD_ID", "genes", "cells")
normalised <- toTable(0, length(colnames), nrow(samples0.filtered), colnames)
normalised$PD_ID <- rownames(samples0.filtered)
rownames(normalised) <- rownames(samples0.filtered)

for (s in 1:nrow(samples0.filtered)) {
	  # Initialize the Seurat object with the raw (non-normalized data)
	  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
	  data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0.filtered$V1[s], "filtered_feature_bc_matrix"))
	  so <- CreateSeuratObject(counts=data, project=samples0.filtered$V3[s], min.cells=3, min.features=200)
	
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
	  ids = c(ids, samples0.filtered$V3[s])
	
  	if (length(genes) != 0) {
	    	genes <- intersect(genes, rownames(so))
	  } else {
		    genes <- rownames(so)
	  }
}
writeTable(normalised, file.path(wd.de.data, "ssc_filtered_normalised.txt"), colnames=T, rownames=F, sep="\t")
save(filtered, normalised, samples0, samples0.filtered, so.list, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))

# Merge Based on Normalized Data
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
so.merged <- merge(x=so.list[[1]], y=so.list[-1], add.cell.ids=ids, project="SSC", merge.data=T)

ids <- c()
for (s in 1:nrow(samples0.filtered)) {
	  ids <- c(ids, rep(samples0.filtered$V3[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

ages <- c()
for (s in 1:nrow(samples0.filtered)) {
	  ages <- c(ages, rep(samples0.filtered$V4[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

n2s <- c()
for (s in 1:nrow(samples0.filtered)) {
	  n2s <- c(n2s, rep(samples0.filtered$V8[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

batches <- c()
for (s in 1:nrow(samples0.filtered)) {
	  batches <- c(batches, rep(samples0.filtered$V9[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

so.merged@meta.data$sample.id <- ids
so.merged@meta.data$age <- ages
so.merged@meta.data$age <- factor(so.merged@meta.data$age, levels = c("25","37","48","57","60","71"))
so.merged@meta.data$n2 <- n2s
so.merged@meta.data$batch <- batches
head(so.merged@meta.data)

save(filtered, normalised, samples0, samples0.filtered, so.merged, ids, ages, n2s, batches, file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))
