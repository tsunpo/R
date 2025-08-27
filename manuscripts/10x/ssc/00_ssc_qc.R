# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 08/08/25
# =============================================================================

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch125/casm/staging/team294/ty2"   ## ty2@farm22
#wd <- "/Users/ty2/Work/sanger/ty2"                   ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE)
wd.rna.raw <- file.path(wd.rna, "ngs", "10x")

wd.de    <- file.path(wd.rna, "analysis", paste0(base, ""))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
dir.create(file.path(wd.de.data, "00_QC"), showWarnings = FALSE)
dir.create(file.path(wd.de.plots), showWarnings = FALSE)

samples0 <- read.table(file.path(wd.rna.raw, "GRCh38-2020-A_mp29.list"), header=F, sep="\t")
samples1 <- read.table(file.path(wd.rna.raw, "GRCh38-2020-A_lm26_ST.list"), header=F, sep="\t")

samples <- rbind(samples0, samples1)
rownames(samples) <- samples$V3

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 00_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#library(sctransform)

for (s in 1:nrow(samples)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
	  data <- Read10X_h5(file.path(wd.rna.raw, "GRCh38-2020-A", samples$V3[s], "filtered_feature_bc_matrix.h5"))
	  so <- CreateSeuratObject(counts=data, project="SSC", min.cells=3, min.features=200)
	  
	  # QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  
	  file.name <- file.path(wd.de.data, "00_QC", paste0(samples$V3[s], "_VlnPlot"))
	  pdf(paste0(file.name, ".pdf"), width=10)
	  VlnPlot(so, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
	  dev.off()
	  
	  plot1 <- FeatureScatter(so, feature1="nCount_RNA", feature2="percent.mt")
	  plot2 <- FeatureScatter(so, feature1="nCount_RNA", feature2="nFeature_RNA")
	  file.name <- file.path(wd.de.data, "00_QC", paste0(samples$V3[s], "_VlnPlot_plot1+plot2"))
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
colnames <- c("PD_ID", "genes", "cells", "cells_QC")
filtered <- data.frame(matrix(ncol = length(colnames), 
                              nrow = nrow(samples)))
colnames(filtered) <- colnames
filtered$PD_ID <- rownames(samples)
rownames(filtered) <- rownames(samples)

for (s in 1:nrow(samples)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
   data <- Read10X_h5(file.path(wd.rna.raw, "GRCh38-2020-A", samples$V3[s], "filtered_feature_bc_matrix.h5"))
   so <- CreateSeuratObject(counts=data, project=samples$V3[s], min.cells=3, min.features=200)
   filtered[s, 2] <- nrow(so)
   filtered[s, 3] <- ncol(so)
   
  	# QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  so <- subset(so, subset=nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA < 50000 & percent.mt < 5)
	  filtered[s, 4] <- ncol(so)
}
write.table(filtered, 
            file.path(wd.de.data, "ssc_filtered.txt"), 
            col.names = TRUE, 
            row.names = FALSE, 
            sep = "\t")
save(samples, filtered, file=file.path(wd.de.data, "ssc_filtered.RData"))

# -----------------------------------------------------------------------------
# Standard Seurat pre-processing workflow (SCT)
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered.RData"))
samples.filtered <- samples[subset(filtered, cells > 200)$PD_ID,]   # Remove PD53622b_2N and PD53622b_M1
samples.filtered$V8 <- mapply(x = 1:nrow(samples.filtered), function(x) unlist(strsplit(samples.filtered$V3[x], "_"))[2])
samples.filtered$V9 <- 1
samples.filtered$V9[grep("M", samples.filtered$V8)] <- 2
samples.filtered$V9[grep("ST", samples.filtered$V8)] <- 3

so.list <- c()
ids = c()
colnames <- c("PD_ID", "genes", "cells")
normalised <- data.frame(matrix(ncol = length(colnames), 
                                nrow = nrow(samples.filtered)))
colnames(normalised) <- colnames
normalised$PD_ID <- rownames(samples.filtered)
rownames(normalised) <- rownames(samples.filtered)

for (s in 1:nrow(samples.filtered)) {
	  # Initialize the Seurat object with the raw (non-normalized data)
	  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
   data <- Read10X_h5(file.path(wd.rna.raw, "GRCh38-2020-A", samples.filtered$V3[s], "filtered_feature_bc_matrix.h5"))
   so <- CreateSeuratObject(counts=data, project=samples.filtered$V3[s], min.cells=3, min.features=200)
	
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
	  ids = c(ids, samples.filtered$V3[s])
}
write.table(normalised, 
            file.path(wd.de.data, "ssc_filtered_normalised.txt"), 
            col.names = TRUE, 
            row.names = FALSE, 
            sep = "\t")
save(filtered, normalised, samples, samples.filtered, so.list, ids, file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))

# Merge Based on Normalized Data
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
so.merged <- merge(x=so.list[[1]], y=so.list[-1], add.cell.ids=ids, project="SSC", merge.data=T)

ids <- c()
for (s in 1:nrow(samples.filtered)) {
	  ids <- c(ids, rep(samples.filtered$V3[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

ages <- c()
for (s in 1:nrow(samples.filtered)) {
	  ages <- c(ages, rep(samples.filtered$V4[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

n2s <- c()
for (s in 1:nrow(samples.filtered)) {
	  n2s <- c(n2s, rep(samples.filtered$V8[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

batches <- c()
for (s in 1:nrow(samples.filtered)) {
	  batches <- c(batches, rep(samples.filtered$V9[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

so.merged@meta.data$sample.id <- ids
so.merged@meta.data$age <- ages
#so.merged@meta.data$age <- factor(so.merged@meta.data$age, levels = c("25","37","48","57","60","71"))
so.merged@meta.data$n2 <- n2s
so.merged@meta.data$batch <- batches
head(so.merged@meta.data)

save(filtered, normalised, samples, samples.filtered, so.merged, ids, ages, n2s, batches, file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))

# -----------------------------------------------------------------------------
# QC and selecting cells for further analysis
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#setup-the-seurat-object
# -----------------------------------------------------------------------------
library(qs)

samples2 <- read.table(file.path(wd.rna.raw, "GRCh38-2020-A_arussell.list"), header=F, sep="\t")
so <- qread(file.path(wd.rna.raw, "GRCh38-2020-A", "05_SEURATS_240523_VL00297_289_AAFJN3FM5_SI-TT-B1_seurat.qs"))
#so <- CreateSeuratObject(counts=data@assays$RNA@layers$counts, project=samples0$V3[1], min.cells=3, min.features=200)

#so[["pca"]] <- NULL
#so[["umap"]] <- NULL
so[["spatial"]] <- NULL  # If 'spatial' is an assay or similar slot
#so[["RNA"]]@layers$scale.data <- NULL
#so[["RNA"]]@layers$data <- NULL

# QC and selecting cells for further analysis
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
nrow(so)
# [1] 36601
ncol(so)
# [1] 17583
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA < 50000 & percent.mt < 5)

# Normalizing the data
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#normalizing-the-data
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
ncol(so)
# [1] 15890

samples.2 <- samples.2[, -c(5,6,7)]
rownames(samples.2) <- samples.2$V3
samples.2$V8 <- "A"
samples.2$V9 <- 4
so.list.2 <- c(so)
ids.2 <- c("VL00297")
save(samples.2, so.list.2, ids.2, file=file.path(wd.de.data, "ssc_filtered_normalised.2.RData"))
