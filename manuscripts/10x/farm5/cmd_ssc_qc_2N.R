#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
s <- as.numeric(args[1])

install.packages('BiocManager')
BiocManager::install('glmGamPoi')

# install.packages("devtools")
devtools::install_github("lazappi/clustree")

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
handbooks  <- c("Commons.R", "Graphics.R")
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
wd.de.plots <- file.path(wd.de, "plots/01_QC")

samples0 <- readTable(file.path(wd.rna.raw, "scRNA_GRCh38-2020.list"), header=F, rownames=3, sep="\t")
samples1 <- readTable(file.path(wd.rna.raw, "scRNA_homemade_ref.list"), header=F, rownames=3, sep="\t")
samples1 <- samples1[rownames(samples0),]

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(clustree)

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
# -----------------------------------------------------------------------------
#library(glmGamPoi)

load(file=file.path(wd.de.data, "ssc_filtered.RData"))
load(file=file.path(wd.de.data, "ssc_filtered_normalised_2N+4N.RData"))
load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_2N+4N.RData"))

samples0$V8 <- mapply(x = 1:nrow(samples0), function(x) unlist(strsplit(samples0$V3[x], "_"))[2])
samples0.filtered <- samples0[subset(filtered, cells > 1000)$PD_ID,]
#samples0.filtered.2n <- subset(samples0.filtered, V8 == "2N")
samples0.filtered.2n <- subset(samples0.filtered, V8 != "M")
samples0.filtered.2n <- subset(samples0.filtered.2n, V8 != "M1")
samples0.filtered.2n <- subset(samples0.filtered.2n, V8 != "M2")
so.list <- c()
ids = c()
genes <- c()

colnames <- c("PD_ID", "genes", "cells")
normalised <- toTable(0, length(colnames), nrow(samples0.filtered.2n), colnames)
normalised$PD_ID <- rownames(samples0.filtered.2n)
rownames(normalised) <- rownames(samples0.filtered.2n)

for (s in 1:nrow(samples0.filtered.2n)) {
	  # Initialize the Seurat object with the raw (non-normalized data)
  	# https://satijalab.org/seurat/articles/pbmc3k_tutorial
	  data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0.filtered.2n$V1[s], "filtered_feature_bc_matrix"))
	  so <- CreateSeuratObject(counts=data, project=samples0.filtered.2n$V3[s], min.cells=3, min.features=200)
	
	  # QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  so <- subset(so, subset=nFeature_RNA > 1000 & nFeature_RNA < 10000 & nCount_RNA > 2000 & nCount_RNA < 50000 & percent.mt < 5)
	  
	  # Apply sctransform normalization
	  # https://satijalab.org/seurat/articles/sctransform_vignette.html
	  so <- SCTransform(so, vars.to.regress="percent.mt", verbose=F)
	  normalised[s, 2] <- nrow(so)
	  normalised[s, 3] <- ncol(so)
	  
	  so.list <- c(so.list, so)
	  ids = c(ids, samples0.filtered.2n$V3[s])
	  
	  if (length(genes) != 0) {
	  	  genes <- intersect(genes, rownames(so))
	  } else {
	  	  genes <- rownames(so)
	  }
}
#writeTable(normalised, file.path(wd.de.data, "ssc_filtered_normalised.txt"), colnames=T, rownames=F, sep="\t")
save(filtered, normalised, samples0, samples0.filtered, samples0.filtered.2n, so.list, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised_2N+4N.RData"))

# Merge Based on Normalized Data
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
so.merged.2n <- merge(x=so.list[[1]], y=so.list[-1], add.cell.ids=ids, project="SSC", merge.data=T)
save(filtered, normalised, samples0, samples0.filtered, samples0.filtered.2n, so.merged.2n, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_2N+4N.RData"))

# -----------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
# 02_UMAP
# https://satijalab.org/seurat/articles/multimodal_vignette
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))

ids <- c()
for (s in 1:nrow(samples0.filtered.2n)) {
	  ids <- c(ids, rep(samples0.filtered.2n$V3[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

ages <- c()
for (s in 1:nrow(samples0.filtered.2n)) {
	  ages <- c(ages, rep(samples0.filtered.2n$V4[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

so.merged.2n@meta.data$sample.id <- ids
so.merged.2n@meta.data$age       <- ages
so.merged.2n@meta.data$age <- factor(so.merged.2n@meta.data$age, levels = c("25","37","48","57","60"))
head(so.merged.2n@meta.data)

# Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA
DefaultAssay(so.merged.2n) <- "RNA"
DefaultAssay(so.merged.2n)
## [1] "RNA"

# perform visualization and clustering steps
so.merged.2n <- NormalizeData(so.merged.2n)
so.merged.2n <- FindVariableFeatures(so.merged.2n)
so.merged.2n
# An object of class Seurat 
# 68104 features across 46987 samples within 2 assays 
# Active assay: RNA (34615 features, 2000 variable features)
# 26 layers present: counts.PD53621b_2N, counts.PD53623b_2N, counts.PD53623b_4N, counts.PD53624b_2N, counts.PD53625b_2N, counts.PD53626b_2N, counts.PD53621b_M, counts.PD53623b_M, counts.PD53624b_M, counts.PD53625b_M, counts.PD53626b_M, counts.PD40746e_M1, counts.PD40746e_M2, data.PD53621b_2N, data.PD53623b_2N, data.PD53623b_4N, data.PD53624b_2N, data.PD53625b_2N, data.PD53626b_2N, data.PD53621b_M, data.PD53623b_M, data.PD53624b_M, data.PD53625b_M, data.PD53626b_M, data.PD40746e_M1, data.PD40746e_M2
# 1 other assay present: SCT

so.merged.2n <- ScaleData(so.merged.2n)
# Centering and scaling data matrix
# |======================================================================| 100%
so.merged.2n <- RunPCA(so.merged.2n, verbose = FALSE)

pdf('human_adult_firstPass.elbow_plot_RNA.pdf')
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.merged.2n, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.merged.2n[["pca"]]@stdev / sum(so.merged.2n[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp,file='human_adult_firstPass.elbow_PC_RNA_2N+4N.txt',row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
resolution.range <- seq(from = 0, to = 0.5, by = 0.05)

so.merged.2n <- FindNeighbors(so.merged.2n, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.merged.2n <- FindClusters(so.merged.2n, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged.2n <- RunUMAP(so.merged.2n, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)

pdf("DimPlot_UMAP_RNA_dim=16_2N+4N_RNA.pdf")
DimPlot(so.merged.2n, label = TRUE)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=16_grouped_by_sampleID_2N+4N_RNA.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=16_grouped_by_age_2N+4N_RNA.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()






# print the Clustree plot, so that we may later check which cluster resolution is appropriate
# in the next script, we will re-calculate the clusters according to our manually-chosen resolution
pdf('human_adult_firstPass.clustree_plot.pdf', paper="a4r")
clustree_fig<-clustree(so.merged.2n, prefix="RNA_snn_res.")
print(clustree_fig)
dev.off()




# print metadata
# note that testis.combined@meta.data$seurat_clusters == testis.combined@meta.data$integrated_snn_res.0.15, so we don't need to export the latter
so.merged.2n@meta.data$integrated_snn_res.0.15<-NULL
write.table(so.merged.2n@meta.data, file='human_adult.metadata.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# determine cluster markers
# we will later run enrichment analyses on these gene lists to assign a predicted function to each cluster
# IMPORTANT: we must set the DefaultAssay to RNA before running differential expression: "we don't recommend using the integrated matrix for differential expression" and "As a general rule, we always recommend performing DE on originally measured values - instead of on batch-corrected, imputed, etc. values. This ensures that the measurements that enter the DE test are indeed independent from each other, which is a requirement of any statistical DE test." (https://github.com/satijalab/seurat/issues/1057, https://github.com/satijalab/seurat/issues/1256 and https://github.com/satijalab/seurat/issues/2136)
DefaultAssay(so.merged.2n) <- 'RNA'

testis.markers <- FindAllMarkers(so.merged.2n, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(testis.markers,file='human_adult.cluster_markers.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

# save the Seurat object for later use
saveRDS(so.merged.2n, file = "so.merged.2n+4n.rds")

# convert the Seurat object to an h5ad object for visualisation with cellxgene
so.merged.2n <- NormalizeData(so.merged.2n)
so.merged.2n <- FindVariableFeatures(so.merged.2n, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(so.merged.2n)
testis <- ScaleData(so.merged.2n, features = all.genes)
sceasy::convertFormat(so.merged.2n, from="seurat", to="anndata", outFile='so.merged.2n+4n.h5ad')







DefaultAssay(so.merged.2n) <- "SCT"
# Error in PrepDR; told to run FindVariableFeatures after SCTransform?
# https://github.com/satijalab/seurat/issues/2852
# https://github.com/satijalab/seurat/issues/4145
obj.features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = 2000)
VariableFeatures(so.merged.2n[["SCT"]]) <- obj.features

so.merged.2n <- RunPCA(so.merged.2n, verbose = FALSE)

pdf('human_adult_firstPass.elbow_plot_SCT_2N+4N.pdf')
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.merged.2n, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.merged.2n[["pca"]]@stdev / sum(so.merged.2n[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp,file='human_adult_firstPass.elbow_PC_SCT_2N+4N.txt',row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."

resolution.range <- seq(from = 0, to = 0.5, by = 0.05)

so.merged.2n <- FindNeighbors(so.merged.2n, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.merged.2n <- FindClusters(so.merged.2n, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged.2n <- RunUMAP(so.merged.2n, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)

pdf("DimPlot_UMAP_RNA_dim=21_2N+4N_SCT.pdf")
DimPlot(so.merged.2n, label = TRUE)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=21_grouped_by_sampleID_2N+4N_SCT.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=21_grouped_by_age_2N+4N_SCT.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()









DefaultAssay(so.merged.2n) <- "SCT"
# Error in PrepDR; told to run FindVariableFeatures after SCTransform?
# https://github.com/satijalab/seurat/issues/2852
# https://github.com/satijalab/seurat/issues/4145
obj.features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = 2000)
VariableFeatures(so.merged.2n[["SCT"]]) <- obj.features

so.merged.2n <- RunPCA(so.merged.2n, verbose = FALSE)
so.merged.2n <- RunUMAP(so.merged.2n, dims = 1:20, verbose = FALSE)

so.merged.2n <- FindNeighbors(so.merged.2n, dims = 1:20, verbose = FALSE)
so.merged.2n <- FindClusters(so.merged.2n, verbose = FALSE)

pdf("DimPlot_UMAP_RNA_dim=20_2N_SCT.pdf")
DimPlot(so.merged.2n, label = TRUE)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=20_grouped_by_sampleID_2N_SCT.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=20_grouped_by_age_2N_SCT.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()



so.merged.2n <- RunPCA(so.merged.2n, verbose = FALSE)














so.merged.2n <- FindNeighbors(so.merged.2n, dims = 1:20)
# Computing nearest neighbor graph
# Computing SNN

so.merged.2n <- FindClusters(so.merged.2n, resolution = 0.8, verbose = FALSE)
so.merged.2n <- RunUMAP(so.merged.2n, dims = 1:20)
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# 15:15:04 UMAP embedding parameters a = 0.9922 b = 1.112
# 15:15:04 Read 46987 rows and found 30 numeric columns
# 15:15:04 Using Annoy for neighbor search, n_neighbors = 30
# 15:15:04 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 15:15:15 Writing NN index file to temp file /tmp/RtmpFKq9At/file1e6e714df04e1
# 15:15:15 Searching Annoy index using 1 thread, search_k = 3000
# 15:15:46 Annoy recall = 100%
# 15:15:46 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 15:15:49 Initializing from normalized Laplacian + noise (using RSpectra)
# 15:15:58 Commencing optimization for 200 epochs, with 2067094 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 15:16:24 Optimization finished
#save(so.merged.2n, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_RNA_PCA_UMAP_2N.RData"))

#pdf('DimPlot_UMAP_RNA_dim=30.pdf')
#DimPlot(so.merged, label = TRUE)
#dev.off()

pdf("DimPlot_UMAP_RNA_dim=17_2N.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=17_grouped_by_sampleID_2N.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=17_grouped_by_age_2N.pdf")
tplot = DimPlot(so.merged.2n, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()




# -----------------------------------------------------------------------------
# Perform dimensionality reduction by PCA and UMAP embedding
# 02_UMAP
# https://satijalab.org/seurat/articles/sctransform_vignette.html
# -----------------------------------------------------------------------------
Assays(so.merged)
# [1] "RNA" "SCT"
DefaultAssay(so.merged)
# [1] "SCT"

#features <- FindVariableFeatures(so.merged)     # Error: SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-
features <- SelectIntegrationFeatures(so.list)   # https://github.com/satijalab/seurat/issues/5761
so.merged <- RunPCA(so.merged, features = features, verbose = FALSE)
#so.merged <- RunUMAP(so.merged, dims = 1:20)
#DimPlot(so.merged, label = TRUE)

pdf('human_adult_firstPass.elbow_plot.pdf')
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.merged, ndims = 50)
dev.off()

# https://hbctraining.github.io/scRNA-seq_online/lessons/03_SC_quality_control-setup.html
head(so.merged@meta.data)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_SCT_PCA_UMAP.RData"))

pct <- so.merged[["pca"]]@stdev / sum(so.merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file='human_adult_firstPass.elbow_PC.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
resolution.range <- seq(from = 0, to = 0.5, by = 0.05)

so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)

pdf('DimPlot_UMAP.pdf')
DimPlot(so.merged, label = TRUE)
dev.off()









so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged))









# -----------------------------------------------------------------------------
# Perform dimensionality reduction by PCA and UMAP embedding
# 02_UMAP
# https://satijalab.org/seurat/articles/sctransform_vignette.html
# -----------------------------------------------------------------------------
features <- SelectIntegrationFeatures(object.list = so.list, normalization.method = "SCT", nfeatures = 5000)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

so.merged <- RunPCA(so.merged, features = features)



so.anchors <- FindIntegrationAnchors(object.list = so.list, anchor.features = features, reduction = "rpca") # we create as our standard reference the Di Persio, Sohni and Zhao samples as these contain the greatest number of cells
so.combined <- IntegrateData(anchorset = so.anchors, normalization.method = "SCT") # see https://github.com/satijalab/seurat/issues/3930 for discussion of k.weight









features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = 5000)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)
testis.anchors <- FindIntegrationAnchors(object.list = so.list, anchor.features = features, normalization.method = "SCT", reference=c(1,2,3,7,8,10,11,12), reduction = "rpca") # we create as our standard reference the Di Persio, Sohni and Zhao samples as these contain the greatest number of cells
testis.combined <- IntegrateData(anchorset = testis.anchors, normalization.method = "SCT") # see https://github.com/satijalab/seurat/issues/3930 for discussion of k.weight



# Centering and scaling data matrix
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#identification-of-highly-variable-features-feature-selection
#all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = genes)   ## Only using 23246 genes shared across all samples
so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged))
so.merged <- RunUMAP(so.merged, dims = 1:30, verbose = FALSE)

save(filtered, normalised, samples0, samples0.filtered, so.merged, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_scaled_23246.RData"))



# These are now standard steps in the Seurat workflow for visualization and clustering
so.merged <- RunPCA(so.merged, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:30, verbose = FALSE)

so.merged <- FindNeighbors(so.merged, dims = 1:30, verbose = FALSE)
so.merged <- FindClusters(so.merged, verbose = FALSE)
DimPlot(so.merged, label = TRUE)





for (s in 1:nrow(samples0.filtered)) {
	  nrow(so.list[[1]])
}





	  # Normalizing the data
	  so <- NormalizeData(so)
	  so <- GetAssayData(object=so, layer="counts")
	  writeTable(as.data.frame(so), gzfile(file.path(wd.de.data, "GRCh38-2020-A", paste0(samples0$V3[s], ".txt.gz"))), colnames=T, rownames=F, sep="\t")
	
	  # Scaling the data (required for PCA)
	  all.genes <- rownames(so)
	  so <- ScaleData(so, features=all.genes)
	
	 # Perform linear dimensional reduction (required for UMAP)
	 #ssc0.pca <- RunPCA(ssc0, features=VariableFeatures(object=ssc0))
	
	# Run non-linear dimensional reduction (UMAP/tSNE)
	#ssc0.umap <- RunUMAP(ssc0.pca, dims=1:10)
	
	#file.name <- file.path(wd.de.plots, "02_UMAP", paste0(samples0$V3[s], "_DimPlot"))
	#pdf(paste0(file.name, ".pdf"))
	#DimPlot(ssc0.umap, reduction="umap")
	#dev.off()
}
filtered <- raw