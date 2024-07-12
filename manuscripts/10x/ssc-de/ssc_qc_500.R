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
wd.de.data  <- file.path(wd.de, "data_500")
wd.de.plots <- file.path(wd.de, "plots_500")

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

# -----------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
# 02_UMAP
# https://satijalab.org/seurat/articles/multimodal_vignette
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))

# Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA
DefaultAssay(so.merged) <- "RNA"

# perform visualization and clustering steps
so.merged <- NormalizeData(so.merged)
so.merged <- FindVariableFeatures(so.merged, selection.method = "vst", nfeatures = 5000)
so.merged
# An object of class Seurat 
# 34615 features across 57031 samples within 1 assay 
# Active assay: RNA (34615 features, 5000 variable features)
# 26 layers present: counts.PD53621b_2N, counts.PD53623b_2N, counts.PD53623b_4N, counts.PD53624b_2N, counts.PD53625b_2N, counts.PD53626b_2N, counts.PD53621b_M, counts.PD53623b_M, counts.PD53624b_M, counts.PD53625b_M, counts.PD53626b_M, counts.PD40746e_M1, counts.PD40746e_M2, data.PD53621b_2N, data.PD53623b_2N, data.PD53623b_4N, data.PD53624b_2N, data.PD53625b_2N, data.PD53626b_2N, data.PD53621b_M, data.PD53623b_M, data.PD53624b_M, data.PD53625b_M, data.PD53626b_M, data.PD40746e_M1, data.PD40746e_M2

all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = all.genes)
#so.merged <- ScaleData(so.merged)
# Centering and scaling data matrix
# |======================================================================| 100%
so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged))
#so.merged <- RunPCA(so.merged, verbose = FALSE)

pdf(file.path(wd.de.data, "ssc_filtered_normalised_merged_ElbowPlot.pdf"))
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.merged, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.merged[["pca"]]@stdev / sum(so.merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(filtered, normalised, samples0, samples0.filtered, so.merged, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)

so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
save(filtered, normalised, samples0, samples0.filtered, so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP.RData"))

# Find neighbors and clusters
load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))

so.merged <- FindNeighbors(so.merged, dims = 1:prin_comp, k.param = 20)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = 0.3)

# Optional: Run UMAP for visualization
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20)
save(filtered, normalised, samples0, samples0.filtered, so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_resolution=0.3.RData"))

##
pdf(file=file.path(wd.de.plots, "UMAP_dims=18_resolution=0.3.pdf"))
DimPlot(so.merged, label = TRUE)
dev.off()

pdf(file=file.path(wd.de.plots, "UMAP_dims=18_resolution=0.3_SampleID.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, "UMAP_dims=18_resolution=0.3_Age.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, "UMAP_dims=18_resolution=0.3_2N.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="n2")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, "UMAP_dims=18_resolution=0.3_Batch.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="batch")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

##
feature_plot <- FeaturePlot(so.merged, reduction = "umap", features = c("VIM", "CD14", "CD163", "C1QA", "CXCR4", "VWF", "PECAM1", "NOTCH4", "ACTA2", "DLK1"),	ncol = 5)
pdf(file = file.path(wd.de.plots, "UMAP_dims=18_resolution=0.25_feature_plot_Marcrophase_Endothelial.pdf"), width = 14, height = 5)
print(feature_plot)
dev.off()

feature_plot <- FeaturePlot(so.merged, reduction = "umap", features = c("DAZL", "MAGEA4", "UTF1", "ID4", "FGFR3", "KIT", "DMRT1", "DMRTB1", "STRA8"),	ncol = 5)
pdf(file = file.path(wd.de.plots, "UMAP_dims=18_resolution=0.25_feature_plot_SSC_differentiating.pdf"), width = 14, height = 5)
print(feature_plot)
dev.off()

feature_plot <- FeaturePlot(so.merged, reduction = "umap", features = c("SYCP3", "SPO11", "MLH3", "SPAG6", "CAMK4", "ZPBP", "CREM", "TNP1", "PRM2"),	ncol = 5)
pdf(file = file.path(wd.de.plots, "UMAP_dims=18_resolution=0.25_feature_plot_Meiosis_Spermatid structure proteins_Nuclear condensation.pdf"), width = 14, height = 5)
print(feature_plot)
dev.off()

##
feature_plot <- FeaturePlot(so.merged, reduction = "umap", features = c("C19orf84", "EGR4", "MAGEA4", "PIWIL4", "TSPAN33", "UTF1", "FGFR3", "NANOS2"),	ncol = 5)
pdf(file = file.path(wd.de.plots, "UMAP_dims=18_resolution=0.25_feature_plot_ST4_SSC_State0.pdf"), width = 14, height = 5)
print(feature_plot)
dev.off()

feature_plot <- FeaturePlot(so.merged, reduction = "umap", features = c("GFRA1", "NANOS3", "DMRT1", "DNMT1", "KIT", "MKI67", "SOHLH2", "MAGE4", "REC8", "STRA8"),	ncol = 5)
pdf(file = file.path(wd.de.plots, "UMAP_dims=18_resolution=0.25_feature_plot_ST4_SSC_Stage1+2+3.pdf"), width = 14, height = 5)
print(feature_plot)
dev.off()

# -----------------------------------------------------------------------------
# DotPlot
# -----------------------------------------------------------------------------
genes_of_interest <- c("AMH", "WT1", "CD34", "CD163", "CSF1R", "CYP1B1", "ACTA2", "MYH11", "DLK1", "INHBA", "UTF1", "EGR4", "PIWIL4", "TSPAN33", "FGFR3", "NANOS2", "NANOS3", "GFRA1", "DMRT1", "MAGEA4", "KIT", "MKI67", "DPEP3", "GINS2", "MEIOB", "SCML1", "SYCP2", "SYCP3", "TEX101", "SPO11", "STRA8", "SPATA8", "OVOL2", "CLDND2", "FAM24A", "SPACA1", "CCDC168", "SIRT2", "TEX29", "HOOK1", "PRM1", "PRM2")

dot_plot <- DotPlot(so.merged, features = genes_of_interest)  +
	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, "DotPlot_dims=18_resolution=0.3.pdf"), width = 12, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("17", "19", "8", "16", "10", "15", "7", "4", "2", "18", "9", "12", "13", "5", "6", "11", "3", "0", "14", "1")
Idents(so.merged) <- factor(Idents(so.merged), levels = new_order)

dot_plot <- DotPlot(so.merged, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
   scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, "DotPlot_dims=18_resolution=0.3_ordered.pdf"), width = 12, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.merged))
cluster_to_celltype <- c('1' = 'Undiff. SPG', '14' = 'Undiff. SPG', '0' = 'Undiff. SPG', 
																									'3' = 'Diff. SPG',
																									'11' = 'Meiosis', '6' = 'Meiosis',
																									'5' = 'Spermatocyte', '13' = 'Spermatocyte', '12' = 'Spermatocyte', 
																									'9' = 'Early spermatid (1)',
																									'4' = 'Early spermatid (2)',
																									'7' = 'Late spermatid',
																									'10' = 'Myoid & Leydig',
																									'16' = 'Sertoli',	'8' = 'Sertoli',
																									'19' = 'Endothelial',
																									'17' = 'Macrophage')

# Update the identities using this mapping
new_order <- c("17", "19", "8", "16", "10", "15", "2", "18", "7", "4", "9", "12", "13", "5", "6", "11", "3", "0", "14", "1")
Idents(so.merged) <- factor(Idents(so.merged), levels = new_order)

Idents(so.merged) <- plyr::mapvalues(x = Idents(so.merged), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.merged, features = genes_of_interest)  +
	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, "DotPlot_dims=18_resolution=0.3_ordered_annotated.pdf"), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.merged, label = TRUE)
ggsave(file.path(wd.de.plots, "UMAP_dims=18_resolution=0.3_ordered_annotated.png"), plot = dim_plot, width = 10, height = 8, dpi = 300)









# -----------------------------------------------------------------------------
# To Andy
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "ssc_filtered.RData"))
samples0   ## From file scRNA_GRCh38-2020.list

so.list <- c()
for (s in 1:nrow(samples0)) {
	  # Initialize the Seurat object with the raw (non-normalized data)
	  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
	  data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0$V1[s], "filtered_feature_bc_matrix"))
	  so <- CreateSeuratObject(counts=data, project=samples0$V3[s], min.cells=3, min.features=200)
	
	  # QC and selecting cells for further analysis
	  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	  so.list <- c(so.list, so)
}
save(samples0, so.list, file=file.path(wd.de.data, "ty2_so.list_To-Andy.RData"))

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

so.merged@meta.data$sample.id <- ids
so.merged@meta.data$age <- ages
so.merged@meta.data$age <- factor(so.merged@meta.data$age, levels = c("25","37","48","57","60","71"))
so.merged@meta.data$n2 <- n2s
head(so.merged@meta.data)

save(so.merged, ids, ages, n2s, file=file.path(wd.de.data, "ssc_so.merged_To-Andy.RData"))









# Find markers for all clusters
so.merged <- JoinLayers(so.merged, features = "RNA")
#markers <- FindAllMarkers(so.merged)
markers <- FindAllMarkers(so.merged, assay="RNA", slot="data", test.use = "wilcox", only.pos = TRUE,	min.pct = 0.25, logfc.threshold = 0.25)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
																															"CD8A"))


options(repr.plot.width=16, repr.plot.height=5)

feature_plot <- FeaturePlot(so.merged, reduction = "umap", features = c("FGFR3", "KIT", "DMRT1", "SYCP3", "SPO11", "MLH3", "SPAG6", "ZPBP", "CAMK4", "CREM", "TNP1", "PRM2"),	ncol = 3)
pdf(file = file.path(wd.de.plots, "PCA_UMAP_dims=10_resolution=0.25_feature_plot.pdf"), width = 10, height = 8)
print(feature_plot)
dev.off()

# -----------------------------------------------------------------------------
# 5.integrate_all_whole_testes_datasets.R
# 03_PCA
# -----------------------------------------------------------------------------
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(sctransform)
library(sceasy)
library(reticulate)
library(clustree)

load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_JoinLayers_nfeatures=2000_PCA_UMAP_k=20_n=20.RData"))
load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_JoinLayers_nfeatures=2000_PCA_UMAP_k=20_n=20_markers.RData"))

so.merged@meta.data$RNA_snn_res.0 <- NULL
#so.merged@meta.data$RNA_snn_res.0.2<-NULL   ## In our data, only these three have values
#so.merged@meta.data$RNA_snn_res.0.4<-NULL   ## In our data, only these three have values
#so.merged@meta.data$RNA_snn_res.0.6<-NULL
#so.merged@meta.data$RNA_snn_res.0.8<-NULL
#so.merged@meta.data$RNA_snn_res.1<-NULL
#so.merged@meta.data$RNA_snn_res.1.2<-NULL
#so.merged@meta.data$RNA_snn_res.1.4<-NULL
#so.merged@meta.data$RNA_snn_res.1.6<-NULL
#so.merged@meta.data$RNA_snn_res.1.8<-NULL
#so.merged@meta.data$RNA_snn_res.2<-NULL
#so.merged@meta.data$orig.ident<-NULL       ## In our data, only these three have values
#so.merged@meta.data$old.ident<-NULL
#saveRDS(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20.rds"))

# convert the Seurat object to an h5ad object for visualisation with cellxgene
DefaultAssay(so.merged) <- 'RNA'
so.merged <- JoinLayers(so.merged)               # https://www.biostars.org/p/9581468/

sceasy::convertFormat(so.merged, from="seurat", to="anndata", outFile=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20.h5ad"))

all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = all.genes)

# https://github.com/satijalab/seurat/issues/8304
sceasy::convertFormat(so.merged, assay="RNA", from="seurat", to="anndata", outFile=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20.h5ad"))

# determine cluster markers
# we will later run enrichment analyses on these gene lists to assign a predicted function to each cluster
# IMPORTANT: we must set the DefaultAssay to RNA before running differential expression: 
# "we don't recommend using the integrated matrix for differential expression" and "As a general rule, we always recommend performing DE on originally measured values - instead of on batch-corrected, imputed, etc. values.
# This ensures that the measurements that enter the DE test are indeed independent from each other, which is a requirement of any statistical DE test."
# https://github.com/satijalab/seurat/issues/1057
# https://github.com/satijalab/seurat/issues/1256
# https://github.com/satijalab/seurat/issues/2136
DefaultAssay(so.merged) <- 'RNA'
so.merged <- JoinLayers(so.merged)   # https://www.biostars.org/p/9581468/

markers <- FindAllMarkers(so.merged, test.use="MAST", only.pos=T, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_JoinLayers_nfeatures=2000_PCA_UMAP_k=20_n=20_markers.txt"), row.names=F, col.names=T, quote=F, sep='\t')
save(markers, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_JoinLayers_nfeatures=2000_PCA_UMAP_k=20_n=20_markers.RData"))

top.markers <- subset(markers, cluster == "0")
lfcs <- getMarkerEffects(top.markers)

###
##
DefaultAssay(so.merged) <- 'RNA'
so.merged <- JoinLayers(so.merged)   # https://www.biostars.org/p/9581468/

so.merged <- FindVariableFeatures(so.merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = all.genes)

markers.0 <- FindMarkers(so.merged, ident.1=0)
head(markers.0, n=5)

markers <- FindAllMarkers(so.merged, only.pos=T)
write.table(markers, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20_markers_.txt"), row.names=F, col.names=T, quote=F, sep='\t')
save(markers, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20_markers_.RData"))

markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

VlnPlot(pbmc, features = c("PRM1", "PRM2"))



# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)





##
pdf("DimPlot_UMAP_RNA_dim=18_k=10.pdf")
DimPlot(so.merged, label = TRUE)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=18_k=10_sampleID.pdf")
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=18_k=10_age.pdf")
tplot = DimPlot(so.merged, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=18_k=10_2n.pdf")
tplot = DimPlot(so.merged, reduction = "umap", group.by="n2")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()




























# -----------------------------------------------------------------------------
# Standard Seurat pre-processing workflow (SCT)
# 01_QC
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "ssc_filtered.RData"))
library(glmGamPoi)

samples0.filtered <- samples0[subset(filtered, cells > 1000)$PD_ID,]
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
	so <- subset(so, subset=nFeature_RNA > 1000 & nFeature_RNA < 10000 & nCount_RNA > 2000 & nCount_RNA < 50000 & percent.mt < 5)
	
	# Apply sctransform normalization
	# https://satijalab.org/seurat/articles/sctransform_vignette.html
	so <- SCTransform(so, vars.to.regress="percent.mt", verbose=F)
	
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
save(filtered, normalised, samples0, samples0.filtered, so.merged, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))


so.merged <- FindNeighbors(so.merged, dims = 1:30)
# Computing nearest neighbor graph
# Computing SNN

so.merged <- FindClusters(so.merged, resolution = 0.8, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:30)
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
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_RNA_PCA_UMAP.RData"))

#pdf('DimPlot_UMAP_RNA_dim=30.pdf')
#DimPlot(so.merged, label = TRUE)
#dev.off()

pdf('DimPlot_UMAP_RNA_dim=30_.pdf')
tplot = DimPlot(so.merged, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('DimPlot_UMAP_RNA_dim=30_grouped_by_sampleID.pdf')
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('DimPlot_UMAP_RNA_dim=30_grouped_by_age.pdf')
tplot = DimPlot(so.merged, reduction = "umap", group.by="age")
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