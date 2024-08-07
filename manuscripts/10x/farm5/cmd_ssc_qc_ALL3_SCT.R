#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
nfeatures <- as.numeric(args[1])

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
wd.de.data  <- file.path(wd.de, "data_ALL3")
wd.de.plots <- file.path(wd.de, "plots_ALL3")

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

load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data",      "ssc_filtered_normalised.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26", "ssc_filtered_normalised.1.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26", "ssc_filtered_normalised.2.RData"))

so.list <- c(so.list, so.list.1, so.list.2)
ids <- c(ids, ids.1, ids.2)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.1)
samples0.filtered <- rbind(samples0.filtered, samples0.2)

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
## Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
so.list <- lapply(X = so.list, FUN = SCTransform)

## As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
## Run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

## When running FindIntegrationAnchors(), and IntegrateData(), set the normalization.method parameter to the value SCT.
## When running sctransform-based workflows, including integration, do not run the ScaleData() function
anchors <- FindIntegrationAnchors(object.list = so.list, normalization.method = "SCT", anchor.features = features, dims = 1:30)
so.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30, k.weight = 90)

save(samples0.filtered, so.list, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_SCT_", nfeatures, ".RData")))

# -----------------------------------------------------------------------------
# Perform CCA integration (Seurat 4.3.0; Running CCA)
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
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

so.integrated@meta.data$sample.id <- ids
so.integrated@meta.data$age <- ages
so.integrated@meta.data$age <- factor(so.integrated@meta.data$age, levels = c("25","27","37","40","48","57","60","71"))
so.integrated@meta.data$n2 <- n2s
so.integrated@meta.data$batch <- batches
head(so.integrated@meta.data)

save(samples0.filtered, so.integrated, ids, ages, n2s, batches, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_", nfeatures, ".RData")))

# -----------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
# 02_UMAP
# https://satijalab.org/seurat/articles/multimodal_vignette
# -----------------------------------------------------------------------------
so.integrated <- RunPCA(so.integrated, verbose = F)

pdf(file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_ElbowPlot_SCT_", nfeatures, ".pdf")))
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.integrated, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.integrated[["pca"]]@stdev / sum(so.integrated[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_3000.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(so.integrated, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)

so.integrated <- FindNeighbors(so.integrated, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, ".RData")))

# Find neighbors and clusters
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))

so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = 0.3)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, ".RData")))

##
file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.3")
	
pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
DimPlot(so.integrated, label = TRUE)
dev.off()

pdf(file=file.path(wd.de.plots, paste0(file.name, "_SampleID.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, paste0(file.name, "_Age.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, paste0(file.name, "_2N.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="n2")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots,  paste0(file.name, "_Batch.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="batch")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
DefaultAssay(so.integrated) <- "RNA"

so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()
#missing_genes <- setdiff(top10$gene, rownames(so.integrated))

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3.png")),
							plot = heatmap, width = 10, height = 12, dpi = 300)

# -----------------------------------------------------------------------------
# Jaccard index for overlap of gene sets
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_markers.RData")))

jaccard_matrix <- getJaccardIndex(markers)
# Visualize the Jaccard index matrix
pheatmap_plot <- pheatmap(jaccard_matrix, main = "Jaccard Index Between Clusters")

# Save the plot as a PNG file
png(file = file.path(wd.de.plots, "heatmap_jaccard_markers_SCT_res=0.3.png"), width = 7, height = 7, units = "in", res = 300)
print(pheatmap_plot)
dev.off()
