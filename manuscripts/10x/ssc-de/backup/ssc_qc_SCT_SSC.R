#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
nfeatures <- 5000

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
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

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
#load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500", "ssc_filtered_normalised.1.RData"))
#load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500", "ssc_filtered_normalised.2.RData"))

#so.list <- c(so.list, so.list.1, so.list.2)
#ids <- c(ids, ids.1, ids.2)
#samples0.filtered <- rbind(samples0.filtered, samples0.filtered.1)
#samples0.filtered <- rbind(samples0.filtered, samples0.2)

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
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_UMAP_", nfeatures, ".RData")))

# Find neighbors and clusters
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))

so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = 0.4)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, ".RData")))

##
file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.4")
	
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
# Di Persio et al; DotPlot (resolution = 0.4)
# -----------------------------------------------------------------------------
nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, ".RData")))

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("6", "16", "14", "18", "17", "3", "13", "9", "12", "1", "0", "10", "7", "11", "8", "2", "15", "5", "4")
new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4_ordered_integrated.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
cluster_to_celltype <- c('5' = 'Undiff. SPG', '15' = 'Undiff. SPG', '2' = 'Undiff. SPG',
																									'8' = 'Diff. SPG', '11' = 'Diff. SPG',
																									'7' = 'Leptotene',
																									'10' = 'Zygotene', '0' = 'Zygotene',
																									'1' = 'Pachytene', '12' = 'Pachytene',
																									'9' = 'Diplotene',
																									'13' = 'Meiotic division',
																									'3' = 'Early spermatid',
																									'17' = 'Late spermatid',
																									'18' = 'Macrophage',
																									'14' = 'Endothelial',
																									'6' = 'Sertoli', '16' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4_ordered_annotated_integrated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.4_ordered_annotated_integrated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.4; Six SPG states)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, ".RData")))

genes_of_interest <- c("TAF6", "ST3GAL4", "SH2B2", "MSL3", "PHGDH", "C19orf84", "LIN7B", "FSD1", "TSPAN33", "EGR4", "PIWIL4", "CELF4", "UTF1", "FGFR3", "A2M", "ENO3", "SERPINE2", "SRRT", "BAG6", "DND1", "PELP1", "NANOS2", "C1QBP", "NANOS3", "GFRA2", "GFRA1", "ID2", "ASB9", "L1TD1", "ID4", "MKI67", "PDPN", "KIT", "DMRT1", "DNMT1", "CALR", "SYCP3", "STRA8")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("6", "16", "14", "18", "17", "3", "13", "9", "12", "1", "0", "10", "7", "11", "8", "2", "15", "5", "4")
new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
write.table(table(Idents(so.integrated)), file=file.path(wd.de.data, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4_ordered_annotated.txt")), row.names=T, col.names=T, quote=F, sep='\t')

cluster_to_celltype <- c('4' = 'Stage 0',
	                        '5' = 'Stage 0A', '15' = 'Stage 0B', '2' = 'Stage 1',
																									'8' = 'Stage 2', '11' = 'Stage 3',
																									'7' = 'Leptotene',
																									'10' = 'Zygotene', '0' = 'Zygotene',
																									'1' = 'Pachytene', '12' = 'Pachytene',
																									'9' = 'Diplotene',
																									'13' = 'Meiotic division',
																									'3' = 'Early spermatid',
																									'17' = 'Late spermatid',
																									'18' = 'Macrophage',
																									'14' = 'Endothelial',
																									'6' = 'Sertoli', '16' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.4_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SPG_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.4_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# DefaultAssay(so.integrated) <- "SCT"
# Error in UseMethod(generic = "JoinLayers", object = object) : 
#  	no applicable method for 'JoinLayers' applied to an object of class "c('SCTAssay', 'Assay', 'KeyMixin')"
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, ".RData")))
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

save(so.integrated, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_markers.RData")))

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_markers.RData")))

top10 <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()
#missing_genes <- setdiff(top10$gene, rownames(so.integrated))

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
# Centering and scaling data matrix
# Warning: Different features in new layer data than already exists for scale.data
# Warning: Different cells in new layer data than already exists for scale.data

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.4.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# Set the new order
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_markers.RData")))

new_order <- c("6", "16", "14", "18", "17", "3", "13", "9", "12", "1", "0", "10", "7", "11", "8", "2", "15", "5", "4")
new_cluster_order <- rev(new_order)
so.integrated <- SetIdent(so.integrated, value = factor(Idents(so.integrated), levels = new_cluster_order))

top10 <- markers %>%
	  mutate(cluster = factor(cluster, levels = new_cluster_order)) %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.4_new_order_cluster.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and cluster identities
read_counts <- FetchData(so.integrated, vars = c("nCount_RNA", "seurat_clusters"))

# Rename columns for clarity
colnames(read_counts) <- c("read_count", "cluster")

# Convert cluster to factor for ordered plotting
read_counts$cluster <- factor(read_counts$cluster)

# Calculate the median read count for each cluster
medians_read <- read_counts %>%
	  group_by(cluster) %>%
	  summarize(median_read_count = median(read_count))

# Reorder the cluster factor levels based on median read counts
read_counts$cluster <- factor(read_counts$cluster, levels = medians_read$cluster[order(medians_read$median_read_count)])

# Generate the box plot for read counts
box_plot_reads <- ggplot(read_counts, aes(x = cluster, y = read_count)) +
	  geom_boxplot() +
	  theme_minimal() +
	  labs(title = "Read Counts",
			  			x = "Cluster",
				  		y = "Number of Reads") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1),
	  						plot.title = element_text(hjust = 0.5))  # Center the title

# Part 2: Generate the Box Plot for Global Expression Levels
# Extract the total expression level (sum of counts) and cluster identities
expression_data <- FetchData(so.integrated, vars = c("nFeature_RNA", "seurat_clusters"))

# Rename columns for clarity
colnames(expression_data) <- c("total_expression_level", "cluster")

# Convert cluster to factor for ordered plotting
expression_data$cluster <- factor(expression_data$cluster)

# Calculate the median total expression level for each cluster
medians_expression <- expression_data %>%
	  group_by(cluster) %>%
	  summarize(median_total_expression_level = median(total_expression_level))

# Reorder the cluster factor levels based on median total expression levels
expression_data$cluster <- factor(expression_data$cluster, levels = medians_expression$cluster[order(medians_expression$median_total_expression_level)])

# Generate the box plot for global expression levels
box_plot_expression <- ggplot(expression_data, aes(x = cluster, y = total_expression_level)) +
	  geom_boxplot() +
	  theme_minimal() +
	  labs(title = "Global Expression Levels",
			  			x = "Cluster",
				  		y = "Total Expression Level") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1),
	  						plot.title = element_text(hjust = 0.5))  # Center the title

# Part 3: Combine the Plots
combined_plot <- box_plot_reads / box_plot_expression

# Print the combined plot
print(combined_plot)

# -----------------------------------------------------------------------------
# QC (Sampke ID)
# -----------------------------------------------------------------------------
# Extract sample IDs and cluster identities
sample_cluster_data <- FetchData(so.integrated, vars = c("sample.id", "seurat_clusters"))

# Calculate the number of cells in each cluster
total_cells_per_cluster <- sample_cluster_data %>%
	  group_by(seurat_clusters) %>%
	  summarize(total_cells = n())

# Calculate the number of cells for each sample ID within each cluster
cells_per_sample_per_cluster <- sample_cluster_data %>%
	  group_by(seurat_clusters, sample.id) %>%
	  summarize(cells = n())

# Merge the data frames to calculate proportions
proportion_data <- merge(cells_per_sample_per_cluster, total_cells_per_cluster, by = "seurat_clusters")

# Calculate the proportion of each sample ID within each cluster
proportion_data <- proportion_data %>%
	  mutate(proportion = cells / total_cells) %>%
	  select(seurat_clusters, sample.id, proportion)

# Rename columns for clarity
colnames(proportion_data) <- c("Cluster", "Sample_ID", "Proportion")

# Calculate the proportion of PD40746e_M2 in each cluster
pd40746e_m2_proportion <- proportion_data %>%
	  filter(Sample_ID == "PD40746e_M2") %>%
	  arrange(desc(Proportion)) %>%
	  pull(Cluster)

# Set factor levels for Cluster based on the proportion of PD40746e_M2
proportion_data$Cluster <- factor(proportion_data$Cluster, levels = pd40746e_m2_proportion)

# Create a vector of unique sample IDs
unique_samples <- unique(proportion_data$Sample_ID)

# Create a vector of random colors
set.seed(42)  # For reproducibility
random_colors <- grDevices::colors()[sample(1:length(grDevices::colors()), length(unique_samples))]

# Create a named vector of colors, assigning specific colors to the highlighted samples
color_vector <- setNames(random_colors, unique_samples)
color_vector["PD40746e_M2"] <- "red"
color_vector["PD40746e_M1"] <- "pink"

# Assign colors to the data
proportion_data$color <- color_vector[proportion_data$Sample_ID]

# Create a stacked bar chart with custom colors
bar_chart <- ggplot(proportion_data, aes(x = as.factor(Cluster), y = Proportion, fill = color)) +
	  geom_bar(stat = "identity", color = "black") +  # Add border for better visibility
	  scale_fill_identity() +
	  theme_minimal() +
	  labs(title = "Proportion of Sample IDs in Each Cluster",
			  			x = "Cluster",
				  		y = "Proportion") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the bar chart
print(bar_chart)

# -----------------------------------------------------------------------------
# Jaccard index for overlap of gene sets
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_markers.RData")))

jaccard_matrix <- getJaccardIndex(markers)
# Visualize the Jaccard index matrix
pheatmap_plot <- pheatmap(jaccard_matrix, main = "Jaccard Index Between Clusters")

# Save the plot as a PNG file
png(file = file.path(wd.de.plots, "heatmap_jaccard_markers_SCT_res=0.4.png"), width = 7, height = 7, units = "in", res = 300)
print(pheatmap_plot)
dev.off()

# -----------------------------------------------------------------------------
# GSEA
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers.RData")))

# Load MSigDB gene sets (Hallmark gene sets as an example)
msigdb <- msigdbr(species = "Homo sapiens", category = "H")
msigdb_list <- split(msigdb$gene_symbol, msigdb$gs_name)

msigdb_df <- bind_rows(
	  lapply(names(msigdb_list), function(term) {
		    data.frame(term = term, gene = msigdb_list[[term]], stringsAsFactors = FALSE)
	  })
)

# Run GSEA for each cluster and store the results in a list
clusters <- unique(markers$cluster)
#clusters <- as.numeric(as.vector(clusters))
gsea_results <- lapply(clusters, run_gsea)
names(gsea_results) <- clusters

gsea_results <- gsea_results[!sapply(gsea_results, is.null)]
save(markers, gsea_results, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, "_markers_gsea_results_H.RData")))

# Apply to all GSEA results
significant_terms_list <- lapply(gsea_results, extract_significant_terms)
significant_terms_list <- significant_terms_list[!sapply(significant_terms_list, is.null)]

# Create a matrix of enrichment scores (NES) for significant terms
term_names <- unique(unlist(lapply(significant_terms_list, function(x) x$Description)))
cluster_names <- names(significant_terms_list)

nes_matrix <- matrix(NA, nrow = length(term_names), ncol = length(cluster_names))
rownames(nes_matrix) <- term_names
colnames(nes_matrix) <- cluster_names

for (i in seq_along(significant_terms_list)) {
	  terms <- significant_terms_list[[i]]$Description
	  nes_scores <- significant_terms_list[[i]]$NES
	  nes_matrix[terms, names(significant_terms_list)[i]] <- nes_scores
}

# Replace NA values in nes_matrix with a suitable value (e.g., 0 for non-significant)
nes_matrix[is.na(nes_matrix)] <- 0

# Generate heatmap of NES
ht_nes <- Heatmap(nes_matrix, 
																		name = "NES", 
																		col = colorRampPalette(c(blue, "white", red))(50),
																		show_row_names = TRUE, 
																		show_column_names = TRUE,
																		cluster_rows = TRUE, 
																		cluster_columns = TRUE,
																		heatmap_legend_param = list(title = "Normalized Enrichment Score"),
																		column_names_gp = gpar(fontsize = 12),  # Set font size for column names
																		column_names_rot = 45,  # Rotate column names by 45 degrees
                  column_names_side = "bottom")  # Place y-axis labels on the left side

# Save the NES heatmap to a file
filename <- file.path(wd.de.plots, paste0("GSEA_H_NES_SCT_", nfeatures, "_UMAP_resolution=0.4_p<0.01.pdf"))
pdf(file = filename, width = 10, height = 8)
draw(ht_nes, annotation_legend_side = "bottom",  # Move legend to the left
					padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
dev.off()

# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression data in each cluster
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.4_", nfeatures, ".RData")))

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Ensure age is numeric
so.integrated@meta.data$age <- as.numeric(so.integrated@meta.data$age)

# Extract age and cluster information from the metadata
age_cluster_data <- so.integrated@meta.data %>%
	  select(age, seurat_clusters)

# Calculate correlations for each cluster and store results
cluster_ids <- levels(so.integrated@active.ident)
correlation_results <- lapply(cluster_ids, calculate_spearman, age_cluster_data = age_cluster_data, so.integrated = so.integrated)

# Combine results from all clusters
all_significant_genes <- bind_rows(correlation_results, .id = "cluster")

# View the results
print(all_significant_genes)
# Optionally, save the results to a CSV file
write.table(all_significant_genes, file = "significant_genes_age_correlation_cluster.txt", row.names = FALSE)


