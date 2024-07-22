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
wd.de.data  <- file.path(wd.de, "data_ALL3_500")
wd.de.plots <- file.path(wd.de, "plots_ALL3_500")

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

load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_500",      "ssc_filtered_normalised.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500", "ssc_filtered_normalised.1.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500", "ssc_filtered_normalised.2.RData"))

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
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_UMAP_", nfeatures, ".RData")))

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
# DotPlot (resolution = 0.25)
# -----------------------------------------------------------------------------
nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, ".RData")))

genes_of_interest <- c("AMH", "WT1", "CD34", "CD163", "CSF1R", "CYP1B1", "ACTA2", "MYH11", "DLK1", "INHBA", "UTF1", "EGR4", "PIWIL4", "TSPAN33", "FGFR3", "NANOS2", "NANOS3", "GFRA1", "DMRT1", "MAGEA4", "KIT", "MKI67", "DPEP3", "GINS2", "MEIOB", "SCML1", "SYCP2", "SYCP3", "TEX101", "SPO11", "STRA8", "MEIOSIN", "SPATA8", "OVOL2", "CLDND2", "FAM24A", "SPACA1", "CCDC168", "SIRT2", "TEX29", "HOOK1", "PRM1", "PRM2")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Bush_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("15", "16", "17", "1", "7", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Bush_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3_ordered.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
cluster_to_celltype <- c('0' = 'Undiff. SPG', '2' = 'Undiff. SPG',
																									'10' = 'Diff. SPG',
																									'9' = 'Meiosis', '11' = 'Meiosis', '5' = 'Meiosis',
																									'12' = 'Spermatocyte', '8' = 'Spermatocyte',
																									'13' = 'Early spermatid (1)', '6' = 'Early spermatid (1)',
																									'4' = 'Late spermatid (1)',
																									'7' = 'Myoid & Leydig',
																									'1' = 'Sertoli',	'17' = 'Sertoli',
																									'16' = 'Endothelial & Macrophage')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
  	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Bush_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Bush_SCT_", nfeatures, "_UMAP_dims=13_resolution=0.3_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.25)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, ".RData")))

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")
new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
  scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
cluster_to_celltype <- c('2' = 'Undiff. SPG',
																									'10' = 'Diff. SPG',
																									'9' = 'Leptotene',
																									'11' = 'Zygotene', '5' = 'Zygotene',
																									'12' = 'Pachytene', '8' = 'Pachytene',
																									'13' = 'Diplotene',
																									'6' = 'Meiotic division',
																									'4' = 'Late spermatid',
																									'16' = 'Endothelial & Macrophage',
																									'7' = 'PMC',
																									'15' = 'Leydig',
																									'1' = 'Sertoli',	'17' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SCT_", nfeatures, "_UMAP_dims=14_resolution=0.3_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_resolution=0.3.RData"))
length <- length(unique(so.integrated$seurat_clusters))

##
metadata <- so.integrated@meta.data
metadata$cluster <- so.integrated$seurat_clusters
batch_cluster_counts <- metadata %>%
	  group_by(cluster, batch) %>%
	  summarise(count = n()) %>%
	  ungroup()

facs <- toTable(0, 4, length, c("Cluster", "B1", "B2", "B3"))
facs$Cluster <- as.numeric(rownames(facs)) - 1
for (c in 0:length-1) {
	  subsets <- subset(batch_cluster_counts, cluster == c)
	  facs$B1[c+1] <- subset(subsets, batch == 1)$count + subset(subsets, batch == 2)$count
	  facs$B2[c+1] <- subset(subsets, batch == 3)$count
	  facs$B3[c+1] <- subset(subsets, batch == 4)$count
	
	  #total <- facs$B1[c+1] + facs$B2[c+1] + facs$B3[c+1]
	  #facs$B1[c+1] <- facs$B1[c+1] / total * 100
	  #facs$B2[c+1] <- facs$B2[c+1] / total * 100
	  #facs$B3[c+1] <- facs$B3[c+1] / total * 100
}
write.table(facs, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_counts.txt"), row.names=F, col.names=F, quote=F, sep='\t')

# -----------------------------------------------------------------------------
# 
# Last Modified: 10/10/23
# -----------------------------------------------------------------------------
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")

#facs <- readTable(file.path(wd.de.data, "batch_cluster_SCT_counts.txt"), header=T, rownames=T, sep="")
facs2 <- t(facs[,-1])
colnames(facs2) <- facs$Cluster
facs2 <- facs2[, new_order]
write.table(facs2, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_counts2.txt"), row.names=T, col.names=T, quote=F, sep='\t')

max <- 0
for (c in 1:length) {
	  if (sum(facs2[, c]) > max)
		    max <- sum(facs2[, c])
}

## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.de.plots, "batch_cluster_counts_SCT_0.3_new_order_numbers")
main.text <- c("", "")
xlab.text <- ""
ylab.text <- "Number of cells"
cols <- c(blue, red, yellow)
sequence <- seq(from = 0.8, by = 1.2, length.out = length)

pdf(paste0(file.name, ".pdf"), width = 12, height = 6)
par(mar=c(5.1, 4.6, 4.1, 9), xpd=TRUE)
barplot(facs2, col=cols, ylim=c(0, 10424), ylab=ylab.text, xaxt="n", main="", cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=new_order, x=sequence, y=par("usr")[3] - 360, srt=45, adj=0.965, xpd=NA, cex=1.9)

legend("right", c("B1", "B2", "B3"), text.col="black", pch=c(15, 15, 15), col=cols, pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
dev.off()

##
facs3 <- toTable(0, 4, length, c("Cluster", "B1", "B2", "B3"))
facs3$Cluster <- as.numeric(rownames(facs3)) - 1
for (c in 0:length-1) {
	  total <- facs$B1[c+1] + facs$B2[c+1] + facs$B3[c+1]
	  facs3$B1[c+1] <- facs$B1[c+1] / total * 100
	  facs3$B2[c+1] <- facs$B2[c+1] / total * 100
	  facs3$B3[c+1] <- facs$B3[c+1] / total * 100
}
facs3 <- t(facs3[,-1])
colnames(facs3) <- facs$Cluster
facs3 <- facs3[, new_order]
write.table(facs3, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_proportion.txt"), row.names=T, col.names=T, quote=F, sep='\t')

## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.de.plots, "batch_cluster_counts_SCT_0.3_new_order")
main.text <- c("", "")
xlab.text <- ""
ylab.text <- "Fraction of cells"
#blue  <- "blue"   ## adjustcolor("#619CFF", alpha.f=0.9)
#red   <- "red"   ## adjustcolor("#F8766D", alpha.f=0.9)
#green <- "darkgray"   ## adjustcolor("#00BA38", alpha.f=0.9)
#cols <- c(blue, red, green)   ## #59a523 (Alcro wasabi)
cols <- c(blue, red, yellow)
sequence <- seq(from = 0.8, by = 1.2, length.out = length)

pdf(paste0(file.name, ".pdf"), width = 12, height = 6)
par(mar=c(5.1, 4.6, 4.1, 9), xpd=TRUE)
barplot(facs3, col=cols, ylim=c(0, 100), ylab=ylab.text, xaxt="n", main="", cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=new_order, x=sequence, y=par("usr")[3] - 4, srt=45, adj=0.965, xpd=NA, cex=1.9)

legend("right", c("B1", "B2", "B3"), text.col="black", pch=c(15, 15, 15), col=cols, pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
#mtext(ylab.text, side=2, line=2.75, cex=1.8)
dev.off()

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# DefaultAssay(so.integrated) <- "SCT"
# Error in UseMethod(generic = "JoinLayers", object = object) : 
#  	no applicable method for 'JoinLayers' applied to an object of class "c('SCTAssay', 'Assay', 'KeyMixin')"
DefaultAssay(so.integrated) <- "RNA"

so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

save(so.integrated, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers.RData")))

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
# Centering and scaling data matrix
# Warning: Different features in new layer data than already exists for scale.data
# Warning: Different cells in new layer data than already exists for scale.data

# Create the heatmap
#heatmap <- DoHeatmap(so.integrated, features = top10$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend()
#ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3_group.by.png")),
#							plot = heatmap, width = 15, height = 15, dpi = 300)

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# Set the new order
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")
new_cluster_order <- rev(new_order)
so.integrated <- SetIdent(so.integrated, value = factor(Idents(so.integrated), levels = new_cluster_order))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3_new_order.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)








# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

markers.cluster0 <- FindMarkers(so.integrated, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


# -----------------------------------------------------------------------------
# Monocle
# -----------------------------------------------------------------------------
library(Monocle2)  # Ensure Monocle 3 is installed

# Checking for a valid number of clusters in Seurat object
DefaultAssay(so.integrated) <- "RNA"

num_clusters <- length(unique(Idents(so.integrated)))
assertthat::assert_that(num_clusters > 0)

# Extract necessary components from Seurat object
expression_data <- GetAssayData(so.integrated, slot = "counts")
cell_metadata <- so.integrated@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_data), row.names = rownames(expression_data))

# Create Monocle 3 cell_data_set
cds <- new_cell_data_set(expression_data,
																									cell_metadata = cell_metadata,
																									gene_metadata = gene_metadata)
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_Monocle3_1.RData")))

# Preprocess the Data in Monocle 3
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds)
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_Monocle3_2.RData")))

# Cluster the cells (this step might be optional depending on your use case)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_Monocle3.RData")))

# Plot cells colored by pseudotime
pseudotimeplot <- plot_cells(cds, color_cells_by = "pseudotime",
											label_groups_by_cluster = FALSE,
											label_leaves = FALSE,
											label_branch_points = FALSE)

ggsave(filename = file.path(wd.de.plots, paste0("Monocle3_pseudotime_SCT_", nfeatures, "_UMAP_resolution=0.25.png")),
							plot = pseudotimeplot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
library(AUCell)

# Assuming 'seurat_object' is your Seurat object
expression_matrix <- GetAssayData(so.integrated, assay = "SCT", layer = "data")


library(SingleR)

vignette("ExperimentHub", package = "ExperimentHub")
reference_dataset <- celldex::HumanPrimaryCellAtlasData()  # Example reference

annotations <- SingleR(test = so.integrated@assays$SCT@data, ref = reference_dataset, labels = reference_dataset$label.main)

