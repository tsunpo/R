#!/usr/bin/env Rscript
args     <- commandArgs(TRUE)
cell_type <- as.character(args[1])

# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 25/07/24; 14/03/24
# =============================================================================

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
res <- 0.5
dims <- 20
k.weight <- 100

load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_clean.RData")))
so.integrated$patient.id <- sub("_.*", "", so.integrated@meta.data$sample.id)

# Assign updated metadata back
batch <- read.table(file.path(wd.de.data, paste0("Batch.txt")), header=T)
meta <- so.integrated@meta.data %>%
   mutate(cell_barcode = rownames(.)) %>%
   left_join(batch[, c("Sample", "Batch")], by = c("sample.id" = "Sample")) %>%
   column_to_rownames("cell_barcode")
so.integrated@meta.data <- meta

# Ensure age groups are defined in your metadata
so.integrated@meta.data <- so.integrated@meta.data %>%
   mutate(age_group = cut(
      age,
      breaks = c(-Inf, 50, Inf),
      labels = c("20-50", "50-80"),
      right = FALSE
   ))

cell_type <- "Stage 0"
so.integrated <- subset(so.integrated, idents = gsub("_", " ", cell_type))
wd.de.plots <- file.path(wd.de.plots, paste0(cell_type, "_RNA"))
dir.create(wd.de.plots, showWarnings = FALSE)

# Step 1: Switch to the uncorrected data
DefaultAssay(so.integrated) <- "RNA"
# Step 2: Normalize and scale (log-normalization)
so.integrated <- NormalizeData(so.integrated)
so.integrated <- FindVariableFeatures(so.integrated)
so.integrated <- ScaleData(so.integrated)

# -----------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
# 02_UMAP
# https://satijalab.org/seurat/articles/multimodal_vignette
# -----------------------------------------------------------------------------
so.integrated <- RunPCA(so.integrated, verbose = F)

pdf(file.path(wd.de.plots, paste0("ssc_filtered_normalised_integrated_ElbowPlot_SCT_", nfeatures, ".pdf")))
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
#write.table(prin_comp, file=file.path(wd.de.plots, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_", nfeatures, "_", dims, "_", k.weight, ".txt")), row.names=F, col.names=F, quote=F, sep='\t')
#save(so.integrated, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_", nfeatures, "_", dims, "_", k.weight, ".RData")))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
resolution.range <- seq(from = 0.05, to = 1, by = 0.05)

so.integrated <- FindNeighbors(so.integrated, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
#save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_", nfeatures, "_", dims, "_", k.weight, ".RData")))

# -----------------------------------------------------------------------------
# Define resolution
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
# Find neighbors and clusters
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_", nfeatures, "_", dims, "_", k.weight, ".RData")))

so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = res)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
#save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_", dims, "_", k.weight, ".RData")))

##
file.name <- paste0("SCT_", nfeatures, "_", 20, "_", 100, "_UMAP_dims=", prin_comp, "_res=0.5")

pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
DimPlot(so.integrated, label = TRUE) + NoLegend()
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

pdf(file=file.path(wd.de.plots, paste0(file.name, "_Age_Group.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="age_group")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots,  paste0(file.name, "_Batch.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="Batch")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots,  paste0(file.name, "_Phase.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="Phase")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots,  paste0(file.name, "_Patient.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="patient.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# -----------------------------------------------------------------------------
# Define resolution
# Last Modified: 19/05/25
# -----------------------------------------------------------------------------
library(GGally)
library(ggplot2)

# Extract PCA embeddings PC1–PC5
pcs <- Embeddings(so.integrated, reduction = "pca")[, 1:5]
pcs_df <- as.data.frame(pcs)
pcs_df$Batch <- so.integrated$Batch  # for coloring

# Rename columns for readability
colnames(pcs_df)[1:5] <- paste0("PC", 1:5)

# Plot only lower triangle of 5 PCs (4 x 4 grid of comparisons)
pdf(file = file.path(wd.de.plots, paste0(file.name, "_PC1to5_LowerTriangle_10x10", "_Batch", ".pdf")), width = 10, height = 10)

ggpairs(
   pcs_df,
   columns = 1:5,
   mapping = aes(color = Batch, alpha = 0.6),
   upper = list(continuous = wrap("blank")),       # hide upper triangle
   diag = list(continuous = wrap("blankDiag")),    # hide diagonal
   lower = list(continuous = wrap("points", size = 0.5))
   ) +
   theme_minimal() +
      theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "bottom"
      )

dev.off()