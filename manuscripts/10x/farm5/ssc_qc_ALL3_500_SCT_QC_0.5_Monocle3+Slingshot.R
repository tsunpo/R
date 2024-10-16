# =============================================================================
# Manuscript   :
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 11/09/24
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
#wd <- "/Users/ty2/Work/sanger/ty2"             ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/scRNA")
wd.rna.raw <- file.path(wd.rna, "10x")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-de"))
wd.de.data  <- file.path(wd.de, "data_ALL3_500_QC")
wd.de.plots <- file.path(wd.de, "plots_ALL3_500_QC")

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(pheatmap)

nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
DefaultAssay(so.integrated) <- "RNA"
# > dim(so.integrated)
# [1] 36601 70738

# -----------------------------------------------------------------------------
# Monocle 3
# -----------------------------------------------------------------------------
# Cluster cells using UMAP embeddings
cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = F)
#cds <- order_cells(cds)

# Identify root cells
root_cells <- WhichCells(so.integrated, idents = "Stage 0")
cds <- order_cells(cds, root_cells = root_cells)

# Plot pseudotime
monocle3 <- plot_cells(cds, color_cells_by = "pseudotime",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=0)
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_pseudotime.png")), plot = monocle3, dpi = 300)

# -----------------------------------------------------------------------------
# Slingshot
# -----------------------------------------------------------------------------
library(SingleCellExperiment)
library(slingshot)

DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Extract data for SCE conversion
normalized_matrix <- GetAssayData(so.integrated[["RNA"]], layer = "data")
umap_embeddings <- Embeddings(so.integrated, reduction = "umap")  # UMAP embeddings
cell_metadata <- so.integrated@meta.data  # Metadata

# Create a SCE object
sce <- SingleCellExperiment(
  	assays = list(data = normalized_matrix),
	  colData = cell_metadata,
	  reducedDims = list(UMAP = umap_embeddings)
)
sce$seurat_clusters <- Idents(so.integrated)

# Running Slingshot using UMAP embeddings and cluster information from Seurat
sce <- slingshot(sce, clusterLabels = sce$seurat_clusters, reducedDim = 'UMAP', start.clus = "Stage 0", end.clus = "Late spermatid")

# View pseudotime values
pseudotime_values <- sce$slingPseudotime_1
head(pseudotime_values)

# Convert UMAP embeddings to a data frame for ggplot
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
umap_df$pseudotime <- pseudotime_values

# Plot the UMAP embeddings with pseudotime overlay
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
	  geom_point() +
	  scale_color_viridis_c() +
	  labs(title = "Slingshot Pseudotime Trajectory on UMAP", x = "UMAP 1", y = "UMAP 2") +
	  theme_minimal()
