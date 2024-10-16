#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 25/07/24; 14/03/24
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
wd.de.data  <- file.path(wd.de, "data_ALL3_500_QC")
wd.de.plots <- file.path(wd.de, "plots_ALL3_500_QC")

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
library(pheatmap)

# -----------------------------------------------------------------------------
# To identify six SPG states
# -----------------------------------------------------------------------------
nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
# > dim(so.integrated)
# [1] 5000 70738

# Set the desired DefaultAssay in Seurat (e.g., "RNA")
DefaultAssay(so.integrated) <- "RNA"

# Cluster cells in Monocle 3 using UMAP embeddings
cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)
# Perform trajectory analysis using precomputed UMAP
cds <- cluster_cells(cds, reduction_method = "UMAP")
# Learn the principal graph
cds <- learn_graph(cds, use_partition = F)

# Order cells in pseudotime
#cds <- order_cells(cds)
# Identify root cells based on some known marker or cluster
root_cells <- WhichCells(so.integrated, idents = "Stage 0")  # Replace "early_cluster_name" with your specific cluster
# Ensure these root cells are present in `cds`
root_cells_in_cds <- root_cells[root_cells %in% colnames(cds)]
# Order cells using the root cells
cds <- order_cells(cds, root_cells = root_cells_in_cds)

# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "pseudotime",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=0)
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime.png")), plot = monocle3, dpi = 300)

# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#tscan
library(ggbeeswarm)
library(colorRamps)
library(viridisLite)

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$idents <- Idents(so.integrated)

#pseudotime_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
#pseudotime_palette <- inferno(100)
pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)

ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = pseudotime_monocle3)) +  # Use pseudotime as the color
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_gradientn(colors = pseudotime_palette) +  # Apply the pseudotime gradient
	  theme_classic() +
	  xlab("pseudotime") + 
	  ylab("") +
	  ggtitle("") +
	  theme(legend.position = "none")  # Equivalent to NoLegend()
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered.png")), plot = ordered, dpi = 300)
