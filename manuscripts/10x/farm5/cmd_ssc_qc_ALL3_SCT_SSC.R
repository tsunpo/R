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
wd.de.data  <- file.path(wd.de, "data_ALL3_500_SSC")
wd.de.plots <- file.path(wd.de, "plots_ALL3_500_SSC")

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

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_5000_SSC_>100.RData"))

## As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
## Run the PrepSCTIntegration() function prior to identifying anchors
so.list <- SplitObject(so.integrated, split.by = "orig.ident")

## Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
so.list <- lapply(X = so.list, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

## When running FindIntegrationAnchors(), and IntegrateData(), set the normalization.method parameter to the value SCT.
## When running sctransform-based workflows, including integration, do not run the ScaleData() function
anchors <- FindIntegrationAnchors(object.list = so.list, normalization.method = "SCT", anchor.features = features, dims = 1:30)
# Error in FindIntegrationAnchors(object.list = so.list, normalization.method = "SCT",  : 
# 	Max dimension too large: objects 8 contain fewer than 20 cells. 
#		Please specify a maximum dimensions that is less than the number of cells in any object (19).

# Warning in irlba(A = mat3, nv = num.cc) :
# 	You're computing too large a percentage of total singular values, use a standard svd instead.
# Merging objects
# Finding neighborhoods
# Error in idx[i, ] <- res[[i]][[1]] : 
#   number of items to replace is not a multiple of replacement length
rm(so.integrated)
so.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30, k.weight = 90)

save(so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_", nfeatures, ".RData")))

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
write.table(prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".txt")),row.names=F, col.names=F,quote=F,sep='\t')
save(so.integrated, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
