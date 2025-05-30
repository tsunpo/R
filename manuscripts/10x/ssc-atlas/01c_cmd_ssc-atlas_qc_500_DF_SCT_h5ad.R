# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 09/05/25
# =============================================================================

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/10x")
wd.rna.raw <- file.path(wd.rna, "atlas")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-atlas"))
wd.de.data  <- file.path(wd.de, "data_500")
wd.de.plots <- file.path(wd.de, "plots_500")

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_DF_SCT_", nfeatures, "_so.list.RData")))

# -----------------------------------------------------------------------------
# Convert SCT-normalized, unintegrated to .h5ad
# -----------------------------------------------------------------------------
# Merge all objects (add cell IDs to avoid collisions)
so.merged <- merge(
   x = so.list[[1]],
   y = so.list[-1],
   add.cell.ids = names(so.list),
   project = "SSC"
)

# (Optional) Clean metadata: remove any list-columns that SeuratDisk can’t export
so.merged@meta.data <- so.merged@meta.data %>%
   dplyr::select(where(~ !is.list(.)))

# Step 4: Save as .h5Seurat, then convert to .h5ad
filename <- file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_DF_SCT_", nfeatures, "_so.merged.h5Seurat"))

SaveH5Seurat(so.merged, filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE)
