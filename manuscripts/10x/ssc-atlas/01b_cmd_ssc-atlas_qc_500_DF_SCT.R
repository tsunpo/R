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
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(pheatmap)

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.RData")))

# Initialize a data frame to store the results
gene_cell_table <- data.frame(
	  Sample = character(),
	  Age = numeric(),
	  Genes = numeric(),
	  Cells = numeric(),
	  Cells_DF = numeric()
)

# Loop through the Seurat objects in so.list
for (i in seq_along(so.list)) {
	  # Get the Seurat object
	  sample <- so.list[[i]]
	
	  # Extract sample name (assume rownames(samples0.filtered) match so.list identifiers)
	  sample_name <- samples.filtered$Sample[i]
	
	  # Add the information to the table
	  gene_cell_table <- rbind(gene_cell_table, data.frame(
		    Sample = sample_name,
		    Age = as.numeric(samples.filtered$Age[i]),
		    Genes = nrow(sample),  # Number of features
		    Cells = ncol(sample),  # Number of cells
		    Cells_DF = ncol(subset(sample, idents = "Singlet"))
	  ))
}

# Display the table
print(gene_cell_table)
write.table(gene_cell_table, file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.txt")), col.names=T, row.names=F, sep="\t")

for (i in seq_along(so.list)) {
	  sample <- so.list[[i]]
	  
	  # Subset to remove doublets
	  sample <- subset(sample, idents = "Singlet")

	  # Update the sample in so.list
	  so.list[[i]] <- sample
}

# -----------------------------------------------------------------------------
# Remove PD53622b_2N and PD53622b_M1
# -----------------------------------------------------------------------------
# Remove items 2 and 8 from the list (where there are only 130 and 92 cells in the samples)
#so.list <- so.list[-c(2, 8)]
# Remove rows 2 and 8 from the data frame
#samples0.filtered <- samples0.filtered[-c(2, 8), ]
# Remove elements 2 and 8 from the array
#ids <- ids[-c(2, 8)]

#save(samples0.filtered, ids, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF_-2-8.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF_-2-8.RData")))

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
#res <- 0.5

## Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
so.list <- lapply(X = so.list, FUN = SCTransform)

## As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
## Run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)
save(samples.filtered, ids, features, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_DF_SCT_", nfeatures, "_so.list.RData")))

# -----------------------------------------------------------------------------
# Convert SCT-normalized, unintegrated to .h5ad
# -----------------------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)

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

