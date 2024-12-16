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
wd.de.data  <- file.path(wd.de, paste0("data_ALL4_500"))
wd.de.plots <- file.path(wd.de, paste0("plots_ALL4_500"))

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

load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_500",           "ssc_filtered_normalised.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500",      "ssc_filtered_normalised.1.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500",      "ssc_filtered_normalised.2.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_June_500", "ssc_filtered_normalised.3.RData"))

so.list <- c(so.list, so.list.1, so.list.2, so.list.3)
ids <- c(ids, ids.1, ids.2, ids.3)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.1)
samples0.filtered <- rbind(samples0.filtered, samples0.2)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.3)

#save(samples0.filtered, ids, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list.RData")))

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
library(DoubletFinder)

for (i in seq_along(so.list)) {
	  sample <- so.list[[i]]
	  raw_counts <- GetAssayData(sample, layer = "counts", assay = "RNA")
	  sample[["RNA"]] <- CreateAssayObject(counts = raw_counts)
	  DefaultAssay(sample) <- "RNA"
	  
	  # Preprocessing
	  sample <- NormalizeData(sample, assay = "RNA")
	  sample <- FindVariableFeatures(sample)
	  sample <- ScaleData(sample)
	  sample <- RunPCA(sample)
	
	  # Perform clustering to generate meaningful identities for modelHomotypic
	  sample <- FindNeighbors(sample, dims = 1:10)
	  sample <- FindClusters(sample, resolution = 0.5)  # Adjust resolution as needed
	  
	  # Find optimal pK
	  sweep_res <- paramSweep_v3(sample, PCs = 1:10, sct = FALSE)
	  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
	  # Save the plot
	  png(file.path(wd.de.data, paste0("pK_optimization_", samples0.filtered$V3[i], ".png")), width = 800, height = 600)
	  optimal_pK <- find.pK(sweep_stats)$pK
	  dev.off()
	  optimal_pK <- as.numeric(as.character(optimal_pK))
	  
	  # Estimate homotypic doublet proportion
	  homotypic_prop <- modelHomotypic(Idents(sample))
	  nExp <- round(0.075 * ncol(sample))  # Adjust expected doublet rate
	  nExp_adj <- round(nExp * (1 - homotypic_prop))
	
	  # Run DoubletFinder
	  sample <- doubletFinder_v3(
		    sample,
		    PCs = 1:10,
		    pN = 0.25,
		    pK = optimal_pK,
		    nExp = nExp_adj,
		    reuse.pANN = FALSE,
		    sct = FALSE
	  )
	  
	  # Extract pK corresponding to the maximum BCmetric
	  pK_stats <- find.pK(sweep_stats)
	  optimal_pK <- pK_stats$pK[which.max(pK_stats$BCmetric)]
	  classification_col <- paste0("DF.classifications_0.25_", optimal_pK, "_", nExp_adj)
	  Idents(sample) <- sample@meta.data[[classification_col]]
	  
	  so.list[[i]] <- sample
}

save(samples0.filtered, ids, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.RData")))

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.RData")))

# Create a mapping of sample names to age
age_mapping <- setNames(samples0.filtered$V4, rownames(samples0.filtered))

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
	  sample_name <- samples0.filtered$V3[i]
	
	  # Add the information to the table
	  gene_cell_table <- rbind(gene_cell_table, data.frame(
		    Sample = sample_name,
		    Age = as.numeric(age_mapping[sample_name]),
		    Genes = nrow(sample),  # Number of features
		    Cells = ncol(sample),  # Number of cells
		    Cells_DF = ncol(subset(sample, idents = "Singlet"))
	  ))
}

# Display the table
print(gene_cell_table)


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
so.list <- so.list[-c(2, 8)]
# Remove rows 2 and 8 from the data frame
samples0.filtered <- samples0.filtered[-c(2, 8), ]
# Remove elements 2 and 8 from the array
ids <- ids[-c(2, 8)]

so.list <- c(so.list, so.list.1, so.list.2, so.list.3)
ids <- c(ids, ids.1, ids.2, ids.3)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.1)
samples0.filtered <- rbind(samples0.filtered, samples0.2)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.3)

save(samples0.filtered, ids, features, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_-2-8.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_-2-8.RData")))

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
#res <- 0.6

## Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
#so.list <- lapply(X = so.list, FUN = SCTransform)

## As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
## Run the PrepSCTIntegration() function prior to identifying anchors
#features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
#so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)
#save(samples0.filtered, ids, features, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_SCT_", nfeatures, "_so.list.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_SCT_", nfeatures, "_so.list.RData")))












# Extract raw counts from so.list and merge them
# Step 1: Find the union of all gene names
all_genes <- Reduce(union, lapply(so.list, function(x) rownames(GetAssayData(x, layer = "counts"))))

# Step 2: Standardize the row names of all datasets
raw_counts_list <- lapply(so.list, function(x) {
	  counts <- GetAssayData(x, layer = "counts")
	  missing_genes <- setdiff(all_genes, rownames(counts))
	  if (length(missing_genes) > 0) {
		    zero_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(length(missing_genes), ncol(counts)))
		    rownames(zero_matrix) <- missing_genes
		    counts <- rbind(counts, zero_matrix)
	  }
	  counts[all_genes, ]  # Reorder rows
})

# Step 3: Combine all standardized matrices
combined_raw_counts <- do.call(cbind, raw_counts_list)

# Step 4: Assign raw counts to so.integrated
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- subset(so.integrated, cells = colnames(combined_raw_counts))
so.integrated[["RNA"]] <- CreateAssayObject(counts = combined_raw_counts)

# Step 5: Verify
slotNames(GetAssay(so.integrated, assay = "RNA"))

# From so.integrated
integrated_samples <- colnames(GetAssayData(so.integrated, layer = "counts"))
# From so.list
list_samples <- lapply(so.list, function(x) colnames(GetAssayData(x, layer = "counts")))
combined_list_samples <- unlist(list_samples)  # Flatten the list of samples from so.list

# Check if the sample orders are identical
identical(integrated_samples, combined_list_samples)

# -----------------------------------------------------------------------------
# 
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
library(DoubletFinder)

# Use clustering information to estimate the proportion of homotypic doublets.
homotypic_prop <- modelHomotypic(Idents(so.integrated))

# Step 1: Find the optimal pK value (can skip if known)
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- UpdateSeuratObject(so.integrated)

sweep_res <- paramSweep_v3(so.integrated, PCs = 1:10, sct = FALSE)  # Adjust PCs
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
optimal_pK <- find.pK(sweep_stats)





# Identify the clusters you want to keep
clusters_to_remove <- c(1)
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)

ssc0 <- as.data.frame(table(so.integrated@meta.data$sample.id))
ssc <- ssc0
write.table(ssc0, file=file.path(wd.de.data, "ssc_samples_-1_>0.txt"), row.names=T, col.names=T, quote=F, sep='\t')
save(ssc0, ssc, so.integrated, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_res=0.6_5000_-1_>0.RData"))

# Access the metadata
metadata <- so.integrated@meta.data
# Count the number of cells per sample
cell_counts <- metadata %>%
	  group_by(sample.id) %>%
	  summarise(cell_count = n()) %>%
	  ungroup()

# Filter samples with more than 100 cells
samples_to_keep <- cell_counts %>%
	  filter(cell_count > 100) %>%
	  pull(sample.id)

# Subset the Seurat object to keep only the samples with more than 100 cells
so.integrated <- subset(so.integrated, cells = rownames(metadata[metadata$sample.id %in% samples_to_keep, ]))

# Check the result
ssc <- as.data.frame(table(so.integrated@meta.data$sample.id))
write.table(ssc, file=file.path(wd.de.data, "ssc_samples_-1_>100.txt"), row.names=T, col.names=T, quote=F, sep='\t')
save(ssc0, ssc, so.integrated, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_res=0.6_5000_-1_>100.RData"))

