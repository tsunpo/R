# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 20/02/25; 25/07/24; 14/03/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

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
wd.de.data  <- file.path(wd.de, "data_ALL4_500_QC")
wd.de.plots <- file.path(wd.de, "plots_ALL4_500_QC")

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
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_clean.RData")))

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: 18/11/24
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)

so.integrated@meta.data$cell_type <- Idents(so.integrated)

# Ensure age groups are defined in your metadata
so.integrated@meta.data <- so.integrated@meta.data %>%
	  mutate(age_group = cut(
		    age,
		    breaks = c(-Inf, 50, Inf),
		    labels = c("20-50", "50-80"),
		    right = FALSE
	   ))

# Summarize counts for each cell type within each age group
cell_type_counts <- so.integrated@meta.data %>%
	  group_by(age_group, cell_type) %>%
	  summarize(
		    cell_count = n(),  # Count of cells for each cell type
		    .groups = "drop"
	  )

# Calculate total cells for each age group
age_group_totals <- cell_type_counts %>%
  	group_by(age_group) %>%
	  summarize(
		    total_cells = sum(cell_count),
		    .groups = "drop"
  	)

# Proportions sum to one, which transforms the data from the Euclidean space into a simplex space.
# Canonical statistical tests such as t-test may not be reliable, as these are not designed for relative proportion data
# But Spearman’s correlation measures the monotonic relationship between two variables by ranking the values rather than using their absolute differences. 
# This makes it more robust to non-normality and scale transformations
proportions <- cell_type_counts %>%
	  left_join(age_group_totals, by = "age_group") %>%
  	mutate(proportion = cell_count / total_cells)

# Define the median age for each age group
age_group_to_median <- c("20-50" = 35, "50-80" = 65)
# Add the `age_median` column
proportions$age_median <- as.numeric(sapply(proportions$age_group, function(x) age_group_to_median[x]))

# -----------------------------------------------------------------------------
# Embed genes from co-expression networks of human spermatogenesis (across young vs. old age groups and multiple developmental stages) 
# into a vector space, and analyze how each gene’s network context shifts with age.
#
# We will use Node2Vec (graph node embeddings via random walks + skip-gram) to generate gene embeddings for each condition, 
# then compare embeddings of the same gene between young and old using cosine similarity.
#
# Finally, we rank genes by embedding “drift” and visualize the results (UMAP plots and heatmaps) to highlight aging-associated transcriptome context shifts.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2. Construct gene co-expression networks (GNs) for each group
#
# First, ensure you have gene co-expression networks (GNs) for each combination of age group and cell type (developmental stage).
#
# Each GN should be represented as a graph with genes as nodes and co-expression links as edges. Typically, these can be derived by 
# computing gene–gene correlations or mutual information from single-cell RNA-seq data and thresholding to define edges.
#
# Last Modified: 05/04/25
# -----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(igraph)

# Set age groups and cell types from your object
age_groups <- as.character(unique(so.integrated$age_group))
cell_types <- as.character(unique(so.integrated$cell_type))

# Create output directory for GNs
output_dir <- file.path(wd.de.data, "coexpression_networks")
dir.create(output_dir, showWarnings = FALSE)

Idents(so.integrated) <- "cell_type"
# Loop over each (age group × cell type) pair
for (age in age_groups) {
	  for (ctype in cell_types) {
	  	  age_chr <- as.character(age)
	  	  ctype_chr <- as.character(ctype)
		    message(paste("Processing:", age_chr, "-", ctype_chr))
		
		    # Subset cells for this age group and cell type
	  	  cells_subset <- WhichCells(so.integrated, 
	  	  	    																						expression = age_group == age_chr & cell_type == ctype_chr)
	  	  
		    if (length(cells_subset) < 50) {
			      warning(paste("Skipping", age, "-", ctype, ": fewer than 50 cells"))
			      next
		    }
	    	subset_obj <- subset(so.integrated, cells = cells_subset)
		
		    # Get log-normalized expression matrix
		    expr_mat <- GetAssayData(subset_obj, slot = "data")
		    # Filter out lowly expressed genes
		    expr_mat <- expr_mat[rowSums(expr_mat > 0) >= 10, ]
		    # Compute Pearson correlation between genes
		    gene_cor <- cor(t(as.matrix(expr_mat)), method = "pearson")
		
		    # Threshold weak correlations
		    cor_threshold <- 0.3
		    gene_cor[abs(gene_cor) < cor_threshold] <- 0
		    diag(gene_cor) <- 0
		
		    # Convert to igraph object
		    gene_graph <- graph.adjacency(gene_cor, mode = "undirected", weighted = TRUE)
		
		    # Save the network
		    graph_path <- file.path(
			      output_dir, 
			      paste0("GN_", age, "_", gsub(" ", "_", ctype), ".rds")
		    )
		    saveRDS(gene_graph, graph_path)
	  }
}
