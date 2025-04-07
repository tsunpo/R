# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 20/02/25; 25/07/24; 14/03/24
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

# -----------------------------------------------------------------------------
# 4.1 Load & Stack All Embeddings
# -----------------------------------------------------------------------------
library(tidyverse)
library(uwot)

# Path to your saved embeddings
embedding_dir <- file.path(wd.de.data, "gene_embeddings")

# Get all embedding files
embedding_files <- list.files(embedding_dir, pattern = "^embedding_.*\\.rds$", full.names = TRUE)

# Read and combine
embedding_df <- purrr::map_dfr(embedding_files, function(f) {
	  emb <- readRDS(f)
	  if (!is.null(emb)) {
		    emb_df <- as.data.frame(emb)
		    emb_df$gene <- rownames(emb)
		    emb_df$embedding_key <- gsub("embedding_|\\.rds", "", basename(f))
		    return(emb_df)
	  } else {
		    return(NULL)
	  }
})

# -----------------------------------------------------------------------------
# 4.2 UMAP on Each Embedding
# -----------------------------------------------------------------------------
# Run UMAP per (age_group × cell_type)
umap_list <- embedding_df %>%
	  group_split(embedding_key) %>%
	  map(function(df) {
		    umap_res <- umap(df %>% select(-gene, -embedding_key), n_neighbors = 15, min_dist = 0.1)
		    tibble(
			   gene = df$gene,
			   embedding_key = df$embedding_key[1],
			   UMAP1 = umap_res[,1],
			   UMAP2 = umap_res[,2]
		    )
	  }) %>%
	  bind_rows()

# Plot all UMAPs in facets
plot <- ggplot(umap_list, aes(x = UMAP1, y = UMAP2)) +
   geom_point(size = 0.5, alpha = 0.7) +
   facet_wrap(~embedding_key, scales = "free") +
   theme_minimal() +
   ggtitle("UMAP of Gene Embeddings by Age Group and Cell Type")

# Save the plot
ggsave(file.path(wd.de.plots, paste0("UMAP_embeddings_age-type.png")), plot = plot, dpi = 300, width = 6, height = 6)
