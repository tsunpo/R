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
# 
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(pheatmap)

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
embedding_df <- purrr::map(embedding_files, function(f) {
   emb <- readRDS(f)
 
   if (is.null(emb)) return(NULL)
 
   # Confirm rownames are intact and character-type
   if (!is.null(rownames(emb)) && is.character(rownames(emb))) {
      emb_df <- as.data.frame(emb)
      emb_df$gene <- rownames(emb)
   } else {
      stop("Missing or malformed rownames in: ", f)
   }
 
   emb_df$embedding_key <- gsub("embedding_|\\.rds", "", basename(f))
   return(emb_df)
}) %>% bind_rows()

dplyr::count(embedding_df, gene, sort = TRUE) %>% head(10)
embedding_df %>% filter(grepl("^MT-", gene))

# -----------------------------------------------------------------------------
# 4.2 UMAP on Each Embedding
# -----------------------------------------------------------------------------
umap_df <- embedding_df %>%
  dplyr::group_split(embedding_key) %>%
  purrr::map_dfr(function(df) {
    mat <- df %>%
      dplyr::select(-gene, -embedding_key) %>%
      as.matrix()

    # Skip groups with fewer rows than n_neighbors
    if (nrow(mat) < 15) return(NULL)

    umap_res <- umap::umap(mat, n_neighbors = 15, min_dist = 0.1)

    tibble::tibble(
      gene = df$gene,
      embedding_key = df$embedding_key[1],
      UMAP1 = umap_res$layout[, 1],
      UMAP2 = umap_res$layout[, 2]
    )
  })

# -----------------------------------------------------------------------------
# Step 5: Cosine Drift Between Age Groups
# -----------------------------------------------------------------------------
library(stringr)

drift_results <- embedding_df %>%
   rowwise() %>%
   mutate(
      vec = list(c_across(starts_with("V"))),  # your embedding dimensions
      age_group = str_extract(embedding_key, "20-50|50-80"),
      cell_type = str_replace(embedding_key, "^(20-50|50-80)_", "")
   ) %>%
   ungroup() %>%
   select(gene, age_group, cell_type, vec) %>%
   pivot_wider(names_from = age_group, values_from = vec) %>%
   filter(!is.na(`20-50`), !is.na(`50-80`)) %>%
   mutate(
      cosine_similarity = purrr::map2_dbl(`20-50`, `50-80`, ~ sum(.x * .y) / (sqrt(sum(.x^2)) * sqrt(sum(.y^2)))),
      drift = 1 - cosine_similarity
   ) %>%
   select(gene, cell_type, cosine_similarity, drift)

drift_results %>%
   arrange(desc(drift)) %>%
   head(20)

# Summarise drift to ensure uniqueness per gene × cell_type
drift_summarised <- drift_results %>%
   group_by(gene, cell_type) %>%
   summarise(drift = mean(drift, na.rm = TRUE), .groups = "drop")

# Pivot to wide matrix
drift_mat <- drift_summarised %>%
   pivot_wider(names_from = cell_type, values_from = drift) %>%
   column_to_rownames("gene") %>%
   as.matrix()

# -----------------------------------------------------------------------------
# Step 6: Drift Heatmap
# -----------------------------------------------------------------------------
library(pheatmap)

drift_mat <- drift_results %>%
   pivot_wider(names_from = cell_type, values_from = drift) %>%
   column_to_rownames("gene") %>%
   as.matrix()

# Filter for genes with sufficient coverage
drift_mat_filtered <- drift_mat[rowSums(!is.na(drift_mat)) >= 5, , drop = FALSE]
drift_mat_clean <- drift_mat_filtered[complete.cases(drift_mat_filtered), , drop = FALSE]

# Plot if at least 2 genes remain
if (nrow(drift_mat_clean) >= 2) {
   pheatmap::pheatmap(
      drift_mat_clean,
      scale = "none",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      main = "Gene Drift Heatmap (1 - Cosine Similarity)"
   )
} else {
   message("Not enough genes for clustering.")
}









drift_vectors <- umap_annotated %>%
  select(gene, cell_type, age_group, UMAP1, UMAP2) %>%
  pivot_wider(
    names_from = age_group,
    values_from = c(UMAP1, UMAP2),
    names_sep = "_"
  ) %>%
  filter(!is.na(`UMAP1_20-50`), !is.na(`UMAP1_50-80`)) %>%
  mutate(
    dx = `UMAP1_50-80` - `UMAP1_20-50`,
    dy = `UMAP2_50-80` - `UMAP2_20-50`
  )

geom_segment(
  aes(xend = `UMAP1_50-80`, yend = `UMAP2_50-80`, color = sqrt(dx^2 + dy^2)),
  arrow = arrow(length = unit(0.1, "inches")),
  linewidth = 0.5, alpha = 0.6
)

library(ggplot2)

# Create one plot per cell type
unique_cell_types <- unique(drift_vectors$cell_type)

for (ct in unique_cell_types) {
  p <- drift_vectors %>%
    filter(cell_type == ct) %>%
    ggplot(aes(x = `UMAP1_20-50`, y = `UMAP2_20-50`)) +
    geom_segment(
      aes(xend = `UMAP1_50-80`, yend = `UMAP2_50-80`, color = sqrt(dx^2 + dy^2)),
      arrow = arrow(length = unit(0.1, "inches")),
      size = 0.5, alpha = 0.6
    ) +
    scale_color_viridis_c(name = "Drift magnitude") +
    labs(
      title = paste("Gene Drift Vectors:", ct),
      subtitle = "From Age Group 20-50 → 50-80",
      x = "UMAP1", y = "UMAP2"
    ) +
    theme_minimal(base_size = 12)

  # Save plot
  ggsave(filename = file.path(output_dir, paste0("drift_vectors_", gsub(" ", "_", ct), ".png")),
         plot = p, width = 6, height = 5)
}






# -------------------------------------------------------------------------
# Step 6: Visualization & Summary
# -------------------------------------------------------------------------
library(pheatmap)
library(tidyr)

# Recalculate drift if needed
drift_results <- map_dfr(unique(embedding_df_fixed$embedding_key), function(key) {
  age <- stringr::str_extract(key, "20-50|50-80")
  ctype <- stringr::str_replace(key, "^(20-50|50-80)_", "")
  e <- embedding_df_fixed %>% filter(embedding_key == key)
  tibble(
    gene = e$gene,
    age_group = age,
    cell_type = ctype,
    vec = split(e[, -(1:2)], seq_len(nrow(e)))
  )
}) %>%
  pivot_wider(names_from = age_group, values_from = vec) %>%
  filter(!is.na(`20-50`), !is.na(`50-80`)) %>%
  mutate(
    cosine_similarity = map2_dbl(`20-50`, `50-80`, ~ sum(.x * .y) / (sqrt(sum(.x^2)) * sqrt(sum(.y^2)))),
    drift = 1 - cosine_similarity
  ) %>%
  select(gene, cell_type, cosine_similarity, drift)

# Create matrix: genes x cell types
drift_mat <- drift_results %>%
  select(gene, cell_type, drift) %>%
  pivot_wider(names_from = cell_type, values_from = drift) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Clean up NA rows
drift_mat_clean <- drift_mat[complete.cases(drift_mat) & !apply(drift_mat, 1, anyNA), ]

# Plot heatmap
pheatmap(drift_mat_clean,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Gene Drift Heatmap (1 - Cosine Similarity)")












library(tidyverse)
library(pheatmap)
library(tidyr)

# ---- 1. Top 10 Genes with Highest Drift per Cell Type ----
top_drifting_genes <- drift_results %>%
  group_by(cell_type) %>%
  slice_max(order_by = drift, n = 10) %>%
  arrange(cell_type, desc(drift))

# Optional: View it
print(top_drifting_genes)

# ---- 2. Heatmap of Drift Values (Genes × Cell Types) ----
drift_mat <- drift_results %>%
  select(gene, cell_type, drift) %>%
  pivot_wider(names_from = cell_type, values_from = drift) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Remove rows with any NA, NaN, or Inf values
drift_mat_clean <- drift_mat[complete.cases(drift_mat) & 
                             !apply(drift_mat, 1, function(x) any(is.infinite(x))), ]

# Plot the heatmap (optional: adjust scale = "row" to normalize across genes)
pheatmap(
  drift_mat_clean, 
  scale = "none", 
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Gene Drift Heatmap (1 - Cosine Similarity)"
)

# ---- 3. Histogram of Drift per Cell Type ----
ggplot(drift_results, aes(x = drift, fill = cell_type)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Gene Drift Distribution by Cell Type",
    x = "1 - Cosine Similarity (Drift)",
    y = "Gene Count"
  )

# ---- 4. Save Results to CSV ----
output_path <- file.path(wd.de.data, "cosine_drift_by_gene.csv")
write_csv(drift_results, output_path)
message("✅ Drift results saved to: ", output_path)
