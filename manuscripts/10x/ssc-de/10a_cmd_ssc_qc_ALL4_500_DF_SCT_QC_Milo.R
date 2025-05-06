args <- commandArgs(trailingOnly = TRUE)
ctype <- as.character(args[2])

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
# Pe-processing
# -----------------------------------------------------------------------------
library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(stringr)
library(BiocParallel)
library(tidyr)
library(lme4)
library(tibble)

nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_clean.RData")))
so.integrated@meta.data$cell_type <- Idents(so.integrated)

# -----------------------------------------------------------------------------
# miloR
# -----------------------------------------------------------------------------
# 1. Subset Seurat object to cell type
ct <- "Stage 0"
seurat_ct <- subset(so.integrated, subset = cell_type == ct)

# 2. Convert to Milo object
sce_ct <- as.SingleCellExperiment(seurat_ct)
milo_ct <- Milo(sce_ct)

# 3. Build graph and neighborhoods
milo_ct <- buildGraph(milo_ct, k = 30, d = 30)
milo_ct <- makeNhoods(milo_ct, prop = 0.1, k = 30, d = 30)

# 4. Count cells in each neighborhood per sample
milo_ct <- countCells(milo_ct, meta.data = colData(milo_ct), sample = "orig.ident")

# 5. Get the neighborhood count matrix
nhood_counts <- nhoodCounts(milo_ct)  # rows = nhoods, cols = samples
nhood_counts_t <- t(nhood_counts)

# Rename neighborhood columns
colnames(nhood_counts_t) <- paste0("Nhood.", seq_len(ncol(nhood_counts_t)))

# 6. Prepare sample metadata
meta <- unique(as.data.frame(colData(milo_ct))[, c("orig.ident", "age")])
meta$patient_id <- str_extract(meta$orig.ident, "^[^_]+")
meta$age <- as.numeric(meta$age)
meta$patient_id <- factor(meta$patient_id)
rownames(meta) <- meta$orig.ident

# 7. Transpose count matrix: long format
nhood_df <- as.data.frame(nhood_counts_t) %>%
   rownames_to_column("orig.ident") %>%
   pivot_longer(
      cols = starts_with("Nhood"),
      names_to = "nhood",
      values_to = "count"
   ) %>%
   left_join(meta, by = "orig.ident")

# 8. Loop over neighborhoods and fit GLMM
fit_results <- nhood_df %>%
   group_by(nhood) %>%
   filter(sum(count) > 0) %>%  # remove neighborhoods with no counts
   group_split() %>%
   map_dfr(function(df) {
      tryCatch({
         model <- glmer(count ~ age + (1 | patient_id), data = df, family = poisson())
         summary_model <- summary(model)
         coefs <- coef(summary_model)
         tibble(
            nhood = unique(df$nhood),
            estimate = coefs["age", "Estimate"],
            std_error = coefs["age", "Std. Error"],
            p_value = coefs["age", "Pr(>|z|)"]
         )
      }, error = function(e) {
         tibble(nhood = unique(df$nhood), estimate = NA, std_error = NA, p_value = NA)
      })
   })

# Prepare results
fit_results$fdr <- p.adjust(fit_results$p_value, method = "BH")
fit_results_df <- as.data.frame(fit_results)
rownames(fit_results_df) <- fit_results_df$nhood

# Build neighborhood graph (required for plotting)
milo_ct <- buildNhoodGraph(milo_ct)

# Plot DA neighborhoods — finally!
plotNhoodGraphDA(milo_ct, milo_res = fit_results_df, alpha = 0.1)


