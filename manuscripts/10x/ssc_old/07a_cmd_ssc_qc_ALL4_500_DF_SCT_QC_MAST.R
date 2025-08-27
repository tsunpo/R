args <- commandArgs(trailingOnly = TRUE)
cell_type <- as.character(args[1])

# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 09/04/25
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
library(dplyr)
library(Seurat)

nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_clean.RData")))

# Ensure age groups are defined in your metadata
so.integrated@meta.data <- so.integrated@meta.data %>%
   mutate(age_group = cut(
      age,
      breaks = c(-Inf, 50, Inf),
      labels = c("20-50", "50-80"),
      right = FALSE
   ))

# Set age groups and cell types from your object
#cell_types <- as.character(unique(so.integrated$cell_type))
#for (ctype in cell_types) {
   # Subset by a specific cell type 
   so.integrated <- subset(so.integrated, idents = gsub("_", " ", cell_type))
   Idents(so.integrated) <- so.integrated$age_group
 
   # Switch to the SCT assay before doing MAST
   DefaultAssay(so.integrated) <- "SCT"
   so.integrated <- PrepSCTFindMarkers(so.integrated)
   
   # Perform differential expression analysis, controlling for batch
   # Replace 'batch' with the actual name of your batch variable in the metadata
   markers <- FindMarkers(
      so.integrated,
      ident.1 = "50-80",
      ident.2 = "20-50",
      group.by = "age_group",
      latent.vars = c("batch", "nFeature_RNA"),
      test.use = "MAST"  # if using SCTransform
   )
 
   # 5. View top markers
   # Set path to co-expression networks and output for embeddings
   wd.de.data <- file.path(wd.de.data, "FindMarkers_young_vs_old_MAST")
   wd.de.plots <- file.path(wd.de.plots, "FindMarkers_young_vs_old_MAST")
   dir.create(wd.de.data, showWarnings = FALSE)
   dir.create(wd.de.plots, showWarnings = FALSE)
   
   head(markers)
   save(markers, file=file.path(wd.de.data, paste0("markers_SCT_MAST_batch_nFeature_", cell_type, ".RData")))
   
   # Assume results_age is the DESeq2 results object for the effect of age
   # Convert to a data frame for ggplot2
   result_df <- as.data.frame(markers)
   result_df$gene <- rownames(result_df)
   result_df$color <- ifelse(result_df$avg_log2FC < -0.05, "blue",
                             ifelse(result_df$avg_log2FC > 0.05, "red", "grey"))
   
   # Create a volcano plot
   text.Log10.P    <- expression("-log" * ""[10] * " Adjusted " * italic("P"))
   text.Log2.FC       <- expression("Log" * ""[2] * " Fold Change")
   
   volcano_plot <- ggplot(result_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = color)) +
      scale_color_manual(values = c(blue, "grey", red)) +
      theme_minimal() +
      labs(title = paste0(cell_type),
           x = text.Log2.FC,
           y = text.Log10.P) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = red) +
      geom_vline(xintercept = -0.05, linetype = "dashed", color = blue) +
      geom_text_repel(data = subset(result_df, abs(avg_log2FC) > 0.05),
                      aes(label = gene), size = 3, max.overlaps = 20)
   
     # Print the volcano plot
     pdf(file = file.path(wd.de.plots, paste0("volcano_markers_SCT_MAST_batch_nFeature_", cell_type, ".pdf")), width = 5, height = 5)
     print(volcano_plot)
     dev.off()
#}
