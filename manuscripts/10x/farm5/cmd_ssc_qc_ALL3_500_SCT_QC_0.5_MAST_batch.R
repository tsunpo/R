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

# -----------------------------------------------------------------------------
# Methods: DE using MAST
# Last Modified: 17/10/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
# > dim(so.integrated)
# [1] 5000 70738

DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
so.integrated$cell_type <- Idents(so.integrated)

# Create age groups ("young" and "old")
#median_age <- median(so.integrated$age)  # Define the threshold for age
#so.integrated$age_group <- ifelse(so.integrated$age <= median_age, "young", "old")

# List of unique cell types
cell_types <- unique(so.integrated$cell_type)  # Adjust if your metadata column for cell type is named differently
# Initialize a list to store DE results
de_results_list <- list()

# Loop over each cell type
for (ct in cell_types) {
	  # Subset Seurat object to include only the cells of the current cell type
	  cells_of_ct <- WhichCells(so.integrated, expression = cell_type == ct)
	  so_ct <- subset(so.integrated, cells = cells_of_ct)
	
	  # Perform differential expression analysis between "old" and "young" using MAST
	  de_results <- FindMarkers(
		    so_ct,
		    #ident.1 = "old",
		    #ident.2 = "young",
		    #group.by = "age_group",  # Use age group for ident.1 and ident.2
		    ident.1 = NULL,  # No group comparison, just association with age
		    group.by = "age",  # Using age as a continuous variable
		    test.use = "MAST",  # MAST for zero-inflation modeling
		    latent.vars = "batch"  # Control for batch effect
	  )
	
	  # Store the results for this cell type
  	de_results_list[[ct]] <- de_results
  	
  	# Define color categories based on q_value and log2_fold_change conditions
  	de_results$color <- "gray"  # Default color
  	de_results$color[de_results$p_val_adj < 1E-9 & de_results$avg_log2FC > 1] <- red
  	de_results$color[de_results$p_val_adj < 1E-9 & de_results$avg_log2FC < -1] <- blue
  	
  	# Create a new column for labeling significant genes
  	de_results$label <- NA
  	de_results$label[de_results$p_val_adj < 1E-9 & abs(de_results$avg_log2FC) > 1] <- rownames(de_results)
  	
  	# Plot the volcano plot
  	volcano_plot <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  		  geom_point(aes(color = color), alpha = 0.7) +
  		  scale_color_identity() +  # Use the color column for point colors
  		  geom_text_repel(aes(label = label), max.overlaps = 20, size = 3) +  # Add gene labels
  		  geom_vline(xintercept = 1, linetype = "dashed", color = red) +
  		  geom_vline(xintercept = -1, linetype = "dashed", color = blue) +
  		  theme_minimal() +
  		  xlab("Log2 FC") +
  		  ylab("-log10(P)") +
  		  ggtitle(ct) +
  		  theme(legend.position = "none",  # Hide the legend for custom color  theme(legend.position = "none",            # Hide the legend
  							plot.title = element_text(hjust = 0.5))  # Center the title
  	
  	# Print the volcano plot
  	pdf(file = file.path(wd.de.plots, paste0("MAST_volcano_", ct, "_RNA_batch.pdf")), width = 5, height = 5)
  	print(volcano_plot)
  	dev.off()
}

save(de_results_list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_MAST_RNA_batch.RData")))
