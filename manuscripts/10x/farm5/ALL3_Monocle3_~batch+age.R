#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

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
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_cds_36601.RData")))

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 09/09/24
# -----------------------------------------------------------------------------
cell_types <- unique(colData(cds)$celltype_column)
results_list <- list()

for (cell_type in cell_types) {
	  # Subset for the current cell type
	  cds_subset <- cds[, colData(cds)$celltype_column == cell_type]
	  
	  # Run the model to test for association with age, correcting for batch effects
	  fit_models_age <- fit_models(cds_subset, 
																														model_formula_str = "~age + batch",  # Include batch as a covariate
																														expression_family = "negbinomial")
	  # Extract the model coefficients
	  coef_table <- coefficient_table(fit_models_age)
	  # Filter significant genes based on p-value or q-value
	  significant_genes <- subset(coef_table, term == "age" & q_value < 1)
	
	  # Store results for this cell type
	  results_list[[cell_type]] <- significant_genes
	  
	  # Convert to a data frame for ggplot2
	  result_df <- as.data.frame(significant_genes)
	  result_df$gene <- result_df$gene_short_name
	  result_df$color <- ifelse(result_df$estimate < -0.05, "blue",
	  																										ifelse(result_df$estimate > 0.05, "red", "grey"))
	  
	  # Create a volcano plot
	  volcano_plot <- ggplot(result_df, aes(x = estimate, y = -log10(q_value))) +
	  	  geom_point(aes(color = color)) +
	  	  scale_color_manual(values = c(blue, "grey", red)) +
	  	  theme_minimal() +
	  	  labs(title = paste0(cell_type),
	  			  			x = "Log2 Fold Change",
	  				  		y = "-log10 Adjusted P-value") +
	  	  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
	  	  geom_vline(xintercept = 0.05, linetype = "dashed", color = red) +
	    	geom_vline(xintercept = -0.05, linetype = "dashed", color = blue) +
	    	geom_text_repel(data = subset(result_df, abs(estimate) > 0.05),
	  		  															aes(label = gene), size = 3, max.overlaps = 20)
	  
	    # Print the volcano plot
	    pdf(file = file.path(wd.de.plots, paste0("volcano_Monocle3_~batch+age_", cell_type, ".pdf")), width = 5, height = 5)
	    print(volcano_plot)
	    dev.off()
}

save(results_list, file=file.path(wd.de.data, paste0("Monocle3_~batch+age.RData")))
