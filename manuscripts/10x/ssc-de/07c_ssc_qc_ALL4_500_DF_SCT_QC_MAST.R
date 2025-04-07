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
# MAST
# -----------------------------------------------------------------------------
nfeatures <- 5000
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_markers_MAST_SSC.RData")))

filtered_markers <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
de_gene_list <- unique(filtered_markers$gene)
length(de_gene_list)
# [1] 2706

# Count the number of marker genes in each cluster
marker_gene_counts <- filtered_markers %>%
	  dplyr::count(cluster)
print(marker_gene_counts)
# 1 Stage 0   1069
# 2 Stage 0A   307
# 3 Stage 0B   359
# 4 Stage 1    371
# 5 Stage 2    808
# 6 Stage 3    501

# -----------------------------------------------------------------------------
# K-means
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_markers_MAST_SSC_SSC_wss_k=5.RData")))

#plot(1:15, wss, type = "b", pch = 19, frame = FALSE, 
#					xlab = "Number of clusters K", 
#					ylab = "Total within-clusters sum of squares")


# Split the gene names by their cluster assignments
genes_by_cluster <- split(names(gene_clusters), gene_clusters)

# Print the list of genes in each cluster
for (i in 1:k) {
	  cat(paste("Cluster", i, "genes:\n"))
	  length(genes_by_cluster[[i]])
	  cat("\n")  # Add a line break for better readability
}

# Assuming gene_clusters contains cluster assignments
go_summary <- list()  # Create a list to store GO summaries for each cluster

# Assuming gene_clusters contains cluster assignments
for (i in unique(gene_clusters)) {
	  # Get genes in the current cluster
	  genes_in_cluster <- names(gene_clusters[gene_clusters == i])
	  
	  # Perform GO enrichment
	  go_enrichment <- enrichGO(gene = genes_in_cluster,
			  																								OrgDb = org.Hs.eg.db, 
				  																							keyType = "SYMBOL", 
					  																						ont = "BP",  # Use "BP" for Biological Process, or "MF", "CC" for other categories
						  																					pvalueCutoff = 0.05,
	  																										qvalueCutoff = 0.05)
	  
	  # Simplify GO terms to group similar terms together
	  go_enrichment_simplified <- simplify(
	  	  go_enrichment, 
	  	  cutoff = 0.7,  # Similarity cutoff for simplification
	  	  by = "p.adjust",  # Use adjusted p-value to rank terms
	  	  select_fun = min  # Keep the term with the minimum p.adjust
	  )
	  
	  # Use simplify to group GO terms under broader umbrella terms
	  bp_results_simplified <- simplify(go_enrichment, cutoff = 0.7, by = "p.adjust", measure = "Wang")
	  
	  # Summarize the umbrella terms and count their occurrences
	  umbrella_summary <- bp_results_simplified %>%
	  	  group_by(Description) %>%
	  	  summarize(n_BP_terms = n())  # Count number of significant BP terms under each umbrella term
	  
	  # Capitalize the first letter of each umbrella term
	  #umbrella_summary$Description <- str_to_title(umbrella_summary$Description)
	  
	  # Sort by the most common umbrella terms and take the top 3
	  top_umbrella_summary <- umbrella_summary %>%
	  	  arrange(desc(n_BP_terms)) %>%
	  	  head(3)
	 
	  # Add the summary to the list
	  go_summary[[paste0("Cluster_", i)]] <- top_umbrella_summary
}

# Convert the summary list into a data frame for easier viewing
go_summary_df <- do.call(rbind, lapply(names(go_summary), function(cluster) {
	  data.frame(Cluster = cluster, go_summary[[cluster]])
}))

# Print the GO summary table
print(go_summary_df)
write.csv(go_summary_df, file=file.path(wd.de.data, "GO_BiologicalProcess_Cluster_Umbrella_Summary_Top3.csv"), row.names = FALSE)






# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------
# Load necessary libraries
library(pheatmap)
library(Seurat)

# Step 1: Prepare the expression matrix
# Assuming `so.integrated` is your Seurat object with RNA expression data
# and `gene_clusters` contains the K-means clustering results for the genes

# Extract the expression matrix
#expression_matrix <- GetAssayData(so.integrated, assay = "RNA", slot = "data")

# Step 1: Create a factor with the desired order
cluster_order <- factor(gene_clusters, levels = c(1, 2, 3, 4, 5))

# Step 2: Subset the expression matrix to only include genes in your K-means clusters
ordered_genes <- names(gene_clusters[order(cluster_order)])  # Order genes by cluster
expression_matrix_ordered <- expression_matrix[ordered_genes, ]

# Create the heatmap
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

heatmap <- DoHeatmap(so.integrated, features = ordered_genes) + NoLegend() + theme(axis.text.y = element_blank())  # Remove gene names from the y-axis
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_", nfeatures, "_UMAP_resolution=0.5_MAST_kmeans_n=2706_pos=T_1,2,3,4,5.png")),
							plot = heatmap, width = 10, height = 20, dpi = 300)

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(cowplot)  # Make sure cowplot is loaded

# Ensure that gene_clusters is a factor and matches the ordered genes
gene_clusters_factor <- factor(gene_clusters[ordered_genes])

# Create the heatmap using DoHeatmap, without showing the y-axis gene names
heatmap <- DoHeatmap(
	  object = so.integrated, 
	  features = ordered_genes,  # Use the ordered genes
	  group.by = "ident",  # Use identity to order cells
	  slot = "data"  # Use normalized data (or "scale.data" if scaling is preferred)
) + NoLegend() + theme(axis.text.y = element_blank())  # Remove gene names from the y-axis

# Create a plot for the k-means cluster labels
cluster_labels <- data.frame(
	  y = seq_along(ordered_genes), 
	  cluster = as.numeric(gene_clusters_factor)
)

label_plot <- ggplot(cluster_labels, aes(x = 1, y = y, label = cluster)) +
	  geom_text(size = 3) +
	  theme_void() +
	  theme(axis.text = element_blank(),  # Remove axis text
			  				axis.ticks = element_blank(), # Remove axis ticks
				  			panel.grid = element_blank(), # Remove grid
						  	plot.margin = margin(0, 0, 0, 0)
	)

# Combine the heatmap and the cluster label plot
combined_plot <- cowplot::plot_grid(
	  heatmap, label_plot, 
	  ncol = 2, rel_widths = c(0.95, 0.05)  # Adjust the width ratio
)

# Save the combined plot
ggsave(
	  filename = file.path(wd.de.plots, paste0("heatmap_with_clusters_", nfeatures, "_UMAP_resolution=0.5_MAST_kmeans_DoHeatmap.png")),
	  plot = combined_plot, width = 12, height = 20, dpi = 300
)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Create the heatmap
# Ensure that gene_clusters is a factor so it can be used to group the rows
gene_clusters_factor <- factor(gene_clusters[ordered_genes])

# Create a heatmap with k-means cluster annotations on the right
heatmap <- DoHeatmap(
	  object = so.integrated, 
	  features = ordered_genes,  # Use the ordered genes
	  group.by = "ident",  # Use identity to order cells
	  slot = "data"  # Use normalized data (or "scale.data" if scaling is preferred)
) + NoLegend()

# Add k-means cluster labels as row annotations
heatmap <- heatmap + 
	  geom_text(aes(y = 1:nrow(expression_matrix_ordered), 
			  												label = as.numeric(gene_clusters_factor)), 
					  						size = 3, hjust = -0.5, vjust = 0.5)

# Save the heatmap with k-means cluster labels
ggsave(
	  filename = file.path(wd.de.plots, paste0("heatmap_", nfeatures, "_UMAP_resolution=0.5_MAST_kmeans_n=2544_with_clusters.png")),
	  plot = heatmap, width = 10, height = 20, dpi = 300
)



# Create a heatmap using DoHeatmap and split the rows by the k-means clusters
heatmap <- DoHeatmap(
	  object = so.integrated, 
	  features = ordered_genes,  # Use the ordered genes
	  group.by = "ident",  # Use identity for ordering cells
	  slot = "data",  # Use normalized data (or "scale.data" if scaling is preferred)
	  split.by = gene_clusters_factor  # Split rows by the k-means clusters
) + NoLegend() + theme(axis.text.y = element_blank())  # Remove gene names from the y-axis

# Save the heatmap
ggsave(
	  filename = file.path(wd.de.plots, paste0("heatmap_", nfeatures, "_UMAP_resolution=0.5_MAST_kmeans_n=2544_kmeans.png")),
	  plot = heatmap, 
	  width = 10, 
	  height = 20, 
	  dpi = 300
)





# Step 3: Generate the heatmap
# Using pheatmap to visualize the ordered expression matrix
pheatmap(
	  expression_matrix_ordered,
	  show_rownames = FALSE,  # Hide gene names for cleaner visualization
	  show_colnames = FALSE,  # Hide cell names for cleaner visualization
	  cluster_rows = FALSE,  # Do not cluster rows (genes are already ordered by K-means clusters)
	  cluster_cols = FALSE,   # Cluster cells (x-axis)
	  color = colorRampPalette(c("blue", "white", "red"))(100),  # Color gradient for expression levels
	  main = "Gene Expression Heatmap by K-means Clustering"
)

# Alternatively, use ComplexHeatmap for more control over annotations
library(ComplexHeatmap)

heatmap <- Heatmap(
	  expression_matrix_ordered, 
	  name = "Expression", 
	  show_row_names = FALSE, 
	  show_column_names = FALSE, 
	  cluster_rows = FALSE,  # Do not cluster rows
	  cluster_columns = FALSE,  # Cluster columns
	  col = colorRamp2(c(-2, 0, 2), c(blue, "white", red)),  # Set color gradient
	  heatmap_legend_param = list(title = "Z-score"), use_raster = TRUE 
)
# Draw the heatmap
draw(heatmap)







# -----------------------------------------------------------------------------
# K-means Clustering of Genes
# -----------------------------------------------------------------------------
# Order cells by their assigned cell types
cell_metadata <- so.integrated@meta.data
cell_metadata$cell_type <- Idents(so.integrated)
ordered_cells <- cell_metadata[order(cell_metadata$cell_type), ]
# Reorder the expression matrix by cell types
ordered_expression_matrix <- expression_matrix_de[, rownames(ordered_cells)]
ordered_expression_matrix <- as.matrix(ordered_expression_matrix)
	
ordered_gene_clusters <- gene_clusters[rownames(ordered_expression_matrix)]

# Load necessary libraries for heatmap
library(Seurat)
library(ggplot2)
# Create a data frame of ordered_gene_clusters to use as a custom feature annotation
ordered_gene_clusters_df <- data.frame(
	  gene = rownames(ordered_expression_matrix),
	  cluster = as.factor(ordered_gene_clusters)
)
# Set color gradient
col_fun <- scale_fill_gradientn(colors = c(blue, "white", red), limits = c(-2, 2))
# Create heatmap using DoHeatmap
heatmap <- DoHeatmap(
	  so.integrated, 
	  features = rownames(ordered_expression_matrix),  # Specify the ordered genes
	  group.by = "ident",  # Use cell identity or group to order cells
	  slot = "data"  # Use the normalized expression values
) + col_fun + NoLegend()
# Add labels to rows based on gene clusters (manually annotate or adjust DoHeatmap for your needs)
heatmap <- heatmap + facet_grid(rows = vars(ordered_gene_clusters_df$cluster))
# Save the heatmap using ggsave
#ggsave(
#	  filename = file.path(wd.de.plots, "DoHeatmap_res=0.5_-12-15-17_clustering_SSC_fc=2.png"), 
#	  plot = heatmap, 
#	  width = 10, 
#	  height = 10, 
#	  dpi = 300
#)
save(heatmap, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering_SSC_MAST_heatmap.RData")))

# -----------------------------------------------------------------------------
# K-means Clustering of Genes
# -----------------------------------------------------------------------------

# Convert Seurat markers to a list of gene symbols for each cluster
markers_list <- split(markers$gene, markers$cluster)

# Perform GO enrichment for each cluster
enrich_results <- lapply(markers_list, function(genes) {
	  enrichGO(gene         = genes,
			  							OrgDb        = org.Hs.eg.db,
					  					keyType      = "SYMBOL",
							  			ont          = "BP",         # Ontology: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
								  		pAdjustMethod = "BH",
						  				pvalueCutoff = 0.01,
						  				qvalueCutoff = 0.05)
})






