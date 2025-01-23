#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
nfeatures <- as.numeric(args[1])

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
# K-means Clustering of Genes
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
# > dim(so.integrated)
# [1] 5000 70738

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "Endothelial and Macrophage", "Leydig")
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)
dim(so.integrated)
# [1]  5000 56396

# Create a new column in the metadata for spermatogonia type
so.integrated@meta.data$stage_column <- as.character(Idents(so.integrated))
table(so.integrated@meta.data$stage_column)
# Reassign the stages
so.integrated@meta.data$stage_column[so.integrated@meta.data$stage_column %in% c("Stage 0", "Stage 0A", "Stage 0B", "Stage 1")] <- "Undiff. spermatogonia"
so.integrated@meta.data$stage_column[so.integrated@meta.data$stage_column %in% c("Stage 2", "Stage 3")] <- "Diff. spermatogonia"
# Set the desired order of cell types
desired_order <- c(
	  "Undiff. spermatogonia", 
	  "Diff. spermatogonia", 
	  "Leptotene", 
	  "Zygotene", 
	  "Pachytene", 
	  "Diplotene", 
	  "Meiotic division", 
	  "Early spermatid", 
	  "Late spermatid"
)
# Convert stage_column to a factor with the specified order
so.integrated@meta.data$stage_column <- factor(so.integrated@meta.data$stage_column, levels = desired_order)
# Update Idents to reflect the new ordered factor
Idents(so.integrated) <- so.integrated@meta.data$stage_column
table(Idents(so.integrated))

# Extract expression data from Seurat object
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
expression_matrix <- GetAssayData(so.integrated[["RNA"]], layer = "data")
expression_matrix <- as.matrix(expression_matrix)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
filtered_markers <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
# Count the number of marker genes in each cluster
marker_gene_counts <- filtered_markers %>%
	  dplyr::count(cluster)
print(marker_gene_counts)
# 1 Undiff. spermatogonia  1369
# 2 Diff. spermatogonia     465
# 3 Leptotene               351
# 4 Zygotene                682
# 5 Pachytene              2378
# 6 Diplotene               627
# 7 Meiotic division        895
# 8 Early spermatid        1054
# 9 Late spermatid          892

save(so.integrated, expression_matrix, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering0_MAST_markers.RData")))

de_gene_list <- unique(filtered_markers$gene)
length(de_gene_list)
# [1] 6478
# [1] 

expression_matrix_de <- expression_matrix[de_gene_list, ]

# Determine optimal number of clusters (Elbow Method)
wss <- sapply(1:15, function(k){
	  kmeans(expression_matrix_de, centers = k, nstart = 10)$tot.withinss
})
#plot(1:15, wss, type = "b", pch = 19, frame = FALSE, 
#					xlab = "Number of clusters K", 
#					ylab = "Total within-clusters sum of squares")
save(so.integrated, expression_matrix, markers, expression_matrix_de, wss, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering0_MAST_fc=1.RData")))

# -----------------------------------------------------------------------------
# K-means Clustering of Genes
# -----------------------------------------------------------------------------
library(stringr)

# Set the number of clusters (k)
k <- 5  # You can adjust this based on your needs
# Perform K-means clustering on genes (rows)
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(expression_matrix_de, centers = k)
# Assign cluster labels to the genes
gene_clusters <- kmeans_result$cluster

save(so.integrated, expression_matrix, markers, expression_matrix_de, wss, kmeans_result, gene_clusters, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering_MAST_only.pos=T_k=5.RData")))

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
write.csv(go_summary_df, file=file.path(wd.de.data, "GO_BiologicalProcess_Cluster_Umbrella_Summary_Top3_MAST_k=5.csv"), row.names = FALSE)

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
ordered_genes <- names(gene_clusters)  # Order genes by cluster
expression_matrix_ordered <- expression_matrix[ordered_genes, ]

# Create the heatmap
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

heatmap <- DoHeatmap(so.integrated, features = ordered_genes) + NoLegend() + theme(axis.text.y = element_blank())  # Remove gene names from the y-axis
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_", nfeatures, "_UMAP_resolution=0.5_MAST_kmeans_n=6478_only.pos=T_k=5.png")),
							plot = heatmap, width = 10, height = 20, dpi = 300)




















# Order cells by their assigned cell types
cell_metadata <- so.integrated@meta.data
cell_metadata$cell_type <- Idents(so.integrated)
ordered_cells <- cell_metadata[order(cell_metadata$cell_type), ]
# Reorder the expression matrix by cell types
ordered_expression_matrix <- expression_matrix[, rownames(ordered_cells)]
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
#	  filename = file.path(wd.de.plots, "DoHeatmap_res=0.5_-12-15-17_clustering.png"), 
#	  plot = heatmap, 
#	  width = 10, 
#	  height = 10, 
#	  dpi = 300
#)

save(heatmap, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering_MAST_heatmap.RData")))
