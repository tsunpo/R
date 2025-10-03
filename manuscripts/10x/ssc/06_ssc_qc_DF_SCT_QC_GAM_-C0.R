# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 11/09/25
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch125/casm/staging/team294/ty2"   ## ty2@farm22
#wd <- "/Users/ty2/Work/sanger/ty2"                   ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE)
wd.rna.raw <- file.path(wd.rna, "ngs", "10x")

wd.de    <- file.path(wd.rna, "analysis", paste0(base, ""))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase.RData")))

# -----------------------------------------------------------------------------
# Heat map for positively selected genes
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(pheatmap)

#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_joined.RData")))
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "ARHGAP35", "CBL", "SSX1", "RBM12", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
#so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Step 2: Extract the scaled data for the genes in mn7
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]
#save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))

# Step 3: Calculate the average expression for each gene in each cluster
average_expression_per_cluster <- sapply(levels(so.integrated), function(cluster) {
	  # Ensure no NAs by filling missing values
	  rowMeans(scaled_data[, Idents(so.integrated) == cluster, drop = FALSE], na.rm = TRUE)
})

# Step 4: Sort the genes sequentially by expression across clusters
ordered_genes <- rownames(scaled_data)[do.call(order, as.data.frame(-average_expression_per_cluster))]
scaled_data <- scaled_data[ordered_genes,]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + theme(axis.text.y = element_text(size = 12))
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_mn7_20x12_monocle3+phase_2_scaled_data.png")),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_scaled.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Last Modified: 23/02/25
# -----------------------------------------------------------------------------
# Step 2: Extract the scaled data for the genes in mn7
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]

# Step 4: Extract Pseudotime Values from Monocle3
pseudotime_values <- cds@principal_graph_aux$UMAP$pseudotime  # Extract Monocle3 pseudotime
pseudotime_values <- pseudotime_values[match(colnames(scaled_data), names(pseudotime_values))]  # Match order

# Step 5: Compute Correlation Between Gene Expression and Pseudotime
correlation_values <- apply(scaled_data, 1, function(gene_expression) {
	  cor(gene_expression, pseudotime_values, method = "spearman", use = "pairwise.complete.obs")
})

# Step 6: Sort Genes by Pseudotime Correlation (High Positive Correlation First)
ordered_genes <- names(sort(correlation_values, decreasing = TRUE))
# Step 7: Apply New Sorting to Scaled Data
scaled_data <- scaled_data[ordered_genes, ]

# Step 8: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident",	disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))
# Step 9: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_pseudotime_TRUE.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
#so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
  	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_pseudotime_TRUE.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Smooth Expression Trends (on pseudotime)
# Last Modified: 24/02/25
# -----------------------------------------------------------------------------
library(mgcv)  # For GAM-based smoothing
library(reshape2)

# Step 2: Extract scaled data for genes of interest
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
gene_expression <- GetAssayData(so.integrated, assay = "RNA", layer = "data")[matched_genes, ]  # Use log-normalized expression
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]

# Step 3: Extract Pseudotime Values from Monocle3
pseudotime_values <- cds@principal_graph_aux$UMAP$pseudotime  # Extract pseudotime
pseudotime_values <- pseudotime_values[match(colnames(scaled_data), names(pseudotime_values))]  # Match order

# Step 4: Fit a Smooth Trend Line for Each Gene Using Generalized Additive Model (GAM)
smoothed_trends <- sapply(rownames(gene_expression), function(gene) {
	  fit <- gam(gene_expression[gene, ] ~ s(pseudotime_values, bs="cs"), method="REML")  # Use cubic smoothing spline
	  predict(fit)  # Get predicted smooth values
})

# Convert to matrix format
smoothed_trends <- as.matrix(smoothed_trends)
# Step 5: Compute Peak Pseudotime for Each Gene (First Max Expression in Smoothed Data)
peak_pseudotime <- apply(smoothed_trends, 2, function(x) which.max(x))
# Step 6: Order Genes by Peak Activation Time in Pseudotime
ordered_genes <- names(sort(peak_pseudotime, decreasing = FALSE))  # Earlier peak = higher rank
# Step 7: Apply New Sorting to Scaled Data
scaled_data <- scaled_data[ordered_genes, ]

# Step 8: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident",	
																					disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))

# Step 9: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_smoothed_pseudotime_.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
#so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_smoothed_pseudotime.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Smooth Expression Trends (on clusters)
# Last Modified: 24/02/25
# -----------------------------------------------------------------------------
# Step 2: Extract cluster identities
cluster_ids <- Idents(so.integrated)  # Get cluster identities per cell
# Ensure cluster_ids follow the correct order
desired_order <- c("Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3", "Leptotene", "Zygotene", "Pachytene", "Diplotene", "Meiotic division", "Early spermatid", "Late spermatid")
cluster_ids <- factor(cluster_ids, levels = desired_order)

# Step 3: Fit GAM smoothing per gene per cluster
smoothed_trends <- sapply(rownames(gene_expression), function(gene) {
	  cluster_means <- tapply(gene_expression[gene, ], cluster_ids, mean, na.rm = TRUE)  # Mean per cluster
	  #cluster_order <- as.numeric(as.factor(names(cluster_means)))  # Convert cluster names to numeric order
	  cluster_order <- match(names(cluster_means), desired_order)  # Correct biological order
	  
	  # Fit GAM model for smoothing within clusters
	  fit <- gam(cluster_means ~ s(cluster_order, bs="cs"), method="REML")
	  predict(fit)  # Get smoothed expression values per cluster
})

# Convert to matrix format
smoothed_trends <- as.matrix(smoothed_trends)
# Step 4: Compute Peak Cluster (Where Gene Expression Peaks)
peak_cluster <- apply(smoothed_trends, 2, function(x) which.max(x))  # Find the cluster with highest expression
# Step 5: Order Genes by Peak Expression in Clusters
ordered_genes <- names(sort(peak_cluster, decreasing = FALSE))  # Sort genes by earliest peak expression
# Step 6: Apply New Sorting to Scaled Data for Heatmap Visualization
scaled_data <- scaled_data[ordered_genes, ]

# Step 7: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", 
																					disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))

# Step 8: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_smoothed_clusters_desired_order.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_smoothed_clusters_desired_order.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Smooth Expression Trends (on clusters; using smoothed trends)
# Last Modified: 26/02/25
# -----------------------------------------------------------------------------
# Step 7: Generate the Heatmap Using `smoothed_trends`
heatmap <- DoHeatmap(
	  object = NULL,  # Remove Seurat object to use a custom matrix
	  features = rownames(smoothed_trends),  # Use smoothed expression trends
	  assay = NULL,  # Not required since we're using custom data
	  slot = NULL  # Not needed as `smoothed_trends` is precomputed
) + 
	  scale_fill_gradientn(colors = c("blue", "white", "red")) +  # Adjust heatmap colors
	  theme(axis.text.y = element_text(size = 12)) +
	  NoLegend()

# Save Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_smoothed_clusters_desired_order_smoothed.png"), plot = heatmap, width = 10, height = 8, dpi = 300)
