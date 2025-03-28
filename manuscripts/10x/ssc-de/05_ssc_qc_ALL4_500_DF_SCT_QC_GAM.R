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
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))

# -----------------------------------------------------------------------------
# Heat map for positively selected genes
# -----------------------------------------------------------------------------
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
#so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_scaled_rev.pdf")), width = 14, height = 6)
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
ordered_genes <- names(sort(correlation_values, decreasing = FALSE))
# Step 7: Apply New Sorting to Scaled Data
scaled_data <- scaled_data[ordered_genes, ]

# Step 8: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident",	disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))
# Step 9: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_pseudotime_FALSE.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
#so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
  	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_pseudotime_FALSE.pdf")), width = 14, height = 6)
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
so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

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

# -----------------------------------------------------------------------------
# Output an .h5ad file from a Seurat object 
# Last Modified: 28/03/25
# -----------------------------------------------------------------------------
#install.packages("remotes")
#remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)

h5Seurat <- file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/ngs/10x/", "ssc_annotated_monocle3+phase.h5Seurat")
SaveH5Seurat(so.integrated, filename = h5Seurat)
Convert(h5Seurat, dest = "h5ad")

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 20/01/25; 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
# > dim(so.integrated)
# [1] 35742 84822

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "PMC", "Leydig", "Endothelial", "Macrophage", "19", "6", "23")
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
#clusters_to_keep <- c("Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3", "Unknown")
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)
# > dim(so.integrated)
# [1]  5000 58098

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
	  )
ggsave(file.path(wd.de.plots, paste0("DimPlot_no-title_6*6_new-color_16_Spermatogenesis_19-6-23.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

##
# Step 1: Set the default assay to "integrated"
DefaultAssay(so.integrated) <- "integrated"
# Step 2: Extract the integrated expression data from Seurat
expression_matrix <- GetAssayData(so.integrated, assay = "integrated", slot = "data")
# Step 3: Create Monocle3 CDS object manually with the integrated data
cds <- new_cell_data_set(
	  expression_matrix,
	  cell_metadata = so.integrated@meta.data,
	  gene_metadata = data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
)
# Step 4: Transfer the UMAP embeddings if available in Seurat
reducedDims(cds)$UMAP <- Embeddings(so.integrated, reduction = "umap")

cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)
# Perform trajectory analysis using precomputed UMAP
cds <- cluster_cells(cds, reduction_method = "UMAP")
# Learn the principal graph
cds <- learn_graph(cds, use_partition = F)

# Order cells in pseudotime
#cds <- order_cells(cds)
# Identify root cells based on some known marker or cluster
root_cells <- WhichCells(so.integrated, idents = "Stage 0")  # Replace "early_cluster_name" with your specific cluster
# Ensure these root cells are present in `cds`
root_cells_in_cds <- root_cells[root_cells %in% colnames(cds)]
# Order cells using the root cells
cds <- order_cells(cds, root_cells = root_cells_in_cds)

plot_cells(cds, show_trajectory_graph = F)

# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "pseudotime",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=0) +
	  guides(color = guide_colorbar(title = "Pseudotime")) +  # Force rename legend title
	  theme_classic() +
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("pseudotime_23_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_integrated_19-23-6_6.5x6.png")), plot = monocle3, dpi = 300, width = 6.5, height = 6)

# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#tscan
library(ggbeeswarm)
library(colorRamps)
library(viridisLite)

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$idents <- Idents(so.integrated)
pdata_cds$Phase <- colData(cds)$Phase  # Ensure Phase is included

#pseudotime_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
#pseudotime_palette <- inferno(100)
pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)

ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = pseudotime_monocle3)) +  # Use pseudotime as the color
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_gradientn(colors = pseudotime_palette) +  # Apply the pseudotime gradient
	  scale_y_discrete(limits = rev(levels(pdata_cds$idents))) +  # Reverse y-axis order
	  theme_classic() +
  	xlab("Pseudotime") +
   ylab("") +
	  theme(plot.title = element_blank(),  # Remove title
			  				axis.text = element_text(size = 16),  # Increase axis text size
				  			axis.title = element_text(size = 18), # Increase axis title size
						  	legend.position = "none")  # Remove legend
ggsave(file.path(wd.de.plots, paste0("pseudotime_23_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered_integrated_19-23-6_7x6_rev.png")), plot = ordered, dpi = 300, width = 7, height = 6)

##
ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = Phase)) +  # Use pseudotime as the color
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply the same Phase colours as DimPlot
	  scale_y_discrete(limits = rev(levels(pdata_cds$idents))) +  # Reverse y-axis order
	  theme_classic() +
	  xlab("Pseudotime") +
	  ylab("") +
	  guides(color = guide_legend(override.aes = list(size = 5))) +  # Increase legend dot size
	  theme(plot.title = element_blank(),  # Remove title
			  				axis.text = element_text(size = 16),  # Increase axis text size
						  	axis.title = element_text(size = 18), # Increase axis title size
	  						legend.position = c(1, 1),  # Move legend to top-right
	  						legend.justification = c(1, 1),  # Align legend correctly
	  						legend.text = element_text(size = 16),  # Increase legend text size
	  						legend.title = element_text(size = 16))
ggsave(file.path(wd.de.plots, paste0("pseudotime_23_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered_integrated_19-23-6_7x6_rev_phase_legend.png")), plot = ordered, dpi = 300, width = 7, height = 6)

# -----------------------------------------------------------------------------
# Cell cycle analysis
# Last Modified: 05/02/25
# -----------------------------------------------------------------------------
nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))

cc.genes <- Seurat::cc.genes.updated.2019  # Updated gene lists
s.genes <- cc.genes$s.genes  # S phase markers
g2m.genes <- cc.genes$g2m.genes  # G2/M phase markers

so.integrated <- CellCycleScoring(so.integrated, 
																																		s.features = s.genes, 
																																		g2m.features = g2m.genes, 
																																		set.ident = F)  # Updates the cell identities (Idents(seurat_obj)) with the assigned cell cycle phases if set.ident = TRUE.

# Ensure Phase is a factor and ordered correctly
so.integrated@meta.data$Phase <- factor(so.integrated@meta.data$Phase, levels = c("G1", "S", "G2M"))
# Define consistent Phase colors (matching your figure)
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")
# Now plot the UMAP with the correct legend order and custom colors
dim_plot <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap") + 
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply custom colors & legend order
  	labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("Cell-cycle_annotated_19-6-23_DimPlot_no-title_8*8_new-color_16.png")), plot = dim_plot, width = 8, height = 8, dpi = 300)

VlnPlot(so.integrated, group.by = "Phase", cols= phase_colors, features = c("S.Score", "G2M.Score"), pt.size = 0)

# -----------------------------------------------------------------------------
# After running Monocle3
# -----------------------------------------------------------------------------
# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "Phase",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=0) +
	  theme_classic() +
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) + # Apply new color
		 theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("Cell-cycle_Phase_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_integrated_new-color_6.5x6_16.png")), plot = monocle3, width = 6.5, height = 6, dpi = 300)

library(ggbeeswarm)
library(colorRamps)
library(viridisLite)
library(ggplot2)

# Extract metadata
pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$Phase <- factor(so.integrated@meta.data$Phase, levels = c("G2M", "S", "G1"))  # Use cell cycle phase

# Define the Seurat default cell cycle colours (modify if needed)
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Plot pseudotime ordered by Phase
ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = Phase, 
																						colour = Phase)) +  # Use Phase for colour
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply the same Phase colours as DimPlot
	  scale_y_discrete(labels = c("G1" = "G1", "S" = "S", "G2M" = "G2/M")) +  # Rename y-axis labels
	  theme_classic() +
	  xlab("Pseudotime") + 
	  ylab("Phase") +
	  theme(plot.title = element_blank(),  # Remove title
	  	   axis.text = element_text(size = 16),  # Increase axis text size
							axis.title = element_text(size = 18), # Increase axis title size
							legend.position = "none")  # Remove legend

# Save the plot
ggsave(file.path(wd.de.plots, paste0("cell_cycle_phase_pseudotime_by_phase_new-color_16_6x6.png")),	plot = ordered, dpi = 300, width = 6, height = 6)

save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
