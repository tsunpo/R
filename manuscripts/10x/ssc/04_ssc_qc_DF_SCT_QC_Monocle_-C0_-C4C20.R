# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 19/08/25
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

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_-C4C20_ordered_annotated.RData")))

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "PMC", "Leydig", "Endothelial and Macrophage")
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
#clusters_to_keep <- c("Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3", "Unknown")

so.integrated$idents <- Idents(so.integrated) 
so.integrated <- subset(so.integrated, idents = clusters_to_keep)

# Rename specific cluster identity
Idents(so.integrated) <- plyr::mapvalues(
   Idents(so.integrated),
   from = "Meiotic division",
   to   = "Meiotic\ndivision"
)
Idents(so.integrated) <- plyr::mapvalues(
   Idents(so.integrated),
   from = "Early spermatid",
   to   = "Early\nspermatid"
)

dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
   labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
   theme(
      plot.title = element_blank(),  # Remove title
      axis.title = element_text(size = 18),  # Increase axis title size
      axis.text = element_text(size = 16),  # Increase axis tick label size
   )
ggsave(file.path(wd.de.plots, paste0("DimPlot_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_-C4C20_6*6_SSC.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

Idents(so.integrated) <- plyr::mapvalues(
   Idents(so.integrated),
   from = "Meiotic\ndivision",
   to   = "Meiotic division"
)
Idents(so.integrated) <- plyr::mapvalues(
   Idents(so.integrated),
   from = "Early\nspermatid",
   to   = "Early spermatid"
)

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 20/01/25; 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
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
      legend.title = element_text(size = 16),  # Increase legend title size
      legend.position = c(0.95, 0.05),       # top-right corner
      legend.justification = c("right", "bottom")  # anchor to that corner
   )
ggsave(file.path(wd.de.plots, paste0("pseudotime_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_C4C20_6x6.png")), plot = monocle3, dpi = 300, width = 6, height = 6)

# -----------------------------------------------------------------------------
# Methods: ggbeeswarm
# Last Modified: 20/01/25; 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
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
ggsave(file.path(wd.de.plots, paste0("pseudotime_ggbeeswarm_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_C4C20.png")), plot = ordered, dpi = 300, width = 7, height = 6)

## After running cell cycle analysis
ordered <- ggplot(as.data.frame(pdata_cds), 
                  aes(x = pseudotime_monocle3, 
                      y = idents, colour = phase_colors)) +  # Use pseudotime as the color
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
         legend.title = element_text(size = 16),
         legend.position = c(0.95, 0.05),       # top-right corner
         legend.justification = c("right", "bottom"))  # anchor to that corner)
ggsave(file.path(wd.de.plots, paste0("pseudotime_ggbeeswarm_cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_-C4C20.png")), plot = ordered, dpi = 300, width = 7, height = 6)

# -----------------------------------------------------------------------------
# Cell cycle analysis
# Last Modified: 05/02/25
# -----------------------------------------------------------------------------
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
# Plot the UMAP with Phase colors
dim_plot <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap", label = FALSE) + 
   scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply custom colors & legend order
   labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
   theme(
      plot.title = element_blank(),  # Remove title
      axis.title = element_text(size = 18),  # Increase axis title size
      axis.text = element_text(size = 16),  # Increase axis tick label size
      legend.text = element_text(size = 16),  # Increase legend text size
      legend.title = element_text(size = 16),  # Increase legend title size
      legend.position = c(0.95, 0.05),        # bottom-right corner
      legend.justification = c("right", "boottom")  # anchor to that corner
   ))
# **Extract cluster positions separately**
cluster_positions <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap", label = TRUE, repel = TRUE)
# **Overlay cluster labels on the same plot**
dim_plot <- LabelClusters(plot = dim_plot, id = "Phase", text = cluster_positions$data$Phase, repel = TRUE, size = 6, color = "black")
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_C4C20_6x6.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

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
ggsave(file.path(wd.de.plots, paste0("Cycle_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_integrated_new-color_6.5x6_16.png")), plot = monocle3, width = 6.5, height = 6, dpi = 300)

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

save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_monocle3+phase.RData")))












# -----------------------------------------------------------------------------
# Cell cycle analysis
# Last Modified: 10/02/25
# -----------------------------------------------------------------------------
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
# Plot the UMAP with Phase colors
dim_plot <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap", label = FALSE) + 
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply custom colors & legend order
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16),  # Increase legend title size
		    legend.position = c(0.95, 0.05),        # bottom-right corner
		    legend.justification = c("right", "bottom")  # anchor to that corner
	  )
# **Extract cluster positions separately**
cluster_positions <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap", label = TRUE, repel = TRUE)
# **Overlay cluster labels on the same plot**
dim_plot <- LabelClusters(plot = dim_plot, id = "Phase", text = cluster_positions$data$Phase, repel = TRUE, size = 6, color = "black")
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_C4C20_6x6.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

VlnPlot(so.integrated, group.by = "Phase", cols= phase_colors, features = c("S.Score", "G2M.Score"), pt.size = 0)
VlnPlot(so.integrated, group.by = "seurat_clusters", features = c("S.Score"), pt.size = 0) + NoLegend()
VlnPlot(so.integrated, group.by = "seurat_clusters", features = c("G2M.Score"), pt.size = 0) + NoLegend()

##
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggrepel)

# Ensure `seurat_clusters` exists
so.integrated$seurat_clusters <- Idents(so.integrated)

# Define consistent Phase colors
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Extract UMAP embeddings
umap_data <- as.data.frame(Embeddings(so.integrated, reduction = "umap"))
umap_data$seurat_clusters <- as.factor(so.integrated$seurat_clusters)  # Ensure clusters are factors

# Compute cluster centroids
cluster_centers <- umap_data %>%
	  group_by(seurat_clusters) %>%
	  summarise(UMAP_1 = median(umap_1), UMAP_2 = median(umap_2))  # Use median to find label position

# Create UMAP plot colored by Phase
dim_plot <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap", label = FALSE) + 
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  
	  labs(x = "UMAP 1", y = "UMAP 2") +  
	  theme(
		    plot.title = element_blank(),  
		    axis.title = element_text(size = 18),  
		    axis.text = element_text(size = 16),  
		    legend.text = element_text(size = 16),  
		    legend.title = element_text(size = 16),
		    legend.position = c(0.95, 0.05),        # bottom-right corner
		    legend.justification = c("right", "bottom")  # anchor to that corner
	  )
# Overlay cluster labels at computed centroids using `geom_text_repel()`
dim_plot <- dim_plot + 
	  geom_text_repel(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters), 
				  													size = 4, color = "black")
# Save the final plot
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_-C4C20_6x6.png")), 
							plot = dim_plot, width = 6, height = 6, dpi = 300)





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
ggsave(file.path(wd.de.plots, paste0("Cycle_pseudotime_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_C4C20_6.5x6.png")), plot = monocle3, width = 6.5, height = 6, dpi = 300)





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
ggsave(file.path(wd.de.plots, paste0("Cycle_pseudotime_new-color_-C0_-C4C20_6x6.png")),	plot = ordered, dpi = 300, width = 6, height = 6)

save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))






# -----------------------------------------------------------------------------
# Top sample per cluster
# Last Modified: 03/03/25
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))

# Identify columns to remove
columns_to_remove <- grep("^pANN_0.25|^DF.classifications_0.25|^integrated_snn_res", colnames(so.integrated@meta.data), value = TRUE)

# Remove the identified columns
so.integrated@meta.data <- so.integrated@meta.data[, !colnames(so.integrated@meta.data) %in% columns_to_remove]
print(colnames(so.integrated@meta.data))

# Add pseudotime to Seurat object's metadata
pseudotime_values <- monocle3::pseudotime(cds)
pseudotime_values <- pseudotime_values[colnames(so.integrated)]
so.integrated$pseudotime <- pseudotime_values

# Add cell types
so.integrated@meta.data$cell_type <- Idents(so.integrated)

# Fix incorrect patient ID in meta.data
so.integrated@meta.data$orig.ident <- ifelse(
   so.integrated@meta.data$orig.ident == "SeuratProject",
   "VL00297",
   so.integrated@meta.data$orig.ident
)

save(cds, so.integrated, phase_colors, plot_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_clean.RData")))
