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

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated.RData")))

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "PMC", "Leydig", "Endothelial and Macrophage", "4", "20")
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
ggsave(file.path(wd.de.plots, paste0("DimPlot_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6*6_SSC.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

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
so.integrated <- ScaleData(so.integrated, verbose = FALSE)

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
      legend.position = c(0.95, 0.95),       # top-right corner
      legend.justification = c("right", "top")  # anchor to that corner
   )
ggsave(file.path(wd.de.plots, paste0("pseudotime_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6.png")), plot = monocle3, dpi = 300, width = 6, height = 6)

# -----------------------------------------------------------------------------
# Re-root on the end cells (then reverse later)
# -----------------------------------------------------------------------------
# choose your end cells (example: a terminal cluster)
end_cells <- WhichCells(so.integrated, idents = "Late spermatid")  # Replace "early_cluster_name" with your specific cluster
# Ensure these root cells are present in `cds`
end_cells_in_cds <- end_cells[end_cells %in% colnames(cds)]

# re-root (this will make them START at 0)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = end_cells_in_cds)
# reverse to make them the END (max pseudotime)
pt <- monocle3::pseudotime(cds)
pt_rev <- max(pt, na.rm = TRUE) - pt
# Cannot overwrite the default pseudotime
colData(cds)$pseudotime_rev <- pt_rev

# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "pseudotime_rev",
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
      legend.position = c(0.95, 0.95),       # top-right corner
      legend.justification = c("right", "top")  # anchor to that corner
   )
ggsave(file.path(wd.de.plots, paste0("pseudotime_rev_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6.png")), plot = monocle3, dpi = 300, width = 6, height = 6)

# -----------------------------------------------------------------------------
# Methods: ggbeeswarm
# Last Modified: 20/01/25; 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#tscan
library(ggbeeswarm)
library(colorRamps)
library(viridisLite)

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- colData(cds)$pseudotime_rev  #monocle3::pseudotime(cds)
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
ggsave(file.path(wd.de.plots, paste0("pseudotime_rev_ggbeeswarm_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0.png")), plot = ordered, dpi = 300, width = 7, height = 6)

## After running cell cycle analysis
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
ggsave(file.path(wd.de.plots, paste0("pseudotime_rev_ggbeeswarm_cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0.png")), plot = ordered, dpi = 300, width = 7, height = 6)

# -----------------------------------------------------------------------------
# Cell cycle analysis
# Last Modified: 20/10/25; 10/02/25
# -----------------------------------------------------------------------------
DefaultAssay(so.integrated) <- "RNA" # Very important!

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
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G0/G1", "S", "G2/M")) +  # Apply custom colors & legend order
	  labs(x = "UMAP 1", y = "UMAP 2", color = "Phase") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16),  # Increase legend title size
		    legend.position = c(0.95, 0.95),        # bottom-right corner
		    legend.justification = c("right", "top")  # anchor to that corner
	  )
# **Extract cluster positions separately**
cluster_positions <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap", label = TRUE, repel = TRUE)
# **Overlay cluster labels on the same plot**
dim_plot <- LabelClusters(plot = dim_plot, id = "Phase", text = cluster_positions$data$Phase, repel = TRUE, size = 6, color = "black")
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6_RNA_00.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

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
umap_data$Phase <- so.integrated$Phase  # Add this line!
umap_data$seurat_clusters <- as.factor(so.integrated$seurat_clusters)  # Ensure clusters are factors

# Compute cluster centroids
cluster_centers <- umap_data %>%
	  group_by(seurat_clusters) %>%
	  summarise(UMAP_1 = median(umap_1), UMAP_2 = median(umap_2))  # Use median to find label position

# Create UMAP plot with manual ggplot (same style as your loop code)
dim_plot <- ggplot() +
 # Plot all cells with their respective phase colors
 geom_point(data = umap_data, 
            aes(x = umap_1, y = umap_2, color = Phase), 
            size = 0.5, alpha = 0.8) +
 # Add cluster labels
 geom_text_repel(data = cluster_centers, 
                 aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters), 
                 size = 5, color = "black") +
 # Manual color scale
 scale_color_manual(values = phase_colors, 
                    breaks = c("G1", "S", "G2M"), 
                    labels = c("G0/G1", "S", "G2/M"),
                    name = "Phase") +
 labs(x = "UMAP 1", y = "UMAP 2") +
 theme_classic() +
 theme(
  plot.title = element_blank(),  
  axis.title = element_text(size = 18),  
  axis.text = element_text(size = 16),  
  legend.text = element_text(size = 16),  
  legend.title = element_text(size = 16),
  legend.position = c(0.95, 0.95),
  legend.justification = c("right", "top")
 ) +
 guides(color = guide_legend(override.aes = list(size = 3)))

# Save the final plot
ggsave(file.path(wd.de.plots, paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6_RNA.png")), 
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
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G0/G1", "S", "G2/M")) + # Apply new color
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("Cycle_pseudotime_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_6x6_RNA.png")), plot = monocle3, width = 6, height = 6, dpi = 300)

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
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G0/G1", "S", "G2/M")) +  # Apply the same Phase colours as DimPlot
	  scale_y_discrete(labels = c("G1" = "G1", "S" = "S", "G2M" = "G2/M")) +  # Rename y-axis labels
	  theme_classic() +
	  xlab("Pseudotime") + 
	  ylab("Phase") +
	  theme(plot.title = element_blank(),  # Remove title
			  				axis.text = element_text(size = 16),  # Increase axis text size
					  		axis.title = element_text(size = 18), # Increase axis title size
						  	legend.position = "none")  # Remove legend

# Save the plot
ggsave(file.path(wd.de.plots, paste0("Cycle_pseudotime_new-color_-C0_6x6_RNA.png")),	plot = ordered, dpi = 300, width = 6, height = 6)

#save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_monocle3+phase.RData")))

# -----------------------------------------------------------------------------
# Top sample per cluster
# Last Modified: 03/03/25
# -----------------------------------------------------------------------------
# Identify columns to remove
#columns_to_remove <- grep("^pANN_0.25|^DF.classifications_0.25|^integrated_snn_res", colnames(so.integrated@meta.data), value = TRUE)
# Remove the identified columns
#so.integrated@meta.data <- so.integrated@meta.data[, !colnames(so.integrated@meta.data) %in% columns_to_remove]
#print(colnames(so.integrated@meta.data))

# Add pseudotime to Seurat object's metadata
pseudotime_values <- colData(cds)$pseudotime_rev
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

save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase_RNA.RData")))

# -----------------------------------------------------------------------------
# Plot phase-by-phase - Each phase highlighted separately in a loop
# Last Modified: 15/09/25
# -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggrepel)

# Ensure `seurat_clusters` exists
so.integrated$seurat_clusters <- Idents(so.integrated)

# Extract UMAP embeddings and phase information
umap_data <- as.data.frame(Embeddings(so.integrated, reduction = "umap"))
umap_data$Phase <- so.integrated$Phase
umap_data$seurat_clusters <- as.factor(so.integrated$seurat_clusters)

# Compute cluster centroids (same for all plots)
cluster_centers <- umap_data %>%
 group_by(seurat_clusters) %>%
 summarise(UMAP_1 = median(umap_1), UMAP_2 = median(umap_2))

# Define phases to plot and their colors
phases_to_plot <- c("G1", "S", "G2M")
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")
phase_labels <- c("G1" = "G0/G1", "S" = "S", "G2M" = "G2/M")

# Loop through each phase
for (phase in phases_to_plot) {
  # Separate target phase cells from others
  target_cells <- umap_data %>% filter(Phase == phase)
  other_cells <- umap_data %>% filter(Phase != phase)
 
  # Create plot with target phase highlighted
  dim_plot <- ggplot() +
  # Plot other cells first (background) in light gray
  geom_point(data = other_cells, 
             aes(x = umap_1, y = umap_2), 
             color = "#D3D3D3", size = 0.5, alpha = 0.8) +
  # Plot target phase cells on top (highlighted)
  geom_point(data = target_cells, 
             aes(x = umap_1, y = umap_2, color = phase_labels[[phase]]), 
             size = 0.5, alpha = 0.8) +
  # Add cluster labels
  geom_text_repel(data = cluster_centers, 
                  aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters), 
                  size = 4, color = "black") +
  # Manual color scale for target phase only
  scale_color_manual(values = setNames(phase_colors[[phase]], phase_labels[[phase]]),
                     name = "Phase") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  theme(
   plot.title = element_blank(),  
   axis.title = element_text(size = 18),  
   axis.text = element_text(size = 16),  
   legend.text = element_text(size = 16),  
   legend.title = element_text(size = 16),
   legend.position = c(0.95, 0.95),
   legend.justification = c("right", "top")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
 
 # Create filename for each phase
 filename <- paste0("Cycle_SSC_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_", 
                    phase, "_highlighted_6x6_RNA.png")
 
 # Save the plot
 ggsave(file.path(wd.de.plots, filename), 
        plot = dim_plot, width = 6, height = 6, dpi = 300)
 
 # Optional: print progress
 cat("Saved plot for", phase_labels[[phase]], "phase\n")
}

# -----------------------------------------------------------------------------
# Bar chart: Cell cycle phase proportions by cell type
# Cell types on y-axis, proportions on x-axis
# Phase order: G1, S, G2/M from left to right
# No percentage labels on bars
# Last Modified: 15/09/25
# -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)

# Extract cell type and phase information
cell_data <- data.frame(
 cell_type = so.integrated$seurat_clusters,
 phase = so.integrated$Phase
)

# Calculate proportions of each phase within each cell type
phase_proportions <- cell_data %>%
 group_by(cell_type, phase) %>%
 summarise(count = n(), .groups = "drop") %>%
 group_by(cell_type) %>%
 mutate(
  total = sum(count),
  proportion = count / total
 ) %>%
 ungroup()

# Reorder levels so "Stage 0" appears at the top of the y-axis
phase_proportions$cell_type <- factor(
 phase_proportions$cell_type,
 levels = rev(levels(factor(phase_proportions$cell_type)))
)

# Reorder phase levels to show G1, S, G2/M from left to right
# Need to specify in reverse order for stacked bars
phase_proportions$phase <- factor(
 phase_proportions$phase,
 levels = c("G2M", "S", "G1")
)

# Define phase colors
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Create stacked bar chart without percentage labels
bar_plot <- ggplot(phase_proportions, aes(x = proportion, y = cell_type, fill = phase)) +
 geom_bar(stat = "identity", width = 0.7) +
 scale_fill_manual(values = phase_colors, 
                   breaks = c("G1", "S", "G2M"), 
                   labels = c("G0/G1", "S", "G2/M")) +
 scale_x_continuous(labels = scales::percent_format(), 
                    expand = c(0, 0)) +
 labs(
  x = "Proportion", 
  y = "Cell type", 
  fill = "Phase"
 ) +
 theme_classic() +
 theme(
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  legend.position = "right",
  panel.grid.major.x = element_line(color = "grey90", size = 0.5),
  panel.grid.minor.x = element_line(color = "grey95", size = 0.25)
 )

# Save the plot
ggsave(file.path(wd.de.plots, "cell_cycle_proportions_by_cell_typ_6x14_RNA.png"), 
       plot = bar_plot, width = 14, height = 6, dpi = 300)

# Print the proportions table
print(phase_proportions)
