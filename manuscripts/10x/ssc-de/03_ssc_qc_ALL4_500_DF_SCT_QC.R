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
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
res <- 0.5
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_", nfeatures, "_-2_", 25, "_", 100, ".RData")))

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
library(dplyr)

# Extract metadata from so.integrated object
metadata <- so.integrated@meta.data

# Summarize the data by sample ID
summary_table <- metadata %>%
	  group_by(sample.id, age) %>%
	  summarise(
		    Reads = sum(nFeature_RNA),         # Total reads
		    Cells = n()                      # Total number of cells
	  ) %>%
	  ungroup()

# Convert to a data frame and print
summary_table <- as.data.frame(summary_table)
print(summary_table)

# Optionally save to a CSV file
write.csv(summary_table, "so_integrated_summary.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Define resolution
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, ".RData")))
resolution.range <- seq(from = 0.5, to = 1, by = 0.05)

for (r in 1:length(resolution.range)) {
	  res <- resolution.range[r]
	  resolution_name <- paste0("integrated_snn_res.", res)
	  #Idents(so.integrated) <- so.integrated[[resolution_name]][[1]]
	  so.integrated <- SetIdent(so.integrated, value = resolution_name)
	
	  # Extract the UMAP embeddings
	  #umap_embeddings <- Embeddings(so.integrated, reduction = "umap")
	  # Flip the UMAP 2 (second dimension)
	  #umap_embeddings[, 2] <- umap_embeddings[, 2] * -1
	  # Re-insert the modified embeddings back into the Seurat object
	  #so.integrated[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAP_", assay = DefaultAssay(so.integrated))
	  
	  ##
	  file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_res=", res, "_flip")
	  pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
	  DimPlot(so.integrated, label = TRUE)
	  dev.off()
}

# Find neighbors and clusters
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))

so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = 0.6)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
save(prin_comp, samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, ".RData")))

##
dims <- 25
k.weight <- 100
file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_-2_res=0.5_", dims, "_", k.weight)

pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
DimPlot(so.integrated, label = TRUE) + NoLegend()
dev.off()

pdf(file=file.path(wd.de.plots, paste0(file.name, "_SampleID.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, paste0(file.name, "_Age.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, paste0(file.name, "_2N.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="n2")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots,  paste0(file.name, "_Batch.pdf")))
tplot = DimPlot(so.integrated, reduction = "umap", group.by="batch")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and cluster identities
read_counts <- get_cluster_factor_levels(so.integrated, "nCount_RNA")
box_plot_reads <- get_box_plot(read_counts, "Read counts", "Number")

# Part 2: Generate the Box Plot for Global Expression Levels
# Extract the total expression level (sum of counts) and cluster identities
expression_data <- get_cluster_factor_levels(so.integrated, "nFeature_RNA")
box_plot_expression <- get_box_plot(expression_data, "Global expression levels", "Expression")

# Part 3: Generate the Box Plot for Mitochondrial Contents
# Extract the mitochondrial gene percentage  (percent.mt) and cluster identities
mt_data <- get_cluster_factor_levels(so.integrated, "percent.mt")
box_plot_mt <- get_box_plot(mt_data, "Mitochondrial contents", "Percentage")

# Part 4: Combine the Plots
combined_plot <- box_plot_reads / box_plot_expression / box_plot_mt

# Print the combined plot
print(combined_plot)

# -----------------------------------------------------------------------------
# Sampke ID
# -----------------------------------------------------------------------------
so.integrated$seurat_clusters <- Idents(so.integrated)

# Extract sample IDs and cluster identities
sample_cluster_data <- FetchData(so.integrated, vars = c("sample.id", "seurat_clusters"))

# Calculate the number of cells in each cluster
total_cells_per_cluster <- sample_cluster_data %>%
	  group_by(seurat_clusters) %>%
	  summarize(total_cells = n())

# Calculate the number of cells for each sample ID within each cluster
cells_per_sample_per_cluster <- sample_cluster_data %>%
	  group_by(seurat_clusters, sample.id) %>%
	  summarize(cells = n())

# Merge the data frames to calculate proportions
proportion_data <- merge(cells_per_sample_per_cluster, total_cells_per_cluster, by = "seurat_clusters")

# Calculate the proportion of each sample ID within each cluster
proportion_data <- proportion_data %>%
	  mutate(proportion = cells / total_cells) %>%
	  dplyr::select(seurat_clusters, sample.id, proportion)

# Rename columns for clarity
colnames(proportion_data) <- c("Cluster", "Sample_ID", "Proportion")

# Calculate the proportion of PD40746e_M2 in each cluster
pd40746e_m2_proportion <- proportion_data %>%
	  filter(Sample_ID == "PD40746e_M2") %>%
	  arrange(desc(Proportion)) %>%
	  pull(Cluster)

# Set factor levels for Cluster based on the proportion of PD40746e_M2
proportion_data$Cluster <- factor(proportion_data$Cluster, levels = pd40746e_m2_proportion)

# Create a vector of unique sample IDs
unique_samples <- unique(proportion_data$Sample_ID)

# Create a vector of random colors
set.seed(42)  # For reproducibility
random_colors <- grDevices::colors()[sample(1:length(grDevices::colors()), length(unique_samples))]

# Create a named vector of colors, assigning specific colors to the highlighted samples
color_vector <- setNames(random_colors, unique_samples)
color_vector["PD40746e_M2"] <- "red"
color_vector["PD40746e_M1"] <- "pink"
color_vector["VL00297"] <- "yellow"

# Assign colors to the data
proportion_data$color <- color_vector[proportion_data$Sample_ID]

# Create a stacked bar chart with custom colors
bar_chart <- ggplot(proportion_data, aes(x = as.factor(Cluster), y = Proportion, fill = color)) +
	  geom_bar(stat = "identity", color = "black") +  # Add border for better visibility
  	scale_fill_identity() +
  	theme_minimal() +
  	labs(title = "Proportion of Sample IDs in Each Cluster",
		  				x = "Cluster",
		  				y = "Proportion") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the bar chart
print(bar_chart)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.25)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, ".RData")))

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=0.5_14x6.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("2", "5", "22", "1", "27", "21", "16", "19", "25", "18", "10", "9", "3", "6", "7", "14", "24", "30", "20", "17", "29", "13", "4", "8", "23", "31", "32", "12", "11", "26", "28", "15", "0")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
   scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=0.5_ordered_14x7.pdf")), width = 14, height = 7)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
	  )
ggsave(file.path(wd.de.plots, paste0("DimPlot_no-title_6*6_new-color_16.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.4; Six SPG states)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.45_", nfeatures, ".RData")))

genes_of_interest <- c("TAF6", "ST3GAL4", "SH2B2", "MSL3", "PHGDH", "C19orf84", "LIN7B", "FSD1", "TSPAN33", "EGR4", "PIWIL4", "CELF4", "UTF1", "FGFR3", "A2M", "ENO3", "SERPINE2", "SRRT", "BAG6", "DND1", "PELP1", "NANOS2", "C1QBP", "NANOS3", "GFRA2", "GFRA1", "ID2", "ASB9", "L1TD1", "ID4", "MKI67", "PDPN", "KIT", "DMRT1", "DNMT1", "CALR", "SYCP3", "STRA8")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=", res, ".pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
#new_order <- c("8", "2", "21", "20", "11", "25", "19", "4", "10", "22", "15", "23", "17", "14", "3", "6", "9", "12", "13", "24", "7", "18", "0", "16", "5", "1")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=", res, "_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
write.table(table(Idents(so.integrated)), file=file.path(wd.de.data, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=", res, "_ordered_annotated.txt")), row.names=T, col.names=T, quote=F, sep='\t')

cluster_to_celltype <- c('2' = 'Stage 0', '5' = 'Stage 0A', '22' = 'Stage 0B', '1' = 'Stage 1',
																									'27' = 'Stage 2', '21' = 'Stage 2', '16' = 'Stage 3', #'19' = '',
																									'25' = 'Leptotene', '18' = 'Leptotene',
																									'10' = 'Zygotene', '9' = 'Zygotene', '3' = 'Zygotene',
																								 '7' = 'Pachytene', '14' = 'Pachytene',  #'6' = 'Pachytene',
																									'24' = 'Diplotene', '30' = 'Diplotene',
																									'20' = 'Meiotic division', '17' = 'Meiotic division', '29' = 'Meiotic division',
																									'13' = 'Early spermatid',
																									'4' = 'Late spermatid', '8' = 'Late spermatid', #'23' = '',
																									'31' = 'Macrophage',
																									'32' = 'Endothelial', 
																									'12' = 'PMC',
																									'11' = 'Fibrotic PMC',
																									'26' = 'Leydig',
																									'28' = 'Sertoli', '15' = 'Sertoli', '0' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=", res, "_ordered_annotated_19-6-23.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend()
ggsave(file.path(wd.de.plots, paste0("Di Persio_SPG_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_res=", res, "_ordered_annotated_NoLegend_8x8_19-6-23.png")), plot = dim_plot, width = 8, height = 8, dpi = 300)

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))

# -----------------------------------------------------------------------------
# Matt's postitively selected genes
# -----------------------------------------------------------------------------
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "SEMG1", "ARHGAP35", "CBL", "SSX1", "RBM12", "TRERF1", "FAT1", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Ensure age is numeric
so.integrated@meta.data$age <- as.numeric(so.integrated@meta.data$age)

# Extract age and cluster information from the metadata
age_cluster_data <- so.integrated@meta.data %>%
  	select(age, seurat_clusters)

# Calculate correlations for each cluster and store results
cluster_ids <- levels(so.integrated@active.ident)
correlation_results <- lapply(cluster_ids, calculate_spearman_mn7, age_cluster_data = age_cluster_data, so.integrated = so.integrated)

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, ".RData")))

cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)

# Cluster cells in Monocle 3 using UMAP embeddings
cds <- cluster_cells(cds, reduction_method = "UMAP")
# Perform trajectory analysis using precomputed UMAP
cds <- learn_graph(cds, use_partition = F)
# Order cells in pseudotime
cds <- order_cells(cds)

# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "cluster",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=1.5)
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_partition.png")), plot = monocle3, dpi = 300)

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
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
# **Extract cluster positions separately**
cluster_positions <- DimPlot(so.integrated, group.by = "seurat_clusters", reduction = "umap", label = TRUE, repel = TRUE)
# **Overlay cluster labels on the same plot**
dim_plot <- LabelClusters(plot = dim_plot, id = "seurat_clusters", text = cluster_positions$data$seurat_clusters, repel = TRUE, size = 6, color = "black")
ggsave(file.path(wd.de.plots, paste0("Cell-cycle_annotated_19-6-23-29_DimPlot_no-title_7.1*6_new-color_16.png")), plot = dim_plot, width = 7.1, height = 6, dpi = 300)

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
		    legend.title = element_text(size = 16)  
	  )
# Overlay cluster labels at computed centroids using `geom_text_repel()`
dim_plot <- dim_plot + 
	  geom_text_repel(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters), 
				  													size = 4, color = "black")
# Save the final plot
ggsave(file.path(wd.de.plots, paste0("Cell-cycle_annotated_UMAP_with_clusters_7.1x6.png")), 
							plot = dim_plot, width = 7.1, height = 6, dpi = 300)
