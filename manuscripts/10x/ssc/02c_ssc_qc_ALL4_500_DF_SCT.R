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
wd.de.data  <- file.path(wd.de, "data_ALL4_500")
wd.de.plots <- file.path(wd.de, "plots_ALL4_500")

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
res <- 0.5
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_20_100.RData")))

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
save(prin_comp, samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, ".RData")))

##
file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.6_NoLegend")

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
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, ".RData")))

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.6.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("8", "9", "10", "2", "7", "5", "4", "6", "11", "12", "13", "19", "16", "15", "17", "14", "18", "1", "21", "3", "23", "20", "22", "0")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
   scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.6_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

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

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=", res, ".pdf")), width = 14, height = 5)
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

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=", res, "_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
write.table(table(Idents(so.integrated)), file=file.path(wd.de.data, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=", res, "_ordered_annotated.txt")), row.names=T, col.names=T, quote=F, sep='\t')

cluster_to_celltype <- c('4' = 'Stage 0', '1' = 'Stage 0A', '19' = 'Stage 0B', '2' = 'Stage 1',
																									'16' = 'Stage 2', '6' = 'Stage 3',
																									'9' = 'Leptotene',
																									'7' = 'Zygotene', '5' = 'Zygotene', '8' = 'Zygotene',
																									'18' = 'Pachytene', '10' = 'Pachytene', '24' = 'Pachytene',# '15' = 'Pachytene',
																									'22' = 'Diplotene', 	
																									'21' = 'Meiotic division',
																									'11' = 'Early spermatid',
																									'13' = 'Late spermatid', '14' = 'Late spermatid',
																									'25' = 'Endothelial and Macrophage',
																									'3' = 'Fibrotic PMC',
																									'23' = 'Leydig',
																									'0' = 'Sertoli', '26' = 'Sertoli')#, '17' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=", res, "_ordered_annotated_15_17.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend()
ggsave(file.path(wd.de.plots, paste0("Di Persio_SPG_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_res=", res, "_ordered_annotated_NoLegend.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_annotated.RData")))

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
