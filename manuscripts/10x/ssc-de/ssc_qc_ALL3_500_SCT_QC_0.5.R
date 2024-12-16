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

wd.de.data0 <- "/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_ALL3_500"
load(file=file.path(wd.de.data0, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data0, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, ".RData")))

# Identify the clusters you want to keep
clusters_to_remove <- c(12, 17)
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)
7
ssc0 <- as.data.frame(table(so.integrated@meta.data$sample.id))
write.table(ssc0, file=file.path(wd.de.data, "ssc_samples_-12-15-17_>0.txt"), row.names=T, col.names=T, quote=F, sep='\t')

# Access the metadata
metadata <- so.integrated@meta.data
# Count the number of cells per sample
cell_counts <- metadata %>%
  	group_by(sample.id) %>%
	  summarise(cell_count = n()) %>%
	  ungroup()

# Filter samples with more than 100 cells
samples_to_keep <- cell_counts %>%
	  filter(cell_count > 100) %>%
	  pull(sample.id)

# Subset the Seurat object to keep only the samples with more than 100 cells
so.integrated <- subset(so.integrated, cells = rownames(metadata[metadata$sample.id %in% samples_to_keep, ]))

# Check the result
ssc <- as.data.frame(table(so.integrated@meta.data$sample.id))
write.table(ssc, file=file.path(wd.de.data, "ssc_samples_-12-15-17_>100.txt"), row.names=T, col.names=T, quote=F, sep='\t')

save(ssc0, ssc, so.integrated, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_5000_-12-17_>100.RData"))

# -----------------------------------------------------------------------------
# Define resolution
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, "_-12-15-17.RData")))
resolution.range <- seq(from = 0.4, to = 0.6, by = 0.05)

for (r in 1:length(resolution.range)) {
	  res <- resolution.range[r]
	  resolution_name <- paste0("integrated_snn_res.", res)
	  #Idents(so.integrated) <- so.integrated[[resolution_name]][[1]]
	  so.integrated <- SetIdent(so.integrated, value = resolution_name)
	
	  # Extract the UMAP embeddings
	  umap_embeddings <- Embeddings(so.integrated, reduction = "umap")
	  # Flip the UMAP 2 (second dimension)
	  umap_embeddings[, 2] <- umap_embeddings[, 2] * -1
	  # Re-insert the modified embeddings back into the Seurat object
	  so.integrated[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAP_", assay = DefaultAssay(so.integrated))
	
	  ##
	  file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_-12-15-17_res=", res, "")
	  pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
	  DimPlot(so.integrated, label = TRUE)
	  dev.off()
}

# Find neighbors and clusters
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, "_-12-15-17.RData")))

so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = 0.5)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17.RData")))

# Extract the UMAP embeddings
umap_embeddings <- Embeddings(so.integrated, reduction = "umap")
# Flip the UMAP 2 (second dimension)
umap_embeddings[, 1] <- umap_embeddings[, 1] * -1
umap_embeddings[, 2] <- umap_embeddings[, 2] * -1
# Re-insert the modified embeddings back into the Seurat object
so.integrated[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAP_", assay = DefaultAssay(so.integrated))

##
file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_-12-15-17_resolution=0.5")

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
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, "_-12-15-17.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17.RData")))

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
# Sample ID
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
	  select(seurat_clusters, sample.id, proportion)

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
color_vector["AMSBIO"] <- "yellow"

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
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17.RData")))

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_-12-15-17_resolution=0.5.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("0", "21", "18", "4", "20", "13", "9", "11", "15", "17", "19", "14", "7", "2", "3", "8", "12", "10", "5", "16", "6", "1")
new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_-12-15-17_resolution=0.5_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
#table(Idents(so.integrated))
#cluster_to_celltype <- c('2' = 'Undiff. SPG',
#	 																							'10' = 'Diff. SPG',
#																									'9' = 'Leptotene',
#																									'11' = 'Zygotene', '5' = 'Zygotene',
#																									'12' = 'Pachytene', '8' = 'Pachytene',
#																									'13' = 'Diplotene',
#																									'6' = 'Meiotic division',
#																									'4' = 'Late spermatid',
#																									'16' = 'Endothelial & Macrophage',
#																									'7' = 'PMC',
#																									'15' = 'Leydig',
#																									'1' = 'Sertoli',	'17' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
#Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)
#
#Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)
#
#dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
#	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
#	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
#scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting
#
#pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.6_ordered_annotated.pdf")), width = 12, height = 5)
#print(dot_plot)
#dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SCT_", nfeatures, "_UMAP_dims=14_resolution=0.6_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.4; Six SPG states)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.45_", nfeatures, ".RData")))

genes_of_interest <- c("TAF6", "ST3GAL4", "SH2B2", "MSL3", "PHGDH", "C19orf84", "LIN7B", "FSD1", "TSPAN33", "EGR4", "PIWIL4", "CELF4", "UTF1", "FGFR3", "A2M", "ENO3", "SERPINE2", "SRRT", "BAG6", "DND1", "PELP1", "NANOS2", "C1QBP", "NANOS3", "GFRA2", "GFRA1", "ID2", "ASB9", "L1TD1", "ID4", "MKI67", "PDPN", "KIT", "DMRT1", "DNMT1", "CALR", "SYCP3", "STRA8")

#DefaultAssay(so.integrated) <- "SCT"
#dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
#	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
#	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
#
#pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=", res, ".pdf")), width = 14, height = 5)
#print(dot_plot)
#dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
#new_order <- c("8", "2", "21", "20", "11", "25", "19", "4", "10", "22", "15", "23", "17", "14", "3", "6", "9", "12", "13", "24", "7", "18", "0", "16", "5", "1")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_-12-15-17_resolution=", res, "_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
write.table(table(Idents(so.integrated)), file=file.path(wd.de.data, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_-12-15-17_resolution=", res, "_ordered_annotated.txt")), row.names=T, col.names=T, quote=F, sep='\t')

cluster_to_celltype <- c('1' = 'Stage 0', '6' = 'Stage 0A', '16' = 'Stage 0B', '5' = 'Stage 1',
																									'10' = 'Stage 2', '12' = 'Stage 3',
																									'8' = 'Leptotene',
																									'3' = 'Zygotene', '2' = 'Zygotene',
																									'7' = 'Pachytene', '14' = 'Pachytene',# '15' = 'Pachytene',
																									'19' = 'Diplotene', 	
																									'17' = 'Meiotic division', '15' = 'Meiotic division',
																									'11' = 'Early spermatid',
																									'9' = 'Late spermatid', '13' = 'Late spermatid',
																									'20' = 'Endothelial and Macrophage',
																									'4' = 'Fibrotic PMC',
																									'18' = 'Leydig',
																									'0' = 'Sertoli', '21' = 'Sertoli')#, '17' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_-12-15-17_resolution=", res, "_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend()
ggsave(file.path(wd.de.plots, paste0("Di Persio_SPG_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_-12-15-17_res=", res, "_ordered_annotated_NoLegend_10x10_PMC.png")), plot = dim_plot, width = 10, height = 10, dpi = 300)

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
DefaultAssay(so.integrated) <- "RNA"

so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)
de_gene_list <- unique(markers$gene)

save(so.integrated, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_markers_-12-15-17.RData")))
save(markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_markers_markers_-12-15-17.RData")))

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()
#missing_genes <- setdiff(top10$gene, rownames(so.integrated))

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.6_group.by_-12-15-17.png")),
							plot = heatmap, width = 15, height = 17, dpi = 300)

# -----------------------------------------------------------------------------
# Jaccard index for overlap of gene sets
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_markers_markers.RData")))

jaccard_matrix <- getJaccardIndex(markers)
# Visualize the Jaccard index matrix
pheatmap_plot <- pheatmap(jaccard_matrix, main = "Jaccard Index Between Clusters")

# Save the plot as a PNG file
png(file = file.path(wd.de.plots, "heatmap_jaccard_markers_SCT_res=0.6_-12-15-17.png"), width = 7, height = 7, units = "in", res = 300)
print(pheatmap_plot)
dev.off()

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
# > dim(so.integrated)
# [1]  5000 56396

# Extract expression data from Seurat object
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
expression_matrix <- GetAssayData(so.integrated[["RNA"]], layer = "data")
expression_matrix <- as.matrix(expression_matrix)

# Set the number of clusters (k)
k <- 5  # You can adjust this based on your needs
# Perform K-means clustering on genes (rows)
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(t(expression_matrix), centers = k)
# Assign cluster labels to the genes
gene_clusters <- kmeans_result$cluster

save(gene_clusters, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering0.RData")))

# Order cells by their assigned cell types
cell_metadata <- so.integrated@meta.data
ordered_cells <- cell_metadata[order(cell_metadata$cell_type), ]
# Reorder the expression matrix by cell types
ordered_expression_matrix <- expression_matrix[, rownames(ordered_cells)]

# Load necessary libraries for heatmap
library(ComplexHeatmap)
library(circlize)
# Define color palette
col_fun <- colorRamp2(c(-2, 0, 2), c(blue, "white", red))
# Create heatmap with clustered genes and ordered cells
Heatmap(ordered_expression_matrix, 
								name = "Expression", 
								show_row_names = FALSE, 
								show_column_names = FALSE, 
								cluster_rows = TRUE, 
								cluster_columns = FALSE,  # Keep columns ordered by cell type
								row_split = gene_clusters, # Split rows by k-means clusters
								col = col_fun, 
								heatmap_legend_param = list(title = "Z-score"))

# Determine optimal number of clusters (Elbow Method)
wss <- sapply(1:15, function(k){
	  kmeans(t(expression_matrix), centers = k, nstart = 10)$tot.withinss
})

plot(1:15, wss, type = "b", pch = 19, frame = FALSE, 
					xlab = "Number of clusters K", 
					ylab = "Total within-clusters sum of squares")

save(kmeans_result, ordered_expression_matrix, wss, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_clustering.RData")))

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
median_age <- median(so.integrated$age)  # Define the threshold for age
so.integrated$age_group <- ifelse(so.integrated$age <= median_age, "young", "old")

# List of unique cell types
cell_types <- unique(so.integrated$cell_type)  # Adjust if your metadata column for cell type is named differently
# Initialize a list to store DE results
de_results_list <- list()

# Loop over each cell type
for (ct in cell_types) {
	  # Subset Seurat object to include only the cells of the current cell type
	  cells_of_ct <- WhichCells(so.integrated, expression = cell_type == ct)
	  so_ct <- subset(so.integrated, cells = cells_of_ct)
	
	  # Perform differential expression analysis between "young" and "old" using MAST
	  de_results <- FindMarkers(
		    so_ct,
		    ident.1 = "old",
		    ident.2 = "young",
		    group.by = "age_group",  # Use age group for ident.1 and ident.2
		    test.use = "MAST"
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
  		  ggtitle("Expression vs. Age") +
  		  theme(legend.position = "none",  # Hide the legend for custom color  theme(legend.position = "none",            # Hide the legend
         plot.title = element_text(hjust = 0.5))  # Center the title
  	
  	# Print the volcano plot
  	pdf(file = file.path(wd.de.plots, paste0("MAST_volcano_", ct, "_RNA.pdf")), width = 5, height = 5)
  	print(volcano_plot)
  	dev.off()
}

save(de_results_list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_MAST_RNA.RData")))

# -----------------------------------------------------------------------------
# Methods: Age group
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

age_info <- so.integrated@meta.data %>%
   dplyr::mutate(sample_id_clean = sub("_.*", "", sample.id)) %>%
   dplyr::select(sample_id_clean, age) %>%
   dplyr::distinct()

# Define age groups
age_info$age_group <- cut(
	  age_info$age,
	  breaks = c(20, 30, 40, 50, 60, 70, 80),
	  labels = c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"),
	  right = FALSE
)

# Assign median age for each group
age_info$age_median = dplyr::case_when(
	  age_info$age_group == "20-30" ~ 25,
	  age_info$age_group == "30-40" ~ 35,
	  age_info$age_group == "40-50" ~ 45,
	  age_info$age_group == "50-60" ~ 55,
	  age_info$age_group == "60-70" ~ 65,
	  age_info$age_group == "70-80" ~ 75,
	  TRUE ~ NA_real_
)

# Count the number of samples in each age group
age_group_counts <- as.data.frame(table(age_info$age_group))
colnames(age_group_counts) <- c("Age_Group", "Count")

# Plot the bar chart
ggplot(age_group_counts, aes(x = Age_Group, y = Count)) +
	  geom_bar(stat = "identity", fill = "white", color = "black") +
	  labs(title = "Number of individuals", x = "Age group", y = "Number") +
	  theme_minimal() +
	  theme(
		    axis.text.x = element_text(angle = 45, hjust = 1),
		    plot.title = element_text(hjust = 0.5)  # Center-align the title
	  )  +
	  scale_y_continuous(breaks = c(0, 1, 2))

# -----------------------------------------------------------------------------
# Methods: Age group (order by sample ID)
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
# Create a summary table with the data
summary_table <- so.integrated@meta.data %>%
	  group_by(sample.id) %>%
	  summarize(
		    Age = unique(age),  # Assuming `age` is a column with age data
		    Genes = sum(nFeature_RNA),  # `nFeature_RNA` is typically gene count per cell
		    Cells = n()  # Number of cells per sample
	  ) %>%
	  arrange(sample.id)

# Calculate total cells
total_cells <- sum(summary_table$Cells)

# Add the total row
summary_table <- bind_rows(summary_table, data.frame(
	  sample.id = "Total", Age = "", Genes = "", Cells = total_cells
))
writeTable(summary_table, file.path(wd.de.data, "summary_table.txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Methods: Age group (order by age)
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
# Create a summary table with the data
summary_table <- so.integrated@meta.data %>%
	  group_by(sample.id) %>%
	  summarize(
		    Age = unique(age),  # Assuming `age` is a column with age data
		    Genes = sum(nFeature_RNA),  # `nFeature_RNA` is typically gene count per cell
		    Cells = n()  # Number of cells per sample
	  ) %>%
	  arrange(Age)

# Calculate total cells
total_cells <- sum(summary_table$Cells)

# Add the total row
writeTable(summary_table, file.path(wd.de.data, "summary_table_age.txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Methods: Age group and cell types
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# Assuming `so.integrated` is your Seurat object
cell_type_colors <- Idents(so.integrated)  # Extract the colors used for each cell type

# Check the colors assigned to each identity
DimPlot(so.integrated)  # This shows the DimPlot with the current colors

# Extract the colors for each identity (cell type) programmatically
plot_colors <- scales::hue_pal()(length(unique(cell_type_colors)))
names(plot_colors) <- levels(cell_type_colors)

# Display the color assignments for reference
plot_colors

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: 18/11/24
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)

so.integrated@meta.data$cell_type <- Idents(so.integrated)

# Ensure age groups are defined in your metadata
so.integrated@meta.data <- so.integrated@meta.data %>%
	  mutate(age_group = cut(
		    age,
		    breaks = c(20, 30, 40, 50, 60, 70, 80),
		    labels = c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"),
		    right = FALSE
	   ))

# Summarize counts for each cell type within each age group
cell_type_counts <- so.integrated@meta.data %>%
	  group_by(age_group, cell_type) %>%
	  summarize(
		    cell_count = n(),  # Count of cells for each cell type
		    .groups = "drop"
	  )

# Calculate total cells for each age group
age_group_totals <- cell_type_counts %>%
  	group_by(age_group) %>%
	  summarize(
		    total_cells = sum(cell_count),
		    .groups = "drop"
  	)

# Merge totals and calculate proportions
proportions <- cell_type_counts %>%
	  left_join(age_group_totals, by = "age_group") %>%
  	mutate(proportion = cell_count / total_cells)

# Define the median age for each age group
age_group_to_median <- c("20-30" = 25, "30-40" = 35, "40-50" = 45, 
																									"50-60" = 55, "60-70" = 65, "70-80" = 75)
# Add the `age_median` column
proportions$age_median <- as.numeric(sapply(proportions$age_group, function(x) age_group_to_median[x]))

# Display the resulting table
print(proportions)
# Optionally, save the table to a file
write.table(proportions, "cell_type_proportions_by_age_group.txt", row.names = FALSE)

# Plot with age groups on x-axis and cell type proportions
plot <- ggplot(proportions, aes(x = age_group, y = proportion, fill = cell_type)) +
	  geom_bar(stat = "identity", position = "fill", color = "black") +
	  scale_fill_manual(values = plot_colors, name = "Cell Type") +
	  labs(title = "", x = "Age group", y = "Proportion of cell type") +
	  theme_minimal(base_size = 14) +
	  theme(
		    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
		    axis.text.y = element_text(size = 11),
		    axis.title.x = element_text(size = 13),
		    axis.title.y = element_text(size = 13),
		    plot.title = element_text(hjust = 0.5),
		    legend.text = element_text(size = 10),   # Increase legend text size
		    legend.title = element_text(size = 12)  # Increase legend title size
	  )

# Print the plot
png(file = file.path(wd.de.plots, "barplot_age_group_cell_type_DimPlot.png"), width = 6, height = 7, , units = "in", res = 300)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Methods: Correlation 1, 2, 3, 4, 5, 6 OR 25, 35, 45, 55, 65, 75
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
# Step 2: Perform Spearman's correlation and get p-value for each cell type
correlation_results <- proportions %>%
	  dplyr::group_by(cell_type) %>%
	  dplyr::summarize(
		    p_value = cor.test(as.numeric(as.factor(age_group)), proportion, method = "spearman")$p.value,
				  spearman_correlation = cor(as.numeric(as.factor(age_group)), proportion, method = "spearman"),
    		.groups = "drop"
  	)

# Step 2: Perform Spearman's correlation and get p-value for each cell type
correlation_results <- proportions %>%
	  dplyr::group_by(cell_type) %>%
	  dplyr::summarize(
		    p_value = cor.test(age_median, proportion, method = "spearman")$p.value,
		    spearman_correlation = cor(age_median, proportion, method = "spearman"),
		    .groups = "drop"
	  )

# Step 3: Display the results, sorted by correlation
correlation_results <- correlation_results %>%
	  dplyr::arrange(desc(spearman_correlation))

# Print the results
print(correlation_results)
writeTable(correlation_results, file.path(wd.de.data, "summary_table_correlation_results_age_median.txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Methods: Age group and cell types
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
library(ggplot2)

# Assuming correlation_results from the previous step contains the correlation and p-value for each cell type
# Add a column to indicate significance level based on p-value (optional)
correlation_results <- correlation_results %>%
	  dplyr::mutate(significance = dplyr::case_when(
		    p_value < 0.001 ~ "***",
		    p_value < 0.01 ~ "**",
		    p_value < 0.05 ~ "*",
		    TRUE ~ ""
	  ))

# Plot the bar plot
plot <- ggplot(correlation_results, aes(x = reorder(cell_type, spearman_correlation), y = spearman_correlation, fill = spearman_correlation)) +
	  geom_bar(stat = "identity") +
	  coord_flip() +
	  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, name = "Correlation") +
	  geom_text(aes(label = significance, hjust = ifelse(spearman_correlation < 0, 1.5, -0.5)),  # Adjust hjust for more space
			  								vjust = 0.5, color = "black", size = 6) + 
	  labs(
	    	title = "",
	    	x = "",
		    y = "Spearmans' rho"
	  ) +
	  theme_minimal() +
	  theme(
		    axis.text.x = element_text(size = 11),   # Adjust y-axis text size
		    axis.text.y = element_text(size = 11),
		    axis.title.x = element_text(size = 12, margin = margin(t = 10)),  # Adjust x-axis title size and position
		    plot.title = element_text(hjust = 0.5),
		    legend.text = element_text(size = 11),   # Increase legend text size
		    legend.title = element_text(size = 12)  # Increase legend title size
	  ) +
	  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))  # Add margins to both sides of the y-axis

# Print the plot
pdf(file = file.path(wd.de.plots, "correlation_results_age_median.pdf"), width = 6, height = 6)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Positively seleced genes
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_joined.RData")))
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "SEMG1", "ARHGAP35", "CBL", "SSX1", "RBM12", "TRERF1", "FAT1", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
# Step 2: Extract the scaled data for the genes in mn7
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]

# Step 3: Calculate the average expression for each gene in each cluster
average_expression_per_cluster <- sapply(levels(so.integrated), function(cluster) {
	  # Ensure no NAs by filling missing values
	  rowMeans(scaled_data[, Idents(so.integrated) == cluster, drop = FALSE], na.rm = TRUE)
})

# Step 4: Sort the genes sequentially by expression across clusters
ordered_genes <- rownames(scaled_data)[do.call(order, as.data.frame(-average_expression_per_cluster))]
scaled_data <- scaled_data[ordered_genes,]

# Step 4: Order the genes based on their expression levels from cluster 1 to the last cluster
# Sort each cluster's expression and then combine orders without dropping any genes
#gene_order <- order(rowMeans(average_expression_per_cluster, na.rm = TRUE), decreasing = TRUE)
#ordered_genes <- rownames(scaled_data)[gene_order]
#scaled_data <- scaled_data[gene_order,]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + theme(axis.text.y = element_text(size = 12))
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_mn7.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# -----------------------------------------------------------------------------
# Postitively selected genes
# -----------------------------------------------------------------------------
library(dplyr)

mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "SEMG1", "ARHGAP35", "CBL", "SSX1", "RBM12", "TRERF1", "FAT1", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Extract age and cluster information from the metadata
so.integrated@meta.data$cell_types <- so.integrated@active.ident
age_cluster_data <- so.integrated@meta.data %>%
	  dplyr::select(age, cell_types)

# Calculate correlations for each cluster and store results
#significant_genes <- calculate_DESeq2_pseudobulk("Stage 0", age_cluster_data, so.integrated, mn7)
cell_types <- levels(so.integrated@active.ident)
correlation_results <- lapply(cell_types, calculate_DESeq2_pseudobulk, age_cluster_data, so.integrated, mn7)









# -----------------------------------------------------------------------------
# Doublet Finder
# -----------------------------------------------------------------------------
library(DoubletFinder)

if (!"counts" %in% names(so.integrated[["RNA"]])) {
	  so.integrated[["RNA"]] <- CreateAssayObject(counts = GetAssayData(so.integrated, slot = "counts"))
}












# -----------------------------------------------------------------------------
# Methods: DE using Monocle 3
# Last Modified: 17/10/24
# -----------------------------------------------------------------------------
# Set the default assay to "integrated"
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Cluster cells in Monocle 3 using UMAP embeddings
#cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)
# Extract expression data from Seurat object
expression_matrix <- GetAssayData(so.integrated[["RNA"]], slot = "data")

# Create Monocle3 CDS manually
cds <- new_cell_data_set(
	  expression_data = expression_matrix,
	  cell_metadata = as.data.frame(so.integrated@meta.data),
	  gene_metadata = data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
)
cds$idents <- Idents(so.integrated)

# Optionally add UMAP embeddings if available
reducedDim(cds, "UMAP") <- Embeddings(so.integrated, "umap")
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_RNA_cds.RData")))

# Get unique cell types
cell_types <- unique(cds$idents)
# Loop through each cell type and perform differential expression analysis
for (ct in cell_types) {
	  # Subset CDS for each cell type
	  cds_subset <- cds[, cds$idents == ct]
	
	  # Fit model using age as a covariate
	  diff_test_res <- fit_models(cds_subset, model_formula_str = "~ age + batch", cores = 4)
	
	  # Extract and save results for each cell type
	  coefficient_table <- coefficient_table(diff_test_res)
	  coefficient_table_age <- subset(coefficient_table, term == "age")
	  significant_genes <- subset(coefficient_table_age, q_value < 1E-3)
	
	  # Ensure coefficient_table_age is a tibble
	  coefficient_table_age <- as.data.frame(coefficient_table_age)
	  rownames(coefficient_table_age) <- coefficient_table_age$gene_id
	  
	  # Add a binary age group: "young" or "old"
	  median_age <- median(cds$age)  # Adjust to your specific threshold
	  cds$age_group <- ifelse(cds$age <= median_age, "young", "old")
	  
	  # Extract the expression matrix from the CDS object
	  expression_matrix <- exprs(cds)
	  
	  # Split the cells into two groups based on age
	  young_cells <- colnames(cds)[cds$age_group == "young"]
	  old_cells <- colnames(cds)[cds$age_group == "old"]
	  
	  # Calculate average expression for each gene in young and old groups
	  avg_expr_young <- rowMeans(expression_matrix[, young_cells])
	  avg_expr_old <- rowMeans(expression_matrix[, old_cells])
	  
	  # Calculate fold change (log2 fold change between old and young)
	  log2_fold_change <- log2(avg_expr_old / avg_expr_young)
	  
	  # Create a data frame with gene names and log2 fold change
	  fold_change_data <- data.frame(
	  	  gene_short_name = rownames(expression_matrix),
	  	  log2_fold_change = log2_fold_change
	  )
	  
	  # Combine fold change data with p-values from Monocle
	  volcano_data <- merge(fold_change_data, coefficient_table_age[, c("gene_short_name", "p_value", "q_value")], by = "gene_short_name")
	  
	  # Calculate -log10(p-value) for y-axis
	  volcano_data <- volcano_data %>%
	  	  mutate(log10_pvalue = -log10(p_value))
	  
	  ggplot(volcano_data, aes(x = log2_fold_change, y = log10_pvalue)) +
	  	  geom_point(aes(color = p_value < 1E-3), alpha = 0.7) +  # Highlight significant genes
	  	  theme_minimal() +
	  	  xlab("Log2 Fold Change") +
	  	  ylab("-log10(p-value)") +
	  	ggtitle("Volcano Plot: Differential Expression vs. Age")
}














# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression data in each cluster
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Ensure age is numeric
#so.integrated@meta.data$age <- as.numeric(so.integrated@meta.data$age)

# Extract age and cluster information from the metadata
age_cluster_data <- so.integrated@meta.data %>%
	  select(age, seurat_clusters)

# Calculate correlations for each cluster and store results
cluster_ids <- levels(so.integrated@active.ident)
correlation_results <- lapply(cluster_ids, calculate_spearman, age_cluster_data = age_cluster_data, so.integrated = so.integrated)

# Combine results from all clusters
all_significant_genes <- bind_rows(correlation_results, .id = "cluster")

# View the results
print(all_significant_genes)
# Optionally, save the results to a CSV file
write.table(all_significant_genes, file = "significant_genes_age_correlation_cluster_-12-15-17.txt", row.names = FALSE)

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_joined.RData")))

# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression data in each cluster
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

mt <- c("MT-ND5", "MT-CO1", "MT-ND4L", "MT-ND2", "MT-ND4", "MT-CYB", "MT-CO3", "MT-ATP6", "MT-ND1", "AC090403.1", "MT-CO2", "MT-ND1")
mp29 <- c("AC090403.1", "AC026116.1", "TPRG1", "AC009271.1", "LINC02074", "DANT2", "WSB1", "GNGT1", "TENT2")
#mp29 <- c("HMGB1", "GCNT2", "KDM2A", "NKAIN2", "TPTE", "SDK1", "CACNA2D3", "FBXL20", "AGAP1", "RUNX1", "LINC01206", "LINC01508", "ANKRD30A", "ANKRD30B", "ANKRD36C", "AP000550.1")
#mp29 <- intersect(mp29, rownames(so.integrated))

# Extract the expression data for these genes
expression_data <- GetAssayData(so.integrated, assay = "RNA", layer = "data")[mt, ]
# Convert to a data frame and transpose for easier manipulation
expression_data_df <- as.data.frame(t(expression_data))
# Add batch information from metadata
expression_data_df$batch <- so.integrated@meta.data$batch
# Reshape the data for ggplot2
expression_data_long <- reshape2::melt(expression_data_df, id.vars = "batch", variable.name = "gene", value.name = "expression")

# Plot using ggplot2
boxplot_expression <- ggplot(expression_data_long, aes(x = factor(batch), y = expression, fill = gene)) +
	  geom_boxplot(fill = blue) +
	  facet_wrap(~ gene, scales = "free_y") +  # Create a separate plot for each gene
	  labs(title = "", x = "Batch", y = "Expression") +
	  theme_minimal() +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend() 

# Print the plot
print(boxplot_expression)

# -----------------------------------------------------------------------------
# Postitively selected genes
# -----------------------------------------------------------------------------
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "SEMG1", "ARHGAP35", "CBL", "SSX1", "RBM12", "TRERF1", "FAT1", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Ensure age is numeric
#so.integrated@meta.data$age <- as.numeric(so.integrated@meta.data$age)

# Extract age and cluster information from the metadata
age_cluster_data <- so.integrated@meta.data %>%
	  select(age, seurat_clusters)

# Calculate correlations for each cluster and store results
cluster_ids <- levels(so.integrated@active.ident)
correlation_results <- lapply(cluster_ids, calculate_spearman_mn7, age_cluster_data = age_cluster_data, so.integrated = so.integrated)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_joined.RData")))
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "SEMG1", "ARHGAP35", "CBL", "SSX1", "RBM12", "TRERF1", "FAT1", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
# Step 2: Extract the scaled data for the genes in mn7
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[mn7, ]

# Step 3: Calculate the average expression for each gene in each cluster
average_expression_per_cluster <- sapply(levels(so.integrated), function(cluster) {
	  # Ensure no NAs by filling missing values
	  rowMeans(scaled_data[, Idents(so.integrated) == cluster, drop = FALSE], na.rm = TRUE)
})

# Step 4: Sort the genes sequentially by expression across clusters
ordered_genes <- rownames(scaled_data)[do.call(order, as.data.frame(-average_expression_per_cluster))]
scaled_data <- scaled_data[ordered_genes,]

# Step 4: Order the genes based on their expression levels from cluster 1 to the last cluster
# Sort each cluster's expression and then combine orders without dropping any genes
#gene_order <- order(rowMeans(average_expression_per_cluster, na.rm = TRUE), decreasing = TRUE)
#ordered_genes <- rownames(scaled_data)[gene_order]
#scaled_data <- scaled_data[gene_order,]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + theme(axis.text.y = element_text(size = 12))
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.5_group.by_-12-15-17_mn7_ordered_12.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
# > dim(so.integrated)
# [1] 5000 70738

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "Endothelial and Macrophage")
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)
# > dim(so.integrated)
# [1]  5000 58098

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
																							graph_label_size=0)
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_integrated.png")), plot = monocle3, dpi = 300)

# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#tscan
library(ggbeeswarm)
library(colorRamps)
library(viridisLite)

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$idents <- Idents(so.integrated)

#pseudotime_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
#pseudotime_palette <- inferno(100)
pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)

ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = pseudotime_monocle3)) +  # Use pseudotime as the color
  	geom_quasirandom(groupOnX = FALSE) +
	  scale_color_gradientn(colors = pseudotime_palette) +  # Apply the pseudotime gradient
	  theme_classic() +
	  xlab("pseudotime") + 
	  ylab("") +
	  ggtitle("") +
	  theme(legend.position = "none")  # Equivalent to NoLegend()
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered_integrated.png")), plot = ordered, dpi = 300)






# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 09/09/24
# -----------------------------------------------------------------------------
# Extract expression data from Seurat object
expression_matrix <- GetAssayData(so.integrated[["RNA"]], slot = "data")

# Create Monocle3 CDS manually
cds <- new_cell_data_set(
	  expression_data = expression_matrix,
	  cell_metadata = as.data.frame(so.integrated@meta.data),
	  gene_metadata = data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
)
# Optionally add UMAP embeddings if available
reducedDim(cds, "UMAP") <- Embeddings(so.integrated, "umap")

colData(cds)$celltype_column <- Idents(so.integrated)
cds$age <- as.numeric(colData(cds)$age)  # Ensure age is numeric
cds$batch <- as.factor(colData(cds)$batch)  # Ensure batch is a factor
save(cds, expression_matrix, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_cds_RNA.data_36601.RData")))

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

# Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds)
#plot_cells(cds)
plot_cells(cds, color_cells_by="celltype_column")

# Find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, group_cells_by="celltype_column",	reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
	  filter(fraction_expressing >= 0.10) %>%
	  group_by(cell_group) %>%
	  top_n(10, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
																				top_specific_marker_ids,
																				group_cells_by="celltype_column",
																				ordering_type="maximal_on_diag",
																				max.size=3)

# Learn the trajectory graph
# Perform trajectory analysis using precomputed UMAP
cds <- cluster_cells(cds, reduction_method = "UMAP")
# Learn the principal graph
cds <- learn_graph(cds, use_partition = F)




plot_cells(cds,
											color_cells_by = "celltype_column",
											label_groups_by_cluster=FALSE,
											label_leaves=FALSE,
											label_branch_points=FALSE)

plot_cells(cds,
											color_cells_by = "celltype_column",
											label_cell_groups=FALSE,
											label_leaves=TRUE,
											label_branch_points=TRUE,
											graph_label_size=1.5)







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
	  significant_genes <- subset(coef_table, term == "age" & q_value < 0.05)
	  
	  # Store results for this cell type
	  results_list[[cell_type]] <- significant_genes
}

# -----------------------------------------------------------------------------
# DESeq2
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
genes <- rownames(so_subset2)
for (gene in genes) {
	  # Extract expression data for a gene of interest (e.g., "GENE_NAME")
	  metadata_cols <- c("age", "batch", "sample.id")
	  gene_expression <- FetchData(so_subset, vars = c(gene, metadata_cols))
	  pseudobulk_counts <- gene_expression_filtered %>%
	                       group_by(sample.id) %>%
	  	                    summarize(total_expression = sum(!!sym(gene))) %>%
	  	                    as.data.frame()
	  
	  pseudobulk_metadata <- gene_expression %>%
	  	  select(sample.id, age, batch) %>%
	  	  distinct()
	  pseudobulk_metadata <- pseudobulk_metadata[match(pseudobulk_counts$sample.id, pseudobulk_metadata$sample.id),]
	  
	  pseudobulk_counts <- cbind(pseudobulk_counts, pseudobulk_metadata[, 2-3])
	  
	  pseudobulk_counts
	  
	  colnames(gene_expression)[1] <- "expression"
	  # Convert to a data frame for ggplot2
	  gene_expression_df <- as.data.frame(gene_expression)
	  # Filter cells with non-zero expression for the gene of interest
	  gene_expression_df_filtered <- subset(gene_expression_df, expression > 0)
	
	  # Define the custom color palette for pseudotime
	  pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)
	
	  # Proceed with plotting
	  ggplot(gene_expression_df_filtered, aes(x = age, y = expression, color = pseudotime)) +
		    geom_point(size = 3) +
		    scale_color_gradientn(colors = pseudotime_palette, limits = c(0, 8)) +
		    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add linear fit
		    theme_minimal() +
		    labs(title = gene, 
					    		x = "Age", 
							    y = "Expression", 
							    color = "pseudotime") +
		    theme(plot.title = element_text(hjust = 0.5))
}

# Define the custom color palette for pseudotime
age_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)

# Proceed with plotting
ggplot(gene_expression_df_filtered, aes(x = age, y = expression, color = pseudotime)) +
	geom_point(size = 3) +
	scale_color_gradientn(colors = pseudotime_palette, limits = c(0, 8)) +
	geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add linear fit
	theme_minimal() +
	labs(title = gene, 
						x = "Age", 
						y = "Expression", 
						color = "pseudotime") +
	theme(plot.title = element_text(hjust = 0.5))







# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
# > dim(so.integrated)
# [1] 5000 70738

# Set the desired DefaultAssay in Seurat (e.g., "RNA")
DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Cluster cells in Monocle 3 using UMAP embeddings
#cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)
# Extract expression data from Seurat object
expression_matrix <- GetAssayData(so.integrated[["RNA"]], slot = "data")

# Create Monocle3 CDS manually
cds <- new_cell_data_set(
	  expression_data = expression_matrix,
	  cell_metadata = as.data.frame(so.integrated@meta.data),
	  gene_metadata = data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
)
# Optionally add UMAP embeddings if available
reducedDim(cds, "UMAP") <- Embeddings(so.integrated, "umap")

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

# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "pseudotime",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=0)
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_RNA.data_36601_UMAP.png")), plot = monocle3, dpi = 300)

save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_cds_RNA.RData")))

# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#tscan
library(ggbeeswarm)
library(colorRamps)
library(viridisLite)

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$idents <- Idents(so.integrated)

#pseudotime_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
#pseudotime_palette <- inferno(100)
pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)

ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = pseudotime_monocle3)) +  # Use pseudotime as the color
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_gradientn(colors = pseudotime_palette) +  # Apply the pseudotime gradient
	  theme_classic() +
	  xlab("pseudotime") + 
	  ylab("") +
	  ggtitle("") +
	  theme(legend.position = "none")  # Equivalent to NoLegend()
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered_RNA_36601.png")), plot = ordered, dpi = 300)

# -----------------------------------------------------------------------------
# Find temporally expressed genes associated with age
# Last Modified: 06/09/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_cds.RData")))
load(file=file.path(wd.de.data, paste0("DESeq_~batch+age.RData")))

so.integrated@meta.data$pseudotime <- monocle3::pseudotime(cds)

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Get the list of unique cell types
cell_types <- levels(Idents(so.integrated))
# Initialize a list to store results for each cell type
results_list2 <- list()

for (cell_type in cell_types) {
	  # Step 0: Subset the Seurat object to include only the current cell type
	  so_subset <- subset(so.integrated, idents = cell_type)
	
	  results_age <- results_list[[cell_type]]
	  result_df <- as.data.frame(results_age)
	  result_df$gene <- rownames(result_df)
	  result_df$color <- ifelse(result_df$log2FoldChange < -0.05, "blue",
	  																										ifelse(result_df$log2FoldChange > 0.05, "red", "grey"))
	  result_df <- subset(result_df, abs(log2FoldChange) > 0.05)
	  
	  so_subset <- so_subset[rownames(result_df),]
	  
	  # Step 1: Extract count data and metadata
	  counts <- GetAssayData(so_subset, layer = "counts")
	  # Filter out genes with no counts in any samples
	  non_zero_genes <- rowSums(counts > 0) > 0
	  counts_filtered <- counts[non_zero_genes, ]
	
	  metadata <- so_subset@meta.data
	  # Ensure the metadata columns of interest are present
	  # Assuming 'age' and 'batch' are stored in metadata
	  metadata$pseudotime <- as.numeric(metadata$pseudotime)  # Convert age to numeric if it isn't already
	  metadata$batch <- as.factor(metadata$batch)  # Ensure batch is a factor
	  # Center and scale the age variable
	  #metadata$age_scaled <- scale(metadata$age)
	
	  results_list2 <- list()
	  for (gene in rownames(counts_filtered)) {
	  	  # Subset the data to keep only cells with non-zero expression for the current gene
	  	  non_zero_cells <- counts_filtered[gene, ] > 0
	  	  counts_filtered_gene <- counts_filtered[gene, non_zero_cells]
	  	  # Subset metadata to match the filtered cells
	  	  metadata_filtered <- metadata[non_zero_cells, ]

	  	  counts_filtered_gene <- as.matrix(t(counts_filtered_gene))
	  	  rownames(counts_filtered_gene) <- gene
	  	  
	  	  # Step 2: Create a DESeqDataSet for the current gene
	  	  dds <- DESeqDataSetFromMatrix(countData = counts_filtered_gene,
	  	  																														colData = metadata_filtered,
	  	  																														design = ~ batch + pseudotime)
	  	  # Run DESeq2
	  	  #dds <- DESeq(dds)
	  	  # Step 1: Estimate size factors
	  	  dds <- estimateSizeFactors(dds)
	  	  # Step 2: Estimate dispersions using gene-wise estimates instead of the default curve fitting
	  	  dds <- estimateDispersionsGeneEst(dds)
	  	  # Step 3: Assign the gene-wise dispersion estimates as final dispersions
	  	  dispersions(dds) <- mcols(dds)$dispGeneEst
	  	  # Step 4: Continue with testing using nbinomWaldTest (default test)
	  	  dds <- nbinomWaldTest(dds)
	  	  
	  	  # You can now extract results for the current gene
	  	  results_pseudotime <- results(dds, name = "pseudotime")
	  	  # Step 6: Store the results for this gene
	  	  results_list1[[gene]] <- results_pseudotime
	  }
	  
	  # Step 2: Create a DESeqDataSet
	  dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
					  																										colData = metadata,
			  																												design = ~ batch + pseudotime)
	  # If necessary, filter out samples with very low counts (e.g., below 10)
	  #dds <- dds[, sample_sums > 10]
	  # Use the poscounts method for size factor estimation
	  dds <- estimateSizeFactors(dds, type = "poscounts")
	
	  # Step 4: Run DESeq2 with the LRT to test for the effect of age
	  dds <- DESeq(dds)
	  # Step 5: Extract the results for the age effect
	  results_pseudotime <- results(dds, name = "pseudotime")
	  
	  results_pseudotime_df <- as.data.frame(results_pseudotime)
	  results_pseudotime_sorted <- results_pseudotime_df[order(results_pseudotime_df$padj), ]
	  
	  # View the top results
	  head(results_pseudotime[order(results_pseudotime$pvalue), ])
	  # Store the results in the list with the cell type as the name
	  results_list2[[cell_type]] <- results_pseudotime
	
	  # Assume results_age is the DESeq2 results object for the effect of age
	  # Convert to a data frame for ggplot2
	  result_df <- as.data.frame(results_pseudotime)
	  result_df$gene <- rownames(result_df)
	  result_df$color <- ifelse(result_df$log2FoldChange < -0.05, "blue", 
	  																										ifelse(result_df$log2FoldChange > 0.05, "red", "grey"))
	  
	  pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)
	  
	  
}

plot_temporal_expression <- funciton(so_subset, results_pseudotime, p) {
	  results_pseudotime_df <- as.data.frame(results_pseudotime)
	  results_pseudotime_sorted <- results_pseudotime_df[order(results_pseudotime_df$padj), ]
	  results_pseudotime_sorted_padj <- subset(results_pseudotime_sorted, padj < p)
	
	  so_subset2 <- so_subset[rownames(results_pseudotime_sorted_padj),]
	  genes <- rownames(so_subset2)
	  for (gene in genes) {
	  	  # Extract expression data for a gene of interest (e.g., "GENE_NAME")
	  	  gene_expression <- FetchData(so_subset2, vars = c(gene, "age", "pseudotime"))
	  	  colnames(gene_expression)[1] <- "expression"
	  	  # Convert to a data frame for ggplot2
	  	  gene_expression_df <- as.data.frame(gene_expression)
	  	  # Filter cells with non-zero expression for the gene of interest
	  	  gene_expression_df_filtered <- subset(gene_expression_df, expression > 0)
	  	  
	  	  # Define the custom color palette for pseudotime
	  	  pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)
	  	  
	  	  # Proceed with plotting
	  	  ggplot(gene_expression_df, aes(x = age, y = expression, color = pseudotime)) +
	  	  	  geom_point(size = 3) +
	  	  	  scale_color_gradientn(colors = pseudotime_palette) +
	  	    	geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add linear fit
	  	  	  theme_minimal() +
	  	  	  labs(title = gene, 
	  	  			  			x = "Age", 
	  	  					  	y = "Expression", 
	  	  						  color = "pseudotime") +
	  	  	  theme(plot.title = element_text(hjust = 0.5))
	  }
}










# -----------------------------------------------------------------------------
# Methods: Slingshot
# Last Modified: 10/09/24
# -----------------------------------------------------------------------------
library(SingleCellExperiment)
library(slingshot)

load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Extract necessary data for SCE conversion
normalized_matrix <- GetAssayData(so.integrated, assay = "RNA", layer = "data")  # For raw counts
umap_embeddings <- Embeddings(so.integrated, reduction = "umap")  # UMAP embeddings
cell_metadata <- so.integrated@meta.data  # Metadata

# Create a SingleCellExperiment object
sce <- SingleCellExperiment(
	  assays = list(data = normalized_matrix),
	  colData = cell_metadata,
	  reducedDims = list(UMAP = umap_embeddings)
)
sce$seurat_clusters <- Idents(so.integrated)

# Running Slingshot using UMAP embeddings and cluster information from Seurat
#sce <- slingshot(sce, clusterLabels = sce$seurat_clusters, reducedDim = 'UMAP')
# Add clustering information (from Seurat) to the SCE object as a factor
sce <- slingshot(sce, clusterLabels = sce$seurat_clusters, reducedDim = 'UMAP', start.clus = "Stage 0", end.clus = "Late spermatid")  # '0' is the start cluster

# View pseudotime values
pseudotime_values <- sce$slingPseudotime_1
head(pseudotime_values)

# Convert UMAP embeddings to a data frame for ggplot
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
umap_df$pseudotime <- pseudotime_values

# Plot the UMAP embeddings with pseudotime overlay
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
	  geom_point() +
	  scale_color_viridis_c() +
	  labs(title = "Slingshot Pseudotime Trajectory on UMAP", x = "UMAP 1", y = "UMAP 2") +
	  theme_minimal()

# -----------------------------------------------------------------------------
# DESeq2
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("DESeq_~batch+age.RData")))

cell_type <- "Sertoli"
results_age <- results_list[[cell_type]]

# Assume results_age is the DESeq2 results object for the effect of age
# Convert to a data frame for ggplot2
result_df <- as.data.frame(results_age)
result_df$gene <- rownames(result_df)
result_df$color <- ifelse(result_df$log2FoldChange < -0.05, "blue",
																										ifelse(result_df$log2FoldChange > 0.05, "red", "grey"))

# Create a volcano plot
volcano_plot <- ggplot(result_df, aes(x = log2FoldChange, y = -log10(padj))) +
	geom_point(aes(color = color)) +
	scale_color_manual(values = c(blue, "grey", red)) +
	theme_minimal() +
	labs(title = paste0(cell_type),
						x = "Log2 Fold Change",
						y = "-log10 Adjusted P-value") +
	theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
	geom_vline(xintercept = 0.05, linetype = "dashed", color = red) +
	geom_vline(xintercept = -0.05, linetype = "dashed", color = blue) +
	geom_text_repel(data = subset(result_df, abs(log2FoldChange) > 0.05),
																	aes(label = gene), size = 3, max.overlaps = 20)

# Print the volcano plot
pdf(file = file.path(wd.de.plots, paste0("volcano_DeSeq_~batch+age_", cell_type, ".pdf")), width = 5, height = 5)
print(volcano_plot)
dev.off()

# -----------------------------------------------------------------------------
# DESeq2 (mn7)
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))
cell_types <- levels(Idents(so.integrated))
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "SEMG1", "ARHGAP35", "CBL", "SSX1", "RBM12", "TRERF1", "FAT1", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

for (cell_type in cell_types) {
	results_age <- results_list[[cell_type]]
	
	# Assume results_age is the DESeq2 results object for the effect of age
	# Convert to a data frame for ggplot2
	result_df <- as.data.frame(results_age)
	result_df$gene <- rownames(result_df)
	result_df <- result_df[mn7,]
	result_df$color <- ifelse(result_df$log2FoldChange < 0, "blue",
																											ifelse(result_df$log2FoldChange > 0, "red", "grey"))
	
	# Create a volcano plot
	volcano_plot <- ggplot(result_df, aes(x = log2FoldChange, y = -log10(padj))) +
		geom_point(aes(color = color)) +
		scale_color_manual(values = c(blue, red)) +
		theme_minimal() +
		labs(title = paste0(cell_type),
							x = "Log2 Fold Change",
							y = "-log10 Adjusted P-value") +
		theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
		geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
		geom_text_repel(data = subset(result_df, padj < 1),
																		aes(label = gene), size = 3, max.overlaps = 30)
	
	# Print the volcano plot
	pdf(file = file.path(wd.de.plots, paste0("mn7_volcano_DeSeq_~batch+age_", cell_type, ".pdf")), width = 5, height = 5)
	print(volcano_plot)
	dev.off()
}



# -----------------------------------------------------------------------------
# Find temporally expressed genes associated with age
# Last Modified: 04/09/24
# -----------------------------------------------------------------------------
# 1. Fit a Model to Pseudotime:
# Fit models to gene expression along pseudotime
gene_fits <- fit_models(cds, model_formula_str = "~ pseudotime")

# 2. Test for Temporal Expression:
# Perform the graph test to find temporally expressed genes
temporal_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
# Filter for significant genes
temporal_genes_filtered <- temporal_genes %>% 
  	filter(q_value < 0.05)  # Adjust significance threshold as needed

# 3. Correlate Gene Expression with Age:
# Fit models that incorporate age as a covariate
age_fits <- fit_models(cds, model_formula_str = "~ pseudotime + age")

# 4. Subset by Cell Type:
# Identify the different cell types
cds@colData$cell_type <- Idents(so.integrated)
cell_types <- unique(cds@colData$cell_type)

# Loop over each cell type to perform the analysis
for (cell_type in cell_types) {
	  # Subset cds by the current cell type
	  cds_subset <- cds[, cds@colData$cell_type == cell_type]
	  # Fit models for pseudotime and age
	  age_fits <- fit_models(cds_subset, model_formula_str = "~ age")
	  # Extract significant genes
	  temporal_genes <- graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 4)
	
	  # Filter for significant genes
	  significant_genes <- temporal_genes %>% 
		    filter(q_value < 0.01)  # Adjust threshold as needed
	
	  # Store or visualize results for the current cell type
	  print(paste("Significant genes for cell type:", cell_type))
	  print(significant_genes)
	  
	  #  5. Visualize Results:
	  # Plot gene expression over pseudotime for a specific gene
	  plot_cells(cds_subset, genes = "CFAP61", color_cells_by = "pseudotime")
	  # Plot gene expression with respect to age
	  mn7_5000 <- intersect(rownames(cds_subset), mn7)
	  plot_genes_in_pseudotime(cds_subset[mn7_5000,], color_cells_by = "age", min_expr = 0.1)
}










# -----------------------------------------------------------------------------
# Find temporally expressed genes associated with age
# Last Modified: 04/09/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_cds.RData")))
cds@colData$cell_type <- Idents(so.integrated)

for (cell_type in cell_types) {
	# Subset cds by the current cell type
	cds_subset <- cds[, cds@colData$cell_type == cell_type]
	
	# 1. Fit a Model with Pseudotime and Age
	# Fit models to gene expression along pseudotime
	model_fits <- fit_models(cds_subset, model_formula_str = "~ pseudotime + age")
	
	# 2. Test for Temporal Expression Related to Pseudotime
	# Test for genes whose expression is associated with pseudotime
	temporal_genes <- graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 4)
	# Filter for significant genes based on pseudotime
	temporal_genes_filtered <- temporal_genes %>% 
		filter(q_value < 0.05)  # Adjust significance threshold as needed
	
	# 3. Test for Correlation with Age
	# Extract genes significantly associated with age
	age_genes <- as.data.frame(model_fits) %>%
		filter(p_adj_value_age < 0.05)  # Adjust the p-value threshold as needed
	
	# 4. Find Genes Related to Both Age and Pseudotime
	# Find genes related to both pseudotime and age
	genes_both <- intersect(temporal_genes_filtered$gene_id, age_genes$gene_id)
	# Subset data to get only these genes
	genes_both_df <- temporal_genes_filtered %>%
		filter(gene_id %in% genes_both)
	
	# 5. Visualize the Results
	# Plot expression of genes over pseudotime
	plot_genes_in_pseudotime(cds, genes = genes_both_df$gene_id, color_cells_by = "age")
	# Alternatively, plot pseudotime and age association
	plot_cells(cds, genes = genes_both_df$gene_id, color_cells_by = "pseudotime")
}







# -----------------------------------------------------------------------------
# DESeq2
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
so.integrated <- JoinLayers(so.integrated, assays = "RNA")

# Get the list of unique cell types
cell_types <- levels(Idents(so.integrated))
# Initialize a list to store results for each cell type
results_list <- list()

for (cell_type in cell_types) {
	  # Step 0: Subset the Seurat object to include only the current cell type
	  so_subset <- subset(so.integrated, idents = cell_type)
	
	  # Step 1: Extract count data and metadata
	  counts <- GetAssayData(so_subset, layer = "counts")
	  # Filter out genes with no counts in any samples
	  non_zero_genes <- rowSums(counts > 0) > 0
	  counts_filtered <- counts[non_zero_genes, ]
	  
	  metadata <- so_subset@meta.data
	  # Ensure the metadata columns of interest are present
	  # Assuming 'age' and 'batch' are stored in metadata
	  metadata$age <- as.numeric(metadata$age)  # Convert age to numeric if it isn't already
	  #metadata$batch <- as.factor(metadata$batch)  # Ensure batch is a factor
	  # Center and scale the age variable
	  #metadata$age_scaled <- scale(metadata$age)
	  
	  # Step 2: Create a DESeqDataSet
	  dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
					  																										colData = metadata,
			  																												design = ~ age)
	  # Use the poscounts method for size factor estimation
	  dds <- estimateSizeFactors(dds, type = "poscounts")
	  
  	# Step 4: Run DESeq2 with the LRT to test for the effect of age
  	dds <- DESeq(dds)
	  # Step 5: Extract the results for the age effect
	  results_age <- results(dds, name = "age")
	  # View the top results
	  head(results_age[order(results_age$pvalue), ])
	  # Store the results in the list with the cell type as the name
	  results_list[[cell_type]] <- results_age
	  
	  # Assume results_age is the DESeq2 results object for the effect of age
	  # Convert to a data frame for ggplot2
	  result_df <- as.data.frame(results_age)
	  result_df$gene <- rownames(result_df)
	  result_df$color <- ifelse(result_df$log2FoldChange < -0.05, "blue",
	  																										ifelse(result_df$log2FoldChange > 0.05, "red", "grey"))
	  
	  # Create a volcano plot
	  volcano_plot <- ggplot(result_df, aes(x = log2FoldChange, y = -log10(padj))) +
	  	  geom_point(aes(color = color)) +
	  	  scale_color_manual(values = c(blue, "grey", red)) +
	  	  theme_minimal() +
	  	  labs(title = paste0(cell_type),
	  			  			x = "Log2 Fold Change",
	  				  		y = "-log10 Adjusted P-value") +
	  	  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
	    	geom_vline(xintercept = 0.05, linetype = "dashed", color = red) +
	    	geom_vline(xintercept = -0.05, linetype = "dashed", color = blue) +
	    	geom_text_repel(data = subset(significant_genes, abs(log2FoldChange) > 0.05),
	  														  			aes(label = gene), size = 3, max.overlaps = 30)
	  
	  # Print the volcano plot
	  pdf(file = file.path(wd.de.plots, paste0("volcano_DeSeq_~age_", cell_type, ".pdf")), width = 5, height = 5)
	  print(volcano_plot)
	  dev.off()
}




















# Ensure the gene list is in your Seurat object
genes_of_interest <- mn7[mn7 %in% rownames(so.integrated)]
# Subset the scaled data for these genes
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[genes_of_interest, ]

library(pheatmap)

summary(is.na(scaled_data))
summary(is.nan(scaled_data))
summary(is.infinite(scaled_data))

# Plot the heatmap
pheatmap(scaled_data, 
									cluster_rows = TRUE, 
									cluster_cols = TRUE, 
									show_rownames = TRUE, 
									show_colnames = FALSE,
									scale = "row", 
									color = colorRampPalette(c("blue", "white", "red"))(50))





# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
















load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_500",      "ssc_filtered_normalised.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500", "ssc_filtered_normalised.1.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500", "ssc_filtered_normalised.2.RData"))

so.list <- c(so.list, so.list.1, so.list.2)
ids <- c(ids, ids.1, ids.2)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.1)
samples0.filtered <- rbind(samples0.filtered, samples0.2)

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000

## Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
so.list <- lapply(X = so.list, FUN = SCTransform)

## As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
## Run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

## When running FindIntegrationAnchors(), and IntegrateData(), set the normalization.method parameter to the value SCT.
## When running sctransform-based workflows, including integration, do not run the ScaleData() function
anchors <- FindIntegrationAnchors(object.list = so.list, normalization.method = "SCT", anchor.features = features, dims = 1:30)
so.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30, k.weight = 90)

save(samples0.filtered, so.list, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_SCT_", nfeatures, ".RData")))

# -----------------------------------------------------------------------------
# Perform CCA integration (Seurat 4.3.0; Running CCA)
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
ids <- c()
for (s in 1:nrow(samples0.filtered)) {
	  ids <- c(ids, rep(samples0.filtered$V3[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

ages <- c()
for (s in 1:nrow(samples0.filtered)) {
	  ages <- c(ages, rep(samples0.filtered$V4[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

n2s <- c()
for (s in 1:nrow(samples0.filtered)) {
	  n2s <- c(n2s, rep(samples0.filtered$V8[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

batches <- c()
for (s in 1:nrow(samples0.filtered)) {
	  batches <- c(batches, rep(samples0.filtered$V9[s], ncol(so.list[[s]]@assays$RNA$counts)))
}

so.integrated@meta.data$sample.id <- ids
so.integrated@meta.data$age <- ages
so.integrated@meta.data$age <- factor(so.integrated@meta.data$age, levels = c("25","27","37","40","48","57","60","71"))
so.integrated@meta.data$n2 <- n2s
so.integrated@meta.data$batch <- batches
head(so.integrated@meta.data)

save(samples0.filtered, so.integrated, ids, ages, n2s, batches, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_", nfeatures, ".RData")))

# -----------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
# 02_UMAP
# https://satijalab.org/seurat/articles/multimodal_vignette
# -----------------------------------------------------------------------------
so.integrated <- RunPCA(so.integrated, verbose = F)

pdf(file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_ElbowPlot_SCT_", nfeatures, ".pdf")))
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.integrated, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.integrated[["pca"]]@stdev / sum(so.integrated[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_3000.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(so.integrated, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
resolution.range <- seq(from = 0.05, to = 1, by = 0.05)

so.integrated <- FindNeighbors(so.integrated, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
save(samples0.filtered, so.integrated, prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, ".RData")))

# -----------------------------------------------------------------------------
# Define resolution
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, ".RData")))
resolution.range <- seq(from = 0.5, to = 1, by = 0.05)

for (r in 1:length(resolution.range)) {
	  res <- resolution.range[r]
	  resolution_name <- paste0("integrated_snn_res.", res)
	  Idents(so.integrated) <- so.integrated[[resolution_name]][[1]]
	
	  # Extract the UMAP embeddings
	  umap_embeddings <- Embeddings(so.integrated, reduction = "umap")
	  # Flip the UMAP 2 (second dimension)
	  umap_embeddings[, 2] <- umap_embeddings[, 2] * -1
	  # Re-insert the modified embeddings back into the Seurat object
	  so.integrated[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAP_", assay = DefaultAssay(so.integrated))
	  
	  ##
	  file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_res=", res, "_flip")
	  pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
	  DimPlot(so.integrated, label = TRUE)
	  dev.off()
}

# Find neighbors and clusters
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))

so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = 0.55)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.55_", nfeatures, ".RData")))

##
file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.55")

pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
DimPlot(so.integrated, label = TRUE)
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
	  select(seurat_clusters, sample.id, proportion)

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
color_vector["AMSBIO"] <- "yellow"

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
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, ".RData")))

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.45.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")
new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
  scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.6_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
cluster_to_celltype <- c('2' = 'Undiff. SPG',
																									'10' = 'Diff. SPG',
																									'9' = 'Leptotene',
																									'11' = 'Zygotene', '5' = 'Zygotene',
																									'12' = 'Pachytene', '8' = 'Pachytene',
																									'13' = 'Diplotene',
																									'6' = 'Meiotic division',
																									'4' = 'Late spermatid',
																									'16' = 'Endothelial & Macrophage',
																									'7' = 'PMC',
																									'15' = 'Leydig',
																									'1' = 'Sertoli',	'17' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.3_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SCT_", nfeatures, "_UMAP_dims=14_resolution=0.3_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.4; Six SPG states)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.45_", nfeatures, ".RData")))

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
new_order <- c("8", "2", "21", "20", "11", "25", "19", "4", "10", "22", "15", "23", "17", "14", "3", "6", "9", "12", "13", "24", "7", "18", "0", "16", "5", "1")
new_order <- rev(new_order)
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

cluster_to_celltype <- c('1' = 'Stage 0', '5' = 'Stage 0A', '16' = 'Stage 0B', '0' = 'Stage 1',
																									'18' = 'Stage 2', '7' = 'Stage 2', '24' = 'Stage 3',
																									'13' = 'Leptotene',
																									'12' = 'Zygotene', '9' = 'Zygotene', '6' = 'Zygotene',
																									'3' = 'Pachytene', '14' = 'Pachytene', '17' = 'Pachytene',
																									'23' = 'Diplotene', 	
																									'15' = 'Meiotic division', '22' = 'Meiotic division',
																									'10' = 'Early spermatid',
																									'4' = 'Late spermatid', '19' = 'Late spermatid',
																									'25' = 'Endothelial and Macrophage',
																									'11' = 'Fibrotic PMC',
																									'20' = 'Leydig',
																									'21' = 'Sertoli', '2' = 'Sertoli', '8' = 'Sertoli')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=", res, "_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SPG_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_res=", res, "_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)






# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
DefaultAssay(so.integrated) <- "RNA"

so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
  	group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- markers %>%
  	group_by(cluster) %>%
  	dplyr::filter(avg_log2FC > 1) %>%
  	slice_head(n = 10) %>%
	  ungroup()
#missing_genes <- setdiff(top10$gene, rownames(so.integrated))

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3_group.by.png")),
							plot = heatmap, width = 10, height = 12, dpi = 300)






# -----------------------------------------------------------------------------
# Batch
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_resolution=0.3.RData"))
length <- length(unique(so.integrated$seurat_clusters))

##
metadata <- so.integrated@meta.data
metadata$cluster <- so.integrated$seurat_clusters
batch_cluster_counts <- metadata %>%
	  group_by(cluster, batch) %>%
	  summarise(count = n()) %>%
	  ungroup()

facs <- toTable(0, 4, length, c("Cluster", "B1", "B2", "B3"))
facs$Cluster <- as.numeric(rownames(facs)) - 1
for (c in 0:length-1) {
	  subsets <- subset(batch_cluster_counts, cluster == c)
	  facs$B1[c+1] <- subset(subsets, batch == 1)$count + subset(subsets, batch == 2)$count
	  facs$B2[c+1] <- subset(subsets, batch == 3)$count
	  facs$B3[c+1] <- subset(subsets, batch == 4)$count
	
	  #total <- facs$B1[c+1] + facs$B2[c+1] + facs$B3[c+1]
	  #facs$B1[c+1] <- facs$B1[c+1] / total * 100
	  #facs$B2[c+1] <- facs$B2[c+1] / total * 100
	  #facs$B3[c+1] <- facs$B3[c+1] / total * 100
}
write.table(facs, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_counts.txt"), row.names=F, col.names=F, quote=F, sep='\t')

# -----------------------------------------------------------------------------
# 
# Last Modified: 10/10/23
# -----------------------------------------------------------------------------
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")

#facs <- readTable(file.path(wd.de.data, "batch_cluster_SCT_counts.txt"), header=T, rownames=T, sep="")
facs2 <- t(facs[,-1])
colnames(facs2) <- facs$Cluster
facs2 <- facs2[, new_order]
write.table(facs2, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_counts2.txt"), row.names=T, col.names=T, quote=F, sep='\t')

max <- 0
for (c in 1:length) {
	  if (sum(facs2[, c]) > max)
		    max <- sum(facs2[, c])
}

## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.de.plots, "batch_cluster_counts_SCT_0.3_new_order_numbers")
main.text <- c("", "")
xlab.text <- ""
ylab.text <- "Number of cells"
cols <- c(blue, red, yellow)
sequence <- seq(from = 0.8, by = 1.2, length.out = length)

pdf(paste0(file.name, ".pdf"), width = 12, height = 6)
par(mar=c(5.1, 4.6, 4.1, 9), xpd=TRUE)
barplot(facs2, col=cols, ylim=c(0, 10424), ylab=ylab.text, xaxt="n", main="", cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=new_order, x=sequence, y=par("usr")[3] - 360, srt=45, adj=0.965, xpd=NA, cex=1.9)

legend("right", c("B1", "B2", "B3"), text.col="black", pch=c(15, 15, 15), col=cols, pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
dev.off()

##
facs3 <- toTable(0, 4, length, c("Cluster", "B1", "B2", "B3"))
facs3$Cluster <- as.numeric(rownames(facs3)) - 1
for (c in 0:length-1) {
	  total <- facs$B1[c+1] + facs$B2[c+1] + facs$B3[c+1]
	  facs3$B1[c+1] <- facs$B1[c+1] / total * 100
	  facs3$B2[c+1] <- facs$B2[c+1] / total * 100
	  facs3$B3[c+1] <- facs$B3[c+1] / total * 100
}
facs3 <- t(facs3[,-1])
colnames(facs3) <- facs$Cluster
facs3 <- facs3[, new_order]
write.table(facs3, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_proportion.txt"), row.names=T, col.names=T, quote=F, sep='\t')

## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.de.plots, "batch_cluster_counts_SCT_0.3_new_order")
main.text <- c("", "")
xlab.text <- ""
ylab.text <- "Fraction of cells"
#blue  <- "blue"   ## adjustcolor("#619CFF", alpha.f=0.9)
#red   <- "red"   ## adjustcolor("#F8766D", alpha.f=0.9)
#green <- "darkgray"   ## adjustcolor("#00BA38", alpha.f=0.9)
#cols <- c(blue, red, green)   ## #59a523 (Alcro wasabi)
cols <- c(blue, red, yellow)
sequence <- seq(from = 0.8, by = 1.2, length.out = length)

pdf(paste0(file.name, ".pdf"), width = 12, height = 6)
par(mar=c(5.1, 4.6, 4.1, 9), xpd=TRUE)
barplot(facs3, col=cols, ylim=c(0, 100), ylab=ylab.text, xaxt="n", main="", cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=new_order, x=sequence, y=par("usr")[3] - 4, srt=45, adj=0.965, xpd=NA, cex=1.9)

legend("right", c("B1", "B2", "B3"), text.col="black", pch=c(15, 15, 15), col=cols, pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
#mtext(ylab.text, side=2, line=2.75, cex=1.8)
dev.off()






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
# Extract the total expression level (sum of counts) and cluster identities
mt_data <- get_cluster_factor_levels(so.integrated, "percent.mt")
box_plot_mt <- get_box_plot(mt_data, "Mitochondrial contents", "Percentage")

# Part 4: Combine the Plots
combined_plot <- box_plot_reads / box_plot_expression / box_plot_mt

# Print the combined plot
print(combined_plot)

# -----------------------------------------------------------------------------
# QC (Sampke ID)
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
	  select(seurat_clusters, sample.id, proportion)

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
# Sampke ID
# -----------------------------------------------------------------------------
# Extract sample IDs and cluster identities
sample_cluster_data <- FetchData(so.integrated, vars = c("sample_id", "seurat_clusters"))

# Calculate the number of cells in each cluster
total_cells_per_cluster <- sample_cluster_data %>%
	  group_by(seurat_clusters) %>%
	  summarize(total_cells = n())

# Calculate the number of cells for each sample ID within each cluster
cells_per_sample_per_cluster <- sample_cluster_data %>%
	  group_by(seurat_clusters, sample_id) %>%
	  summarize(cells = n())

# Merge the data frames to calculate proportions
proportion_data <- merge(cells_per_sample_per_cluster, total_cells_per_cluster, by = "seurat_clusters")

# Calculate the proportion of each sample ID within each cluster
proportion_data <- proportion_data %>%
	  mutate(proportion = cells / total_cells) %>%
	  select(seurat_clusters, sample_id, proportion)

# Rename columns for clarity
colnames(proportion_data) <- c("Cluster", "Sample_ID", "Proportion")

# Create a vector of unique sample IDs
unique_samples <- unique(proportion_data$Sample_ID)

# Create a vector of random colors
set.seed(42)  # For reproducibility
random_colors <- grDevices::colors()[sample(1:length(grDevices::colors()), length(unique_samples))]

# Create a named vector of colors, assigning specific colors to the highlighted samples
color_vector <- setNames(random_colors, unique_samples)
color_vector["PD40746e_M2"] <- "red"
color_vector["PD40746e_M1"] <- "orange"
color_vector["AMSBIO"] <- "yellow"

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
	  select(seurat_clusters, sample.id, proportion)

# Rename columns for clarity
colnames(proportion_data) <- c("Cluster", "Sample_ID", "Proportion")

# View the summary
print(proportion_data)






facs <- toTable(0, 1+length(samples), length, c("Cluster", samples))
facs$Cluster <- as.numeric(rownames(facs)) - 1
for (c in 0:length-1) {
	  subsets <- subset(batch_cluster_counts, cluster == c)
	  facs$B1[c+1] <- subset(subsets, batch == 1)$count + subset(subsets, batch == 2)$count
	  facs$B2[c+1] <- subset(subsets, batch == 3)$count
	  facs$B3[c+1] <- subset(subsets, batch == 4)$count
	
	#total <- facs$B1[c+1] + facs$B2[c+1] + facs$B3[c+1]
	#facs$B1[c+1] <- facs$B1[c+1] / total * 100
	#facs$B2[c+1] <- facs$B2[c+1] / total * 100
	#facs$B3[c+1] <- facs$B3[c+1] / total * 100
}
write.table(facs, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_counts.txt"), row.names=F, col.names=F, quote=F, sep='\t')

# -----------------------------------------------------------------------------
# 
# Last Modified: 10/10/23
# -----------------------------------------------------------------------------
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")

#facs <- readTable(file.path(wd.de.data, "batch_cluster_SCT_counts.txt"), header=T, rownames=T, sep="")
facs2 <- t(facs[,-1])
colnames(facs2) <- facs$Cluster
facs2 <- facs2[, new_order]
write.table(facs2, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_counts2.txt"), row.names=T, col.names=T, quote=F, sep='\t')

max <- 0
for (c in 1:length) {
	if (sum(facs2[, c]) > max)
		max <- sum(facs2[, c])
}

## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.de.plots, "batch_cluster_counts_SCT_0.3_new_order_numbers")
main.text <- c("", "")
xlab.text <- ""
ylab.text <- "Number of cells"
cols <- c(blue, red, yellow)
sequence <- seq(from = 0.8, by = 1.2, length.out = length)

pdf(paste0(file.name, ".pdf"), width = 12, height = 6)
par(mar=c(5.1, 4.6, 4.1, 9), xpd=TRUE)
barplot(facs2, col=cols, ylim=c(0, 10424), ylab=ylab.text, xaxt="n", main="", cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=new_order, x=sequence, y=par("usr")[3] - 360, srt=45, adj=0.965, xpd=NA, cex=1.9)

legend("right", c("B1", "B2", "B3"), text.col="black", pch=c(15, 15, 15), col=cols, pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
dev.off()

##
facs3 <- toTable(0, 4, length, c("Cluster", "B1", "B2", "B3"))
facs3$Cluster <- as.numeric(rownames(facs3)) - 1
for (c in 0:length-1) {
	  total <- facs$B1[c+1] + facs$B2[c+1] + facs$B3[c+1]
	  facs3$B1[c+1] <- facs$B1[c+1] / total * 100
	  facs3$B2[c+1] <- facs$B2[c+1] / total * 100
	  facs3$B3[c+1] <- facs$B3[c+1] / total * 100
}
facs3 <- t(facs3[,-1])
colnames(facs3) <- facs$Cluster
facs3 <- facs3[, new_order]
write.table(facs3, file=file.path(wd.de.data, "batch_cluster_SCT_0.3_proportion.txt"), row.names=T, col.names=T, quote=F, sep='\t')

## https://stackoverflow.com/questions/7588020/how-to-write-labels-in-barplot-on-x-axis-with-duplicated-names
file.name <- file.path(wd.de.plots, "batch_cluster_counts_SCT_0.3_new_order")
main.text <- c("", "")
xlab.text <- ""
ylab.text <- "Fraction of cells"
#blue  <- "blue"   ## adjustcolor("#619CFF", alpha.f=0.9)
#red   <- "red"   ## adjustcolor("#F8766D", alpha.f=0.9)
#green <- "darkgray"   ## adjustcolor("#00BA38", alpha.f=0.9)
#cols <- c(blue, red, green)   ## #59a523 (Alcro wasabi)
cols <- c(blue, red, yellow)
sequence <- seq(from = 0.8, by = 1.2, length.out = length)

pdf(paste0(file.name, ".pdf"), width = 12, height = 6)
par(mar=c(5.1, 4.6, 4.1, 9), xpd=TRUE)
barplot(facs3, col=cols, ylim=c(0, 100), ylab=ylab.text, xaxt="n", main="", cex.names=1.7, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=new_order, x=sequence, y=par("usr")[3] - 4, srt=45, adj=0.965, xpd=NA, cex=1.9)

legend("right", c("B1", "B2", "B3"), text.col="black", pch=c(15, 15, 15), col=cols, pt.cex=3, cex=1.9, horiz=F, bty="n", inset=c(-0.16, 0))
#mtext(ylab.text, side=2, line=2.75, cex=1.8)
dev.off()






# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# DefaultAssay(so.integrated) <- "SCT"
# Error in UseMethod(generic = "JoinLayers", object = object) : 
#  	no applicable method for 'JoinLayers' applied to an object of class "c('SCTAssay', 'Assay', 'KeyMixin')"
DefaultAssay(so.integrated) <- "RNA"

so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

save(so.integrated, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers.RData")))
save(markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers_markers.RData")))

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()
#missing_genes <- setdiff(top10$gene, rownames(so.integrated))

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
# Centering and scaling data matrix
# Warning: Different features in new layer data than already exists for scale.data
# Warning: Different cells in new layer data than already exists for scale.data

# Create the heatmap
#heatmap <- DoHeatmap(so.integrated, features = top10$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend()
#ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3_group.by.png")),
#							plot = heatmap, width = 15, height = 15, dpi = 300)

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# Set the new order
new_order <- c("1", "17", "15", "7", "16", "3", "4", "6", "13", "14", "8", "12", "5", "11", "9", "10", "2", "0")
new_cluster_order <- rev(new_order)
so.integrated <- SetIdent(so.integrated, value = factor(Idents(so.integrated), levels = new_cluster_order))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.3_new_order.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and cluster identities
read_counts <- FetchData(so.integrated, vars = c("nCount_RNA", "seurat_clusters"))

# Rename columns for clarity
colnames(read_counts) <- c("read_count", "cluster")

# Convert cluster to factor for ordered plotting
read_counts$cluster <- factor(read_counts$cluster)

# Calculate the median read count for each cluster
medians_read <- read_counts %>%
	  group_by(cluster) %>%
  	summarize(median_read_count = median(read_count))

# Reorder the cluster factor levels based on median read counts
read_counts$cluster <- factor(read_counts$cluster, levels = medians_read$cluster[order(medians_read$median_read_count)])

# Generate the box plot for read counts
box_plot_reads <- ggplot(read_counts, aes(x = cluster, y = read_count)) +
	  geom_boxplot() +
	  theme_minimal() +
	  labs(title = "Read Counts",
			  			x = "Cluster",
				  		y = "Number of Reads") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1),
			  				plot.title = element_text(hjust = 0.5))  # Center the title

# Part 2: Generate the Box Plot for Global Expression Levels
# Extract the total expression level (sum of counts) and cluster identities
expression_data <- FetchData(so.integrated, vars = c("nFeature_RNA", "seurat_clusters"))

# Rename columns for clarity
colnames(expression_data) <- c("total_expression_level", "cluster")

# Convert cluster to factor for ordered plotting
expression_data$cluster <- factor(expression_data$cluster)

# Calculate the median total expression level for each cluster
medians_expression <- expression_data %>%
	  group_by(cluster) %>%
	  summarize(median_total_expression_level = median(total_expression_level))

# Reorder the cluster factor levels based on median total expression levels
expression_data$cluster <- factor(expression_data$cluster, levels = medians_expression$cluster[order(medians_expression$median_total_expression_level)])

# Generate the box plot for global expression levels
box_plot_expression <- ggplot(expression_data, aes(x = cluster, y = total_expression_level)) +
  	geom_boxplot() +
	  theme_minimal() +
	  labs(title = "Global Expression Levels",
			  			x = "Cluster",
				  		y = "Total Expression Level") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1),
			  				plot.title = element_text(hjust = 0.5))  # Center the title

# Part 3: Combine the Plots
combined_plot <- box_plot_reads / box_plot_expression

# Print the combined plot
print(combined_plot)

# -----------------------------------------------------------------------------
# QC (Sampke ID)
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
	  select(seurat_clusters, sample.id, proportion)

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
# scMayoMap
# -----------------------------------------------------------------------------
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
sapply(pkgs, require, character.only = TRUE)
library(scMayoMap)

obj <- scMayoMap(data = markers, database=scMayoMapDatabase, tissue = 'testis')
plt <- scMayoMap.plot(scMayoMap.object = obj)

gns <- obj$markers$genes[obj$markers$cluster==3 & obj$markers$cluster==14 & obj$markers$celltype=='Spermatogonial stem cell']
gns <- strsplit(gns, ',')[[1]]
DotPlot(so.integrated, features = gns)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# Filter results to get marker genes for cluster 3
markers_cluster3 <- markers %>% filter(cluster == 3)

# View the first few rows of the marker genes for cluster 3
head(markers_cluster3)

# Optionally, save the results to a CSV file
write.csv(markers_cluster3, file = "markers_cluster3.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# QC (Number of reads)
# -----------------------------------------------------------------------------
# Extract the number of reads (nCount_RNA) and cluster identities
read_counts <- FetchData(so.integrated, vars = c("nCount_RNA", "seurat_clusters"))

# Rename columns for clarity
colnames(read_counts) <- c("read_count", "cluster")

# Calculate the median read count for each cluster
medians <- read_counts %>%
	  group_by(cluster) %>%
	  summarize(median_read_count = median(read_count))

# Reorder the cluster factor levels based on median read counts
read_counts$cluster <- factor(read_counts$cluster, levels = medians$cluster[order(medians$median_read_count)])

# Generate the box plot
box_plot <- ggplot(read_counts, aes(x = cluster, y = read_count)) +
	  geom_boxplot() +
	  theme_minimal() +
	  labs(title = "Box Plot of Read Counts per Cluster (Sorted by Median Read Count)",
			  			x = "Cluster",
				  		y = "Number of Reads") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(box_plot)

# -----------------------------------------------------------------------------
# QC (Global expression level)
# -----------------------------------------------------------------------------
# Extract the total expression level (sum of counts) and cluster identities
expression_data <- FetchData(so.integrated, vars = c("nFeature_RNA", "seurat_clusters"))

# Rename columns for clarity
colnames(expression_data) <- c("total_expression_level", "cluster")

# Convert cluster to factor for ordered plotting
expression_data$cluster <- factor(expression_data$cluster)

# Calculate the median total expression level for each cluster
medians <- expression_data %>%
	  group_by(cluster) %>%
	  summarize(median_total_expression_level = median(total_expression_level))

# Reorder the cluster factor levels based on median total expression levels
expression_data$cluster <- factor(expression_data$cluster, levels = medians$cluster[order(medians$median_total_expression_level)])

# Generate the box plot
box_plot <- ggplot(expression_data, aes(x = cluster, y = total_expression_level)) +
	  geom_boxplot() +
	  theme_minimal() +
	  labs(title = "Box Plot of Global Expression Levels per Cluster (Sorted by Median Total Expression Level)",
			  			x = "Cluster",
				  		y = "Total Expression Level") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(box_plot)


# -----------------------------------------------------------------------------
# DE
# -----------------------------------------------------------------------------
markers.3.14 <- FindMarkers(so.integrated, ident.1 = 3, ident.2 = 14)
head(markers.3.14, n = 5)

markers.3.7 <- FindMarkers(so.integrated, ident.1 = 3, ident.2 = 7)
# Add gene names as a column for labeling
markers.3.7$gene <- rownames(markers.3.7)
# Create a new column for labeling significant genes
markers.3.7$significant <- ifelse(markers.3.7$p_val_adj < 0.05 & abs(markers.3.7$avg_log2FC) > 0.25, "Significant", "Not Significant")

# Create the volcano plot
volcano_plot <- ggplot(markers.3.7, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
	geom_point(alpha = 0.8) +
	scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
	theme_minimal() +
	labs(title = "Clusters 3 vs 7",
						x = "Log2 Fold Change",
						y = "-Log10 Adjusted P-Value") +
	theme(legend.position = "top")

# Print the plot
print(volcano_plot)


markers.3.1 <- FindMarkers(so.integrated, ident.1 = 3, ident.2 = 1)
# Add gene names as a column for labeling
markers.3.1$gene <- rownames(markers.3.1)
# Create a new column for labeling significant genes
markers.3.1$significant <- ifelse(markers.3.1$p_val_adj < 0.05 & abs(markers.3.1$avg_log2FC) > 0.25, "Significant", "Not Significant")

# Create the volcano plot
volcano_plot <- ggplot(markers.3.1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
	geom_point(alpha = 0.8) +
	scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
	theme_minimal() +
	labs(title = "Clusters 3 vs 1",
						x = "Log2 Fold Change",
						y = "-Log10 Adjusted P-Value") +
	theme(legend.position = "top")

# Print the plot
print(volcano_plot)

# -----------------------------------------------------------------------------
# Jaccard index for overlap of gene sets
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers.RData")))

jaccard_matrix <- getJaccardIndex(markers)
# Visualize the Jaccard index matrix
pheatmap_plot <- pheatmap(jaccard_matrix, main = "Jaccard Index Between Clusters", angle_col = 0)

# Save the plot as a PNG file
png(file = file.path(wd.de.plots, "heatmap_jaccard_markers_SCT_res=0.3_02.png"), width = 7, height = 7, units = "in", res = 300)
print(pheatmap_plot)
dev.off()

# -----------------------------------------------------------------------------
# Functional enrichment analysis
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers.RData")))

go_results <- getGO(markers)
save(go_results, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers_GO.RData")))

# Loop through each cluster and generate plots
for (i in 1:length(go_results)) {
	  cluster_name <- paste0("Cluster ", i-1)
	  p <- plot_custom_barplot(go_results[[i]], cluster_name)
	
	  if (!is.null(p)) {
		    pdf(file = file.path(wd.de.plots, paste0("barplot_markers_SCT_res=0.3_cluster", i-1, ".pdf")), width = 7, height = 2.5)
		    print(p)
		    dev.off()
	  }
}

# -----------------------------------------------------------------------------
# GSEA
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers.RData")))

# Load MSigDB gene sets (Hallmark gene sets as an example)
msigdb <- msigdbr(species = "Homo sapiens", category = "H")
msigdb_list <- split(msigdb$gene_symbol, msigdb$gs_name)

msigdb_df <- bind_rows(
	  lapply(names(msigdb_list), function(term) {
		    data.frame(term = term, gene = msigdb_list[[term]], stringsAsFactors = FALSE)
	  })
)

# Run GSEA for each cluster and store the results in a list
clusters <- unique(markers$cluster)
clusters <- as.numeric(as.vector(clusters))
gsea_results <- lapply(clusters, run_gsea)
names(gsea_results) <- clusters

gsea_results <- gsea_results[!sapply(gsea_results, is.null)]
save(markers, gsea_results, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_markers_gsea_results_H.RData")))

# Apply to all GSEA results
significant_terms_list <- lapply(gsea_results, extract_significant_terms)
significant_terms_list <- significant_terms_list[!sapply(significant_terms_list, is.null)]

# Create a matrix of p-values for significant terms
term_names <- unique(unlist(lapply(significant_terms_list, function(x) x$Description)))
cluster_names <- names(significant_terms_list)



pvalue_matrix <- matrix(NA, nrow = length(term_names), ncol = length(cluster_names))
rownames(pvalue_matrix) <- term_names
colnames(pvalue_matrix) <- cluster_names

for (i in seq_along(significant_terms_list)) {
	  terms <- significant_terms_list[[i]]$Description
	  pvalues <- significant_terms_list[[i]]$p.adjust
	  pvalue_matrix[terms, names(significant_terms_list)[i]] <- pvalues
}

# Replace NA values in pvalue_matrix with a suitable value (e.g., maximum p-value)
pvalue_matrix[is.na(pvalue_matrix)] <- 1

# Generate heatmap of p-values
ht_pvalues <- Heatmap(pvalue_matrix, 
																						name = "p-value", 
																						col = colorRampPalette(c(blue, "white", red))(50),
																						show_row_names = TRUE, 
																						show_column_names = TRUE,
																						cluster_rows = TRUE, 
																						cluster_columns = TRUE,
																						heatmap_legend_param = list(title = "Adjusted P-Value"))

# Save the heatmap to a file
filename <- file.path(wd.de.plots, paste0("GSEA_H_Pvalue_SCT_", nfeatures, "_UMAP_resolution=0.3_cluster", cluster, "_p<1.pdf"))
pdf(file = filename, width = 10, height = 20)
draw(ht_pvalues)
dev.off()



# Create a matrix of enrichment scores (NES) for significant terms
nes_matrix <- matrix(NA, nrow = length(term_names), ncol = length(cluster_names))
rownames(nes_matrix) <- term_names
colnames(nes_matrix) <- cluster_names

for (i in seq_along(significant_terms_list)) {
	  terms <- significant_terms_list[[i]]$Description
	  nes_scores <- significant_terms_list[[i]]$NES
	  nes_matrix[terms, names(significant_terms_list)[i]] <- nes_scores
}

# Replace NA values in nes_matrix with a suitable value (e.g., 0 for non-significant)
nes_matrix[is.na(nes_matrix)] <- 0

# Generate heatmap of NES
ht_nes <- Heatmap(nes_matrix, 
																		name = "NES", 
																		col = colorRampPalette(c("white", red))(50),
																		show_row_names = TRUE, 
																		show_column_names = TRUE,
																		cluster_rows = TRUE, 
																		cluster_columns = TRUE,
																		heatmap_legend_param = list(title = "Normalized Enrichment Score"),
																		column_names_gp = gpar(fontsize = 12),  # Set font size for column names
																		column_names_rot = 0)  # Rotate column names by 45 degrees
                  #row_names_side = "left")  # Place y-axis labels on the left side

# Save the NES heatmap to a file
filename <- file.path(wd.de.plots, paste0("GSEA_H_NES_SCT_", nfeatures, "_UMAP_resolution=0.3_p<0.05_02.pdf"))
pdf(file = filename, width = 10, height = 8)
draw(ht_nes, annotation_legend_side = "bottom",  # Move legend to the left
					padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
dev.off()

# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression data in each cluster
# -----------------------------------------------------------------------------
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
correlation_results <- lapply(cluster_ids, calculate_spearman, age_cluster_data = age_cluster_data, so.integrated = so.integrated)

# Combine results from all clusters
all_significant_genes <- bind_rows(correlation_results, .id = "cluster")

# View the results
print(all_significant_genes)
# Optionally, save the results to a CSV file
write.table(all_significant_genes, file = "significant_genes_age_correlation.txt", row.names = FALSE)









# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

markers.cluster0 <- FindMarkers(so.integrated, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


# -----------------------------------------------------------------------------
# Monocle
# -----------------------------------------------------------------------------
library(Monocle2)  # Ensure Monocle 3 is installed

# Checking for a valid number of clusters in Seurat object
DefaultAssay(so.integrated) <- "RNA"

num_clusters <- length(unique(Idents(so.integrated)))
assertthat::assert_that(num_clusters > 0)

# Extract necessary components from Seurat object
expression_data <- GetAssayData(so.integrated, slot = "counts")
cell_metadata <- so.integrated@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_data), row.names = rownames(expression_data))

# Create Monocle 3 cell_data_set
cds <- new_cell_data_set(expression_data,
																									cell_metadata = cell_metadata,
																									gene_metadata = gene_metadata)
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_Monocle3_1.RData")))

# Preprocess the Data in Monocle 3
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds)
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_Monocle3_2.RData")))

# Cluster the cells (this step might be optional depending on your use case)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
save(cds, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.3_", nfeatures, "_Monocle3.RData")))

# Plot cells colored by pseudotime
pseudotimeplot <- plot_cells(cds, color_cells_by = "pseudotime",
											label_groups_by_cluster = FALSE,
											label_leaves = FALSE,
											label_branch_points = FALSE)

ggsave(filename = file.path(wd.de.plots, paste0("Monocle3_pseudotime_SCT_", nfeatures, "_UMAP_resolution=0.25.png")),
							plot = pseudotimeplot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
library(AUCell)

# Assuming 'seurat_object' is your Seurat object
expression_matrix <- GetAssayData(so.integrated, assay = "SCT", layer = "data")


library(SingleR)

vignette("ExperimentHub", package = "ExperimentHub")
reference_dataset <- celldex::HumanPrimaryCellAtlasData()  # Example reference

annotations <- SingleR(test = so.integrated@assays$SCT@data, ref = reference_dataset, labels = reference_dataset$label.main)

