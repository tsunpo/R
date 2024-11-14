# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 14/03/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Graphics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg38.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "Shosei"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

# -----------------------------------------------------------------------------
# Standard Seurat re-processing workflow
# 01_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#library(sctransform)

count_matrix <- readTable(file.path(wd.de.data, "GRCm38_FeatureCounts.txt"), header=T, rownames=T, sep="\t")
count_matrix <- count_matrix[,-c(1:6)]
so <- CreateSeuratObject(counts = count_matrix)

DefaultAssay(so) <- "RNA"
so <- NormalizeData(so)

metadata <- readTable(file.path(wd.de.data, "GRCm38_FeatureCounts_metadata.txt"), header=T, rownames=T, sep="\t")
so@meta.data$sample    <- metadata$Sample
so@meta.data$sample.id <- metadata$Sample_ID
so@meta.data$cell.type <- metadata$Cell_Type
so@meta.data$set       <- metadata$Set
head(so@meta.data)

# QC and selecting cells for further analysis
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^mt-")
save(so, file=file.path(wd.de.data, "so.RData"))

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and cluster identities
#read_counts <- get_cluster_factor_levels(so, "nCount_RNA", factor = "sample")
#box_plot_reads <- get_box_plot(read_counts, "Read counts", "Number")

# Part 2: Generate the Box Plot for Global Expression Levels
# Extract the total expression level (sum of counts) and cluster identities
#expression_data <- get_cluster_factor_levels(so, "nFeature_RNA", factor = "sample")
#box_plot_expression <- get_box_plot(expression_data, "Global expression levels", "Expression")

# Part 3: Generate the Box Plot for Mitochondrial Contents
# Extract the total expression level (sum of counts) and cluster identities
#mt_data <- get_cluster_factor_levels(so, "percent.mt", factor = "sample")
#box_plot_mt <- get_box_plot(mt_data, "Mitochondrial contents", "Percentage")

# Part 4: Combine the Plots
#combined_plot <- box_plot_reads / box_plot_expression / box_plot_mt

# Print the combined plot
#print(combined_plot)

# -----------------------------------------------------------------------------
# QC (Read counts)
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and sample identities
read_counts <- data.frame(cell_type = so@meta.data$cell.type, nCount_RNA = so@meta.data$nCount_RNA, set = so@meta.data$set)
# Set the order of cell types
read_counts$cell_type <- factor(read_counts$cell_type, levels = c("Liver", "Spleen",	"Spermatogonia_PIvap_positive",
																																																																		"Spermatogonia_CD9_positive", "Spermatogonia_TA", 
																																																																		"Spermatocytes", "Spermatids", "Sperm"))
# Generate a box plot for read counts grouped by sample ID
box_plot_reads <- ggplot(read_counts, aes(x = cell_type, y = nCount_RNA, color = set)) +
	  geom_boxplot(outlier.shape = NA) +  # Hide outliers
	  geom_jitter(width = 0.2, size = 1) +  # Add jitter to show dots
	  scale_color_manual(values = c("A" = red, "B" = blue)) +  # Set colors for A and B
	  theme_minimal() +
	  labs(title = "Total read counts", x = "", y = "Count (nCount_RNA)") + 
	  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Center the title
print(box_plot_reads)

# Save the plot as a PNG file for the current cell type
ggsave(
	  filename = file.path(wd.de.plots, paste0("QC_plot_Reads.png")),
	  plot = box_plot_reads,
	  width = 4,
	  height = 6,
	  dpi = 300
)

# -----------------------------------------------------------------------------
# QC (Expression)
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and sample identities
expression_data <- data.frame(cell_type = so@meta.data$cell.type, nFeature_RNA = so@meta.data$nFeature_RNA, set = so@meta.data$set)
# Set the order of cell types
expression_data$cell_type <- factor(expression_data$cell_type, levels = c("Liver", "Spleen",	"Spermatogonia_PIvap_positive",
																																																																		"Spermatogonia_CD9_positive", "Spermatogonia_TA", 
																																																																		"Spermatocytes", "Spermatids", "Sperm"))
# Generate a box plot for read counts grouped by sample ID
box_plot_reads <- ggplot(expression_data, aes(x = cell_type, y = nFeature_RNA, color = set)) +
	  geom_boxplot(outlier.shape = NA) +  # Hide outliers
	  geom_jitter(width = 0.2, size = 1) +  # Add jitter to show dots
	  scale_color_manual(values = c("A" = red, "B" = blue)) +  # Set colors for A and B
	  theme_minimal() +
	  labs(title = "Number of detected genes", x = "", y = "Number (nFeature_RNA)") + 
	  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Center the title
print(box_plot_reads)

# Save the plot as a PNG file for the current cell type
ggsave(
	  filename = file.path(wd.de.plots, paste0("QC_plot_Expression.png")),
  	plot = box_plot_reads,
	  width = 4,
	  height = 6,
	  dpi = 300
)

# -----------------------------------------------------------------------------
# QC (MT)
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and sample identities
mt_data <- data.frame(cell_type = so@meta.data$cell.type, percent.mt = so@meta.data$percent.mt, set = so@meta.data$set)
# Set the order of cell types
mt_data$cell_type <- factor(mt_data$cell_type, levels = c("Liver", "Spleen",	"Spermatogonia_PIvap_positive",
																																																																		"Spermatogonia_CD9_positive", "Spermatogonia_TA", 
																																																																		"Spermatocytes", "Spermatids", "Sperm"))
# Generate a box plot for read counts grouped by sample ID
box_plot_reads <- ggplot(mt_data, aes(x = cell_type, y = percent.mt, color = set)) +
	  geom_boxplot(outlier.shape = NA) +  # Hide outliers
	  geom_jitter(width = 0.2, size = 1) +  # Add jitter to show dots
	  scale_color_manual(values = c("A" = red, "B" = blue)) +  # Set colors for A and B
	  theme_minimal() +
	  labs(title = "Mitochondrial content", x = "", y = "Percentage (percent.mt)") + 
	  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Center the title
print(box_plot_reads)

# Save the plot as a PNG file for the current cell type
ggsave(
	  filename = file.path(wd.de.plots, paste0("QC_plot_MT.png")),
	  plot = box_plot_reads,
	  width = 4,
	  height = 6,
	  dpi = 300
)

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------
library(ggplot2)
library(patchwork)

# List of unique cell types
cell_types <- unique(so@meta.data$cell.type)

# Loop over each cell type
for (cell_type in cell_types) {
	  # Filter the data for the current cell type
	  cell_data <- so@meta.data[so@meta.data$cell.type == cell_type, ]
	
	  # Part 1: Generate the Box Plot for Read Counts
	  read_counts <- data.frame(sample.id = cell_data$sample.id, nCount_RNA = cell_data$nCount_RNA)
	  box_plot_reads <- ggplot(read_counts, aes(x = sample.id, y = nCount_RNA)) +
		    geom_boxplot() +
		    theme_minimal() +
		    labs(title = "Read counts", x = "", y = "Number") +
		    theme(axis.text.x = element_text(angle = 45, hjust = 1),
						    		plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
	
  	# Part 2: Generate the Box Plot for Global Expression Levels
	  expression_data <- data.frame(sample.id = cell_data$sample.id, nFeature_RNA = cell_data$nFeature_RNA)
	  box_plot_expression <- ggplot(expression_data, aes(x = sample.id, y = nFeature_RNA)) +
		    geom_boxplot() +
		    theme_minimal() +
		    labs(title = "Global expression levels", x = "", y = "Expression") +
		    theme(axis.text.x = element_text(angle = 45, hjust = 1),
						    		plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
	
  	# Part 3: Generate the Box Plot for Mitochondrial Contents
	  mt_data <- data.frame(sample.id = cell_data$sample.id, percent.mt = cell_data$percent.mt)
	  box_plot_mt <- ggplot(mt_data, aes(x = sample.id, y = percent.mt)) +
		    geom_boxplot() +
		    theme_minimal() +
		    labs(title = "Mitochondrial content", x = "Sample ID", y = "Percentage (percent.mt)") +
		    theme(axis.text.x = element_text(angle = 45, hjust = 1),
						    		plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
	
	  # Combine the plots
	  combined_plot <- box_plot_reads / box_plot_expression / box_plot_mt
	
	  # Add the cell type name as a title on top of the combined plot
	  combined_plot <- combined_plot +
		    plot_annotation(title = paste0(gsub("_", " ", cell_type), "\n")) &
		    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Center the title and adjust size
	
	  # Print the combined plot for the current cell type
	  print(combined_plot)
	
	  # Save the plot as a PNG file for the current cell type
	  ggsave(
		    filename = file.path(wd.de.plots, paste0("QC_plot_", gsub(" ", "_", cell_type), ".png")),
		    plot = combined_plot,
		    width = 4,
		    height = 6,
		    dpi = 300
	  )
}

# -----------------------------------------------------------------------------
# All cell types 
# -----------------------------------------------------------------------------
DefaultAssay(so) <- "RNA"

# perform visualization and clustering steps
so.merged <- NormalizeData(so)
so.merged <- FindVariableFeatures(so.merged, selection.method = "vst", nfeatures = 5000)
so.merged
# An object of class Seurat 
# 55367 features across 35 samples within 1 assay 
# Active assay: RNA (55367 features, 5000 variable features)
# 2 layers present: counts, data

all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = all.genes)
#so.merged <- ScaleData(so.merged)
# Centering and scaling data matrix
# |======================================================================| 100%
so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged), npcs = 10)
#so.merged <- RunPCA(so.merged, verbose = FALSE)

pdf(file.path(wd.de.data, "ssc_filtered_normalised_merged_ElbowPlot_ALL.pdf"))
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.merged, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.merged[["pca"]]@stdev / sum(so.merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(so.merged, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_ALL.RData"))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
#resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)
#
#so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 10, verbose = FALSE)
#so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
#so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 10, verbose = FALSE)
#save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP.RData"))

# Find neighbors and clusters
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_ALL.RData"))

so.merged <- FindNeighbors(so.merged, dims = 1:prin_comp, k.param = 20)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = 0.25)

# Optional: Run UMAP for visualization
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_resolution=0.25_ALL.RData"))

##
pdf(file=file.path(wd.de.plots, "ALL_dims=10_resolution=0.25.pdf"))
DimPlot(so.merged, label = TRUE, pt.size = 5)
dev.off()

pdf(file=file.path(wd.de.plots, "ALL_dims=10_resolution=0.25_SampleID.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, "ALL_dims=10_resolution=0.25_CellType.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="cell.type", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

##
library(ggplot2)

# Extract UMAP coordinates and metadata
umap_data <- Embeddings(so.merged, "umap")
metadata <- so.merged@meta.data

# Create a data frame with UMAP coordinates and metadata
umap_df <- data.frame(umap_data, sample.id = metadata$sample.id, cell.type = metadata$cell.type)

# Plot with ggplot2 using DimPlot for the base
pdf(file=file.path(wd.de.plots, "ALL_dims=10_resolution=0.25_CellType_Spleen.pdf"))

# Base DimPlot for cell types
tplot <- DimPlot(so.merged, reduction = "umap", group.by = "cell.type", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5

# Overlay text labels for the "Spermatagonia_TA" cluster
tplot <- tplot + 
	  geom_text(
		    data = subset(umap_df, cell.type == "Spleen"), 
		    aes(x = umap_1, y = umap_2, label = sample.id),
		    size = 3, color = "black", vjust = -0.5, hjust = -0.5
	  )

print(tplot)
dev.off()

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
file <- "Fig.1E_Hermann_undiff_diffA2"

load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_ALL.RData"))
mn7 <- readTable(file.path(wd.de.plots, paste0(file, ".txt")), header=F, rownames=F, sep="\t")
mn7 <- intersect(mn7, rownames(so.merged))
so.integrated <- so.merged

# Define the desired order of cell types
desired_order <- c("Liver", "Spleen", "Spermatogonia_PIvap_positive", 
																			"Spermatogonia_CD9_positive", "Spermatogonia_TA", 
																			"Spermatocytes", "Spermatids", "Sperm")
# Set the order of Idents based on the desired order
Idents(so.integrated) <- factor(so.integrated@meta.data$cell.type, levels = desired_order)

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
writeTable(ordered_genes, file.path(wd.de.plots, paste0(file, ".list")), colnames=F, rownames=F, sep="\t")

# Step 4: Order the genes based on their expression levels from cluster 1 to the last cluster
# Sort each cluster's expression and then combine orders without dropping any genes
#gene_order <- order(rowMeans(average_expression_per_cluster, na.rm = TRUE), decreasing = TRUE)
#ordered_genes <- rownames(scaled_data)[gene_order]
#scaled_data <- scaled_data[gene_order,]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12),
	  						legend.text = element_text(size = 20),  # Increase legend text size
	  						legend.title = element_text(size = 20))  # Increase legend title size
#ggsave(filename = file.path(wd.de.plots, paste0(file, ".pdf")),
#						plot = heatmap, width = 15, height = 45, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
file <- "Fig.1C_Hermann_Mouse Marker Genes"

load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_ALL.RData"))
mn70 <- readTable(file.path(wd.de.plots, paste0(file, ".txt")), header=F, rownames=F, sep="\t")
mn70 <- intersect(mn70, rownames(so.merged))
so.integrated <- so.merged

# Define the desired order of cell types
desired_order <- c("Liver", "Spleen", "Spermatogonia_PIvap_positive", 
																			"Spermatogonia_CD9_positive", "Spermatogonia_TA", 
																			"Spermatocytes", "Spermatids", "Sperm")
# Set the order of Idents based on the desired order
Idents(so.integrated) <- factor(so.integrated@meta.data$cell.type, levels = desired_order)

# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
# Step 2: Extract the scaled data for the genes in mn7
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[mn7, ]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 16))
ggsave(filename = file.path(wd.de.plots, paste0(file, ".pdf")),
						plot = heatmap, width = 10, height = 18, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
file <- "Markers of mouse spermatogenesis_sy2"

load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_ALL.RData"))
mn7 <- readTable(file.path(wd.de.plots, paste0(file, ".txt")), header=F, rownames=F, sep="\t")
mn7 <- intersect(mn7, rownames(so.merged))
so.integrated <- so.merged

# Define the desired order of cell types
desired_order <- c("Liver", "Spleen", "Spermatogonia_PIvap_positive", 
																			"Spermatogonia_CD9_positive", "Spermatogonia_TA", 
																			"Spermatocytes", "Spermatids", "Sperm")
# Set the order of Idents based on the desired order
Idents(so.integrated) <- factor(so.integrated@meta.data$cell.type, levels = desired_order)

# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))
# Step 2: Extract the scaled data for the genes in mn7
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[mn7, ]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 16))
ggsave(filename = file.path(wd.de.plots, paste0(file, ".pdf")),
							plot = heatmap, width = 10, height = 18, dpi = 300)

# -----------------------------------------------------------------------------
# Hierarchical clustering (No dendrogram)
# -----------------------------------------------------------------------------
# Step 1: Perform hierarchical clustering on the genes in mn7
gene_clustering <- hclust(dist(scaled_data))
ordered_genes <- rownames(scaled_data)[gene_clustering$order]  # Get the order of genes from clustering

# Step 2: Reorder the scaled data by the clustered gene order
scaled_data <- scaled_data[ordered_genes, ]

# Step 3: Create the heatmap with the reordered genes
heatmap <- DoHeatmap(so.integrated, features = ordered_genes, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + 
   NoLegend() + 
	  theme(axis.text.y = element_text(size = 17))

# Save the heatmap to a PDF file
ggsave(filename = file.path(wd.de.plots, paste0("Markers of mouse spermatogenesis_sy2_HC_Combined", ".pdf")),
							plot = heatmap, width = 13, height = 25, dpi = 300)

# -----------------------------------------------------------------------------
# Hierarchical clustering (with dendrogram)
# -----------------------------------------------------------------------------
# 1. Re-order genes in `scaled_data` based on hierarchical clustering
ordered_genes <- rownames(scaled_data)[gene_clustering$order]
scaled_data_ordered <- scaled_data[ordered_genes, ]

# 2. Convert dendrogram to data for ggplot using ggdendro
dendrogram_data <- ggdendro::dendro_data(as.dendrogram(gene_clustering))

# Plot the dendrogram with ggdendro
dendrogram_gg <- ggplot() + 
	  geom_segment(data = dendrogram_data$segments, 
			  											aes(x = y, y = x, xend = yend, yend = xend)) +
	  theme_void() + 
  	coord_flip() + 
  	labs(x = NULL, y = NULL)

# 3. Plot the heatmap with reordered genes
heatmap <- DoHeatmap(
	  so.integrated,
	  features = ordered_genes,
	  group.by = "ident",
	  disp.min = -2.5,
	  disp.max = 2.5
) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 16))

# 4. Combine plots
pdf(file.path(wd.de.plots, "Markers of mouse spermatogenesis_sy2_HC_dendrogram.pdf"), width = 12, height = 18)
grid.arrange(dendrogram_gg, heatmap, ncol = 2, widths = c(1, 5))
dev.off()





# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Step 2: Extract the scaled data for the genes in mn7
scaled_data <- GetAssayData(so.integrated, assay = "RNA", slot = "scale.data")[mn7, ]
heatmap_matrix <- as.matrix(scaled_data)

# Perform hierarchical clustering on the genes
gene_clustering <- hclust(dist(scaled_data))

# Create the heatmap with ComplexHeatmap
pdf(file = file.path(wd.de.plots, "Markers of mouse spermatogenesis_sy2_HCluster.pdf"), width = 10, height = 18)
draw(Heatmap(
	  scaled_data,
	     name = "Expression",
	     col = colorRamp2(c(-1, 0, 2), c("magenta", "black", "yellow")),  # Color scale
	     cluster_rows = gene_clustering,
	     cluster_columns = FALSE,
	     show_row_names = TRUE,
	     show_column_names = FALSE,
	     row_names_gp = gpar(fontsize = 10),  # Adjust row name font size as needed
	     heatmap_legend_param = list(title = "Expression", title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 12))
))
dev.off()











# -----------------------------------------------------------------------------
# Germ cells only
# -----------------------------------------------------------------------------
# Identify the clusters you want to keep
Idents(so.merged) <- so.merged@meta.data$cell.type
clusters_to_remove <-c("Liver", "Spleen")
clusters_to_keep <- setdiff(unique(Idents(so.merged)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.merged <- subset(so.merged, idents = clusters_to_keep)
dim(so.merged)

# -----------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
# 02_UMAP
# https://satijalab.org/seurat/articles/multimodal_vignette
# -----------------------------------------------------------------------------
# Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA
DefaultAssay(so.merged) <- "RNA"

# perform visualization and clustering steps
so.merged <- NormalizeData(so.merged)
so.merged <- FindVariableFeatures(so.merged, selection.method = "vst", nfeatures = 5000)
so.merged
# An object of class Seurat 
# 55367 features across 35 samples within 1 assay 
# Active assay: RNA (55367 features, 5000 variable features)
# 2 layers present: counts, data

all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = all.genes)
#so.merged <- ScaleData(so.merged)
# Centering and scaling data matrix
# |======================================================================| 100%
so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged), npcs = 10)
#so.merged <- RunPCA(so.merged, verbose = FALSE)

#pdf(file.path(wd.de.data, "ssc_filtered_normalised_merged_ElbowPlot.pdf"))
#options(repr.plot.width=9, repr.plot.height=6)
#ElbowPlot(so.merged, ndims = 50)
#dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- so.merged[["pca"]]@stdev / sum(so.merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_SSC.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(so.merged, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_SSC.RData"))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
#resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)
#
#so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 10, verbose = FALSE)
#so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
#so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 10, verbose = FALSE)
#save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_SSC_UMAP.RData"))

# Find neighbors and clusters
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_SSC.RData"))

so.merged <- FindNeighbors(so.merged, dims = 1:prin_comp, k.param = 20)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = 0.25)

# Optional: Run UMAP for visualization
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_resolution=0.25_SSC.RData"))

##
pdf(file=file.path(wd.de.plots, "UMAP_dims=10_resolution=0.25.pdf"))
DimPlot(so.merged, label = TRUE, pt.size = 5)
dev.off()

pdf(file=file.path(wd.de.plots, "UMAP_dims=10_resolution=0.25_SampleID.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf(file=file.path(wd.de.plots, "UMAP_dims=10_resolution=0.25_CellType.pdf"))
tplot = DimPlot(so.merged, reduction = "umap", group.by="cell.type", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

##
library(ggplot2)

# Extract UMAP coordinates and metadata
umap_data <- Embeddings(so.merged, "umap")
metadata <- so.merged@meta.data

# Create a data frame with UMAP coordinates and metadata
umap_df <- data.frame(umap_data, sample.id = metadata$sample.id, cell.type = metadata$cell.type)

# Plot with ggplot2 using DimPlot for the base
pdf(file=file.path(wd.de.plots, "UMAP_dims=10_resolution=0.25_CellType_Spermatagonia_TA.pdf"))

# Base DimPlot for cell types
tplot <- DimPlot(so.merged, reduction = "umap", group.by = "cell.type", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5

# Overlay text labels for the "Spermatagonia_TA" cluster
tplot <- tplot + 
	  geom_text(
		    data = subset(umap_df, cell.type == "Spermatogonia_TA"), 
		    aes(x = umap_1, y = umap_2, label = sample.id),
		    size = 3, color = "black", vjust = -0.5, hjust = -0.5
	  )

print(tplot)
dev.off()

# Plot with ggplot2 using DimPlot for the base
pdf(file=file.path(wd.de.plots, "UMAP_dims=10_resolution=0.25_CellType_Spermatogonia_PIvap_positive.pdf"))

# Base DimPlot for cell types
tplot <- DimPlot(so.merged, reduction = "umap", group.by = "cell.type", pt.size = 5)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5

# Overlay text labels for the "Spermatagonia_TA" cluster
tplot <- tplot + 
	geom_text(
		data = subset(umap_df, cell.type == "Spermatogonia_PIvap_positive"), 
		aes(x = umap_1, y = umap_2, label = sample.id),
		size = 3, color = "black", vjust = -0.5, hjust = -0.5
	) +
geom_text(
	data = subset(umap_df, cell.type == "Spermatogonia_CD9_positive"), 
	aes(x = umap_1, y = umap_2, label = sample.id),
	size = 3, color = "black", vjust = -0.5, hjust = -0.5
)

print(tplot)
dev.off()

# -----------------------------------------------------------------------------
# 5.integrate_all_whole_testes_datasets.R
# 03_PCA
# -----------------------------------------------------------------------------
library(Seurat)
library(Matrix)
#library(tidyverse)
library(patchwork)
library(sctransform)
#library(sceasy)
#library(reticulate)
#library(clustree)

# determine cluster markers
# we will later run enrichment analyses on these gene lists to assign a predicted function to each cluster
# IMPORTANT: we must set the DefaultAssay to RNA before running differential expression: 
# "we don't recommend using the integrated matrix for differential expression" and "As a general rule, we always recommend performing DE on originally measured values - instead of on batch-corrected, imputed, etc. values.
# This ensures that the measurements that enter the DE test are indeed independent from each other, which is a requirement of any statistical DE test."
# https://github.com/satijalab/seurat/issues/1057
# https://github.com/satijalab/seurat/issues/1256
# https://github.com/satijalab/seurat/issues/2136
DefaultAssay(so.merged) <- 'RNA'
so.merged <- JoinLayers(so.merged)   # https://www.biostars.org/p/9581468/

Idents(so.merged) <- so.merged@meta.data$cell.type
#so.merged@meta.data$seurat_clusters <- so.merged@meta.data$cell.type

markers <- FindAllMarkers(so.merged, test.use = "wilcox", only.pos=F, min.pct=0.25, logfc.threshold=0.25)
filtered_markers <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
# Count the number of marker genes in each cluster
marker_gene_counts <- filtered_markers %>%
	  dplyr::count(cluster)
print(marker_gene_counts)

write.table(markers, file=file.path(wd.de.data, "so.merged_SSC_markers_wilcox.txt"), row.names=F, col.names=T, quote=F, sep='\t')
save(so.merged, markers, file=file.path(wd.de.data, "so.merged_SSC_markers_wilcox.RData"))

# Reorder the identities of the Seurat object based on your desired order
so.merged$cell.type <- factor(Idents(so.merged), levels = c(
	  "Spermatogonia_PIvap_positive", 
	  "Spermatogonia_CD9_positive", 
	  "Spermatogonia_TA", 
	  "Spermatocytes", 
	  "Spermatids", 
	  "Sperm"
))

# Assign reordered identities to the Seurat object
Idents(so.merged) <- so.merged$cell.type

# Desired order of clusters
desired_order <- c("Spermatogonia_PIvap_positive", 
																			"Spermatogonia_CD9_positive", 
																			"Spermatogonia_TA", 
																			"Spermatocytes", 
																			"Spermatids", 
																			"Sperm")
# Reorder the 'cluster' column in the markers data frame
markers$cluster <- factor(markers$cluster, levels = desired_order)

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
# Filter and select the top genes for each cluster
top10 <- markers %>%
	  group_by(cluster) %>%
	  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
	  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%  # Get top 10 or fewer
	  ungroup()

top10_filtered <- top10 %>%
	  filter(!cluster %in% c("Spermatogonia_PIvap_positive", "Spermatids"))

# Re-scale the data including all genes in top10$gene
so.merged <- ScaleData(so.merged, features = unique(top10_filtered$gene))

# Reassign identities if needed
so.filtered <- subset(so.merged, idents = c("Spermatogonia_CD9_positive", "Spermatogonia_TA", "Spermatocytes", "Sperm"))
#Idents(so.merged) <- so.merged$cell.type  # Ensure this matches top10's cluster column

# Create the heatmap
heatmap <- DoHeatmap(so.filtered, features = top10_filtered$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
  	theme(axis.text.y = element_text(size = 16))
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_filtered_markers_SSC_desired_order.png")),
							plot = heatmap, width = 13, height = 20, dpi = 300)

# -----------------------------------------------------------------------------
# Spermatogonia_PIvap_positive vs. Spermatogonia_TA in Set A
# -----------------------------------------------------------------------------
# Subset for Set A
setA_subset <- subset(so.merged, subset = set == "A" & cell.type %in% c("Spermatogonia_PIvap_positive", "Spermatogonia_TA"))

# Run differential expression
de_results_setA <- FindMarkers(
	  setA_subset,
	  ident.1 = "Spermatogonia_PIvap_positive",
	  ident.2 = "Spermatogonia_TA",
	  test.use = "wilcox"  # Or other test as needed
)

# Add a column to label significant genes for volcano plot
de_results_setA <- de_results_setA %>%
	  mutate(
		    significance = case_when(
			   p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Upregulated",
			   p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Downregulated",
			   TRUE ~ "Not significant"
		 )
	)

# Volcano plot
ggplot(de_results_setA, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
	  geom_point(aes(color = significance), alpha = 0.8) +
	  scale_color_manual(values = c("Upregulated" = red, "Downregulated" = blue, "Not significant" = "gray")) +
	  labs(title = "Volcano Plot: Spermatogonia PIvap+ vs Spermatogonia TA in Set A",
			  			x = "Log2 Fold Change", y = "-log10 Adjusted p-value") +
	  theme_minimal() +
	  theme(plot.title = element_text(hjust = 0.5))

# -----------------------------------------------------------------------------
# Spermatogonia_PIvap_positive vs. Spermatogonia_TA in Set A
# -----------------------------------------------------------------------------
# Subset for Set B
setB_subset <- subset(so.merged, subset = set == "B" & cell.type %in% c("Spermatogonia_CD9_positive", "Spermatogonia_TA"))

# Run differential expression
de_results_setB <- FindMarkers(
	  setB_subset,
	  ident.1 = "Spermatogonia_CD9_positive",
	  ident.2 = "Spermatogonia_TA",
	  test.use = "wilcox"  # Or other test as needed
)

# View significant genes in Set A comparison
head(de_results_setA %>% filter(p_val_adj < 0.05))

# View significant genes in Set B comparison
head(de_results_setB %>% filter(p_val_adj < 0.05))






# -----------------------------------------------------------------------------
# K-means Clustering of Genes
# -----------------------------------------------------------------------------
table(Idents(so.merged))

# Extract expression data from Seurat object
expression_matrix <- GetAssayData(so.merged[["RNA"]], layer = "data")
expression_matrix <- as.matrix(expression_matrix)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
#markers <- FindAllMarkers(so.merged, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
filtered_markers <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
# Count the number of marker genes in each cluster
marker_gene_counts <- filtered_markers %>%
	  dplyr::count(cluster)
print(marker_gene_counts)
# 1 Spermatogonia_PIvap_positive   109
# 2 Spermatogonia_TA               304
# 3 Sperm                         8443
# 4 Spermatogonia_CD9_positive     693
# 5 Spermatocytes                  590
# 6 Spermatids                     123

