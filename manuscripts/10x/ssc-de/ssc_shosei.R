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
read_counts <- get_cluster_factor_levels(so, "nCount_RNA", factor = "sample")
box_plot_reads <- get_box_plot(read_counts, "Read counts", "Number")

# Part 2: Generate the Box Plot for Global Expression Levels
# Extract the total expression level (sum of counts) and cluster identities
expression_data <- get_cluster_factor_levels(so, "nFeature_RNA", factor = "sample")
box_plot_expression <- get_box_plot(expression_data, "Global expression levels", "Expression")

# Part 3: Generate the Box Plot for Mitochondrial Contents
# Extract the total expression level (sum of counts) and cluster identities
mt_data <- get_cluster_factor_levels(so, "percent.mt", factor = "sample")
box_plot_mt <- get_box_plot(mt_data, "Mitochondrial contents", "Percentage")

# Part 4: Combine the Plots
combined_plot <- box_plot_reads / box_plot_expression / box_plot_mt

# Print the combined plot
print(combined_plot)

# -----------------------------------------------------------------------------
# QC (Read counts)
# -----------------------------------------------------------------------------
# Part 1: Generate the Box Plot for Read Counts
# Extract the number of reads (nCount_RNA) and sample identities
read_counts <- data.frame(cell_type = so@meta.data$cell.type, nCount_RNA = so@meta.data$nCount_RNA, set = so@meta.data$set)
# Set the order of cell types
read_counts$cell_type <- factor(read_counts$cell_type, levels = c("Liver", "Spleen", "Spermatogonia_TA", 
																																																																		"Spermatogonia_PIvap_positive", 
																																																																		"Spermatogonia_CD9_positive", 
																																																																		"Spermatocytes", "Spermatids", "Sperm"))
# Generate a box plot for read counts grouped by sample ID
box_plot_reads <- ggplot(read_counts, aes(x = cell_type, y = nCount_RNA, color = set)) +
	  geom_boxplot(outlier.shape = NA) +  # Hide outliers
	  geom_jitter(width = 0.2, size = 1) +  # Add jitter to show dots
	  scale_color_manual(values = c("A" = red, "B" = blue)) +  # Set colors for A and B
	  theme_minimal() +
	  labs(title = "Read counts", x = "", y = "Number") + 
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
expression_data$cell_type <- factor(expression_data$cell_type, levels = c("Liver", "Spleen", "Spermatogonia_TA", 
																																																																		"Spermatogonia_PIvap_positive", 
																																																																		"Spermatogonia_CD9_positive", 
																																																																		"Spermatocytes", "Spermatids", "Sperm"))
# Generate a box plot for read counts grouped by sample ID
box_plot_reads <- ggplot(expression_data, aes(x = cell_type, y = nFeature_RNA, color = set)) +
	  geom_boxplot(outlier.shape = NA) +  # Hide outliers
	  geom_jitter(width = 0.2, size = 1) +  # Add jitter to show dots
	  scale_color_manual(values = c("A" = red, "B" = blue)) +  # Set colors for A and B
	  theme_minimal() +
	  labs(title = "Global expression levels", x = "", y = "Expression") + 
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
mt_data$cell_type <- factor(mt_data$cell_type, levels = c("Liver", "Spleen", "Spermatogonia_TA", 
																																																																										"Spermatogonia_PIvap_positive", 
																																																																										"Spermatogonia_CD9_positive", 
																																																																										"Spermatocytes", "Spermatids", "Sperm"))
# Generate a box plot for read counts grouped by sample ID
box_plot_reads <- ggplot(mt_data, aes(x = cell_type, y = percent.mt, color = set)) +
	  geom_boxplot(outlier.shape = NA) +  # Hide outliers
	  geom_jitter(width = 0.2, size = 1) +  # Add jitter to show dots
	  scale_color_manual(values = c("A" = red, "B" = blue)) +  # Set colors for A and B
	  theme_minimal() +
	  labs(title = "Mitochondrial content", x = "", y = "Percentage") + 
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
		    labs(title = "Mitochondrial content", x = "Sample ID", y = "Percentage") +
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
resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)

so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 10, verbose = FALSE)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 10, verbose = FALSE)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP.RData"))

# Find neighbors and clusters
load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_ALL.RData"))

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
# Germ cells only
# -----------------------------------------------------------------------------
# Identify the clusters you want to keep
Idents(so) <- so@meta.data$cell.type
clusters_to_remove <-c("Liver", "Spleen")
clusters_to_keep <- setdiff(unique(Idents(so)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 1
so.merged <- subset(so, idents = clusters_to_keep)
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
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(so.merged, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)

so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 10, verbose = FALSE)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 10, verbose = FALSE)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP.RData"))

# Find neighbors and clusters
load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))

so.merged <- FindNeighbors(so.merged, dims = 1:prin_comp, k.param = 20)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = 0.25)

# Optional: Run UMAP for visualization
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_resolution=0.25.RData"))

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
library(tidyverse)
library(patchwork)
library(sctransform)
library(sceasy)
library(reticulate)
library(clustree)

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

markers <- FindAllMarkers(so.merged, test.use="MAST", only.pos=T, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, file=file.path(wd.de.data, "so.merged_nfeatures=5000_PCA_UMAP_k=10_n=10_JoinLayers_markers_SSC.txt"), row.names=F, col.names=T, quote=F, sep='\t')
save(so.merged, markers, file=file.path(wd.de.data, "so.merged_nfeatures=5000_PCA_UMAP_k=10_n=10_JoinLayers_markers_SSC.RData"))

# Reorder the identities of the Seurat object based on your desired order
so.merged$cell.type <- factor(Idents(so.merged), levels = c(
	  "Spermatogonia_TA", 
	  "Spermatogonia_PIvap_positive", 
	  "Spermatogonia_CD9_positive", 
	  "Spermatocytes", 
	  "Spermatids", 
	  "Sperm"
))

# Assign reordered identities to the Seurat object
Idents(so.merged) <- so.merged$cell.type

# Desired order of clusters
desired_order <- c("Spermatogonia_TA", 
																			"Spermatogonia_PIvap_positive", 
																			"Spermatogonia_CD9_positive", 
																			"Spermatocytes", 
																			"Spermatids", 
																			"Sperm")
# Reorder the 'cluster' column in the markers data frame
markers$cluster <- factor(markers$cluster, levels = desired_order)

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()
#missing_genes <- setdiff(top10$gene, rownames(so.integrated))

# Re-scale the data including all genes in top10$gene
so.merged <- ScaleData(so.merged, features = rownames(so.merged))

# Create the heatmap
heatmap <- DoHeatmap(so.merged, features = top10$gene, group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
  	theme(axis.text.y = element_text(size = 16))
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_NEW.png")),
							plot = heatmap, width = 13, height = 20, dpi = 300)

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
















top.markers <- subset(markers, cluster == "0")
lfcs <- getMarkerEffects(top.markers)

###
##
DefaultAssay(so.merged) <- 'RNA'
so.merged <- JoinLayers(so.merged)   # https://www.biostars.org/p/9581468/

so.merged <- FindVariableFeatures(so.merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = all.genes)

markers.0 <- FindMarkers(so.merged, ident.1=0)
head(markers.0, n=5)

markers <- FindAllMarkers(so.merged, only.pos=T)
write.table(markers, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20_markers_.txt"), row.names=F, col.names=T, quote=F, sep='\t')
save(markers, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA_UMAP_k=20_n=20_markers_.RData"))

markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

VlnPlot(pbmc, features = c("PRM1", "PRM2"))



# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)





##
pdf("DimPlot_UMAP_RNA_dim=18_k=10.pdf")
DimPlot(so.merged, label = TRUE)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=18_k=10_sampleID.pdf")
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=18_k=10_age.pdf")
tplot = DimPlot(so.merged, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf("DimPlot_UMAP_RNA_dim=18_k=10_2n.pdf")
tplot = DimPlot(so.merged, reduction = "umap", group.by="n2")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()




























# -----------------------------------------------------------------------------
# Standard Seurat pre-processing workflow (SCT)
# 01_QC
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "ssc_filtered.RData"))
library(glmGamPoi)

samples0.filtered <- samples0[subset(filtered, cells > 1000)$PD_ID,]
so.list <- c()
ids = c()
genes <- c()

colnames <- c("PD_ID", "genes", "cells")
normalised <- toTable(0, length(colnames), nrow(samples0.filtered), colnames)
normalised$PD_ID <- rownames(samples0.filtered)
rownames(normalised) <- rownames(samples0.filtered)

for (s in 1:nrow(samples0.filtered)) {
	# Initialize the Seurat object with the raw (non-normalized data)
	# https://satijalab.org/seurat/articles/pbmc3k_tutorial
	data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0.filtered$V1[s], "filtered_feature_bc_matrix"))
	so <- CreateSeuratObject(counts=data, project=samples0.filtered$V3[s], min.cells=3, min.features=200)
	
	# QC and selecting cells for further analysis
	so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^MT-")
	so <- subset(so, subset=nFeature_RNA > 1000 & nFeature_RNA < 10000 & nCount_RNA > 2000 & nCount_RNA < 50000 & percent.mt < 5)
	
	# Apply sctransform normalization
	# https://satijalab.org/seurat/articles/sctransform_vignette.html
	so <- SCTransform(so, vars.to.regress="percent.mt", verbose=F)
	
	normalised[s, 2] <- nrow(so)
	normalised[s, 3] <- ncol(so)
	
	so.list <- c(so.list, so)
	ids = c(ids, samples0.filtered$V3[s])
	
	if (length(genes) != 0) {
		genes <- intersect(genes, rownames(so))
	} else {
		genes <- rownames(so)
	}
}
writeTable(normalised, file.path(wd.de.data, "ssc_filtered_normalised.txt"), colnames=T, rownames=F, sep="\t")
save(filtered, normalised, samples0, samples0.filtered, so.list, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))

# Merge Based on Normalized Data
# https://satijalab.org/seurat/archive/v4.3/merge#:~:text=Merge%20Based%20on%20Normalized%20Data,data%20%3D%20TRUE%20
so.merged <- merge(x=so.list[[1]], y=so.list[-1], add.cell.ids=ids, project="SSC", merge.data=T)
save(filtered, normalised, samples0, samples0.filtered, so.merged, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised_merged.RData"))


so.merged <- FindNeighbors(so.merged, dims = 1:30)
# Computing nearest neighbor graph
# Computing SNN

so.merged <- FindClusters(so.merged, resolution = 0.8, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:30)
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# 15:15:04 UMAP embedding parameters a = 0.9922 b = 1.112
# 15:15:04 Read 46987 rows and found 30 numeric columns
# 15:15:04 Using Annoy for neighbor search, n_neighbors = 30
# 15:15:04 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 15:15:15 Writing NN index file to temp file /tmp/RtmpFKq9At/file1e6e714df04e1
# 15:15:15 Searching Annoy index using 1 thread, search_k = 3000
# 15:15:46 Annoy recall = 100%
# 15:15:46 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 15:15:49 Initializing from normalized Laplacian + noise (using RSpectra)
# 15:15:58 Commencing optimization for 200 epochs, with 2067094 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 15:16:24 Optimization finished
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_RNA_PCA_UMAP.RData"))

#pdf('DimPlot_UMAP_RNA_dim=30.pdf')
#DimPlot(so.merged, label = TRUE)
#dev.off()

pdf('DimPlot_UMAP_RNA_dim=30_.pdf')
tplot = DimPlot(so.merged, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('DimPlot_UMAP_RNA_dim=30_grouped_by_sampleID.pdf')
tplot = DimPlot(so.merged, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('DimPlot_UMAP_RNA_dim=30_grouped_by_age.pdf')
tplot = DimPlot(so.merged, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()




# -----------------------------------------------------------------------------
# Perform dimensionality reduction by PCA and UMAP embedding
# 02_UMAP
# https://satijalab.org/seurat/articles/sctransform_vignette.html
# -----------------------------------------------------------------------------
Assays(so.merged)
# [1] "RNA" "SCT"
DefaultAssay(so.merged)
# [1] "SCT"

#features <- FindVariableFeatures(so.merged)     # Error: SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-
features <- SelectIntegrationFeatures(so.list)   # https://github.com/satijalab/seurat/issues/5761
so.merged <- RunPCA(so.merged, features = features, verbose = FALSE)
#so.merged <- RunUMAP(so.merged, dims = 1:20)
#DimPlot(so.merged, label = TRUE)

pdf('human_adult_firstPass.elbow_plot.pdf')
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(so.merged, ndims = 50)
dev.off()

# https://hbctraining.github.io/scRNA-seq_online/lessons/03_SC_quality_control-setup.html
head(so.merged@meta.data)
save(so.merged, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_SCT_PCA_UMAP.RData"))

pct <- so.merged[["pca"]]@stdev / sum(so.merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data
prin_comp <- min(component1, component2)
write.table(prin_comp, file='human_adult_firstPass.elbow_PC.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
resolution.range <- seq(from = 0, to = 0.5, by = 0.05)

so.merged <- FindNeighbors(so.merged, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.merged <- FindClusters(so.merged, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)

pdf('DimPlot_UMAP.pdf')
DimPlot(so.merged, label = TRUE)
dev.off()









so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged))









# -----------------------------------------------------------------------------
# Perform dimensionality reduction by PCA and UMAP embedding
# 02_UMAP
# https://satijalab.org/seurat/articles/sctransform_vignette.html
# -----------------------------------------------------------------------------
features <- SelectIntegrationFeatures(object.list = so.list, normalization.method = "SCT", nfeatures = 5000)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)

so.merged <- RunPCA(so.merged, features = features)



so.anchors <- FindIntegrationAnchors(object.list = so.list, anchor.features = features, reduction = "rpca") # we create as our standard reference the Di Persio, Sohni and Zhao samples as these contain the greatest number of cells
so.combined <- IntegrateData(anchorset = so.anchors, normalization.method = "SCT") # see https://github.com/satijalab/seurat/issues/3930 for discussion of k.weight









features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = 5000)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)
testis.anchors <- FindIntegrationAnchors(object.list = so.list, anchor.features = features, normalization.method = "SCT", reference=c(1,2,3,7,8,10,11,12), reduction = "rpca") # we create as our standard reference the Di Persio, Sohni and Zhao samples as these contain the greatest number of cells
testis.combined <- IntegrateData(anchorset = testis.anchors, normalization.method = "SCT") # see https://github.com/satijalab/seurat/issues/3930 for discussion of k.weight



# Centering and scaling data matrix
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#identification-of-highly-variable-features-feature-selection
#all.genes <- rownames(so.merged)
so.merged <- ScaleData(so.merged, features = genes)   ## Only using 23246 genes shared across all samples
so.merged <- RunPCA(so.merged, features = VariableFeatures(object = so.merged))
so.merged <- RunUMAP(so.merged, dims = 1:30, verbose = FALSE)

save(filtered, normalised, samples0, samples0.filtered, so.merged, ids, genes, file=file.path(wd.de.data, "ssc_filtered_normalised_merged_scaled_23246.RData"))



# These are now standard steps in the Seurat workflow for visualization and clustering
so.merged <- RunPCA(so.merged, verbose = FALSE)
so.merged <- RunUMAP(so.merged, dims = 1:30, verbose = FALSE)

so.merged <- FindNeighbors(so.merged, dims = 1:30, verbose = FALSE)
so.merged <- FindClusters(so.merged, verbose = FALSE)
DimPlot(so.merged, label = TRUE)





for (s in 1:nrow(samples0.filtered)) {
	  nrow(so.list[[1]])
}





	  # Normalizing the data
	  so <- NormalizeData(so)
	  so <- GetAssayData(object=so, layer="counts")
	  writeTable(as.data.frame(so), gzfile(file.path(wd.de.data, "GRCh38-2020-A", paste0(samples0$V3[s], ".txt.gz"))), colnames=T, rownames=F, sep="\t")
	
	  # Scaling the data (required for PCA)
	  all.genes <- rownames(so)
	  so <- ScaleData(so, features=all.genes)
	
	 # Perform linear dimensional reduction (required for UMAP)
	 #ssc0.pca <- RunPCA(ssc0, features=VariableFeatures(object=ssc0))
	
	# Run non-linear dimensional reduction (UMAP/tSNE)
	#ssc0.umap <- RunUMAP(ssc0.pca, dims=1:10)
	
	#file.name <- file.path(wd.de.plots, "02_UMAP", paste0(samples0$V3[s], "_DimPlot"))
	#pdf(paste0(file.name, ".pdf"))
	#DimPlot(ssc0.umap, reduction="umap")
	#dev.off()
}
filtered <- raw