# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 25/07/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm22
#wd.src <- "/Users/ty2/Work/dev/R"                ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Graphics.R", "SingleCellTranscriptomics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg38.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"    ## ty2@farm22
#wd <- "/Users/ty2/Work/sanger/ty2"              ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/scRNA")
wd.rna.raw <- file.path(wd.rna, "10x")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-de"))
wd.de.data  <- file.path(wd.de, "data_ALL3_500_SPG")
wd.de.plots <- file.path(wd.de, "plots_ALL3_500_SPG")

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
library(monocle3)
library(pheatmap)

# -----------------------------------------------------------------------------
# To identify six SPG states
# -----------------------------------------------------------------------------
nfeatures <- 5000

wd.de.data0 <- "/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_ALL3_500"
#load(file=file.path(wd.de.data0, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data0, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", nfeatures, ".RData")))

# i.e. retain only cells of states 0-4 (predominantly spermatogonial stem cells to differentiating spermatogonia, although the trajectory ends with some primary spermatocytes)
#so.integrated <- subset(so.integrated, (seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters == 10))
so.integrated <- subset(so.integrated, idents = c(0, 2, 10))
length(unique(so.integrated@meta.data$sample.id))

# Access the metadata
metadata <- so.integrated@meta.data

# Count the number of cells per sample
cell_counts <- metadata %>%
	  group_by(sample.id) %>%
	  summarise(cell_count = n()) %>%
	  ungroup()

# Filter samples with more than 1000 cells
samples_to_keep <- cell_counts %>%
	  filter(cell_count > 100) %>%
	  pull(sample.id)

# Subset the Seurat object to keep only the samples with more than 1000 cells
so.integrated <- subset(so.integrated, cells = rownames(metadata[metadata$sample.id %in% samples_to_keep, ]))

# Check the result
ssc <- as.data.frame(table(so.integrated@meta.data$sample.id))
write.table(ssc, file=file.path(wd.de.data, "ssc_samples_>100.txt"), row.names=T, col.names=T, quote=F, sep='\t')

save(so.integrated, ssc, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_5000_SSC_>100.RData"))

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
write.table(prin_comp, file=file.path(wd.de.data, "ssc_filtered_normalised_integrated_SCT_PCA_5000.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
save(so.integrated, pct, cumu, component1, component2, prin_comp, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", 5000, ".RData")))

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_merged_PCA.RData"))
resolution.range <- seq(from = 0.05, to = 0.5, by = 0.05)

so.integrated <- FindNeighbors(so.integrated, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = resolution.range, verbose = FALSE)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
save(so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", 5000, ".RData")))

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
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
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_res=0.25.png")), plot = monocle3, dpi = 300)

# -----------------------------------------------------------------------------
# Define resolution
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
resolution.range <- seq(from = 0.25, to = 0.5, by = 0.05)

for (r in 1:length(resolution.range)) {
	  res <- resolution.range[r]
	  resolution_name <- paste0("integrated_snn_res.", res)
	  Idents(so.integrated) <- so.integrated[[resolution_name]][[1]]
	
	  ##
	  file.name <- paste0("SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_res=", res)
	  pdf(file=file.path(wd.de.plots, paste0(file.name, ".pdf")))
	  DimPlot(so.integrated, label = TRUE)
	  dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", 5000, ".RData")))

# Find neighbors and clusters
so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = 0.25)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)
save(samples0.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", 5000, ".RData")))

##
file.name <- paste0("SCT_", 5000, "_UMAP_dims=", prin_comp, "_resolution=0.25")

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
# Di Persio et al; DotPlot (resolution = 0.4; Six SPG states)
# -----------------------------------------------------------------------------
#nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, ".RData")))

genes_of_interest <- c("TAF6", "ST3GAL4", "SH2B2", "MSL3", "PHGDH", "C19orf84", "LIN7B", "FSD1", "TSPAN33", "EGR4", "PIWIL4", "CELF4", "UTF1", "FGFR3", "A2M", "ENO3", "SERPINE2", "SRRT", "BAG6", "DND1", "PELP1", "NANOS2", "C1QBP", "NANOS3", "GFRA2", "GFRA1", "ID2", "ASB9", "L1TD1", "ID4", "MKI67", "PDPN", "KIT", "DMRT1", "DNMT1", "CALR", "SYCP3", "STRA8")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", 5000, "_SCT_dims=", prin_comp, "_resolution=0.25.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("6", "3", "9", "7", "4", "8", "0", "5", "2", "1")
new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.25_ordered.pdf")), width = 14, height = 5)
print(dot_plot)
dev.off()

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
cluster_to_celltype <- c('1' = 'Stage 0',
																									'2' = 'Stage 0A',
																									'5' = 'Stage 0B',
																									'0' = 'Stage 1', '8' = 'Stage 1',
																									'4' = 'Stage 2', '7' = 'Stage 2',
																									'3' = 'Stage 3', '9' = 'Stage 3')

# Update the identities using this mapping
#new_order <- c("16", "14", "13", "11", "10", "7", "5", "6", "8", "12", "15", "9", "1", "2", "4", "0", "3")
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

Idents(so.integrated) <- plyr::mapvalues(x = Idents(so.integrated), from = names(cluster_to_celltype), to = cluster_to_celltype)

dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
scale_y_discrete(limits = new_order) # Ensure the new order is used in plotting

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.25_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

##
dim_plot <- DimPlot(so.integrated, label = TRUE)
ggsave(file.path(wd.de.plots, paste0("Di Persio_SCT_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.25_ordered_annotated.png")), plot = dim_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Monocle 3
# -----------------------------------------------------------------------------
library(monocle3)

cds <- getMonocle3CDS(so.integrated, umap_embeddings=T)

# Cluster cells in Monocle 3 using UMAP embeddings
cds <- cluster_cells(cds, reduction_method = "UMAP")
# Perform trajectory analysis using precomputed UMAP
cds <- learn_graph(cds, use_partition = T)
# Order cells in pseudotime
cds <- order_cells(cds)

# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "cluster",
																	label_cell_groups=FALSE,
																	label_leaves=FALSE,
																	label_branch_points=FALSE,
																	graph_label_size=1.5)
ggsave(file.path(wd.de.plots, paste0("pseudotime_", nfeatures, "_UMAP_dims=", prin_comp, "_resolution=0.25_use_partition.png")), plot = monocle3, dpi = 300)

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
# DefaultAssay(so.integrated) <- "SCT"
# Error in UseMethod(generic = "JoinLayers", object = object) : 
#  	no applicable method for 'JoinLayers' applied to an object of class "c('SCTAssay', 'Assay', 'KeyMixin')"
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", nfeatures, ".RData")))

DefaultAssay(so.integrated) <- "RNA"
so.integrated <- JoinLayers(so.integrated, assay = DefaultAssay(so.integrated))
#names(so.integrated[[DefaultAssay(so.integrated)]]@data)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(so.integrated, only.pos = TRUE)
markers %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1)

save(so.integrated, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", nfeatures, "_markers.RData")))

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", nfeatures, "_markers.RData")))

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
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.25.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# Set the new order
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", nfeatures, "_markers.RData")))

new_order <- c("6", "3", "9", "7", "4", "8", "0", "5", "2", "1")
new_cluster_order <- rev(new_order)
so.integrated <- SetIdent(so.integrated, value = factor(Idents(so.integrated), levels = new_cluster_order))

top10 <- markers %>%
	  mutate(cluster = factor(cluster, levels = new_cluster_order)) %>%
	  group_by(cluster) %>%
	  dplyr::filter(avg_log2FC > 1) %>%
	  slice_head(n = 10) %>%
	  ungroup()

# Re-scale the data including all genes in top10$gene
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = top10$gene) + NoLegend()
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.25_new_order_cluster_2.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

# -----------------------------------------------------------------------------
# SingleR
# -----------------------------------------------------------------------------
library(SingleR)
library(celldex)

# Load the Human Primary Cell Atlas data
hpca <- celldex::HumanPrimaryCellAtlasData()

# Extract the normalized data matrix from Seurat object
data_matrix <- GetAssayData(so.integrated, slot = "data")

# Run SingleR
singleR_results <- SingleR(test = data_matrix, ref = hpca, labels = hpca$label.main, clusters = Idents(so.integrated))

# Add SingleR annotations to metadata
so.integrated$SingleR_labels <- singleR_results$labels[as.character(Idents(so.integrated))]

# Set cell type as the default identity class
Idents(so.integrated) <- so.integrated$SingleR_labels

# Visualize clusters with cell type annotations
DefaultAssay(so.integrated) <- "SCT"

plot <- DimPlot(so.integrated, group.by = "SingleR_labels", label = TRUE, repel = TRUE)

pdf(file = file.path(wd.de.plots, paste0("SingleR_DimPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.25_ordered.pdf")))
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# GSEA
# -----------------------------------------------------------------------------
library(clusterProfiler)
library(msigdbr)
library(dplyr)

# Load MSigDB gene sets (Hallmark gene sets as an example)
msigdb <- msigdbr(species = "Homo sapiens", category = "C5")
msigdb_list <- split(msigdb$gene_symbol, msigdb$gs_name)

msigdb_df <- bind_rows(
	  lapply(names(msigdb_list), function(term) {
		    data.frame(term = term, gene = msigdb_list[[term]], stringsAsFactors = FALSE)
	  })
)

# Function to run GSEA for a specific cluster
run_gsea <- function(cluster) {
	  # Filter markers for the specific cluster
	  cluster_markers <- markers %>% filter(cluster == !!cluster)
	  cat("Running GSEA for cluster", cluster, "with", nrow(cluster_markers), "markers\n")
	  
	  # Prepare the gene list for GSEA
	  gene_list <- cluster_markers$avg_log2FC
	  names(gene_list) <- cluster_markers$gene
	  
	  # Ensure the gene list has unique names
	  gene_list <- gene_list[!duplicated(names(gene_list))]
	  
	  # Sort the gene list
	  gene_list <- sort(gene_list, decreasing = TRUE)
	
	  # Filter gene_list to include only genes present in msigdb_list
	  all_genes <- unique(unlist(msigdb_list))
	  filtered_gene_list <- gene_list[names(gene_list) %in% all_genes]
	  
	  # Run GSEA
	  gsea_result <- tryCatch(
	  	  GSEA(gene = filtered_gene_list, TERM2GENE = msigdb_df),
	     error = function(e) {
	  	     warning(paste("GSEA failed for cluster", cluster, "with error:", e$message))
	  	     return(NULL)
	     }
	  )
	  
	  return(gsea_result)
}

# Check if the genes in simple_gene_list are present in the simple_gene_set
all_genes <- unique(unlist(msigdb_list))
print(all_genes)

# Filter simple_gene_list to include only genes present in simple_gene_set
filtered_gene_list <- gene_list[names(gene_list) %in% all_genes]

# Print the filtered gene list
#print(filtered_gene_list)

# Run GSEA for each cluster and store the results in a list
clusters <- unique(markers$cluster)
clusters <- as.numeric(as.vector(clusters))
gsea_results <- lapply(clusters, run_gsea)
names(gsea_results) <- clusters

gsea_results <- gsea_results[!sapply(gsea_results, is.null)]
save(markers, gsea_results, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.25_", nfeatures, "_markers_gsea_results_C5.RData")))

## Function to plot GSEA results for a specific cluster
#plot_gsea <- function(gsea_result, cluster) {
#	  # Check if the result contains any significant enrichment
#	  if (length(gsea_result@result) > 0) {
#	    	p1 <- enrichplot::gseaplot2(gsea_result, geneSetID = 1, title = paste("Cluster", cluster))
#	    	#print(p1)
#		    
#		    ggsave(filename = file.path(wd.de.plots, paste0("GSEA_C8_markers_SCT_", nfeatures, "_UMAP_resolution=0.25_cluster", cluster, ".png")),
#		    							plot = p1, dpi = 300)
#	  } else {
#		    message(paste("No significant enrichment found for cluster", cluster))
#	  }
#}

## Inspect the results for a specific cluster
#for (i in seq_along(gsea_results)) {
#	  if (!is.null(gsea_results[[i]])) {
#		    print(paste("Cluster", i))
#	  	  #print(gsea_results[[i]])
#		    plot_gsea(gsea_results[[i]], i)
#	  } else {
#		    print(paste("No GSEA result for cluster", i))
#	  }
#}

# Plot GSEA results for all clusters
# for (cluster in clusters) {
#   plot_gsea(gsea_results[[cluster]], cluster)
# }

# Extract significant terms from GSEA results
extract_significant_terms <- function(gsea_result, pvalue_cutoff = 0.001) {
	  if (is.null(gsea_result)) return(NULL)
	  result <- gsea_result@result
	  significant_terms <- result[result$p.adjust < pvalue_cutoff, ]
	  return(significant_terms)
}

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

# Generate heatmap of p-values
ht_pvalues <- Heatmap(pvalue_matrix, 
																						name = "p-value", 
																						col = colorRamp2(c(0, 0.05, 1), c(red, "white", blue)), 
																						show_row_names = TRUE, 
																						show_column_names = TRUE,
																						cluster_rows = TRUE, 
																						cluster_columns = TRUE,
																						heatmap_legend_param = list(title = "Adjusted P-Value"))

# Save the heatmap to a file
filename <- file.path(wd.de.plots, paste0("GSEA_C5_Pvalue_SCT_", nfeatures, "_UMAP_resolution=0.25_cluster", cluster, "_p<0.001.pdf"))
pdf(file = filename, width = 10, height = 20)
draw(ht_pvalues)
dev.off()

# Generate heatmap of NES
ht_nes <- Heatmap(nes_matrix, 
																		name = "NES", 
																		col = colorRamp2(c(-3, 0, 3), c(blue, "white", red)), 
																		show_row_names = TRUE, 
																		show_column_names = TRUE,
																		cluster_rows = TRUE, 
																		cluster_columns = TRUE,
																		heatmap_legend_param = list(title = "Normalized Enrichment Score"))

# Save the NES heatmap to a file
filename <- file.path(wd.de.plots, paste0("GSEA_C5_NES_SCT_", nfeatures, "_UMAP_resolution=0.25_cluster", cluster, "_p<0.001.pdf"))
pdf(file = filename, width = 10, height = 20)
draw(ht_nes)
dev.off()









# Function to create data frame for each cluster
create_plot_data <- function(cluster_name, significant_terms) {
	  data <- data.frame(
		    term = significant_terms$Description,
		    NES = significant_terms$NES,
		    pvalue = significant_terms$p.adjust
	  )
	  data$cluster <- cluster_name
	  return(data)
}

# Combine data frames for all clusters
plot_data_list <- mapply(create_plot_data, names(significant_terms_list), significant_terms_list, SIMPLIFY = FALSE)
plot_data <- do.call(rbind, plot_data_list)

# Create bar charts for each cluster
create_bar_chart <- function(cluster_data) {
	  cluster_name <- unique(cluster_data$cluster)
	  p <- ggplot(cluster_data, aes(x = reorder(term, NES), y = NES, fill = pvalue)) +
		 geom_bar(stat = "identity") +
		 scale_fill_gradient(low = "red", high = "blue", name = "P-Value") +
		 coord_flip() +
		 theme_minimal() +
		 labs(title = paste("Cluster", cluster_name, "NES and P-Values"),
							x = "Terms",
							y = "NES") +
		 theme(axis.text.x = element_text(angle = 45, hjust = 1),
								axis.title.x = element_text(size = 12),
								axis.title.y = element_text(size = 12),
								plot.title = element_text(size = 14, face = "bold"))
  	return(p)
}

# Apply the function to create bar charts for each cluster
bar_charts <- lapply(split(plot_data, plot_data$cluster), create_bar_chart)

# Save each bar chart to a PDF file
for (i in seq_along(bar_charts)) {
	  chart <- bar_charts[[i]]

	  filename <- file.path(wd.de.plots, paste0("GSEA_H_BAR_NES_SCT_", nfeatures, "_UMAP_resolution=0.25_cluster", i, ".pdf"))
	  pdf(file = filename)
	  draw(chart)
	  dev.off()
}

# -----------------------------------------------------------------------------
# Jaccard Index for Overlap of Gene Sets
# -----------------------------------------------------------------------------
library("pheatmap")

getJaccardIndex(markers)
# Visualize the Jaccard index matrix
pheatmap_plot <- pheatmap(jaccard_matrix, main = "Jaccard Index Between Clusters", angle_col = 0)

# Save the plot as a PNG file
png(file = file.path(wd.de.plots, "heatmap_jaccard_markers_SCT_res=0.5.png"), width = 7, height = 7, units = "in", res = 300)
print(pheatmap_plot)
dev.off()









# Apply to all GSEA results
significant_terms_list <- lapply(gsea_results, extract_significant_terms)
significant_terms_list <- significant_terms_list[!sapply(significant_terms_list, is.null)]

# Create a matrix of enrichment scores (NES) for significant terms
term_names <- unique(unlist(lapply(significant_terms_list, function(x) x$Description)))
cluster_names <- names(significant_terms_list)
enrichment_matrix <- matrix(NA, nrow = length(term_names), ncol = length(cluster_names))
rownames(enrichment_matrix) <- term_names
colnames(enrichment_matrix) <- cluster_names

for (i in seq_along(significant_terms_list)) {
	  terms <- significant_terms_list[[i]]$Description
	  nes_scores <- significant_terms_list[[i]]$NES
	  enrichment_matrix[terms, names(significant_terms_list)[i]] <- nes_scores
}

library(ComplexHeatmap)
library(circlize)

# Plot the heatmap
p2 <- Heatmap(enrichment_matrix, 
								name = "NES", 
								col = colorRamp2(c(-3, 0, 3), c(blue, "white", red)), 
								show_row_names = TRUE, 
								show_column_names = TRUE,
								cluster_rows = TRUE, 
								cluster_columns = TRUE,
								heatmap_legend_param = list(title = "Normalized Enrichment Score"))

pdf(file = file.path(wd.de.plots, paste0("GSEA_H_markers_SCT_", nfeatures, "_UMAP_resolution=0.25_cluster", cluster, ".pdf")))
draw(p2)
dev.off()


ggsave(filename = file.path(wd.de.plots, paste0("GSEA_H_markers_SCT_", nfeatures, "_UMAP_resolution=0.25_cluster", cluster, ".png")),
							plot = p2, dpi = 300)