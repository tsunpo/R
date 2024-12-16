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
wd.de.data  <- file.path(wd.de, "data_ALL4_500_QC")
wd.de.plots <- file.path(wd.de, "plots_ALL4_500_QC")

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

nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_", nfeatures, ".RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_", nfeatures, "_-1_>200_", 20, "_", 100, ".RData")))

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

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_resolution=0.6_ordered_annotated.pdf")), width = 12, height = 5)
print(dot_plot)
dev.off()

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

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, "_annotated.RData")))

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

save(so.integrated, markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, "_markers_15_17.RData")))
save(markers, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.6_", nfeatures, "_markers_markers_15_17.RData")))

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
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_top10_markers_SCT_", nfeatures, "_UMAP_resolution=0.6_group.by_15_17.png")),
							plot = heatmap, width = 15, height = 17, dpi = 300)

# -----------------------------------------------------------------------------
# Jaccard index for overlap of gene sets
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_markers_markers.RData")))

jaccard_matrix <- getJaccardIndex(markers)
# Visualize the Jaccard index matrix
pheatmap_plot <- pheatmap(jaccard_matrix, main = "Jaccard Index Between Clusters")

# Save the plot as a PNG file
png(file = file.path(wd.de.plots, "heatmap_jaccard_markers_SCT_res=0.6_15_17.png"), width = 7, height = 7, units = "in", res = 300)
print(pheatmap_plot)
dev.off()

# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression data in each cluster
# -----------------------------------------------------------------------------
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
write.table(all_significant_genes, file = "significant_genes_age_correlation_cluster.txt", row.names = FALSE)




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

