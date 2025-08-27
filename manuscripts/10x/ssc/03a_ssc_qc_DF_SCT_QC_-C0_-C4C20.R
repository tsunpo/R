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

load(file=file.path(wd.de.data, "ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0.RData"))

# -----------------------------------------------------------------------------
# Remove cluster 0 from the data
# -----------------------------------------------------------------------------
library(dplyr)

load(file=file.path(wd.de.data, "ssc_filtered_normalised.RData"))
nfeatures <- 5000
res <- 0.5

# using seurat_clusters
so.integrated <- subset(so.integrated, subset = seurat_clusters != 0)
so.integrated <- subset(so.integrated, subset = seurat_clusters != 4)
so.integrated <- subset(so.integrated, subset = seurat_clusters != 20)
table(Idents(so.integrated))          # if identities are clusters
table(so.integrated$seurat_clusters)   # if using meta column

# Extract metadata from so.integrated object
metadata <- so.integrated@meta.data

# Summarize the data by sample ID
summary_table <- metadata %>%
   group_by(sample.id, age) %>%
   summarise(
      Reads = sum(nFeature_RNA),       # Total reads
      Cells = n()                      # Total number of cells
   ) %>%
   ungroup()

# Convert to a data frame and print
summary_table <- as.data.frame(summary_table)

rownames(summary_table) <- summary_table$sample.id
qced <- summary_table[rownames(normalised),]
qced <- rbind(qced, summary_table["VL00297",])
print(qced)

# Optionally save to a CSV file
write.table(qced, file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_-C0_-C4C20.txt")), col.names=T, row.names=F, sep="\t")

# -----------------------------------------------------------------------------
# 
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
# work on the integrated assay
DefaultAssay(so.integrated) <- "integrated"

# re-compute after subsetting
so.integrated <- ScaleData(so.integrated, verbose = FALSE)
so.integrated <- RunPCA(so.integrated, npcs = prin_comp, verbose = FALSE)

# then redo graph-based steps & embeddings
so.integrated <- FindNeighbors(so.integrated, dims = 1:prin_comp, k.param = 20)
so.integrated <- FindClusters(so.integrated, algorithm=3, resolution = res)
so.integrated <- RunUMAP(so.integrated, dims = 1:prin_comp, n.neighbors = 20)

# Flip UMAP_1 and UMAP_2
so.integrated@reductions$umap@cell.embeddings[,1] <- -so.integrated@reductions$umap@cell.embeddings[,1]  # flip UMAP_1
so.integrated@reductions$umap@cell.embeddings[,2] <- -so.integrated@reductions$umap@cell.embeddings[,2]  # flip UMAP_2

save(prin_comp, samples.filtered, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_", nfeatures, "_", 25, "_", 100, "_integrated_PCA_UMAP_", prin_comp, "_0.5_-C0_-C4C20.RData")))

##
file.name <- paste0("UMAP_SCT_", nfeatures, "_", 25, "_", 100, "_", prin_comp, "_0.5_-C0_-C4C20")

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
#load(file=file.path(wd.de.data, "ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5-C0.RData"))

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
