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
wd.de.data  <- file.path(wd.de, paste0("data_ALL4_500"))
wd.de.plots <- file.path(wd.de, paste0("plots_ALL4_500"))

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

load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_500",           "ssc_filtered_normalised.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500",      "ssc_filtered_normalised.1.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_500",      "ssc_filtered_normalised.2.RData"))
load(file=file.path("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_lm26_June_500", "ssc_filtered_normalised.3.RData"))

so.list <- c(so.list, so.list.1, so.list.2, so.list.3)
ids <- c(ids, ids.1, ids.2, ids.3)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.1)
samples0.filtered <- rbind(samples0.filtered, samples0.2)
samples0.filtered <- rbind(samples0.filtered, samples0.filtered.3)

#save(samples0.filtered, ids, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list.RData")))

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
library(DoubletFinder)

classification_cols <- list()
for (i in seq_along(so.list)) {
	  sample <- so.list[[i]]
	  raw_counts <- GetAssayData(sample, layer = "counts", assay = "RNA")
	  sample[["RNA"]] <- CreateAssayObject(counts = raw_counts)
	  DefaultAssay(sample) <- "RNA"
	  
	  # Preprocessing
	  sample <- NormalizeData(sample, assay = "RNA")
	  sample <- FindVariableFeatures(sample)
	  sample <- ScaleData(sample)
	  sample <- RunPCA(sample)
	
	  # Perform clustering to generate meaningful identities for modelHomotypic
	  sample <- FindNeighbors(sample, dims = 1:10)
	  sample <- FindClusters(sample, resolution = 0.5)  # Adjust resolution as needed
	  
	  # Find optimal pK
	  sweep_res <- paramSweep_v3(sample, PCs = 1:10, sct = FALSE)
	  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
	  # Save the plot
	  png(file.path(wd.de.data, paste0("pK_optimization_", samples0.filtered$V3[i], ".png")), width = 800, height = 600)
	  optimal_pK <- find.pK(sweep_stats)$pK
	  dev.off()
	  optimal_pK <- as.numeric(as.character(optimal_pK))
	  
	  # Estimate homotypic doublet proportion
	  homotypic_prop <- modelHomotypic(Idents(sample))
	  nExp <- round(0.075 * ncol(sample))  # Adjust expected doublet rate
	  nExp_adj <- round(nExp * (1 - homotypic_prop))
	
	  # Run DoubletFinder
	  sample <- doubletFinder_v3(
		    sample,
		    PCs = 1:10,
		    pN = 0.25,
		    pK = optimal_pK,
		    nExp = nExp_adj,
		    reuse.pANN = FALSE,
		    sct = FALSE
	  )
	  
	  # Extract pK corresponding to the maximum BCmetric
	  pK_stats <- find.pK(sweep_stats)
	  optimal_pK <- pK_stats$pK[which.max(pK_stats$BCmetric)]
	  classification_col <- paste0("DF.classifications_0.25_", optimal_pK, "_", nExp_adj)
	  Idents(sample) <- sample@meta.data[[classification_col]]
	  
	  # Store sample ID and classification column
	  classification_cols[[samples0.filtered$V3[i]]] <- classification_col
	  
	  so.list[[i]] <- sample
}

save(samples0.filtered, ids, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.RData")))

# -----------------------------------------------------------------------------
# DoubletFinder
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.RData")))

# Create a mapping of sample names to age
age_mapping <- setNames(samples0.filtered$V4, rownames(samples0.filtered))

# Initialize a data frame to store the results
gene_cell_table <- data.frame(
	  Sample = character(),
	  Age = numeric(),
	  Genes = numeric(),
	  Cells = numeric(),
	  Cells_DF = numeric()
)

# Loop through the Seurat objects in so.list
for (i in seq_along(so.list)) {
	  # Get the Seurat object
	  sample <- so.list[[i]]
	
	  # Extract sample name (assume rownames(samples0.filtered) match so.list identifiers)
	  sample_name <- samples0.filtered$V3[i]
	
	  # Add the information to the table
	  gene_cell_table <- rbind(gene_cell_table, data.frame(
		    Sample = sample_name,
		    Age = as.numeric(age_mapping[sample_name]),
		    Genes = nrow(sample),  # Number of features
		    Cells = ncol(sample),  # Number of cells
		    Cells_DF = ncol(subset(sample, idents = "Singlet"))
	  ))
}

# Display the table
print(gene_cell_table)
writeTable(gene_cell_table, file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF.txt")), colnames=T, rownames=F, sep="\t")

for (i in seq_along(so.list)) {
	  sample <- so.list[[i]]
	  
	  # Subset to remove doublets
	  sample <- subset(sample, idents = "Singlet")

	  # Update the sample in so.list
	  so.list[[i]] <- sample
}

# -----------------------------------------------------------------------------
# Remove PD53622b_2N and PD53622b_M1
# -----------------------------------------------------------------------------
# Remove items 2 and 8 from the list (where there are only 130 and 92 cells in the samples)
so.list <- so.list[-c(2, 8)]
# Remove rows 2 and 8 from the data frame
samples0.filtered <- samples0.filtered[-c(2, 8), ]
# Remove elements 2 and 8 from the array
ids <- ids[-c(2, 8)]

save(samples0.filtered, ids, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF_-2-8.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_so.list_DF_-2-8.RData")))

# -----------------------------------------------------------------------------
# Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# -----------------------------------------------------------------------------
nfeatures <- 5000
#res <- 0.5

## Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
so.list <- lapply(X = so.list, FUN = SCTransform)

## As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
## Run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = so.list, nfeatures = nfeatures)
so.list <- PrepSCTIntegration(object.list = so.list, anchor.features = features)
save(samples0.filtered, ids, features, so.list, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated0_DF_SCT_", nfeatures, "_so.list.RData")))

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
# Last Modified: 14/02/25
# -----------------------------------------------------------------------------
nfeatures <- 5000
res <- 0.5
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_", nfeatures, "_-2_", 25, "_", 100, ".RData")))

load("/lustre/scratch127/casm/team294rr/ty2/SSC/analysis/expression/ssc-de/data_ALL4_500/ssc_filtered_normalised_integrated0_so.list_DF_classification_cols.RData")
so.list <- so.list[-c(2, 8)]
classification_cols <- classification_cols[-c(2, 8)]

# Check where "PD53621b_2N" appears in sample_ids
#sample_index <- which(sample_ids == "PD53621b_2N")
# Find the corresponding classification column using the same index
#classification_col <- classification_cols[sample_index]
#pANN_col <- gsub("DF.classifications_", "pANN_", classification_col)

# Subset the sample
#sample_cells <- so.integrated@meta.data[so.integrated@meta.data$orig.ident == "PD53621b_2N", ]
# Get classification values for these 10 cells
#original_scores <- sample_cells[1:10, pANN_col]

# Initialize a dataframe to store mapping information
mapping_table <- data.frame(Sample_ID = character(), pANN_Column = character(), stringsAsFactors = FALSE)
sample_ids <- unique(so.integrated@meta.data$sample.id)

# Iterate through each sample
for (i in seq_along(sample_ids)) {
	  sample_id <- sample_ids[i]
	
	  # Get the corresponding classification column
	  classification_col <- classification_cols[i]
	  pANN_col <- gsub("DF.classifications_", "pANN_", classification_col)
	
	  # Append to mapping table
	  mapping_table <- rbind(mapping_table, data.frame(Sample_ID = sample_id, pANN_Column = pANN_col))
}
# Print the mapping table
print(mapping_table)

# -----------------------------------------------------------------------------
# Last Modified: 03/03/25
# -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)

# Step 1: Assign the correct pANN values to each cell in `so.integrated`
so.integrated@meta.data$pANN_score <- NA  # Initialize a column for pANN scores

# Iterate through each sample to assign the correct pANN column values
for (i in seq_along(mapping_table$Sample_ID)) {
	  sample_id <- mapping_table$Sample_ID[i]
	  pANN_col <- mapping_table$pANN_Column[i]
	
	  # Subset cells from the corresponding sample and assign the correct pANN values
	  so.integrated@meta.data$pANN_score[so.integrated@meta.data$sample.id == sample_id] <- 
		    so.integrated@meta.data[so.integrated@meta.data$sample.id == sample_id, pANN_col]
}

# Step 2: Generate a UMAP plot colored by pANN values using FeaturePlot()
umap_plot <- FeaturePlot(so.integrated, features = "pANN_score", reduction = "umap") + NoLegend() +
  	scale_color_viridis_c(option = "magma", name = "pANN score") +  # Use a color scale that emphasizes high values
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.position = c(0.05, 0.95),  # Move legend to top-left corner
		    legend.justification = c(0, 1),  # Anchor legend at top-left
		    legend.key.size = unit(0.5, "cm"),  # Adjust legend size
		    legend.text = element_text(size = 12)  # Adjust legend text size
	  )

# Save the UMAP plot
ggsave(file.path(wd.de.plots, "UMAP_pANN_colored_legend_.png"), plot = umap_plot, width = 6, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Top sample per cluster
# Last Modified: 03/03/25
# -----------------------------------------------------------------------------
library(dplyr)

# Get Seurat metadata containing cluster and sample ID information
metadata <- so.integrated@meta.data

# Count the number of cells per (Cluster, Sample_ID) combination
cluster_sample_counts <- metadata %>%
	  group_by(seurat_clusters, sample.id) %>%
	  summarise(Cell_Count = n(), .groups = "drop")

# Calculate the total number of cells in each cluster
total_cells_per_cluster <- cluster_sample_counts %>%
	  group_by(seurat_clusters) %>%
	  summarise(Total_Cells = sum(Cell_Count))

# Merge the total cell counts into the cluster-sample counts
cluster_sample_counts <- merge(cluster_sample_counts, total_cells_per_cluster, by = "seurat_clusters")

# Calculate the percentage of each sample within the cluster
cluster_sample_counts <- cluster_sample_counts %>%
	  mutate(Percentage = (Cell_Count / Total_Cells) * 100)

# Find the sample with the highest percentage for each cluster
top_sample_per_cluster <- cluster_sample_counts %>%
	  group_by(seurat_clusters) %>%
	  slice_max(order_by = Percentage, n = 1) %>%
	  arrange(desc(Percentage))  # Sort clusters numerically

# Rename columns for clarity
colnames(top_sample_per_cluster) <- c("Cluster", "Top_Sample_ID", "Cell_Count", "Total_Cells", "Percentage")

# Print the summary table
print(top_sample_per_cluster)

# -----------------------------------------------------------------------------
# Top sample per cluster
# Last Modified: 03/03/25
# -----------------------------------------------------------------------------
# Identify columns to remove
columns_to_remove <- grep("^pANN_0.25|^DF.classifications_0.25|^integrated_snn_res", colnames(so.integrated@meta.data), value = TRUE)

# Remove the identified columns
so.integrated@meta.data <- so.integrated@meta.data[, !colnames(so.integrated@meta.data) %in% columns_to_remove]

# Verify removal
print(colnames(so.integrated@meta.data))

# Identify the clusters you want to keep
clusters_to_remove <- c(23)
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
# Subset the Seurat object to exclude cluster 23
so.integrated <- subset(so.integrated, idents = clusters_to_keep)

# Save the updated Seurat object (optional)
#saveRDS(so.integrated, file = file.path(wd.de.data, "ssc_filtered_normalised_integrated0_so.list_cleaned.RDS"))
save(so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_", nfeatures, "_-2_", 25, "_", 100, "-23.RData")))











##
##
library(ggplot2)

# Initialize dataframe to store all scores
density_data <- data.frame()
# Get unique cluster identities
clusters <- unique(so.integrated@meta.data$seurat_clusters)

# Iterate through each cluster
for (cluster in clusters) {
	  # Subset the Seurat object to get cells in this cluster
	  cluster_cells <- subset(so.integrated, idents = cluster)
	  # Get unique sample IDs in this cluster
	  sample_ids_in_cluster <- unique(cluster_cells@meta.data$sample.id)
	
  	# Iterate through each sample in this cluster
	  for (sample_id in sample_ids_in_cluster) {
		    # Lookup pANN column from mapping_table
		    pANN_col <- mapping_table$pANN_Column[mapping_table$Sample_ID == sample_id]
		
		    # Extract cells belonging to this sample
		    sample_cells <- cluster_cells@meta.data[cluster_cells@meta.data$sample.id == sample_id, ]
		
		    # Extract original pANN scores for this sample
		    if (pANN_col %in% colnames(sample_cells)) {  # Ensure column exists
			      original_scores <- sample_cells[[pANN_col]]
			
			      # Store in dataframe
			      density_data <- rbind(density_data, data.frame(
				        Cluster = cluster,
				        Sample = sample_id,
				        Doublet_Score = original_scores
			      ))
		    }
	  }
}

library(ggplot2)

# Define output directory for individual plots
output_dir <- file.path(wd.de.plots, "density_plots_per_cluster")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)  # Create directory if it doesn't exist

# Generate a color palette with enough colors
num_samples <- length(sample_ids)
sample_colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_samples)  # Expand Set1 to required number

# Assign sample colors as a named vector
names(sample_colors) <- sample_ids
# Iterate through each unique cluster and generate a separate plot
for (cluster in unique(density_data$Cluster)) {
	  # Subset data for the current cluster
  	cluster_data <- subset(density_data, Cluster == cluster)
	
  	# Generate density plot
  	p <- ggplot(cluster_data, aes(x = Doublet_Score, color = Sample)) +
		    geom_density(alpha = 1, size = 1.1) +  # Increase line thickness
		    theme_classic() +
		    xlab("Doublet probability score (pANN)") +
		    ylab("Density") +
		    ggtitle(paste("Cluster", cluster)) +
		    scale_color_manual(values = sample_colors) +  # Assign colors to samples
		    theme(
			      legend.title = element_blank(),  # Remove legend title
			      legend.text = element_text(size = 12),
			      legend.position = "right",  # Move legend to the right of the plot
			      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center title
		    ) +
		    guides(color = guide_legend(ncol = 1))  # Force legend to one line
	
	     # Save each plot as a PNG file
	     ggsave(filename = file.path(output_dir, paste0("density_plot_cluster_", cluster, ".png")),
			 					plot = p, width = 6, height = 6, dpi = 300)
}






