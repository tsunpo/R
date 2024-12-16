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

# -----------------------------------------------------------------------------
# Methods: Age group
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated.RData")))

# Subset the Seurat object to exclude "2N" and "4N" samples
so.integrated_no_2N_4N <- subset(
	  so.integrated, 
	  cells = rownames(so.integrated@meta.data[!grepl("2N|4N", so.integrated@meta.data$sample.id), ])
)
# Verify the result
unique(so.integrated_no_2N_4N@meta.data$sample.id)
so.integrated <- so.integrated_no_2N_4N

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
write.table(proportions, "cell_type_proportions_by_age_group_M.txt", row.names = FALSE)

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
png(file = file.path(wd.de.plots, "barplot_age_group_cell_type_DimPlot_M.png"), width = 6, height = 7, , units = "in", res = 300)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Methods: Correlation 1, 2, 3, 4, 5, 6 OR 25, 35, 45, 55, 65, 75
# Last Modified: 14/11/24
# -----------------------------------------------------------------------------
# Step 2: Perform Spearman's correlation and get p-value for each cell type
#correlation_results <- age_info_long %>%
#	  dplyr::group_by(cell_type) %>%
#	  dplyr::summarize(
#		    p_value = cor.test(as.numeric(as.factor(Age_Group)), Proportion, method = "spearman")$p.value,
#				  spearman_correlation = cor(as.numeric(as.factor(Age_Group)), Proportion, method = "spearman"),
#    		.groups = "drop"
#  	)

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
writeTable(correlation_results, file.path(wd.de.data, "summary_table_correlation_results_age_median_M.txt"), colnames=T, rownames=T, sep="\t")

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
pdf(file = file.path(wd.de.plots, "correlation_results_age_median_M.pdf"), width = 6, height = 6)
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
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_mn7_M.png")),
							plot = heatmap, width = 15, height = 15, dpi = 300)

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





# -----------------------------------------------------------------------------
# DESeq2
# Last Modified: 02/09/24
# -----------------------------------------------------------------------------
for (gene in matched_genes) {
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