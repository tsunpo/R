# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 20/02/25; 25/07/24; 14/03/24
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
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))

# -----------------------------------------------------------------------------
# Di Persio et al; DotPlot (resolution = 0.25)
# -----------------------------------------------------------------------------
# Apply this custom order
# Note: Only change levels if 'Idents' are factor. If they are characters, convert them to factors first.
new_order <- c("2", "5", "22", "1", "27", "21", "16", "19", "25", "18", "10", "9", "3", "6", "7", "14", "24", "30", "20", "17", "29", "13", "4", "8", "23", "31", "32", "12", "11", "26", "28", "15", "0")
#new_order <- rev(new_order)
Idents(so.integrated) <- factor(Idents(so.integrated), levels = new_order)

# Create a mapping from cluster ID to cell type name
# This mapping should be adjusted according to your specific dataset and clustering results
table(Idents(so.integrated))
write.table(table(Idents(so.integrated)), file=file.path(wd.de.data, paste0("Di Persio_SPG_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=", res, "_ordered_annotated.txt")), row.names=T, col.names=T, quote=F, sep='\t')

cluster_to_celltype <- c('2' = 'Stage 0', '5' = 'Stage 0A', '22' = 'Stage 0B', '1' = 'Stage 1',
																									'27' = 'Stage 2', '21' = 'Stage 2', '16' = 'Stage 3', '19' = '',
																									'25' = 'Leptotene', '18' = 'Leptotene',
																									'10' = 'Zygotene', '9' = 'Zygotene', '3' = 'Zygotene', '6' = '',
																									'7' = 'Pachytene', '14' = 'Pachytene', 
																									'24' = 'Diplotene', '30' = 'Diplotene',
																									'20' = 'Meiotic division', '17' = 'Meiotic division', '29' = 'Meiotic division',
																									'13' = 'Early spermatid',
																									'4' = 'Late spermatid', '8' = 'Late spermatid', '23' = '',
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

genes_of_interest <- c("DDX4", "MAGEA4", "DMRT1", "SOX4", "ID4", "FGFR3", "TCF3", "GFRA1", "NANOS2", "KIT", "MKI67", "NANOS3", "STRA8", "SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOB", "SCML1", "TEX19", "DPH7", "DMC1", "LY6K", "SELENOT", "TDRG1", "PIWIL1", "POU5F2", "OVOL2", "CCDC112", "AURKA", "CCNA1", "C9orf116", "SLC26A3", "SIRPG", "TEX29", "TNP1", "PRM2", "PRM1", "VIM", "CITED1", "SOX9", "FATE1", "HSD17B3", "STAR", "INSL3", "CLEC3B", "CFD", "MYH11", "ACTA2", "PECAM1", "VWF", "CD68", "LYZ", "C1QA", "CD14")

DefaultAssay(so.integrated) <- "SCT"
dot_plot <- DotPlot(so.integrated, features = genes_of_interest)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=0.5_14x6_19-6-23.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

save(prin_comp, so.integrated, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))

# -----------------------------------------------------------------------------
# Methods: Age group
# Last Modified: 15/01/25; 14/11/24
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))

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
   #arrange(Age)

# Calculate total cells
#total_cells <- sum(summary_table$Cells)

writeTable(summary_table, file.path(wd.de.plots, "summary_table_Germ.txt"), colnames=T, rownames=T, sep="\t")

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
dot_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )# This shows the DimPlot with the current colors

# Extract the colors for each identity (cell type) programmatically
plot_colors <- scales::hue_pal()(length(unique(cell_type_colors)))
names(plot_colors) <- levels(cell_type_colors)

# Display the color assignments for reference
#plot_colors

pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=0.5_6x6_seurat_clusters.pdf")), width = 6, height = 6)
print(dot_plot)
dev.off()

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

# Proportions sum to one, which transforms the data from the Euclidean space into a simplex space.
# Canonical statistical tests such as t-test may not be reliable, as these are not designed for relative proportion data
# But Spearman’s correlation measures the monotonic relationship between two variables by ranking the values rather than using their absolute differences. 
# This makes it more robust to non-normality and scale transformations
proportions <- cell_type_counts %>%
	  left_join(age_group_totals, by = "age_group") %>%
  	mutate(proportion = cell_count / total_cells)

# Define the median age for each age group
age_group_to_median <- c("20-30" = 25, "30-40" = 35, "40-50" = 45, 
																									"50-60" = 55, "60-70" = 65, "70-80" = 75)
# Add the `age_median` column
proportions$age_median <- as.numeric(sapply(proportions$age_group, function(x) age_group_to_median[x]))

# Plot with age groups on x-axis and cell type proportions
plot <- ggplot(proportions, aes(x = age_group, y = proportion, fill = cell_type)) +
	  geom_bar(stat = "identity", position = "fill", color = "black") +
	  scale_fill_manual(values = plot_colors, name = "Cell type") +
	  labs(title = "", x = "Age group", y = "Proportion of cell type") +
	  theme_minimal(base_size = 14) +
	  theme(
		    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, margin = margin(t = -15)),
		    axis.text.y = element_text(size = 11),
		    axis.title.x = element_text(size = 13),
		    axis.title.y = element_text(size = 13),
		    plot.title = element_text(hjust = 0.5),
		    legend.text = element_text(size = 10),   # Increase legend text size
		    legend.title = element_text(size = 12)  # Increase legend title size
	  ) +
  	guides(fill = guide_legend(ncol = 1))       # Force single-row legend

# Print the plot
png(file = file.path(wd.de.plots, "barplot_age_group_cell_type_DimPlot_5.3x7_seurat_clusters.png"), width = 5.3, height = 7, units = "in", res = 300)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Methods: Plot cell counts
# Last Modified: 15/01/25
# -----------------------------------------------------------------------------
# Plot with age groups on x-axis and cell type proportions
plot <- ggplot(proportions, aes(x = age_group, y = cell_count, fill = cell_type)) +
	  geom_bar(stat = "identity", position = "stack", color = "black") +  # Use "stack" for actual counts
	  scale_fill_manual(values = plot_colors, name = "Cell type") +
	  labs(title = "", x = "Age group", y = "Number of cells") +
	  theme_minimal(base_size = 14) +
	  theme(
		    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, margin = margin(t = -15)),
		    axis.text.y = element_text(size = 11),
		    axis.title.x = element_text(size = 13),
		    axis.title.y = element_text(size = 13),
		    plot.title = element_text(hjust = 0.5),
		    legend.text = element_text(size = 10),   # Increase legend text size
		    legend.title = element_text(size = 12)  # Increase legend title size
	  )  +
	  guides(fill = guide_legend(ncol = 1))       # Force single-row legend


# Print the plot
png(file = file.path(wd.de.plots, "barplot_age_group_cell_type_DimPlot_cell_counts_5.3x7_seurat_clusters.png"), width = 5.3, height = 7, units = "in", res = 300)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Methods: Correlation
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

# Step 3: Display the results, sorted by correlation
#correlation_results <- correlation_results %>%
#	  dplyr::arrange(desc(spearman_correlation))

# Print the results
#print(correlation_results)
#writeTable(correlation_results, file.path(wd.de.plots, "summary_table_correlation_results_age_median_19-6-23.txt"), colnames=T, rownames=T, sep="\t")

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

# Reorder levels so "Stage 0" appears at the top of the y-axis
correlation_results$cell_type <- factor(
	  correlation_results$cell_type,
	  levels = rev(correlation_results$cell_type)
)

# Plot the bar plot
#plot <- ggplot(correlation_results, aes(x = reorder(cell_type, spearman_correlation), y = spearman_correlation, fill = spearman_correlation)) +
# Plot the bar plot (without sorting)
plot <- ggplot(correlation_results, aes(x = cell_type, y = spearman_correlation, fill = spearman_correlation)) +
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
pdf(file = file.path(wd.de.plots, "barplot_correlation_results_age_median_unsorted_seurat_clusters.pdf"), width = 6, height = 6)
print(plot)
dev.off()

# -----------------------------------------------------------------------------
# Methods: Centered Log-Ratio (CLR) Transformation
# Last Modified: 20/02/25
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(compositions)

# Reshape data: wide format with cell types as columns
proportions_wide <- proportions %>%
  	select(age_group, cell_type, proportion) %>%
	  pivot_wider(names_from = cell_type, values_from = proportion, values_fill = list(proportion = 0))

# Extract age groups separately to keep track
age_groups <- proportions_wide$age_group

# Apply CLR transformation to the proportion values (excluding age_group)
clr_transformed <- clr(proportions_wide %>% select(-age_group))

# Convert back to a tibble and restore age_group
clr_transformed_df <- as_tibble(clr_transformed)
clr_transformed_df$age_group <- age_groups

# Reshape back to long format
proportions_clr <- clr_transformed_df %>%
	  pivot_longer(cols = -age_group, names_to = "cell_type", values_to = "clr_proportion")

# Merge CLR-transformed data back into original proportions table
proportions <- proportions %>%
	  left_join(proportions_clr, by = c("age_group", "cell_type"))

# Print to check results
print(proportions)

# Step 2: Perform Spearman's correlation and get p-value for each cell type
correlation_results <- proportions %>%
  	dplyr::group_by(cell_type) %>%
	  dplyr::summarize(
	    	p_value = cor.test(as.numeric(as.factor(age_group)), clr_proportion, method = "spearman")$p.value,
		    spearman_correlation = cor(as.numeric(as.factor(age_group)), clr_proportion, method = "spearman"),
		    .groups = "drop"
	  )
correlation_results <- correlation_results %>%
	  dplyr::arrange(match(cell_type, unique(proportions$cell_type)))

# Assuming correlation_results from the previous step contains the correlation and p-value for each cell type
# Add a column to indicate significance level based on p-value (optional)
correlation_results <- correlation_results %>%
	  dplyr::mutate(significance = dplyr::case_when(
	  	  p_value < 0.001 ~ "***",
		    p_value < 0.01 ~ "**",
		    p_value < 0.05 ~ "*",
		    TRUE ~ ""
	  ))

# Reorder levels so "Stage 0" appears at the top of the y-axis
correlation_results$cell_type <- factor(
	  correlation_results$cell_type,
	  levels = rev(correlation_results$cell_type)
)

# Plot the bar plot
#plot <- ggplot(correlation_results, aes(x = reorder(cell_type, spearman_correlation), y = spearman_correlation, fill = spearman_correlation)) +
# Plot the bar plot (without sorting)
plot <- ggplot(correlation_results, aes(x = cell_type, y = spearman_correlation, fill = spearman_correlation)) +
	  geom_bar(stat = "identity") +
	  coord_flip() +
  	scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, name = "Correlation") +
	  geom_text(aes(label = significance, hjust = ifelse(spearman_correlation < 0, 1.5, -0.5)),  # Adjust hjust for more space
			  								vjust = 0.5, color = "black", size = 6) + 
	  labs(
		    title = "",
		    x = "",
		    y = "Spearmans' rho (CLR-transformed)"
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
pdf(file = file.path(wd.de.plots, "barplot_correlation_results_age_median_unsorted_seurat_clusters_CLR.pdf"), width = 6, height = 6)
print(plot)
dev.off()






# -----------------------------------------------------------------------------
# Heat map for positively selected genes
# -----------------------------------------------------------------------------
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_SCT_PCA_UMAP_resolution=0.5_", nfeatures, "_-12-15-17_annotated_joined.RData")))
mn7 <- c("KDM5B", "PTPN11", "NF1", "SMAD6", "CUL3", "MIB1", "RASA2", "PRRC2A", "PTEN", "RIT1", "ROBO1", "DDX3X", "CSNK2B", "KRAS", "FGFR3", "PPM1D", "ARID1A", "BRAF", "HRAS", "KMT2E", "EP300", "SCAF4", "BMPR2", "TCF12", "CCAR2", "DHX9", "NSD1", "LZTR1", "FGFR2", "ARHGAP35", "CBL", "SSX1", "RBM12", "FAM222B", "SMAD4", "AR", "KDM5C", "KMT2D", "CTNNB1", "RAF1")

# Set the default assay to "RNA"
DefaultAssay(so.integrated) <- "RNA"
# Join layers if needed
#so.integrated <- JoinLayers(so.integrated, assays = "RNA")
# Step 1: Re-scale the data including all genes in the dataset
so.integrated <- ScaleData(so.integrated, features = rownames(so.integrated))

# Step 2: Extract the scaled data for the genes in mn7
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]
#save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))

# Step 3: Calculate the average expression for each gene in each cluster
average_expression_per_cluster <- sapply(levels(so.integrated), function(cluster) {
	  # Ensure no NAs by filling missing values
	  rowMeans(scaled_data[, Idents(so.integrated) == cluster, drop = FALSE], na.rm = TRUE)
})

# Step 4: Sort the genes sequentially by expression across clusters
ordered_genes <- rownames(scaled_data)[do.call(order, as.data.frame(-average_expression_per_cluster))]
scaled_data <- scaled_data[ordered_genes,]

# Create the heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", disp.min = -2.5, disp.max = 2.5) + NoLegend() + theme(axis.text.y = element_text(size = 12))
ggsave(filename = file.path(wd.de.plots, paste0("heatmap_mn7_20x12_monocle3+phase_2_scaled_data.png")),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
#so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_scaled_rev.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Last Modified: 23/02/25
# -----------------------------------------------------------------------------
# Step 2: Extract the scaled data for the genes in mn7
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]

# Step 4: Extract Pseudotime Values from Monocle3
pseudotime_values <- cds@principal_graph_aux$UMAP$pseudotime  # Extract Monocle3 pseudotime
pseudotime_values <- pseudotime_values[match(colnames(scaled_data), names(pseudotime_values))]  # Match order

# Step 5: Compute Correlation Between Gene Expression and Pseudotime
correlation_values <- apply(scaled_data, 1, function(gene_expression) {
	  cor(gene_expression, pseudotime_values, method = "spearman", use = "pairwise.complete.obs")
})

# Step 6: Sort Genes by Pseudotime Correlation (High Positive Correlation First)
ordered_genes <- names(sort(correlation_values, decreasing = FALSE))
# Step 7: Apply New Sorting to Scaled Data
scaled_data <- scaled_data[ordered_genes, ]

# Step 8: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident",	disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))
# Step 9: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_pseudotime_FALSE.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
#so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
  	scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_pseudotime_FALSE.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Smooth Expression Trends (on pseudotime)
# Last Modified: 24/02/25
# -----------------------------------------------------------------------------
library(mgcv)  # For GAM-based smoothing
library(reshape2)

# Step 2: Extract scaled data for genes of interest
matched_genes <- mn7[mn7 %in% rownames(so.integrated)]
gene_expression <- GetAssayData(so.integrated, assay = "RNA", layer = "data")[matched_genes, ]  # Use log-normalized expression
scaled_data <- GetAssayData(so.integrated, assay = "RNA", layer = "scale.data")[matched_genes, ]

# Step 3: Extract Pseudotime Values from Monocle3
pseudotime_values <- cds@principal_graph_aux$UMAP$pseudotime  # Extract pseudotime
pseudotime_values <- pseudotime_values[match(colnames(scaled_data), names(pseudotime_values))]  # Match order

# Step 4: Fit a Smooth Trend Line for Each Gene Using Generalized Additive Model (GAM)
smoothed_trends <- sapply(rownames(gene_expression), function(gene) {
	  fit <- gam(gene_expression[gene, ] ~ s(pseudotime_values, bs="cs"), method="REML")  # Use cubic smoothing spline
	  predict(fit)  # Get predicted smooth values
})

# Convert to matrix format
smoothed_trends <- as.matrix(smoothed_trends)
# Step 5: Compute Peak Pseudotime for Each Gene (First Max Expression in Smoothed Data)
peak_pseudotime <- apply(smoothed_trends, 2, function(x) which.max(x))
# Step 6: Order Genes by Peak Activation Time in Pseudotime
ordered_genes <- names(sort(peak_pseudotime, decreasing = FALSE))  # Earlier peak = higher rank
# Step 7: Apply New Sorting to Scaled Data
scaled_data <- scaled_data[ordered_genes, ]

# Step 8: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident",	
																					disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))

# Step 9: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_smoothed_pseudotime_.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_smoothed_pseudotime.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Smooth Expression Trends (on clusters)
# Last Modified: 24/02/25
# -----------------------------------------------------------------------------
# Step 2: Extract cluster identities
cluster_ids <- Idents(so.integrated)  # Get cluster identities per cell
# Ensure cluster_ids follow the correct order
desired_order <- c("Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3", "Leptotene", "Zygotene", "Pachytene", "Diplotene", "Meiotic division", "Early spermatid", "Late spermatid")
cluster_ids <- factor(cluster_ids, levels = desired_order)

# Step 3: Fit GAM smoothing per gene per cluster
smoothed_trends <- sapply(rownames(gene_expression), function(gene) {
	  cluster_means <- tapply(gene_expression[gene, ], cluster_ids, mean, na.rm = TRUE)  # Mean per cluster
	  #cluster_order <- as.numeric(as.factor(names(cluster_means)))  # Convert cluster names to numeric order
	  cluster_order <- match(names(cluster_means), desired_order)  # Correct biological order
	  
	  # Fit GAM model for smoothing within clusters
	  fit <- gam(cluster_means ~ s(cluster_order, bs="cs"), method="REML")
	  predict(fit)  # Get smoothed expression values per cluster
})

# Convert to matrix format
smoothed_trends <- as.matrix(smoothed_trends)
# Step 4: Compute Peak Cluster (Where Gene Expression Peaks)
peak_cluster <- apply(smoothed_trends, 2, function(x) which.max(x))  # Find the cluster with highest expression
# Step 5: Order Genes by Peak Expression in Clusters
ordered_genes <- names(sort(peak_cluster, decreasing = FALSE))  # Sort genes by earliest peak expression
# Step 6: Apply New Sorting to Scaled Data for Heatmap Visualization
scaled_data <- scaled_data[ordered_genes, ]

# Step 7: Generate the Heatmap
heatmap <- DoHeatmap(so.integrated, features = rownames(scaled_data), group.by = "ident", 
																					disp.min = -2.5, disp.max = 2.5) + NoLegend() + 
	  theme(axis.text.y = element_text(size = 12))

# Step 8: Save the Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_smoothed_clusters_desired_order.png"),
							plot = heatmap, width = 20, height = 12, dpi = 300)

# Reorder clusters so Cluster 1 appears at the top
so.integrated@active.ident <- factor(Idents(so.integrated), levels = rev(levels(Idents(so.integrated))))

dot_plot <- DotPlot(so.integrated, features = ordered_genes)  +
	  scale_color_gradientn(colors = c("blue", "white", "red")) + # Change color gradient
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

pdf(file = file.path(wd.de.plots, paste0("DotPlot_SCT_mn7_ordered_genes_smoothed_clusters_desired_order.pdf")), width = 14, height = 6)
print(dot_plot)
dev.off()

# -----------------------------------------------------------------------------
# Smooth Expression Trends (on clusters; using smoothed trends)
# Last Modified: 26/02/25
# -----------------------------------------------------------------------------
# Step 7: Generate the Heatmap Using `smoothed_trends`
heatmap <- DoHeatmap(
	  object = NULL,  # Remove Seurat object to use a custom matrix
	  features = rownames(smoothed_trends),  # Use smoothed expression trends
	  assay = NULL,  # Not required since we're using custom data
	  slot = NULL  # Not needed as `smoothed_trends` is precomputed
) + 
	  scale_fill_gradientn(colors = c("blue", "white", "red")) +  # Adjust heatmap colors
	  theme(axis.text.y = element_text(size = 12)) +
	  NoLegend()

# Save Heatmap
ggsave(filename = file.path(wd.de.plots, "heatmap_sorted_by_smoothed_clusters_desired_order_smoothed.png"), plot = heatmap, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 20/01/25; 30/08/24; 25/07/24
# -----------------------------------------------------------------------------
# > dim(so.integrated)
# [1] 35742 84822

# Identify the clusters you want to keep
clusters_to_remove <- c("Sertoli", "Fibrotic PMC", "PMC", "Leydig", "Endothelial", "Macrophage", "19", "6", "23")
clusters_to_keep <- setdiff(unique(Idents(so.integrated)), clusters_to_remove)
#clusters_to_keep <- c("Stage 0", "Stage 0A", "Stage 0B", "Stage 1", "Stage 2", "Stage 3", "Unknown")
# Subset the Seurat object to exclude cluster 1
so.integrated <- subset(so.integrated, idents = clusters_to_keep)
# > dim(so.integrated)
# [1]  5000 58098

##
dim_plot <- DimPlot(so.integrated, label = TRUE) + NoLegend() +
	  labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
	  )
ggsave(file.path(wd.de.plots, paste0("DimPlot_no-title_6*6_new-color_16_Spermatogenesis_19-6-23.png")), plot = dim_plot, width = 6, height = 6, dpi = 300)

##
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
																							graph_label_size=0) +
	  guides(color = guide_colorbar(title = "Pseudotime")) +  # Force rename legend title
	  theme_classic() +
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("pseudotime_23_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_integrated_19-23-6_6.5x6.png")), plot = monocle3, dpi = 300, width = 6.5, height = 6)

# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#tscan
library(ggbeeswarm)
library(colorRamps)
library(viridisLite)

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$idents <- Idents(so.integrated)
pdata_cds$Phase <- colData(cds)$Phase  # Ensure Phase is included

#pseudotime_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
#pseudotime_palette <- inferno(100)
pseudotime_palette <- colorRampPalette(c("darkblue", "red", "yellow"))(100)

ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = pseudotime_monocle3)) +  # Use pseudotime as the color
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_gradientn(colors = pseudotime_palette) +  # Apply the pseudotime gradient
	  scale_y_discrete(limits = rev(levels(pdata_cds$idents))) +  # Reverse y-axis order
	  theme_classic() +
  	xlab("Pseudotime") +
   ylab("") +
	  theme(plot.title = element_blank(),  # Remove title
			  				axis.text = element_text(size = 16),  # Increase axis text size
				  			axis.title = element_text(size = 18), # Increase axis title size
						  	legend.position = "none")  # Remove legend
ggsave(file.path(wd.de.plots, paste0("pseudotime_23_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered_integrated_19-23-6_7x6_rev.png")), plot = ordered, dpi = 300, width = 7, height = 6)

##
ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = idents, colour = Phase)) +  # Use pseudotime as the color
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply the same Phase colours as DimPlot
	  scale_y_discrete(limits = rev(levels(pdata_cds$idents))) +  # Reverse y-axis order
	  theme_classic() +
	  xlab("Pseudotime") +
	  ylab("") +
	  guides(color = guide_legend(override.aes = list(size = 5))) +  # Increase legend dot size
	  theme(plot.title = element_blank(),  # Remove title
			  				axis.text = element_text(size = 16),  # Increase axis text size
						  	axis.title = element_text(size = 18), # Increase axis title size
	  						legend.position = c(1, 1),  # Move legend to top-right
	  						legend.justification = c(1, 1),  # Align legend correctly
	  						legend.text = element_text(size = 16),  # Increase legend text size
	  						legend.title = element_text(size = 16))
ggsave(file.path(wd.de.plots, paste0("pseudotime_23_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_ordered_integrated_19-23-6_7x6_rev_phase_legend.png")), plot = ordered, dpi = 300, width = 7, height = 6)

# -----------------------------------------------------------------------------
# Cell cycle analysis
# Last Modified: 05/02/25
# -----------------------------------------------------------------------------
nfeatures <- 5000
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))

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
# Now plot the UMAP with the correct legend order and custom colors
dim_plot <- DimPlot(so.integrated, group.by = "Phase", reduction = "umap") + 
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply custom colors & legend order
  	labs(x = "UMAP 1", y = "UMAP 2") +  # Set x-axis and y-axis labels
	  theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("Cell-cycle_annotated_19-6-23_DimPlot_no-title_8*8_new-color_16.png")), plot = dim_plot, width = 8, height = 8, dpi = 300)

VlnPlot(so.integrated, group.by = "Phase", cols= phase_colors, features = c("S.Score", "G2M.Score"), pt.size = 0)

# -----------------------------------------------------------------------------
# After running Monocle3
# -----------------------------------------------------------------------------
# Optionally, plot cells colored by Seurat clusters
monocle3 <- plot_cells(cds, color_cells_by = "Phase",
																							label_cell_groups=FALSE,
																							label_leaves=FALSE,
																							label_branch_points=FALSE,
																							graph_label_size=0) +
	  theme_classic() +
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) + # Apply new color
		 theme(
		    plot.title = element_blank(),  # Remove title
		    axis.title = element_text(size = 18),  # Increase axis title size
		    axis.text = element_text(size = 16),  # Increase axis tick label size
		    legend.text = element_text(size = 16),  # Increase legend text size
		    legend.title = element_text(size = 16)  # Increase legend title size
	  )
ggsave(file.path(wd.de.plots, paste0("Cell-cycle_Phase_", nfeatures, "_UMAP_dims=", prin_comp, "_umap_embeddings_use_pseudotime_integrated_new-color_6.5x6_16.png")), plot = monocle3, width = 6.5, height = 6, dpi = 300)

library(ggbeeswarm)
library(colorRamps)
library(viridisLite)
library(ggplot2)

# Extract metadata
pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdata_cds$Phase <- factor(so.integrated@meta.data$Phase, levels = c("G2M", "S", "G1"))  # Use cell cycle phase

# Define the Seurat default cell cycle colours (modify if needed)
phase_colors <- c("G1" = "#619CFF", "S" = "#F8766D", "G2M" = "#00BA38")

# Plot pseudotime ordered by Phase
ordered <- ggplot(as.data.frame(pdata_cds), 
																		aes(x = pseudotime_monocle3, 
																						y = Phase, 
																						colour = Phase)) +  # Use Phase for colour
	  geom_quasirandom(groupOnX = FALSE) +
	  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M")) +  # Apply the same Phase colours as DimPlot
	  scale_y_discrete(labels = c("G1" = "G1", "S" = "S", "G2M" = "G2/M")) +  # Rename y-axis labels
	  theme_classic() +
	  xlab("Pseudotime") + 
	  ylab("Phase") +
	  theme(plot.title = element_blank(),  # Remove title
	  	   axis.text = element_text(size = 16),  # Increase axis text size
							axis.title = element_text(size = 18), # Increase axis title size
							legend.position = "none")  # Remove legend

# Save the plot
ggsave(file.path(wd.de.plots, paste0("cell_cycle_phase_pseudotime_by_phase_new-color_16_6x6.png")),	plot = ordered, dpi = 300, width = 6, height = 6)

save(cds, so.integrated, phase_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
