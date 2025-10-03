# =============================================================================
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 01/09/25
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

load(file.path(wd.de.data, paste0("ssc_filtered_normalised_DF_SCT_5000_25_100_integrated_PCA_UMAP_23_0.5_-C0_ordered_annotated_monocle3+phase.RData")))

# -----------------------------------------------------------------------------
# Methods: Age group
# Last Modified: 15/01/25; 14/11/24
# -----------------------------------------------------------------------------
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

writeTable(summary_table, file.path(wd.de.plots, "summary_table_germ.txt"), colnames=T, rownames=T, sep="\t")

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

#pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=0.5_6x6_seurat_clusters.pdf")), width = 6, height = 6)
#print(dot_plot)
#dev.off()pdf(file = file.path(wd.de.plots, paste0("Di Persio_DotPlot_SCT_", nfeatures, "_SCT_dims=", prin_comp, "_res=0.5_6x6_seurat_clusters.pdf")), width = 6, height = 6)
#print(dot_plot)
#dev.off()

#save(cds, so.integrated, phase_colors, plot_colors, file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_23_res=0.5_-C0_ordered_annotated_monocle3+phase.RData")))

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
	  labs(title = "", x = "Age group", y = "Proportion of cells") +
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
png(file = file.path(wd.de.plots, "6-bins_barplot_age-type_5.3x7.png"), width = 5.3, height = 7, units = "in", res = 300)
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
png(file = file.path(wd.de.plots, "6-bins_barplot_age-count_5.3x7.png"), width = 5.3, height = 7, units = "in", res = 300)
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
pdf(file = file.path(wd.de.plots, "6-bins_barplot_age-type_correlation_spearman.pdf"), width = 6, height = 6)
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
   dplyr::select(age_group, cell_type, proportion) %>%
   tidyr::pivot_wider(names_from = cell_type, values_from = proportion, values_fill = 0)

# Extract age groups separately to keep track
age_groups <- proportions_wide$age_group

# Apply CLR transformation to the proportion values (excluding age_group)
clr_transformed <- clr(proportions_wide %>% dplyr::select(-age_group))

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
pdf(file = file.path(wd.de.plots, "6-bins_barplot_age-type_correlation_CLR.pdf"), width = 6, height = 6)
print(plot)
dev.off()
