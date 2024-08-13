# =============================================================================
# Library      : Single Cell Transcriptomics
# Name         : handbook-of/SingleCellTranscriptomics.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 26/07/24
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Generate the Box Plot for Read Counts, Global Expression Levels, and Mitochondrial Content
# Last Modified: 02/08/24
# -----------------------------------------------------------------------------
get_cluster_factor_levels <- function(so.integrated, factor) {
	  # Extract the total expression level (sum of counts) and cluster identities
	  expression_data <- FetchData(so.integrated, vars = c(factor, "seurat_clusters"))
	
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
	
	  return(expression_data)
}

get_box_plot <- function(expression_data, title, y) {
	  box_plot_reads <- ggplot(expression_data, aes(x = cluster, y = total_expression_level)) +
	     geom_boxplot() +
	     theme_minimal() +
     	labs(title = title,
					     	x = "Cluster",
						     y = y) +
	     theme(axis.text.x = element_text(angle = 45, hjust = 1),
							     plot.title = element_text(hjust = 0.5))  # Center the title
	  
	  return(box_plot_reads)
}

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
library(monocle3)

getMonocle3CDS <- function(so.integrated, umap_embeddings=T) {
	  # Extract SCT data
	  data <- as(as.matrix(so.integrated@assays$integrated@data), "sparseMatrix")
	
	  # Extract cell metadata
	  cell_metadata <- so.integrated@meta.data
	  # Ensure cell_metadata matches the filtered cells
	  cell_metadata <- cell_metadata[colnames(data), ]
	
	  # Create the gene annotation
	  gene_annotation <- data.frame(gene_short_name = rownames(data))
	  rownames(gene_annotation) <- rownames(data)
	
	  # Create the Monocle 3 Cell Data Set (CDS)
	  cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

	  if (umap_embeddings) {
	  	  # Extract UMAP embeddings from Seurat
	  	  umap_embeddings <- Embeddings(so.integrated, reduction = "umap")
	  	
	  	  # Add UMAP embeddings to Monocle 3 object
	  	  reducedDims(cds)$UMAP <- umap_embeddings
	  }
	  
	  return(cds)
}

# -----------------------------------------------------------------------------
# Jaccard index for overlap of gene sets
# -----------------------------------------------------------------------------
getJaccardIndex <- function(markers) {
	  # Create a list of gene sets for each cluster
	  gene_sets <- markers %>% 
		 group_by(cluster) %>% 
	 	summarise(genes = list(gene))
	
	  # Compute Jaccard index for each pair of clusters
	  jaccard_index <- function(set1, set2) {
		    intersection <- length(intersect(set1, set2))
		    union <- length(union(set1, set2))
		    return(intersection / union)
  	}
	
	  # Create a matrix to store Jaccard indices
	  n_clusters <- length(unique(markers$cluster))
	  jaccard_matrix <- matrix(0, n_clusters, n_clusters, dimnames = list(unique(markers$cluster), unique(markers$cluster)))
	
	  for (i in 1:n_clusters) {
		    for (j in 1:n_clusters) {
			      if (i != j) {
				        jaccard_matrix[i, j] <- jaccard_index(gene_sets$genes[[i]], gene_sets$genes[[j]])
			      }
		    }
	  }
	  
	  return(jaccard_matrix)
}

# -----------------------------------------------------------------------------
# Functional enrichment analysis (ClusterProfiler)
# -----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db) # for human gene annotations; use org.Mm.eg.db for mouse
library(tidyr)

# Perform GO enrichment for each cluster
enrich_result <- function(genes) {
	  enrichGO(gene         = genes,
			  							OrgDb        = org.Hs.eg.db,
					  					keyType      = "SYMBOL",
							  			ont          = "BP",         # Ontology: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
								  		pAdjustMethod = "BH",
									  	pvalueCutoff = 0.01,
										  qvalueCutoff = 0.05)
}

perform_clusterprofiler_enrichment <- function(markers) {
	  # Convert Seurat markers to a list of gene symbols for each cluster
	  markers_list <- split(markers$gene, markers$cluster)
	
	  # Perform GO enrichment for each cluster
	  enrich_results <- lapply(markers_list, function(genes) {
		    enrichGO(gene         = genes,
			 		   						OrgDb        = org.Hs.eg.db,
				 				   			keyType      = "SYMBOL",
						 					   ont          = "BP",         # Ontology: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
					 						   pAdjustMethod = "BH",
							 				   pvalueCutoff = 0.01,
								 			   qvalueCutoff = 0.05)
	  })
	  
	  # Filter out NULL results (clusters with no significant enrichment)
	  enrich_results <- Filter(Negate(is.null), enrich_results)
	  
	  return(enrich_results)
}

# Function to create custom barplot using ggplot2
plot_custom_barplot <- function(enrich_result, cluster_name, col = "gray") {
	  if (is.null(enrich_result)) return(NULL)
	  enrich_data <- as.data.frame(enrich_result)
	  enrich_data <- enrich_data %>%
		    dplyr::mutate(log2P = -log2(p.adjust)) %>%
		    dplyr::arrange(desc(log2P)) %>%
		    dplyr::mutate(Description = factor(Description, levels = rev(Description))) %>%
		    dplyr::slice_head(n = 10)  # Select top 10 terms
	
	  p <- ggplot(enrich_data, aes(x = Description, y = log2P)) +
	  	  geom_bar(stat = "identity", fill = col) +  # Adjust width to make bars narrower
		    coord_flip() +
		    labs(title = paste("", cluster_name),
				       x = "GO:BP Term",
				       y = "-log2[P-value]") +
	     theme_minimal() +
   	  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
	
	  return(p)
}

plot_custom_barplot_flip <- function(enrich_result, cluster_name, col = "gray") {
	  if (is.null(enrich_result)) return(NULL)
	  enrich_data <- as.data.frame(enrich_result)
	  enrich_data <- enrich_data %>%
		    dplyr::mutate(log2P = -log2(p.adjust)) %>%
		    dplyr::arrange(desc(log2P)) %>%
		    dplyr::mutate(Description = factor(Description, levels = rev(Description))) %>%
		    dplyr::slice_head(n = 10)  # Select top 10 terms
	
	  p <- ggplot(enrich_data, aes(x = Description, y = log2P)) +
		    geom_bar(stat = "identity", fill = col) +
		    coord_flip() +
		    scale_y_reverse() +  # Flip the x-axis to have the y-axis on the right
		    labs(title = paste("", cluster_name),
						    	x = "GO:BP Term",
							    y = "-log2[P-value]") +
		 theme_minimal() +
		 theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
	
	  return(p)
}

# Prepare data for heatmap
prepare_heatmap_data <- function(enrich_results) {
	  # Create a data frame with -log2(p.adjust) for heatmap
	  enrich_results_combined <- do.call(rbind, lapply(seq_along(enrich_results), function(i) {
		    result <- enrich_results[[i]]
		    result <- as.data.frame(result)
		    result$Cluster <- names(enrich_results)[i]
		    return(result)
	  }))
	
	  # Calculate -log2(p.adjust) and handle infinite values
	  enrich_results_combined <- enrich_results_combined %>%
	    	mutate(log2P = -log2(p.adjust)) %>%
	  	  dplyr::select(Cluster, Description, log2P) %>%
	    	pivot_wider(names_from = Cluster, values_from = log2P)
	  
	  # Replace NA values with 1 (non-significant)
	  enrich_results_combined <- enrich_results_combined %>%
	  	  mutate(across(everything(), ~ replace_na(.x, 1)))
	  
	  # Convert data frame to matrix for heatmap
	  heatmap_matrix <- as.matrix(enrich_results_combined[,-1])
	  rownames(heatmap_matrix) <- enrich_results_combined$Description
	
	  # Replace NA values with 0
	  heatmap_matrix[is.na(heatmap_matrix)] <- 0
	  
	  return(heatmap_matrix)
}

# Extract p-values and prepare matrix
prepare_pvalue_matrix <- function(enrich_results_list) {
	  # Combine enrichment results into a single data frame
	  enrichment_results_combined <- do.call(rbind, lapply(seq_along(enrich_results_list), function(i) {
		    result <- enrich_results_list[[i]]
		    if (is.null(result) || nrow(result) == 0) return(NULL)
		    result <- as.data.frame(result)
		    result$Cluster <- names(enrich_results_list)[i]
		    return(result)
	  }))
	
	  # Extract p-values and pivot the data frame
	  pvalue_data <- enrichment_results_combined %>%
		    dplyr::select(Cluster, Description, p.adjust) %>%
		    pivot_wider(names_from = Cluster, values_from = p.adjust)
	
	  # Replace NA values with 1 (non-significant)
	  pvalue_data <- pvalue_data %>%
	    	mutate(across(everything(), ~ replace_na(.x, 1)))
	
	  # Convert data frame to matrix for heatmap
	  pvalue_matrix <- as.matrix(pvalue_data[,-1])
	  rownames(pvalue_matrix) <- pvalue_data$Description
	
	  return(pvalue_matrix)
}

# -----------------------------------------------------------------------------
# Functional enrichment analysis (STRINGdb)
# -----------------------------------------------------------------------------
library(STRINGdb)

perform_string_enrichment <- function(gene_list) {
  	# Ensure the gene list is a character vector
	  if (is.null(gene_list) || length(gene_list) == 0) return(NULL)
	  gene_list <- as.character(gene_list)
	
	  # Convert gene list to data frame for mapping
	  gene_df <- data.frame(gene = gene_list)
	
	  # Map gene names to STRING IDs
	  mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
	  if (nrow(mapped_genes) == 0) return(NULL)
	
	  # Perform enrichment analysis
	  enrich_results <- string_db$get_enrichment(mapped_genes$STRING_id, category = "Process")
	  return(enrich_results)
}

get_heatmap_matrix <- function(markers) {
	  # Convert Seurat markers to a list of gene symbols for each cluster
	  markers_list <- split(markers$gene, markers$cluster)
	
	  # Perform enrichment analysis for each cluster
	  enrich_results <- lapply(markers_list, perform_string_enrichment)
	  # Filter out NULL results (clusters with no significant enrichment)
	  enrich_results <- Filter(Negate(is.null), enrich_results)
	
	  # Combine results and focus on GO categories
	  enrich_results_combined <- bind_rows(enrich_results, .id = "Cluster")
	  #enrichment_results_combined <- enrichment_results_combined %>%
	  #	  filter(category == "Process") %>%
	  #	  arrange(FDR)
	
	  # Create a data frame with -log2(p.adjust) for heatmap
	  heatmap_data <- enrich_results_combined %>%
		    filter(category == "Process") %>%
		    mutate(log2P = -log2(`p.adjust`)) %>%
		    select(Cluster, term, log2P) %>%
		    spread(Cluster, log2P)
	  
	  # Replace NA values with 1 (non-significant)
	  heatmap_data <- heatmap_data %>%
	  	  mutate(across(everything(), ~ replace_na(.x, 1)))
	
	  # Convert data frame to matrix for heatmap
	  heatmap_matrix <- as.matrix(heatmap_data[,-1])
	  rownames(heatmap_matrix) <- heatmap_data$term
	
	  return(heatmap_matrix)
}

# -----------------------------------------------------------------------------
# GSEA
# -----------------------------------------------------------------------------
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(grid)

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

# Extract significant terms from GSEA results
extract_significant_terms <- function(gsea_result, pvalue_cutoff = 0.05) {
	  if (is.null(gsea_result)) return(NULL)
	     result <- gsea_result@result
	     significant_terms <- result[result$p.adjust < pvalue_cutoff, ]
	     return(significant_terms)
}

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

# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression for each cluster
# -----------------------------------------------------------------------------
library(ggrepel)

calculate_spearman <- function(cluster, age_cluster_data, so.integrated) {
	  # Subset the data for the specific cluster
	  cells_in_cluster <- WhichCells(so.integrated, idents = cluster)
	  age_data <- age_cluster_data[cells_in_cluster, "age"]
	  gene_expression <- GetAssayData(so.integrated, assay = "RNA", slot = "data")[, cells_in_cluster]
	
	  # Calculate Spearman's correlation for each gene
	  cor_results <- apply(gene_expression, 1, function(gene_expr) {
		    cor.test(gene_expr, age_data, method = "spearman")
	  })
	
	  # Extract correlation coefficients and p-values
	  cor_coefficients <- sapply(cor_results, function(res) res$estimate)
	  p_values <- sapply(cor_results, function(res) res$p.value)
	
	  # Adjust p-values for multiple testing
	  p_adj_values <- p.adjust(p_values, method = "BH")
	
	  # Create a data frame with results
	  result_df <- data.frame(
		    gene = rownames(gene_expression),
		    cor_coefficient = cor_coefficients,
		    p_value = p_values,
		    p_adj_value = p_adj_values
	  )
	
	  # Filter for significant correlations (e.g., p_adj_value < 0.05)
	  significant_genes <- result_df %>%
 		   filter(p_adj_value < 0.05) %>%
	  	  arrange(p_adj_value)  # Order by adjusted p-value
	
	  p_value_threshold <- max(rbind(subset(significant_genes, cor_coefficient >= 0.2), subset(significant_genes, cor_coefficient <= -0.2))$p_adj_value)
 
	  # Create a volcano plot
	  volcano_plot <- ggplot(significant_genes, aes(x = cor_coefficient, y = -log10(p_adj_value))) +
	  	  geom_point(aes(color = ifelse(p_adj_value < p_value_threshold & cor_coefficient < 0, blue,
	  																															ifelse(p_adj_value < p_value_threshold & cor_coefficient > 0, red, "grey")))) +
	  	  scale_color_manual(values = c(blue, red, "grey")) +
	  	  theme_minimal() +
	  	  labs(title = paste0("", cluster),
	  			  			x = "Correlation Coefficient",
	  				  		y = "-log10 Adjusted P-value") +
	    	theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
	  	  geom_vline(xintercept = 0.2, linetype = "dashed", color = red) +
	  	  geom_vline(xintercept = -0.2, linetype = "dashed", color = blue) +
	    	geom_text_repel(data = subset(significant_genes, p_adj_value < p_value_threshold),
	  	  																aes(label = gene), size = 3, max.overlaps = 30)
	  # Print the volcano plot
	  pdf(file = file.path(wd.de.plots, paste0("volcano_0.1_SCT_res=0.5_", cluster, "_volcano.pdf")), width = 5, height = 5)
	  print(volcano_plot)
	  dev.off()
	  
	  #genes_pos <- subset(subset(significant_genes, p_adj_value <= p_value_threshold), cor_coefficient > 0)$gene
	  #genes_neg <- subset(subset(significant_genes, p_adj_value <= p_value_threshold), cor_coefficient < 0)$gene
	  #enrich_result_pos <- enrich_result(genes_pos)
	  #enrich_result_neg <- enrich_result(genes_neg)
	  
	  #p <- plot_custom_barplot(enrich_result_pos, paste0("Cluster ", cluster, " (Up)"), col = red)
	  #if (!is.null(p)) {
	  #	  pdf(file = file.path(wd.de.plots, paste0("volcano_0.1_SCT_res=0.25_cluster", cluster, "_go_pos.pdf")), width = 7, height = 2.5)
	  #  print(p)
	  #	  dev.off()
	  #}

	  #p <- plot_custom_barplot_flip(enrich_result_neg, paste0("", cluster, " (Down)"), col = blue)
	  #if (!is.null(p)) {
	  #	  pdf(file = file.path(wd.de.plots, paste0("volcano_0.1_SCT_res=0.5_", cluster, "_go_neg.pdf")), width = 7, height = 2.5)
	  #	  print(p)
	  #	  dev.off()
	  #}
	  
	  return(significant_genes)
}

calculate_spearman_mn7 <- function(cluster, age_cluster_data, so.integrated) {
	  # Subset the data for the specific cluster
	  cells_in_cluster <- WhichCells(so.integrated, idents = cluster)
	  age_data <- age_cluster_data[cells_in_cluster, "age"]
	  gene_expression <- GetAssayData(so.integrated, assay = "RNA", slot = "data")[, cells_in_cluster]
	
	  # Calculate Spearman's correlation for each gene
	  cor_results <- apply(gene_expression, 1, function(gene_expr) {
	    	cor.test(gene_expr, age_data, method = "spearman")
	  })
	
	  # Extract correlation coefficients and p-values
	  cor_coefficients <- sapply(cor_results, function(res) res$estimate)
	  p_values <- sapply(cor_results, function(res) res$p.value)
	
	  # Adjust p-values for multiple testing
	  p_adj_values <- p.adjust(p_values, method = "BH")
	
	  # Create a data frame with results
	  result_df <- data.frame(
		    gene = rownames(gene_expression),
	 	   cor_coefficient = cor_coefficients,
	 	   p_value = p_values,
	 	   p_adj_value = p_adj_values
	  )
	
	  # Filter for significant correlations (e.g., p_adj_value < 0.05)
	  significant_genes <- result_df %>%
	    	filter(p_adj_value < 1) %>%
		    arrange(p_adj_value)  # Order by adjusted p-value
	
	  p_value_threshold <- max(rbind(subset(significant_genes, cor_coefficient >= 0.1), subset(significant_genes, cor_coefficient <= -0.1))$p_adj_value)
	
	  significant_genes_mn7 <- subset(significant_genes, gene %in% mn7)
	  significant_genes_mn7$color <- ifelse(significant_genes_mn7$p_adj_value < p_value_threshold & significant_genes_mn7$cor_coefficient < -0.1, blue,
	  									ifelse(significant_genes_mn7$p_adj_value < p_value_threshold & significant_genes_mn7$cor_coefficient > 0.1, red, "grey"))
	  
  	# Create a volcano plot
	  volcano_plot <- ggplot(significant_genes_mn7, aes(x = cor_coefficient, y = -log10(p_adj_value))) +
		    geom_point(aes(color = color)) +
	  	  scale_color_manual(values = unique(significant_genes_mn7$color)) +
		    theme_minimal() +
		    labs(title = paste0(cluster),
						    	x = "Correlation Coefficient",
							    y = "-log10 Adjusted P-value") +
		    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
	    	geom_vline(xintercept = 0.1, linetype = "dashed", color = red) +
		    geom_vline(xintercept = -0.1, linetype = "dashed", color = blue) +
		    #geom_text_repel(data = subset(significant_genes, p_adj_value < p_value_threshold),
		    #																aes(label = gene), size = 3, max.overlaps = 30)
		    geom_text_repel(data = subset(significant_genes_mn7, p_adj_value < 1),
						    												aes(label = gene), size = 3, max.overlaps = 30)
  	# Print the volcano plot
	  pdf(file = file.path(wd.de.plots, paste0("mn7_volcano_0.1_SCT_res=0.5_", cluster, "_volcano.pdf")), width = 5, height = 5)
	  print(volcano_plot)
	  dev.off()
	
	  return(significant_genes)
}

# -----------------------------------------------------------------------------
# Spearman's correlation between age and gene expression for each cluster
# -----------------------------------------------------------------------------