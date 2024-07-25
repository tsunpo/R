# =============================================================================
# Library      : Single Cell Transcriptomics
# Name         : handbook-of/SingleCellTranscriptomics.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 25/07/24
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Monocle 3
# Last Modified: 25/07/24
# -----------------------------------------------------------------------------
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
