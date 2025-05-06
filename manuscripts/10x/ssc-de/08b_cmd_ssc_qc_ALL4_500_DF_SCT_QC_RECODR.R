args <- commandArgs(trailingOnly = TRUE)
age <- as.character(args[1])
cell_type <- as.character(args[2])

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
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase.RData")))
#load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_scaled_data.RData")))
load(file=file.path(wd.de.data, paste0("ssc_filtered_normalised_integrated_DF_SCT_PCA_UMAP_res=0.5_", nfeatures, "_annotated_19-6-23_monocle3+phase_clean.RData")))

# -----------------------------------------------------------------------------
# 3. Generate Node2Vec Embeddings for Each Graph
#
# To obtain gene-level embeddings, we apply Node2Vec to each network. Node2Vec performs biased random walks on the graph 
# to generate sequences of nodes, then uses a Word2Vec skip-gram model to learn embeddings such that genes frequently 
# co-occurring on walks have similar vectors​
# 
# This leverages the distributional hypothesis (from NLP) in a biological context: genes with similar “contexts”
# (co-expression partners) will have similar embeddings, analogous to words in similar contexts having similar meanings​
# 
# Last Modified: 07/04/25
# -----------------------------------------------------------------------------
library(node2vec)
library(igraph)

# Set path to co-expression networks and output for embeddings
input_dir <- file.path(wd.de.data, "coexpression_networks")
output_dir <- file.path(wd.de.data, "gene_embeddings")
dir.create(output_dir, showWarnings = FALSE)

# Store embeddings in a list
#embedding_list <- list()

# Loop through all graphs and run Node2Vec
#for (age in age_groups) {
#	  for (cell_type in cell_types) {
	  	  age_chr <- as.character(age)
	  	  cell_type_chr <- as.character(cell_type)
		    key <- paste(age_chr, gsub(" ", "_", cell_type_chr), sep = "_")
		    message(paste("Running Node2Vec for:", key))
		
		    # Load the graph
		    graph_file <- file.path(input_dir, paste0("GN_", key, ".rds"))
		    if (!file.exists(graph_file)) {
			      warning(paste("Graph not found:", graph_file))
			      next
		    }
		    g <- readRDS(graph_file)
		
		    # Extract all nodes
      nodes <- V(g)$name
      # Ensure graph has all nodes in edge list (some libraries drop isolated nodes)
      edges_df <- igraph::as_data_frame(g, what = "edges") %>%
         dplyr::select(from, to, weight)
      # Add isolated nodes explicitly (if needed)
      # Or better: re-build the graph from scratch with full node list
      g2 <- graph_from_data_frame(edges_df, directed = FALSE, vertices = data.frame(name = nodes))
		    
      # Identify isolated nodes
      all_nodes <- V(g2)$name
      connected_nodes <- unique(c(edges_df$from, edges_df$to))
      isolated_nodes <- base::setdiff(all_nodes, connected_nodes)

      # Add self-loops
      if (length(isolated_nodes) > 0) {
         isolated_edges <- data.frame(
            from = isolated_nodes,
            to = isolated_nodes,
            weight = 1e-6  # very small weight
         )
         edges_augmented <- dplyr::bind_rows(edges_df, isolated_edges)
      } else {
         edges_augmented <- edges_df
      }
      # Rebuild graph with full vertex list
      g3 <- igraph::graph_from_data_frame(edges_augmented, directed = FALSE, vertices = data.frame(name = all_nodes))

      # Run Node2Vec again
      embedding <- node2vecR(
         edges_augmented,
         p = 1, q = 1,
         num_walks = 10,
         walk_length = 80,
         dim = 128
      )
      save(embedding, g3, file=file.path(output_dir, paste0("embedding_", key, ".RData")))
      
      # ✅ Ensure correct gene names (not lost during coercion)
      if (nrow(embedding) == length(V(g3))) {
         rownames(embedding) <- V(g3)$name
      } else {
         stop("Mismatch: Number of rows in embedding != number of nodes in graph")
      }
      # Check that result has rownames (node names)
      #if (!is.null(rownames(embedding)) && all(rownames(embedding) %in% V(g3)$name)) {
         #embedding_list[[key]] <- embedding
         saveRDS(embedding, file.path(output_dir, paste0("embedding_", key, ".rds")))
      #} else {
      #   warning("Embedding rownames mismatch for: ", key)
      #}
#	  }
#}
