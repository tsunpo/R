args <- commandArgs(trailingOnly = TRUE)
age <- as.character(args[1])
ctype <- as.character(args[2])

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
embedding_list <- list()

# Loop through all graphs and run Node2Vec
#for (age in age_groups) {
#	  for (ctype in cell_types) {
	  	  age_chr <- as.character(age)
	  	  ctype_chr <- as.character(ctype)
		    key <- paste(age_chr, gsub(" ", "_", ctype_chr), sep = "_")
		    message(paste("Running Node2Vec for:", key))
		
		    # Load the graph
		    graph_file <- file.path(input_dir, paste0("GN_", key, ".rds"))
		    if (!file.exists(graph_file)) {
			      warning(paste("Graph not found:", graph_file))
			      next
		    }
		    g <- readRDS(graph_file)
		
    		# Run Node2Vec (parameters can be adjusted)
		    embedding <- node2vecR(
		    	  igraph::as_data_frame(g, what = "edges")[, c("from", "to", "weight")],
		    	  p = 1, q = 1,
		    	  num_walks = 10,
		    	  walk_length = 80,
		    	  dim = 128
		    )
		
		    # Save embedding to file and to list
		    embed_path <- file.path(output_dir, paste0("embedding_", key, ".rds"))
		    saveRDS(embedding, embed_path)
		    embedding_list[[key]] <- embedding
#	  }
#}
