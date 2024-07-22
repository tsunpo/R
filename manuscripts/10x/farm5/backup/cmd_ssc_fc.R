#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : 
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 14/03/24
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"            ## ty2@farm
#wd.src <- "/Users/tpyang/Work/dev/R"             ## ty2@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Graphics.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg38.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## ty2@farm
#wd <- "/Users/ty2/Work/uni-koeln/tyang2"       ## ty2@localhost
BASE <- "SSC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/scRNA")
wd.rna.raw <- file.path(wd.rna, "10x")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples0 <- readTable(file.path(wd.rna.raw, "scRNA_GRCh38-2020.list"), header=F, rownames=3, sep="\t")
samples1 <- readTable(file.path(wd.rna.raw, "scRNA_homemade_ref.list"), header=F, rownames=3, sep="\t")
samples1 <- samples1[rownames(samples0),]

# -----------------------------------------------------------------------------
# Associating transcripts to gene-level TPM estimates using sleuth (v0.29.0)
# 
# -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)

#colnames <- c("PD_ID", "genes_0", "genes_1", "genes_shared", "cells_0", "cells_1", "cells_shared")
#raw <- toTable(0, length(colnames), nrow(samples0), colnames)
#raw$PD_ID <- rownames(samples0)
#rownames(raw) <- rownames(samples0)

ssc.counts <- c()
for (s in 1:nrow(samples0)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
	  ssc0.data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0$V1[s], "filtered_feature_bc_matrix"))
	  ssc0 <- CreateSeuratObject(counts=ssc0.data, project="ssc", min.cells=3, min.features=200)
	  ssc0 <- NormalizeData(ssc0)
	  ssc0.counts <- GetAssayData(object=ssc0, layer="counts")
	  writeTable(as.data.frame(ssc0.counts), gzfile(file.path(wd.de.data, "GRCh38-2020-A", paste0(samples0$V3[s], ".txt.gz"))), colnames=T, rownames=F, sep="\t")
	  
	  if (s == 1) {
	  	  ssc.counts <- ssc0.counts
	  } else {
	  	  genes <- intersect(rownames(ssc.counts), rownames(ssc0.counts))
	  	  ssc.counts <- cbind(ssc.counts[genes,], ssc0.counts[genes,])
	  }
	  rm(ssc0.counts)
}
