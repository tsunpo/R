#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
s <- as.numeric(args[1])

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
plotCorrelation <- function(file.name, main.text, xlab.text, ylab.text, x, y, pos="bottomleft", cols=c("dimgray", "black"), size=5, pch=1, cex=2, p=12, chr=2) {
	  pdf(paste0(file.name, ".pdf"), height=size, width=5)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
  	plot(y ~ x, ylab=ylab.text, xlab=xlab.text, main=main.text, pch=pch, cex=cex, col=cols[1], cex.axis=1.9, cex.lab=2, cex.main=2.1)
	
	  lm.fit <- lm(y ~ x)
	  abline(lm.fit, lwd=5, col=cols[2])
	
	  cor <- cor.test(y, x, method="spearman", exact=F)
	  legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("P = 1.00E-00")), text.col=c(cols[2], "white"), text.font=1, bty="n", cex=2)
	  legend(pos, c("", expression(italic('P')~"                   ")), text.col=cols[2], text.font=1, bty="n", cex=2)
	  legend(pos, c("", paste0("   = ", scientific(cor[[3]]))), text.col=cols[2], text.font=1, bty="n", cex=2)
	  dev.off()
}

library(dplyr)
library(Seurat)
library(patchwork)

#for (s in 1:nrow(samples0)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
   ssc0.data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0$V1[s], "filtered_feature_bc_matrix"))
	  ssc0 <- CreateSeuratObject(counts=ssc0.data, project="ssc", min.cells=3, min.features=200)
	  ssc0 <- NormalizeData(ssc0)
	  ssc0.counts <- GetAssayData(object=ssc0, layer="counts")
	  
	  ssc1.data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "homemade_ref", samples1$V1[s], "filtered_feature_bc_matrix"))
	  ssc1 <- CreateSeuratObject(counts=ssc1.data, project="ssc", min.cells=3, min.features=200)
	  ssc1 <- NormalizeData(ssc1)
	  ssc1.counts <- GetAssayData(object=ssc1, layer="counts")
	  
   genes <- intersect(rownames(ssc0.counts), rownames(ssc1.counts))
   ssc0.counts <- ssc0.counts[genes,]
   ssc1.counts <- ssc1.counts[genes,]
	  
   sum0 <- mapply(x = 1:nrow(ssc0.counts), function(x) sum(ssc0.counts[x,]))
   sum1 <- mapply(x = 1:nrow(ssc1.counts), function(x) sum(ssc1.counts[x,]))
	  
   file.name <- file.path(wd.de.plots, paste0(rownames(samples1)[s], "_2020-vs-homemade_Normalised"))
   main.text <- c(rownames(samples1)[s], "")
   xlab.text <- paste0("GRCh38-2020-A")
   ylab.text <-  paste0("homemade_ref")
   plotCorrelation(file.name, main.text, xlab.text, ylab.text, x=sum0, y=sum1, pos="bottomright", cols=c("black", "black"), size=5)
#}
