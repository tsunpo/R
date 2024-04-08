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

colnames <- c("PD_ID", "genes_0", "genes_1", "genes_shared", "cells_0", "cells_1", "cells_shared")
raw <- toTable(0, length(colnames), nrow(samples0), colnames)
raw$PD_ID <- rownames(samples0)
rownames(raw) <- rownames(samples0)

for (s in 1:nrow(samples0)) {
	  # Initialize the Seurat object with the raw (non-normalized data).
	  ssc0.data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "GRCh38-2020", samples0$V1[s], "filtered_feature_bc_matrix"))
	  #ssc0 <- CreateSeuratObject(counts=ssc0.data, project="ssc", min.cells=3, min.features=200)
	  #ssc0 <- NormalizeData(ssc0)
	  
	  ssc1.data <- Read10X(data.dir=file.path("/lustre/scratch126/casm/team294rr/mp29/scRNA_10x", "homemade_ref", samples1$V1[s], "filtered_feature_bc_matrix"))
	  #ssc1 <- CreateSeuratObject(counts=ssc1.data, project="ssc", min.cells=3, min.features=200)
	  #ssc1 <- NormalizeData(ssc1)
	  
	  raw[s, 2] <- nrow(ssc0.data)
	  raw[s, 3] <- nrow(ssc1.data)
	  raw[s, 4] <- length(intersect(rownames(ssc0.data), rownames(ssc1.data)))
	  
	  raw[s, 5] <- ncol(ssc0.data)
	  raw[s, 6] <- ncol(ssc1.data)
	  raw[s, 7] <- length(intersect(colnames(ssc0.data), colnames(ssc1.data)))
	  
	  #genes <- intersect(rownames(ssc0.data), rownames(ssc1.data))
	  #ssc0.data <- ssc0.data[genes,]
	  #ssc1.data <- ssc1.data[genes,]
	  #
	  #sum0 <- mapply(x = 1:nrow(ssc0.data), function(x) sum(ssc0.data[x,]))
	  #sum1 <- mapply(x = 1:nrow(ssc1.data), function(x) sum(ssc1.data[x,]))
	  #
	  #file.name <- file.path(wd.de.plots, paste0(raw$PD_ID[s], "_2020-vs-homemade"))
	  #main.text <- c(raw$PD_ID[s], "")
	  #xlab.text <- paste0("GRCh38-2020-A", " (Counts)")
	  #ylab.text <-  paste0("homemade_ref", " (Counts)")
	  #plotCorrelation(file.name, main.text, xlab.text, ylab.text, x=sum0, y=sum1, pos="topright", cols=c("black", "black"), size=5, p=8)
}
filtered <- raw
writeTable(filtered, file.path(wd.de.plots, "ssc_filtered.txt"), colnames=T, rownames=F, sep="\t")







tsv <- file.path(wd.rna.raw, samples)
s2c <- data.frame(path=tsv, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default filter settings
# reading in kallisto results
# dropping unused factor levels
# ......................................................
# normalizing est_counts
# 87172 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 20328 genes passed the filter
# summarizing bootstraps

tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5p47.RData")))

## Gene-level TPMs (All)
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
#writeTable(tpm.gene, gzfile(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.txt.gz"))), colnames=T, rownames=T, sep="\t")
nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs (Detected)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
nrow(tpm.gene)
# [1] 22807

load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene, 1)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median1.RData")))
nrow(tpm.gene)
# [1] 

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
nrow(tpm.gene)
# [1] 18502

## Gene-level TPMs (Detected + Expressed)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.median0 <- tpm.gene
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene.r5p47   <- tpm.gene
overlaps <- intersect(rownames(tpm.gene.median0), rownames(tpm.gene.r5p47))

tpm.gene <- tpm.gene[overlaps,]
#save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.r5p47.RData")))
nrow(tpm.gene)
# [1] 18502

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 06/09/20; 29/05/20
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Ensembl", tpm=NA)

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "expressed", tpm=0)

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median1"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "expressed", tpm=1)

## Consistently expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "expressed", tpm="r5p47")
