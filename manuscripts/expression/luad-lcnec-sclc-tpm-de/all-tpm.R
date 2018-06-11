# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/all-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 30/11/17
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Sleuth installation (v0.29.0)
# -----------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5", "tximport")
install.packages("devtools", "stringr", "cluster", "survival", "tidyr")   ## ADD
install.packages("stringi", "ggplot2", "tidyr", "reshape2")

library("devtools")
devtools::install_github("pachterlab/sleuth")   ## R version 3.3.2 (2016-10-31)

# -----------------------------------------------------------------------------
# Read in transcript TPM estimates/aboundants from kallisto (v0.43.1)
# Based on https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
# -----------------------------------------------------------------------------
getTSV <- function(wd.rna, samples) {
   return(file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias--fusion", samples))
}

wd <- "/ngs/cangen/tyang2"                     ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

wd.rna <- file.path(wd, "LUAD/ngs/RNA")
samples.luad <- readTable(file.path(wd.rna, "luad_rna_n49-1.list"), header=F, rownames=T, sep="")[,1]   ## ADD 28/04/18: S01728 is SCLC
tsv.luad <- getTSV(wd.rna, samples.luad)

wd.rna <- file.path(wd, "LCNEC/ngs/RNA")
samples.lcnec <- readTable(file.path(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=T, sep="")[,1]
tsv.lcnec <- getTSV(wd.rna, samples.lcnec)

wd.rna <- file.path(wd, "SCLC/ngs/RNA")
samples.sclc <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")[,1]
tsv.sclc <- getTSV(wd.rna, samples.sclc)

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ALL"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/RNA")
wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias--fusion")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- c(samples.luad, samples.lcnec, samples.sclc)   ## ALL (n=198)

# -----------------------------------------------------------------------------
# Associating transcripts to gene-level TPM estimates using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# By using sleuth's default filter settings: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.2.2 (on gauss)

tsv <- c(tsv.luad, tsv.lcnec, tsv.sclc)
s2c <- data.frame(path=tsv, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default filter settings
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 88145 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 20593 genes passed the filter
# summarizing bootstraps

tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm_r5_p47.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5_p47.RData")))

###
## Full gene list without any filtering
tpm.gene.patch <- list2Matrix(tpm.norm$tpm, tpm.norm)   ## Gene-level TPMs with patches

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
overlaps <- intersect(rownames(tpm.gene.patch), rownames(ensGene))
tpm.gene <- tpm.gene.patch[overlaps,]                   ## Gene-level TPMs
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
# > nrow(tpm.gene)
# [1] 34908

## Remove not expressed genes (ADD 10/06/18)
tpm.gene <- tpm.gene[getExpressed(tpm.gene),]
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_tpm0.RData")))
# > nrow(tpm.gene)
# [1] 15101
# > 15101 - 15084   ## Different with line 139
# [1] 17   ## Genes with no 0 TPM in any of the samples, but failed at least 5 read in 47% of the samples (i.e. very low-expressed genes)

###
## Gene list after default filtering
tpm.gene.patch <- list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt)   ## Gene-level TPMs with patches
# > nrow(tpm.gene.patch)
# [1] 20593   ## Matched to line 93

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
overlaps <- intersect(rownames(tpm.gene.patch), rownames(ensGene))
tpm.gene <- tpm.gene.patch[overlaps,]                             ## Gene-level TPMs
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
# > nrow(tpm.gene)
# [1] 18898

## Remove not expressed genes (ADD 10/06/18)
tpm.gene <- tpm.gene[getExpressed(tpm.gene),]
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5_p47_tpm0.RData")))
# > nrow(tpm.gene)
# [1] 15084

# =============================================================================
# Density plots
# Last Modified: 11/06/18
# =============================================================================
load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene)
plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.pdf")))

load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_tpm0.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene)
plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_tpm0.pdf")))

load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene)
plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5_p47.pdf")))

load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5_p47_tpm0.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene)
plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5_p47_tpm0.pdf")))
