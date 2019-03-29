# =============================================================================
# Manuscript   :
# Chapter I    : RB1-loss additive effect on gene expression among LUAD, LCNEC and SCLC
# Name         : manuscripts/expression/luad+lcnec+sclc-tpm.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/01/19
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Human Genome
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
# Read in transcript-level estimates from kallisto (v0.43.1)
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
BASE <- "LUAD+LCNEC+SCLC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna <- file.path(wd, BASE, "ngs/RNA")
wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias--fusion")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- c(samples.luad, samples.lcnec, samples.sclc)   ## ALL (n=198)

# -----------------------------------------------------------------------------
# From transcript-level estimates to gene-level TPMs using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# With sleuth's default quality-control filters: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.2.2 (on gauss)

tsv <- c(tsv.luad, tsv.lcnec, tsv.sclc)
s2c <- data.frame(path=tsv, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)   ## Use full Ensembl transcripts with patches/scaffold sequences (*_PATCH) to avoid warning messages
                                    ## Please refer to line 49 (in guide-to-the/hg19.R)
so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default quality-control filters
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

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5p47.RData")))

## Gene-level TPMs without filtering
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
# > nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
# > nrow(tpm.gene)
# [1] 18898

# =============================================================================
# Density plots and histograms (See DifferentialExpression.R)
# Last Modified: 11/06/18
# =============================================================================
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))

tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=0.01)
plotDensity(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_density_pc0.01.pdf"), detected=T, pseudocount=0.01, NULL)
plotHistogram(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_hist_pc0.01.pdf"), detected=T, pseudocount=0.01, NULL)

## Expressed genes (with default filters)
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47"))
load(paste0(file.main, ".RData"))

tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=0.01)
plotDensity(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_density_pc0.01.pdf"), detected=F, pseudocount=0.01, NULL)
plotHistogram(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_hist_pc0.01.pdf"), detected=F, pseudocount=0.01, NULL)
