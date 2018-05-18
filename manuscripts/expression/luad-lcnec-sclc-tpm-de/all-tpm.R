# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/all-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 30/11/17
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))

# -----------------------------------------------------------------------------
# Sleuth installation (v0.29.0)
# -----------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5", "tximport")
#install.packages("devtools", "stringr", "cluster", "survival", "tidyr")   ## ADD
#
#library("devtools")
#devtools::install_github("pachterlab/sleuth")   ## R version 3.3.2 (2016-10-31)

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2/"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/"   ## tyang2@local

wd.all <- paste0(wd, "ALL/")
wd.all.de <- paste0(wd.all, "analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/")
setwd(wd.all)

# -----------------------------------------------------------------------------
# Read in transcript TPM estimates/aboundants from kallisto (v0.43.1)
# Based on https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
# -----------------------------------------------------------------------------
getTSV <- function(wd, samples) {
   return(paste0(wd, "kallisto_hg19.ensembl_quant-b100--bias--fusion/", samples))
}

wd.rna <- paste0(wd, "LUAD/ngs/RNA/")
samples.luad <- readTable(paste0(wd.rna, "luad_rna_n49-1.list"), header=F, rownames=F, sep="\t")   ## ADD 28/04/18: S01728 is SCLC
tsvs.luad    <- getTSV(wd.rna, samples.luad)

wd.rna <- paste0(wd, "LCNEC/ngs/RNA/")
samples.lcnec <- readTable(paste0(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=F, sep="\t")
tsvs.lcnec    <- getTSV(wd.rna, samples.lcnec)

wd.rna <- paste0(wd, "SCLC/ngs/RNA/")
samples.sclc <- readTable(paste0(wd.rna, "sclc_rna_n81.list"), header=F, rownames=F, sep="\t")[,1]
tsvs.sclc    <- getTSV(wd.rna, samples.sclc)

# -----------------------------------------------------------------------------
# Associating transcripts to gene-level TPM estimates using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# By using sleuth's default filter settings: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.2.2 (on gauss)

samples <- c(samples.luad, samples.lcnec, samples.sclc)   ## ALL (n=198)
tsvs    <- c(tsvs.luad, tsvs.lcnec, tsvs.sclc)
s2c <- data.frame(path=tsvs, sample=samples, stringsAsFactors=F)
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
save(tpm.norm,      file=paste0(wd.all.de, "data/all_kallisto_0.43.1_tpm.norm_r5_p47.RData"))
save(tpm.norm.filt, file=paste0(wd.all.de, "data/all_kallisto_0.43.1_tpm.norm.filt_r5_p47.RData"))

tpm.gene.patch <- list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt)
# tpm.gene.patch <- kallisto_table_to_matrix(tpm.norm, min_reads=5, min_prop=0.47)
# > nrow(tpm.gene.patch)   ## Gene-level TPMs with patches
# [1] 20593 

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
overlaps <- intersect(rownames(tpm.gene.patch), rownames(ensGene))
tpm.gene <- tpm.gene.patch[overlaps,]
save(tpm.gene, file=paste0(wd.all.de, "data/all_kallisto_0.43.1_tpm.gene_r5_p47.RData"))
# > nrow(tpm.gene)         ## Gene-level TPMs
# [1] 18898
