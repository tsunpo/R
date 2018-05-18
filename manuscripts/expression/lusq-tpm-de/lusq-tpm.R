# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/lusq-tpm.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 18/05/18
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"                  ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"                 ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd <- "/ngs/cangen/tyang2/LUSQ/"                     ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2/LUSQ/"   ## tyang2@local
wd.rna <- paste0(wd, "ngs/RNA/")
wd.rna.kallisto <- paste0(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias--fusion/")
wd.de <- paste0(wd, "analysis/expression/kallisto/lusq-tpm-de/")

samples <- readTable(paste0(wd.rna, "lusq_rna_n21.list"), header=F, rownames=F, sep="\t")[,1]

# -----------------------------------------------------------------------------
# Associating transcripts to gene-level TPM estimates using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# By using sleuth's default filter settings: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.2.2 (on gauss)

tsvs <- sapply(samples, function(s) file.path(wd.rna.kallisto, s))
s2c <- data.frame(path=tsvs, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default filter settings
# reading in kallisto results
# .....................................................................
# normalizing est_counts
# 80578 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 19995 genes passed the filter
# summarizing bootstraps

tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=paste0(wd.de, "data/lusq_kallisto_0.43.1_tpm.norm_r5_p47.RData"))
save(tpm.norm.filt, file=paste0(wd.de, "data/lusq_kallisto_0.43.1_tpm.norm.filt_r5_p47.RData"))

## Genes with patches
tpm.gene.patch <- list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt)
# tpm.gene.patch <- kallisto_table_to_matrix(tpm.norm, min_reads=5, min_prop=0.47)
# > nrow(tpm.gene.patch)   ## Gene-level TPMs with patches
# [1] 19995

## Remove patches (*_PATCH)
## https://www.ncbi.nlm.nih.gov/grc/help/patches
overlaps <- intersect(rownames(tpm.gene.patch), rownames(ensGene))
tpm.gene <- tpm.gene.patch[overlaps,]
save(tpm.gene, file=paste0(wd.de, "data/lusq_kallisto_0.43.1_tpm.gene_r5_p47.RData"))
# > nrow(tpm.gene)         ## Gene-level TPMs
# [1] 18357
