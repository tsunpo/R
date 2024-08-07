# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/lusc-tpm.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 28/08/18
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "LUSC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                     ## tyang2@gauss
wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/RNA/Normal")
wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0("normal", "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "normal_rna_n41.list"), header=F, rownames=T, sep="")[,1]

# -----------------------------------------------------------------------------
# Associating transcripts to gene-level TPM estimates using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# By using sleuth's default filter settings: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")

tsv <- file.path(wd.rna.raw, samples)
s2c <- data.frame(path=tsv, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default filter settings
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 87040 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 21829 genes passed the filter
# summarizing bootstraps

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5p47.RData")))

###
## Gene-level TPMs without filtering
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
# > nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs (Detected)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene, 0)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
nrow(tpm.gene)
# [1] 23847

load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene, 1)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median1.RData")))
nrow(tpm.gene)
# [1] 15088

###
## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
# > nrow(tpm.gene)
# [1] 19872

tpm.gene.lung.log2 <- log2(tpm.gene + 1)

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










# -----------------------------------------------------------------------------
# Associating transcripts to gene-level TPM estimates using sleuth (v0.30.0)
# sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, gene_mode=T)   ## Different from v0.29.0
# kallisto_table(so, use_filtered=F, normalized=T, include_covariates=T)                                        ## Different from v0.29.0
#
# By using sleuth's default filter settings: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")

tsv <- c(tsv.sclc, tsv.normal)
s2c <- data.frame(path=tsv, sample=c(samples.sclc, samples.normal), stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47, gene_mode=T)   ## Default filter settings
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 89411 targets passed the filter
# normalizing tpm
# merging in metadata
# 21190 genes passed the filter
# summarizing bootstraps

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=T)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=T)
save(tpm.norm,      file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.norm.filt_r5p47.RData")))

###
## Gene-level TPMs without filtering
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.RData")))
nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs (Detected)
load(file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene)
save(tpm.gene, file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.median0.RData")))
nrow(tpm.gene)
# [1] 23737

load(file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene, 1)
save(tpm.gene, file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.median1.RData")))
nrow(tpm.gene)
# [1] 15565

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
nrow(tpm.gene)
# [1] 19411

# =============================================================================
# Density plots
# Last Modified: 11/06/18
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Total Ensembl")

## Expressed genes
file.main <- file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed")

## Expressed genes
file.main <- file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.median1"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed", tpm=1)

## Filtered genes
file.main <- file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.r5p47"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Filtered")






##
load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=1)
plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.pdf")), NULL)

load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=1)
plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.pdf")), NULL)
