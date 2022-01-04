# =============================================================================
# Manuscript   :
# Name         : manuscripts/expression/esad-tpm.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/08/20
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"              ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Mouse Genome
load(file.path(wd.src.ref, "mm10.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "SFB"
base <- tolower(BASE)

#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/projects"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.rna.raw <- file.path(wd.rna, "kallisto_mm10.ensembl_quant-b100--single--bias")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "sfb_3rna_n41.txt"), header=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# From transcript-level estimates to gene-level TPMs using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# With sleuth's default quality-control filters: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.3.2 (on MacBook)

tsv <- file.path(wd.rna.raw, samples$SAMPLE_ID)
s2c <- data.frame(path=tsv, sample=samples$SAMPLE_ID, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)   ## Use full Ensembl transcripts with patches/scaffold sequences (*_PATCH) to avoid warning messages
                                    ## Please refer to line 49 (in guide-to-the/hg19.R)
so <- sleuth_prep(s2c, target_mapping=t2g, gene_mode=T, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default quality-control filters
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 28012 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 15235 genes passed the filter
# summarizing bootstraps
# .........................................
# There were 50 or more warnings (use warnings() to see the first 50)

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5p47.RData")))

## Gene-level TPMs (All)
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
writeTable(tpm.gene, gzfile(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.txt.gz"))), colnames=T, rownames=T, sep="\t")
nrow(tpm.gene)
# [1] 35930

## Gene-level TPMs (Expressed)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
tpm.gene <- removeMedian0(tpm.gene)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
nrow(tpm.gene)
# [1] 18383

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
nrow(tpm.gene)
# [1] 15056

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 13/09/20; 29/05/20
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Total Ensembl")

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed")

## Consistently expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Filtered")
