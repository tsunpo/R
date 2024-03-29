# =============================================================================
# Manuscript   :
# Name         : manuscripts/expression/esad-tpm.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/08/20
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Human Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "ESAD"
base <- tolower(BASE)

#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--single--bias")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt"), header=T, rownames=T, sep="")
samples <- samples[rownames(subset(samples, GROUP_ID == 0)),][,1]

# -----------------------------------------------------------------------------
# From transcript-level estimates to gene-level TPMs using sleuth (v0.29.0)
# Based on https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html
#
# With sleuth's default quality-control filters: minimum 5 reads in at least 47% of the samples
# https://pachterlab.github.io/sleuth/docs/basic_filter.html
# https://groups.google.com/forum/#!topic/kallisto-sleuth-users/QrKxxEEFnE0
# -----------------------------------------------------------------------------
library("sleuth")   ## R version 3.3.2 (on MacBook)

tsv <- file.path(wd.rna.raw, samples)
s2c <- data.frame(path=tsv, sample=samples, stringsAsFactors=F)
t2g <- tx2Ens(ensGene.transcript)   ## Use full Ensembl transcripts with patches/scaffold sequences (*_PATCH) to avoid warning messages
                                    ## Please refer to line 49 (in guide-to-the/hg19.R)
so <- sleuth_prep(s2c, target_mapping=t2g, gene_mode=T, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=5, min_prop=0.47)   ## Default quality-control filters
# reading in kallisto results
# dropping unused factor levels
# .....................................................................
# normalizing est_counts
# 25530 targets passed the filter
# 22896 (N)
# 25362 (TR)
# 29106 (T)
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 14369 genes passed the filter
# 13719 (N)
# 14283 (TR)
# 15247 (T)
# summarizing bootstraps

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm_T.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r5p47_T.RData")))

## Gene-level TPMs (All)
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_T.RData")))
#writeTable(tpm.gene, gzfile(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.txt.gz"))), colnames=T, rownames=T, sep="\t")
nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs (Detected)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_T.RData")))
tpm.gene <- removeMedian0(tpm.gene)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0_T.RData")))
nrow(tpm.gene)
# [1] 15658
# [1] 14838 (N)
# [1] 15515 (TR)
# [1] 16500 (T)


## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47_N.RData")))
nrow(tpm.gene)
# [1] 13569
# [1] 12989 (N)
# [1] 13471 (TR)
# [1] 14367 (T)











## Gene-level TPMs (Detected + Expressed)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.median0 <- tpm.gene
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene.r5p47   <- tpm.gene
overlaps <- intersect(rownames(tpm.gene.median0), rownames(tpm.gene.r5p47))

tpm.gene <- tpm.gene[overlaps,]
#save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.r5p47.RData")))
nrow(tpm.gene)
# [1] 13569

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 29/05/20
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "All")

## Detected genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Detected")

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed")

## Expressed genes
#file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47"))
#load(paste0(file.main, ".RData"))
#plotDensityHistogram(tpm.gene, file.main, "Expressed")
