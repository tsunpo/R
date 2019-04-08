# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/cll-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 23/05/18
# =============================================================================
wd.src <- "/ngs/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "CLL"
base <- tolower(BASE)

wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/RNA_WGS")
wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-wgs-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "EGAD00001002130.list"), header=F, rownames=T, sep="")[,1]

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

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=T, min_reads=30, min_prop=1.00)
# reading in kallisto results
# dropping unused factor levels
# ......................................................
# normalizing est_counts
# 73258 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 29799 genes passed the filter
# summarizing bootstraps

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
tpm.norm      <- kallisto_table(so, use_filtered=F, normalized=T, include_covariates=F)
tpm.norm.filt <- kallisto_table(so, use_filtered=T, normalized=T, include_covariates=F)
save(tpm.norm,      file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.RData")))
save(tpm.norm.filt, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.norm.filt_r30p100.RData")))

## Gene-level TPMs without filtering
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm$tpm, tpm.norm), ensGene)             ## Gene-level TPMs (without filtering)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.RData")))
# > nrow(tpm.gene)
# [1] 34908

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r30p100.RData")))
# > nrow(tpm.gene)
# [1] 28364

save(so, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_so.RData")))

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Last Modified: 25/11/18
# =============================================================================
## Detected genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))

tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=0.01)
plotDensity(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_density_pc0.01.pdf"), detected=T, pseudocount=0.01, NULL)
plotHistogram(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_hist_pc0.01.pdf"), detected=T, pseudocount=0.01, NULL)

## Expressed genes (with default filters)
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r30p100"))
load(paste0(file.main, ".RData"))

tpm.gene.log2 <- getLog2andMedian(tpm.gene, pseudocount=0.01)
plotDensity(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_density_pc0.01.pdf"), detected=F, pseudocount=0.01, NULL)
plotHistogram(tpm.gene.log2$MEDIAN, BASE, paste0(file.main, "_hist_pc0.01.pdf"), detected=F, pseudocount=0.01, NULL)
