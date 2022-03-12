# =============================================================================
# Manuscript   :
# Chapter      :
# Name         : manuscripts/expression/nbl-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/05/20
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"              ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"             ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# https://apple.stackexchange.com/questions/254380/why-am-i-getting-an-invalid-active-developer-path-when-attempting-to-use-git-a
# https://github.com/kharchenkolab/conos/wiki/Installing-Conos-for-Mac-OS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/RNA")
#wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

t2g <- tx2Ens(ensGene.transcript)

tpm <- readTable(file.path(wd.rna, "pcawg.rnaseq.transcript.expr.tpm.tsv.gz"), header=T, rownames=T, sep="")[,-1]
ids <- colnames(tpm)
ids <- gsub("\\.", "-", ids)
ids <- gsub("X", "", ids)

aliquots <- readTable(file.path(wd.rna, "rnaseq.extended.metadata.aliquot_id.V4.txt"), header=T, rownames=T, sep="\t")
overlaps <- intersect(ids, rownames(aliquots))
aliquots <- aliquots[overlaps,]
colnames(tpm) <- aliquots$icgc_specimen_id

genes <- unique(t2g$ens_gene)
tpm.gene <- toTable(NA, ncol(tpm), 0, colnames(tpm))
for (g in 1:length(genes)) {
   transcripts <- subset(t2g, ens_gene== genes[g])
   
   tpm.transcript <- tpm[grep(paste(transcripts$target_id, collapse="|"), rownames(tpm)),]
   if (nrow(tpm.transcript) > 0) {
      tpm.g <- toTable(NA, ncol(tpm), 1, colnames(tpm))
      rownames(tpm.g) <- genes[g]
      
      tpm.g[1,] <- mapply(x = 1:ncol(tpm.transcript), function(x) sum(tpm.transcript[, x]))
      tpm.gene <- rbind(tpm.gene, tpm.g)
   }
}

tpm.1 <- tpm[]






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

so <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column="ens_gene", extra_bootstrap_summary=F, min_reads=5, min_prop=0.47)   ## Default filter settings
# reading in kallisto results
# dropping unused factor levels
# ......................................................
# normalizing est_counts
# 82631 targets passed the filter
# normalizing tpm
# merging in metadata
# aggregating by column: ens_gene
# 20410 genes passed the filter
# summarizing bootstraps

## Transcript-level estimates with patches/scaffold sequences (*_PATCH)   ## See line 30 in guide-to-the/hg19.R
## https://www.ncbi.nlm.nih.gov/grc/help/patches
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
# [1] 22899

## Gene-level TPMs with default filters
tpm.gene <- getGeneTPM(list2Matrix(tpm.norm.filt$tpm, tpm.norm.filt), ensGene)   ## Gene-level TPMs (with default filters)
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
nrow(tpm.gene)
# [1] 18764

## Gene-level TPMs (Detected + Expressed)
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.median0 <- tpm.gene
load(file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene.r5p47   <- tpm.gene
overlaps <- intersect(rownames(tpm.gene.median0), rownames(tpm.gene.r5p47))

tpm.gene <- tpm.gene[overlaps,]
#save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0.r5p47.RData")))
nrow(tpm.gene)
# [1] 18764

# =============================================================================
# Density plot and histograms (See DifferentialExpression.R)
# Figure(s)    : Figure S2 (A and B)
# Last Modified: 06/09/20; 29/05/20
# =============================================================================
## All genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Total Ensembl")

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.median0"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed")

## Expressed genes
file.main <- file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene.r5p47"))
load(paste0(file.main, ".RData"))
plotDensityHistogram(tpm.gene, file.main, "Expressed")
