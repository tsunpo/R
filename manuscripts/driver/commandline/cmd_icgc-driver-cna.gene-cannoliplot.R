#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
H <- as.numeric(args[1])

# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/driver/icgc-driver-cna.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/03/22
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Asymmetry.R", "Mutation.R", "Survival.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.icgc     <- file.path(wd, BASE, "consensus")
wd.icgc.cna <- file.path(wd.icgc, "copy_number_alterations")

wd.meta     <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.driver <- file.path(wd.anlys, "driver", paste0(base, "-driver"))
wd.driver.data  <- file.path(wd.driver, "data")
wd.driver.plots <- file.path(wd.driver, "plots")
wd.rt <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

load(file=file.path(wd.rt.data, paste0("icgc_wgs_samples_n2612.RData")))
samplesSAMPLE_ID <- samples$icgc_specimen_id
nrow(samples)
# [1] 2612

release <- readTable(file.path(wd.meta, "data_release", "release_may2016.v1.4.tsv"), header=T, rownames=F, sep="")
rownames(release) <- release$tumor_wgs_aliquot_id
nrow(release)
# [1] 2834

list <- strsplit0(readTable(file.path(wd.icgc.cna, "copy_number_alterations.list"), header=F, rownames=F, sep=""), "_mutcall_filtered.vcf", 1)
length(list)
# [1] 2777

overlaps <- intersect(rownames(release), list)
length(overlaps)
# [1] 2595
release <- release[overlaps,]
rownames(release) <- release$tumor_wgs_icgc_specimen_id

overlaps <- intersect(rownames(samples), rownames(release))
length(overlaps)
# [1] 2443
release     <- release[overlaps,]
samples.cna <- samples[overlaps,]
samples.cna$tumor_wgs_aliquot_id <- release$tumor_wgs_aliquot_id
#samples.cna <- setProliferation(samples.cna, cor.rho)
nrow(samples.cna)
# [1] 2443
#writeTable(samples.cna[, c("icgc_specimen_id", "tumor_wgs_aliquot_id")], file.path(wd.meta, paste0("copy_number_alterations.txt")), colnames=F, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Assigning copy numbers for each Ensembl gene in each sample
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
## See cmd_icgc-driver-cna.R

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd.rna <- file.path(wd, BASE, "ngs/RNA")
#wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

#load(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.median0.RData")))
#tpm.gene.median0 <- tpm.gene
load(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.RData")))
dim(tpm.gene)
# [1] 20720  1359

overlaps <- intersect(samples.cna$icgc_specimen_id, colnames(tpm.gene))
#overlaps <- intersect(samples.surv$icgc_specimen_id, colnames(tpm.gene))
length(overlaps)
# [1] 902

samples.cna.tpm <- samples.cna[overlaps,]
tpm.gene.cna <- tpm.gene[, overlaps]
dim(samples.cna.tpm)
# [1] 902  14
dim(tpm.gene.cna)
# [1] 18502   902
 
icgc.cna.tpm <- as.data.frame(sort(table(samples.cna.tpm$histology_abbreviation), decreasing=T))
colnames(icgc.cna.tpm) <- c("histology_abbreviation", "N")
# > icgc.cna.tpm
# Var1 Freq
# 1        Lymph-BNHL  105
# 2         Liver-HCC   99
# 3     Ovary-AdenoCA   88
# 4    Breast-AdenoCA   81
# 5        Kidney-RCC   68
# 6  ColoRect-AdenoCA   50
# 7          Lung-SCC   47
# 8       Thy-AdenoCA   47
# 9      Kidney-ChRCC   43
# 10         Head-SCC   42
# 11   Uterus-AdenoCA   42
# 12     Lung-AdenoCA   37
# 13    Skin-Melanoma   36
# 14  Stomach-AdenoCA   29
# 15          CNS-GBM   24
# 16      Bladder-TCC   23
# 17    Prost-AdenoCA   19
# 18  Biliary-AdenoCA   16
# 19      Eso-AdenoCA    6

# -----------------------------------------------------------------------------
# Differential copy number alterations (CNA)
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
#for (h in 1:nrow(icgc.cna.tpm)) {
HIST <- as.vector(icgc.cna.tpm$histology_abbreviation)[H]
samples.hist <- subset(samples.cna.tpm, histology_abbreviation == HIST)
N <- nrow(samples.hist)
#samples.cna.tpm.histSAMPLE_ID <- samples.cna.tpm.hist$icgc_specimen_id

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
wd.rna <- file.path(wd, BASE, "ngs/RNA")
#wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

load(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.RData")))
dim(tpm.gene)

# -----------------------------------------------------------------------------
# Cannoli plot (P < 0.001)
# Last Modified: 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
xlab.text <- "Expression vs. CNA (rho)"
#ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")
ylab.text <- "Expression vs. proliferation (rho)"
pvalue <- 0.001

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes <- toTable(NA, length(colnames), 0, colnames)

## Total
de <- getCannoli(wd.de.data, HIST, N, NULL)
#plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", HIST, "_TPM-CNA-SORTING_P1E03"))
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.de <- paste0(plot.de, ".pdf")
#file.main <- c(paste0(HIST, " (n=", N, ")"), "")
#plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Expressed
total <- intersect(rownames(de), rownames(tpm.gene))
tpm.gene.cna <- tpm.gene[total, rownames(samples.hist)]
expressed <- rownames(removeMedian0(tpm.gene.cna))

de <- getCannoli(wd.de.data, HIST, N, expressed)
plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", HIST, "_TPM-CNA-SORTING_P1E03_MEDIAN0"))
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(HIST, " (n=", N, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Not expressed
de <- getCannoli(wd.de.data, HIST, N, setdiff(total, expressed))
plot.de <- file.path(wd.de.plots,  paste0("cannoliplot_SRC_", HIST, "_TPM-CNA-SORTING_P1E03_MEDIAN0-0"))
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(HIST, " not expressed"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)
