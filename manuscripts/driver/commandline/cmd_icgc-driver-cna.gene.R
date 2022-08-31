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

load(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.median0.RData")))
tpm.gene.median0 <- tpm.gene
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

cna.gene <- getGeneCNV(samples.hist)
cna.gene.nona <- cna.gene[removeMissing(cna.gene),]
save(samples.hist, cna.gene, cna.gene.nona, file=file.path(wd.driver.data, paste0(HIST, "-driver-cna_iCN.seg_NONA.RData")), version=2)
dim(cna.gene)
# [1] 57773     6
dim(cna.gene.nona)
# [1] 57735     6
#}

tpm.gene.hist <- tpm.gene[, rownames(samples.hist)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.hist.log2   <- log2(tpm.gene.hist + 1)

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "M1", "M2", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.hist.log2), colnames)
rownames(src) <- rownames(tpm.gene.hist.log2)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.hist.log2), function(x) cor.test(as.numeric(tpm.gene.hist.log2[x,]), samples.hist$COR, method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.hist.log2), function(x) cor.test(as.numeric(tpm.gene.hist.log2[x,]), samples.hist$COR, method="spearman", exact=F)[[3]])
src <- src[!is.na(src$P),]

## Log2 fold change
#de$M1 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 0)))
#de$M2 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 1)))
#de$LOG2_FC <- de$M2 - de$M1

## FDR
#library(qvalue)
#src$Q <- qvalue(src$P)$qvalue
#src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, paste0("SRC_", HIST, "_tpm-gene_SORTING-vs-TPM_q_n", N, ".txt")), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples.hist, file=file.path(wd.de.data, paste0("SRC_", HIST, "_tpm-gene_SORTING-vs-TPM_q_n", N, ".RData")))
nrow(src.tpm.gene)
# [1] 11708

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and CNAs
# Last Modified: 26/06/22; 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
cna.gene.nona.tpm <- cna.gene.nona[intersect(rownames(cna.gene.nona), rownames(tpm.gene)),]
dim(cna.gene.nona.tpm)
# [1] 20707     6

cna.gene.nona.tpm.src <- cna.gene.nona.tpm[, rownames(samples.hist)]
tpm.gene.log2.cna <- tpm.gene.hist.log2[rownames(cna.gene.nona.tpm.src), rownames(samples.hist)]
#colnames(tpm.gene.log2.cna) <- samples.sclc.tpm$SAMPLE_ID
dim(cna.gene.nona.tpm.src)
# [1] 20707     6
dim(tpm.gene.log2.cna)
# [1] 20707     6

colnames <- c("RHO", "P", "Q", "DEL", "AMP", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2.cna), colnames)
rownames(src) <- rownames(tpm.gene.log2.cna)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) cor.test(as.numeric(tpm.gene.log2.cna[x,]), as.numeric(cna.gene.nona.tpm.src[x,]), method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) cor.test(as.numeric(tpm.gene.log2.cna[x,]), as.numeric(cna.gene.nona.tpm.src[x,]), method="spearman", exact=F)[[3]])

## Log2 fold change
#src$DEL <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) median(as.numeric(cna.gene.nona.tpm.src[x, colnames(cna.gene.nona.tpm.src)[which(cna.gene.nona.tpm.src[x,] < 2)]])))
#src$AMP <- mapply(x = 1:nrow(tpm.gene.log2.cna), function(x) median(as.numeric(cna.gene.nona.tpm.src[x, colnames(cna.gene.nona.tpm.src)[which(cna.gene.nona.tpm.src[x,] > 2)]])))
#src$LOG2_FC <- src$AMP - src$DEL
src <- src[!is.na(src$P),]

## FDR
#library(qvalue)
#src$Q <- qvalue(src$P)$qvalue
#src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!
nrow(src.tpm.gene)
# [1] 11643

writeTable(src.tpm.gene, file.path(wd.de.data, paste0("SRC_", HIST, "_tpm-gene_CNA-vs-TPM_q_n", N, ".txt")), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples.hist, file=file.path(wd.de.data, paste0("SRC_", HIST, "_tpm-gene_CNA-vs-TPM_q_n", N, ".RData")))
#writeRNKformat(src.tpm.gene, wd.de.gsea, "SRC_SCLC_tpm-gene-median0_CNA-vs-TPM_q_n70")   ## GSEA

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
plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", HIST, "_TPM-CNA-SORTING_P1E03"))
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(HIST, " (n=", N, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Expressed
expressed <- intersect(rownames(de), rownames(tpm.gene.median0))

de <- getCannoli(wd.de.data, HIST, N, expressed)
plot.de <- file.path(wd.de.plots, paste0("cannoliplot_SRC_", HIST, "_TPM-CNA-SORTING_P1E03_MEDIAN0"))
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(HIST, " expressed genes"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

## Not expressed
de <- getCannoli(wd.de.data, HIST, N, setdiff(rownames(src.tpm.gene), expressed))
plot.de <- file.path(wd.de.plots,  paste0("cannoliplot_SRC_", HIST, "_TPM-CNA-SORTING_P1E03_MEDIAN0-0"))
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0(HIST, " not expressed genes"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)
