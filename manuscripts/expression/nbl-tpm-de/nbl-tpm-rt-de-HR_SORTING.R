# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/expression/sclc-tpm-rt-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 27/05/20
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R", "TranscriptionReplicationInteraction.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "NBL"
base <- tolower(BASE)

wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "Peifer 2015")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

## Testing for SORTING
samples.nbl.tpm.HR <- subset(samples.tpm.nbl, RISK == "high")

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median1.RData")))
tpm.gene <- tpm.gene[, rownames(samples.nbl.tpm.HR)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
#tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
dim(tpm.gene.log2)
# [1] 22899    37   ## HR

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "G1", "S", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(src) <- rownames(tpm.gene.log2)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.nbl.tpm.HR$COR, method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.nbl.tpm.HR$COR, method="spearman", exact=F)[[3]])
src <- src[!is.na(src$P),]

## Log2 fold change
#src$G1 <- median00(tpm.gene.log2, rownames(subset(samples.nbl.tpm, SORTING == "G1")))
#src$S  <- median00(tpm.gene.log2, rownames(subset(samples.nbl.tpm, SORTING == "S")))
#src$LOG2_FC <- src$S - src$G1

## FDR
library(qvalue)
src$Q <- qvalue(src$P)$qvalue
src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, "2015", "SRC_NBL_tpm-gene_SORTING-vs-TPM_q_n37.txt"), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples, file=file.path(wd.de.data, "2015", "SRC_NBL_tpm-gene_SORTING-vs-TPM_q_n37.RData"))
nrow(src.tpm.gene)
# [1] 22899

# -----------------------------------------------------------------------------
# Cannoli plot (P < 0.001)
# Last Modified: 04/07/22; 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
#ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")
xlab.text <- "Expression vs. SCF index [rho]"
ylab.text <- "Expression vs. CNA [rho]"
pvalue <- 0.001

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("MKI67", "MYCN", "TERT", "NTRK1")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
#genes[4, 2] <- 1

## Total
de <- getCannoli(file.path(wd.de.data, "2015"), BASE, 37, NULL, TEST="SORTING")
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
expressed <- intersect(rownames(de), rownames(tpm.gene))
length(expressed)
# [1] 22886

de <- getCannoli(file.path(wd.de.data, "2015"), BASE, 37, expressed, TEST="SORTING")
plot.de <- file.path(wd.de.plots, "2015", "cannoliplot_SRC_NBL-HR_TPM-CNA-SORTING_P1E03_MEDIAN0")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("All-HR (n=", 37, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

save(samples.nbl.tpm, samples.nbl.tpm.HR, expressed, file=file.path(wd.de.data, "2015", "samples.nbl.tpm_HR_n37.RData"))

# -----------------------------------------------------------------------------
# GSEA
# Last Modified: 15/07/22; 26/06/22; 11/10/17
# -----------------------------------------------------------------------------
de.pos.pos <- subset(subset(de, Effect2 > 0), Effect1 > 0)
de.pos.pos.sig <- subset(subset(de.pos.pos, P1 <= 0.001), P2 <= 0.001)
de.pos.neg <- subset(subset(de, Effect2 > 0), Effect1 < 0)
de.pos.neg.sig <- subset(subset(de.pos.neg, P1 <= 0.001), P2 <= 0.001)

de.neg.neg <- subset(subset(de, Effect2 < 0), Effect1 < 0)
de.neg.neg.sig <- subset(subset(de.neg.neg, P1 <= 0.001), P2 <= 0.001)
de.neg.pos <- subset(subset(de, Effect2 < 0), Effect1 > 0)
de.neg.pos.sig <- subset(subset(de.neg.pos, P1 <= 0.001), P2 <= 0.001)

writeRNKformatCNA(rbind(de.pos.pos, de.pos.neg), wd.de.gsea, "SRC_NBL-HR_tpm-gene-median0_SORTING-CNA-TPM_q_n37_GAIN")   ## GSEA
writeRNKformatCNA(rbind(de.neg.pos, de.neg.neg), wd.de.gsea, "SRC_NBL-HR_tpm-gene-median0_SORTING-CNA-TPM_q_n37_LOSS")   ## GSEA

# -----------------------------------------------------------------------------
# Correlation bwteen TPM and Purities
# Last Modified: 21/09/22
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "G1", "S", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(src) <- rownames(tpm.gene.log2)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.nbl.tpm.HR$purity2, method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.nbl.tpm.HR$purity2, method="spearman", exact=F)[[3]])
src <- src[!is.na(src$P),]

## Log2 fold change
#src$G1 <- median00(tpm.gene.log2, rownames(subset(samples.nbl.tpm, SORTING == "G1")))
#src$S  <- median00(tpm.gene.log2, rownames(subset(samples.nbl.tpm, SORTING == "S")))
#src$LOG2_FC <- src$S - src$G1

## FDR
library(qvalue)
src$Q <- qvalue(src$P)$qvalue
src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, "2015", "SRC_NBL_tpm-gene_Purity-vs-TPM_q_n37.txt"), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples.nbl.tpm.HR, file=file.path(wd.de.data, "2015", "SRC_NBL_tpm-gene_Purity-vs-TPM_q_n37.RData"))
#writeRNKformat(src.tpm.gene, wd.de.gsea, "SRC_NBL_tpm-gene-r5p47_SORTING_q_n54")   ## GSEA
nrow(src.tpm.gene)
# [1] 22899

# -----------------------------------------------------------------------------
# Cannoli plot (P < 0.001)
# Last Modified: 21/09/22
# -----------------------------------------------------------------------------
xlab.text <- "Expression vs. SCF index [rho]"
#ylab.text <- expression("Expression vs."~italic('in silico')~"sorting")
ylab.text <- "Expression vs. Purity [rho]"
pvalue <- 0.001

colnames <- c("GENE", "ADJ_1", "ADJ_2")
genes0 <- c("MKI67", "MYCN", "TERT", "NTRK1")
genes <- toTable(NA, length(colnames), length(genes0), colnames)
genes$GENE <- genes0
#genes[4, 2] <- 1
#genes[5, 2] <- 1

## Expressed
de <- getCannoli(file.path(wd.de.data, "2015"), BASE, 37, NULL, TEST="SORTING", TEST2="Purity")
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
expressed <- intersect(rownames(de), rownames(tpm.gene))

de <- getCannoli(file.path(wd.de.data, "2015"), BASE, 37, expressed, TEST="SORTING", TEST2="Purity", M2="")
plot.de <- file.path(wd.de.plots, "2015", "cannoliplot_SRC_NBL-HR_TPM-Purity-SORTING_P1E03_MEDIAN0_SIMPLES")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(paste0("All-HR (n=", 37, ")"), "")
plotCannoli(de, pvalue, genes, file.de, file.main, xlab.text, ylab.text, "bottomleft", c("", ""), c(red.lighter, blue.lighter), c(red, blue), fold=0, pos="bottomright")

# -----------------------------------------------------------------------------
# GSEA (Purity)
# Last Modified: 21/09/22
# -----------------------------------------------------------------------------
de.pos.pos <- subset(subset(de, Effect2 > 0), Effect1 > 0)
de.pos.pos.sig <- subset(subset(de.pos.pos, P1 <= 0.001), P2 <= 0.001)
de.pos.neg <- subset(subset(de, Effect2 > 0), Effect1 < 0)
de.pos.neg.sig <- subset(subset(de.pos.neg, P1 <= 0.001), P2 <= 0.001)

de.neg.neg <- subset(subset(de, Effect2 < 0), Effect1 < 0)
de.neg.neg.sig <- subset(subset(de.neg.neg, P1 <= 0.001), P2 <= 0.001)
de.neg.pos <- subset(subset(de, Effect2 < 0), Effect1 > 0)
de.neg.pos.sig <- subset(subset(de.neg.pos, P1 <= 0.001), P2 <= 0.001)

writeRNKformatCNA(rbind(de.pos.pos, de.pos.neg), wd.de.gsea, "SRC_NBL_tpm-gene-median0_SORTING-Purity-TPM_q_n54_GAIN")   ## GSEA
writeRNKformatCNA(rbind(de.neg.pos, de.neg.neg), wd.de.gsea, "SRC_NBL_tpm-gene-median0_SORTING-Purity-TPM_q_n54_LOSS")   ## GSEA

nrow(subset(de, Effect2 > 0)) / nrow(de)
# [1] 0.5548714
nrow(subset(de, Effect2 <= 0)) / nrow(de)
# [1] 0.4451286









## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "G1", "HALLMARK", "gsea_report_for_na_pos_1657897584856.tsv"), header=T, rownames=F, sep="\t")
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "HALLMARK_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[5:1, c("NAME", "NES")]
gsea.pos[5, 1] <- "E2F targets"
gsea.pos[4, 1] <- "G2M checkpoint"
gsea.pos[3, 1] <- firstup(gsea.pos[3, 1])
gsea.pos[2, 1] <- firstup(gsea.pos[2, 1])
gsea.pos[1, 1] <- "MTORC1 signaling"

main.text <- c("Hallmark enrichment score                         ", "")
file.de <- file.path(wd.de.gsea, "G1", "HALLMARK_POS.pdf")
pdf(file.de, height=2, width=4)
par(mar=c(2,11.5,2,2))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "G1", "HALLMARK", "gsea_report_for_na_neg_1657897584856.tsv"), header=T, rownames=F, sep="\t")
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[5:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) firstup(gsea.neg$NAME[x]))
gsea.neg[2, 1] <- "KRAS signaling UP"
gsea.neg[3, 1] <- "Interferon alpha res."
gsea.neg[4, 1] <- "Interferon gamma res."

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "G1", "HALLMARK_NEG.pdf")
pdf(file.de, height=1.6, width=4)
par(mar=c(2,5,0,8.5))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

genes <- c("MYCN")
for (g in 1:length(genes)) {
 id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
 plotSRC(file.path(wd.de.plots, "2015"), genes[g], as.numeric(tpm.gene.log2[id, rownames(samples.nbl.tpm)]), samples.nbl.tpm$COR, 1, pos="bottomright", col=blue, size=5)
}

for (g in 1:length(genes)) {
 id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
 plotCNA(file.path(wd.de.plots, "2015"), genes[g], as.numeric(tpm.gene.log2.cna[id, samples.nbl.tpm$SAMPLE_ID]), as.numeric(cna.gene.nona.tpm.src[id, samples.nbl.tpm$SAMPLE_ID]), 1, pos="bottomright", col="black", size=5)
}










## Not expressed
#de <- getCannoli(wd.de.data, BASE, 70, setdiff(rownames(src.tpm.gene), expressed))
#plot.de <- file.path(wd.de.plots, "cannoliplot_SRC_SCLC_TPM-CNA-SORTING_P1E03_MEDIAN0-0")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.de <- paste0(plot.de, ".pdf")
#file.main <- c(paste0(BASE, " not expressed genes"), "")
#plotCannoli(de, pvalue, toTable(NA, length(colnames), 0, colnames), file.de, file.main, xlab.text, ylab.text, "topleft", c("", ""), c(red, blue), c(red, blue), fold=0)

de.neg.pos <- subset(subset(de, Effect1 < 0), Effect2 > 0)
de.neg.pos.sig <- subset(subset(de.neg.pos, P1 <= 0.001), P2 <= 0.001)

de.pos.neg <- subset(subset(de, Effect1 > 0), Effect2 < 0)
de.pos.neg.sig <- subset(subset(de.pos.neg, P1 <= 0.001), P2 <= 0.001)




writeTable(de.negative.sig, file.path(wd.de.data, "SRC_NBL-G1_tpm-gene_SORTING-CNA-TPM_q_n29_neg10.txt"), colnames=T, rownames=T, sep="\t")

