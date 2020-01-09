# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/lcnec-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 08/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
base <- tolower(BASE)

wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata/Tirosh 2016")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
wd.de.gsea  <- file.path(wd.de, "gsea")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
samples$M2 <- as.factor(samples$M2)

##
samples.rna <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")
overlaps <- intersect(rownames(samples.rna), rownames(samples))
samples <- samples[overlaps,]
# > length(which(samples$M2 == 1))
# [1] 35
# > length(which(samples$M2 == 0))
# [1] 35

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]
nrow(tpm.gene)
# [1] 19131
# nrow(tpm.gene)
# [1] 34908
# tpm.gene <- tpm.gene[getExpressed(tpm.gene),]
# nrow(tpm.gene)
# [1] 16573

#overlaps <- intersect(rownames(tpm.gene), rownames(wgs.gene))   ## After sclc-wgs-de-cm2.R
#tpm.gene <- tpm.gene[overlaps,]
tpm.gene.log2 <- log2(tpm.gene + 0.01)
# > nrow(tpm.gene.log2)
# [1] 19131
# > nrow(tpm.gene.log2)
# [1] 34908

# -----------------------------------------------------------------------------
# 
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
tpm.gene.log2 <- tpm.gene.log2[,rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!

colnames <- c("RHO", "P", "Q", "M1", "M2", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$M1 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 0)))
de$M2 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 1)))
de$LOG2_FC <- de$M2 - de$M1

## FDR
library(qvalue)
de$Q <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_n70.RData"))
nrow(de.tpm.gene)
# [1] 19131
# > nrow(de.tpm.gene)
# [1] 34908
# nrow(de.tpm.gene)
# [1] 16573
de.tpm.gene <- de.tpm.gene[!is.na(de.tpm.gene$P),]
writeTable(de.tpm.gene, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_n70.txt"), colnames=T, rownames=F, sep="\t")
# nrow(de.tpm.gene)
# [1] 19131
# nrow(de.tpm.gene)
# [1] 34055
# nrow(de.tpm.gene)
# [1] 16573

# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.sclc   ## MUY MUY IMPORTANTE!!

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene, 0.01)
tpm.gene.log2 <- tpm.gene.log2[rownames(de.tpm.gene),]
nrow(tpm.gene.log2)
# [1] 19131
# nrow(tpm.gene.log2)
# [1] 34055
# nrow(tpm.gene.log2)
# [1] 16573

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene.log2 <- getLog2andMedian(tpm.gene, 0.01)
#overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2))
#tpm.gene.log2 <- tpm.gene.log2[overlaps,]
#nrow(tpm.gene.log2)
# [1] 18123
# > nrow(tpm.gene.log2)
# [1] 16252

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
tpm.gene.log2 <- cbind(annot[rownames(tpm.gene.log2),], tpm.gene.log2[, "MEDIAN"], ensGene.bed[rownames(tpm.gene.log2),]) 
colnames(tpm.gene.log2)[8] <- "MEDIAN"

tpm.gene.log2$TSS_RFD <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$RFD
#tpm.gene.log2$TTS_RFD <- nrds.RT.NRFD[tpm.gene.log2$TTS,]$RFD
tpm.gene.log2$TSS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$SPLINE
#tpm.gene.log2$TTS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2$TTS,]$SPLINE

tpm.gene.log2 <- tpm.gene.log2[!is.na(tpm.gene.log2$TSS_RFD),]
nrow(tpm.gene.log2)
# [1] 18123   ## r5p47
# nrow(tpm.gene.log2)
# [1] 31207
# > nrow(tpm.gene.log2)
# [1] 15818   ## tpm0
# > nrow(tpm.gene.log2)
# [1] 16252   ## r5p47
#tpm.gene.log2 <- subset(tpm.gene.log2, gene_biotype == "protein_coding")
# > nrow(tpm.gene.log2)
# [1] 19007

##
tpm.gene.log2$TRC <- tpm.gene.log2$strand * tpm.gene.log2$TSS_RFD
tpm.gene.log2$TRC[which(tpm.gene.log2$TRC > 0)] <- 1
tpm.gene.log2$TRC[which(tpm.gene.log2$TRC < 0)] <- -1

## CTR (IZ)
tpm.gene.log2$TSS_NRFD <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$NRFD
tpm.gene.log2.ctr <- subset(subset(tpm.gene.log2, TSS_RFD < 0.9), TSS_RFD > -0.9)
tpm.gene.log2.ctr.iz <- subset(tpm.gene.log2.ctr, TSS_NRFD > 0)
tpm.gene.log2.ctr.iz.cd <- subset(tpm.gene.log2.ctr.iz, TRC > 0)
tpm.gene.log2.ctr.iz.ho <- subset(tpm.gene.log2.ctr.iz, TRC < 0)

tpm.gene.log2.ctr.tz <- subset(tpm.gene.log2.ctr, TSS_NRFD < 0)
tpm.gene.log2.ctr.tz.cd <- subset(tpm.gene.log2.ctr.tz, TRC > 0)
tpm.gene.log2.ctr.tz.ho <- subset(tpm.gene.log2.ctr.tz, TRC < 0)

## TTR (r5p47)
tpm.gene.log2.ttr <- rbind(subset(tpm.gene.log2, TSS_RFD >= 0.9), subset(tpm.gene.log2, TSS_RFD <= -0.9))
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr$MEDIAN)
# [1] 0.05385203
testU(tpm.gene.log2.ctr.iz$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.01539443
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.iz$MEDIAN)
# [1] 0.003022463
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.6190699

## TTR (All)
tpm.gene.log2.ttr <- rbind(subset(tpm.gene.log2, TSS_RFD >= 0.9), subset(tpm.gene.log2, TSS_RFD <= -0.9))
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr$MEDIAN)
# [1] 0.001668999
testU(tpm.gene.log2.ctr.iz$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 7.8635e-19
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.iz$MEDIAN)
# [1] 3.268986e-15
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 2.127784e-05

## TTR (tpm0)
tpm.gene.log2.ttr <- rbind(subset(tpm.gene.log2, TSS_RFD >= 0.9), subset(tpm.gene.log2, TSS_RFD <= -0.9))
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr$MEDIAN)
# [1] 0.009274047
testU(tpm.gene.log2.ctr.iz$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.13953
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.iz$MEDIAN)
# [1] 0.002927365
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.5111025

## TTR (PCG)
tpm.gene.log2.ttr <- rbind(subset(tpm.gene.log2, TSS_RFD >= 0.9), subset(tpm.gene.log2, TSS_RFD <= -0.9))
testU(tpm.gene.log2.ctr$MEDIAN, tpm.gene.log2.ttr$MEDIAN)
# [1] 0.02423306
testU(tpm.gene.log2.ctr.iz$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 5.201952e-06
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.iz$MEDIAN)
# [1] 7.280999e-06
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.05546897

## ALL
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2))
de.tpm.gene <- de.tpm.gene[overlaps,]
de.tpm.gene$Q <- qvalue(de.tpm.gene$P)$qvalue
#de.tpm.gene$FDR <- p.adjust(de.tpm.gene$P, method="BH", n=length(de.tpm.gene$P))
writeTable(de.tpm.gene, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.iz))
de.tpm.gene.iz <- de.tpm.gene[overlaps,]
de.tpm.gene.iz$Q <- qvalue(de.tpm.gene.iz$P)$qvalue
writeTable(de.tpm.gene.iz, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_iz_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.iz.ho))
de.tpm.gene.iz.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.iz.ho$Q <- qvalue(de.tpm.gene.iz.ho$P)$qvalue
writeTable(de.tpm.gene.iz.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_iz_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.iz.cd))
de.tpm.gene.iz.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.iz.cd$Q <- qvalue(de.tpm.gene.iz.cd$P)$qvalue
writeTable(de.tpm.gene.iz.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_iz_cd_n70_.txt"), colnames=T, rownames=F, sep="\t")

## TZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.tz))
de.tpm.gene.tz <- de.tpm.gene[overlaps,]
de.tpm.gene.tz$Q <- qvalue(de.tpm.gene.tz$P)$qvalue
writeTable(de.tpm.gene.tz, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_tz_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.tz.cd))
de.tpm.gene.tz.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.tz.cd$Q <- qvalue(de.tpm.gene.tz.cd$P)$qvalue
writeTable(de.tpm.gene.tz.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_tz_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.tz.ho))
de.tpm.gene.tz.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.tz.ho$Q <- qvalue(de.tpm.gene.tz.ho$P)$qvalue
writeTable(de.tpm.gene.tz.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_tz_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Figures
# Last Modified: 06/01/20
# -----------------------------------------------------------------------------
plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
   
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox2 <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names) {
   trait <- rep(0, length(tpm.1))
   trait <- c(trait, rep(1, length(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1, tpm.2))
   ylim <- c(min(expr), max(expr))
   
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

plotBox <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, names=names, ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
  
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

ylim <- c(min(tpm.gene.log2$MEDIAN), 15)
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_TTR-vs-IZ-vs-TZ_")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.ttr, tpm.gene.log2.ctr.iz, tpm.gene.log2.ctr.tz, main="TTR + CTR", names=c("TTR", "IZ", "TZ"), cols=c("black", "red", "blue"), ylim)

###
## IZ vs TZ
# > median(tpm.gene.log2.ctr.iz$MEDIAN)
# [1] 3.43802
# > median(tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 3.205719
# > testU(tpm.gene.log2.ctr.iz$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.01250994
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_IZ-vs-TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.ctr.iz, tpm.gene.log2.ctr.tz, main="CTR", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

## IZ (CD vs. HO)
# > median(tpm.gene.log2.ctr.iz.cd$MEDIAN)
# [1] 3.411179
# > median(tpm.gene.log2.ctr.iz.ho$MEDIAN)
# [1] 3.469866
# > testU(tpm.gene.log2.ctr.iz.cd$MEDIAN, tpm.gene.log2.ctr.iz.ho$MEDIAN)
# [1] 0.8621095
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_IZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.ctr.iz.ho, tpm.gene.log2.ctr.iz.cd, main="IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

## TZ (CD vs. HO)
# > median(tpm.gene.log2.ctr.tz.cd$MEDIAN)
# [1] 2.942382
# > median(tpm.gene.log2.ctr.tz.ho$MEDIAN)
# [1] 3.391783
# > testU(tpm.gene.log2.ctr.tz.cd$MEDIAN, tpm.gene.log2.ctr.tz.ho$MEDIAN)
# [1] 0.01135786
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_TZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.ctr.tz.ho, tpm.gene.log2.ctr.tz.cd, main="TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
plotCYS <- function(gene, cn, snr, pch, col, pos) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "In-silico sorting [rho]"
   ylab.text <- "log2(TPM + 0.01)"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.de.plots, paste0("TPM-vs-SORTING_", genes[g]))
   
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=paste0(gene, " (", id, ")"), col=col, pch=pch, cex=2, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
 
   #text(snr, n, samples, col=col, pos=c(1,3,3,3), cex=1.75)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=3)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col, bty="n", cex=1.75)
   #legend("bottomright", c("A", "B", "C"), col="black", bty="n", pt.cex=1.4, pch=c(1, 2, 0), horiz=T, cex=1.4)
 
   mtext(ylab.text, side=2, line=2.85, cex=1.7)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

##
genes <- c("GTF3C2", "RP11-141C7.3", "SUPT7L", "RAD9A", "E2F3", "ERCC8", "BLM", "POLE")
genes <- c("MYC", "MYCL", "MYCN")
genes <- c("PIF1", "EIF3B", "UBE2I", "BRCA2", "TOR1AIP1")
genes <- c("RAD9A", "PIF1", "BRCA2", "E2F3", "ERCC8", "BLM", "POLE")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 15/08/18
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, Q <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, fdr, genes, file.de, file.main, xlab.text) {
   #adjustcolor.gray <- adjustcolor("lightgray", alpha.f=0.9)  
   pvalue <- fdrToP(fdr, de)
   #fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #xmax <- 1
   ymax <- max(de$log10P)
   #ymax <- 4.5
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=xlab.text, ylab="-log10(p-value)", col="lightgray", main=file.main[1], cex=1.2, cex.axis=1.2, cex.lab=1.2, cex.main=1.25)

   abline(h=c(-log10(pvalue)), lty=5)
   text(xmax*-1 + 2*xmax/13, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)
   #text(xmax*-1 + 2*xmax/11, -log10(pvalue) - ymax/30, paste0("Q=1.00E-16"), cex=1.1)
   
   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="gold", cex=1.2)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col="dodgerblue", cex=1.2)
 
   if (nrow(genes) != 0) {
      for (g in 1:nrow(genes)) {
         gene <- subset(de, external_gene_name == genes[g,]$GENE)
         gene <- cbind(gene, genes[g,])
      
         if (nrow(gene) > 0) {
            points(gene$LOG2_FC, gene$log10P, pch=1, col="black", cex=1.2)
         
            if (!is.na(gene$ADJ_1))
               if (is.na(gene$ADJ_2))
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=1.2)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.2)
            else
               if (gene$LOG2_FC > 0)
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=1.2)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=1.2)
         } else
            print(genes[g])
      }
   }
   
   mtext(file.main[2], cex=1.25, line=0.3)
   legend("topleft", legend=c("Upregulated", "Downregulated"), col=c("gold", "dodgerblue"), pch=19)
   dev.off()
}

##
xlab.text <- "SCLC M2/M1 [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_fdr0.1")

## All genes
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.main <- c("TTR + CTR", paste0("(n=31,207)"))
file.main <- c("TTR + CTR", paste0("(n=18,123)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.14, genes, file.de, file.main, xlab.text)

## IZ
genes <- readTable(paste0(plot.de, "_iz.tab"), header=T, rownames=F, sep="\t")
#file.main <- c("Initiation (IZ)", paste0("(n=4,108)"))
file.main <- c("Initiation (IZ)", paste0("(n=2,598)"))
file.de <- paste0(plot.de, "_iz.pdf")
plotVolcano(de.tpm.gene.iz, 0.07, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# D.E.
# Last Modified: 08/01/20
# -----------------------------------------------------------------------------
#tpm.gene.log2.0      <- tpm.gene.log2[rownames(de.tpm.gene.iz),]
#tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[rownames(de.tpm.gene.iz),]
tpm.gene.log2.0      <- tpm.gene.log2[rownames(de.tpm.gene),]
tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[rownames(tpm.gene.log2.0),]

colnames <- c("P", "FDR", "N", "T", "LOG2_FC")
de2 <- toTable(0, length(colnames), nrow(tpm.gene.log2.0), colnames)
rownames(de2) <- rownames(tpm.gene.log2.0)

## SRC
de2$P <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) testU(as.numeric(tpm.gene.log2.0[x,]), as.numeric(tpm.gene.lung.log2.0[x,])))

## Log2 fold change
de2$N <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) median(as.numeric(tpm.gene.lung.log2.0[x,])))
de2$T <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) median(as.numeric(tpm.gene.log2.0[x,])))
de2$LOG2_FC <- de2$T - de2$N

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de2$FDR <- p.adjust(de2$P, method="BH", n=length(de2$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de2 <- de2[order(de2$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de2.tpm.gene <- cbind(annot[rownames(de2),], de2)   ## BE EXTRA CAREFUL!!

save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41.txt"), colnames=T, rownames=F, sep="\t")

## Volcano
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_q1e-06")

genes <- readTable(paste0(plot.de, "_iz_lung.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", paste0("(n=4,108)"))
file.de <- paste0(plot.de, "_iz_lung.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

de2.tpm.gene <- de2.tpm.gene[rownames(de.tpm.gene),]
genes <- readTable(paste0(plot.de, "_iz_lung_all.tab"), header=T, rownames=F, sep="\t")
file.main <- c("TTR + CTR", paste0("(n=31,207)"))
file.de <- paste0(plot.de, "_iz_lung_all.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
## PIF1
median(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]))
median(as.numeric(tpm.gene.log2["ENSG00000140451",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]), as.numeric(tpm.gene.log2["ENSG00000140451",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_PIF1")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000140451",], tpm.gene.log2["ENSG00000140451",], main="PIF1 (ENSG00000140451)", names=c("Lung", "SCLC"))

## BRCA2
median(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]))
median(as.numeric(tpm.gene.log2["ENSG00000139618",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]), as.numeric(tpm.gene.log2["ENSG00000139618",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRCA2")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139618",], tpm.gene.log2["ENSG00000139618",], main="BRCA2 (ENSG00000139618)", names=c("Lung", "SCLC"))

## BRD9
median(as.numeric(tpm.gene.lung.log2["ENSG00000028310",]))
median(as.numeric(tpm.gene.log2["ENSG00000028310",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000028310",]), as.numeric(tpm.gene.log2["ENSG00000028310",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRD9")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000028310",], tpm.gene.log2["ENSG00000028310",], main="BRD9 (ENSG00000028310)", names=c("Lung", "SCLC"))















## IZ / CD
genes <- readTable(paste0(plot.de, "_iz_cd.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", "Co-directional (CD)")
file.de <- paste0(plot.de, "_iz_cd.pdf")
plotVolcano(de.tpm.gene.iz.cd, 0.1, genes, file.de, file.main, xlab.text)

## IZ / HO
genes <- readTable(paste0(plot.de, "_iz_ho.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", "Head-on (HO)")
file.de <- paste0(plot.de, "_iz_ho.pdf")
plotVolcano(de.tpm.gene.iz.ho, 0.1, genes, file.de, file.main, xlab.text)

## TZ
genes <- readTable(paste0(plot.de, "_tz.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Termination", "")
file.de <- paste0(plot.de, "_tz.pdf")
plotVolcano(de.tpm.gene.tz, 0.1, genes, file.de, file.main, xlab.text)

## TZ / CD
genes <- readTable(paste0(plot.de, "_tz_cd.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Termination (TZ)", "Co-directional (CD)")
file.de <- paste0(plot.de, "_tz_cd.pdf")
plotVolcano(de.tpm.gene.tz.cd, 0.1, genes, file.de, file.main, xlab.text)

## TZ / HO
genes <- readTable(paste0(plot.de, "_tz_ho.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Termination (TZ)", "Head-on (HO)")
file.de <- paste0(plot.de, "_tz_ho.pdf")
plotVolcano(de.tpm.gene.tz.ho, 0.1, genes, file.de, file.main, xlab.text)







##
plot.main <- "Differential expression between M2 and M1 in SCLC"
xlab.text <- "log2FC(SCLC M2/M1)"
plot.de <- file.path(wd.de.plots, "volcanoplot-r5p47-wgs_sclc_p1e-4_m2_")

## E2F3
genes <- readTable(paste0(plot.de, "E2F3.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "")
#file.de <- paste0(plot.de, "_chr2_log2FC1.pdf")
file.de <- paste0(plot.de, "E2F3.pdf")
plotVolcano(de.tpm.gene, 1.00E-04, genes, file.de, file.main, xlab.text)
















# -----------------------------------------------------------------------------
# 
# Last Modified: 
# -----------------------------------------------------------------------------
tpm.gene.log2.rl <- subset(tpm.gene.log2, TSS_RFD > 0)
tpm.gene.log2.ll <- subset(tpm.gene.log2, TSS_RFD < 0)

tpm.gene.log2.rl.cd <- subset(tpm.gene.log2.rl, TRC > 0)
tpm.gene.log2.rl.ho <- subset(tpm.gene.log2.rl, TRC < 0)
tpm.gene.log2.ll.cd <- subset(tpm.gene.log2.ll, TRC > 0)
tpm.gene.log2.ll.ho <- subset(tpm.gene.log2.ll, TRC < 0)

tpm.gene.log2.cd <- subset(tpm.gene.log2, TRC > 0)
tpm.gene.log2.ho <- subset(tpm.gene.log2, TRC < 0)

##
plotTRC <- function(tpm, rfd, main.text, file.name, xlim, ylim, col, pos) {
 xlab.text <- "RFD"
 ylab.text <- "log2(TPM + 0.01)"
 
 pdf(paste0(file.name, ".pdf"), height=6, width=6)
 plot(tpm ~ rfd, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main=main.text[1], col=col[1], pch=1, cex=1, cex.axis=1.7, cex.lab=1.7, cex.main=1.8)
 
 lm.fit <- lm(tpm ~ rfd)
 abline(lm.fit, col=col[2], lwd=3)
 
 cor <- cor.test(tpm, rfd, method="spearman", exact=F)
 legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col[2], bty="n", cex=1.75)
 
 mtext(ylab.text, side=2, line=2.85, cex=1.7)
 mtext(main.text[2], cex=1.2, line=0.3)
 dev.off()
}

xlim <- c(min(tpm.gene.log2$TSS_RFD), max(tpm.gene.log2$TSS_RFD))
ylim <- c(min(tpm.gene.log2$MEDIAN), max(tpm.gene.log2$MEDIAN))

## Right-leading, Co-directional
file.name <- file.path(wd.rt.plots, "TRC_RFD_RL_CD")
plotTRC(tpm.gene.log2.rl.cd$MEDIAN, tpm.gene.log2.rl.cd$TSS_RFD, c("Right-leading", "Co-directional"), file.name, c(0, 1), ylim, col=c("sandybrown", "blue"), "topright")

file.name <- file.path(wd.rt.plots, "TRC_RFD_RL_HO")
plotTRC(tpm.gene.log2.rl.ho$MEDIAN, tpm.gene.log2.rl.ho$TSS_RFD, c("Right-leading", "Head-on"), file.name, c(0, 1), ylim, col=c("sandybrown", "red"), "topright")

## Left-leading, Co-directional
file.name <- file.path(wd.rt.plots, "TRC_RFD_LL_CD")
plotTRC(tpm.gene.log2.ll.cd$MEDIAN, tpm.gene.log2.ll.cd$TSS_RFD, c("Left-leading", "Co-directional"), file.name, c(-1, 0), ylim, col=c("steelblue1", "blue"), "topright")

file.name <- file.path(wd.rt.plots, "TRC_RFD_LL_HO")
plotTRC(tpm.gene.log2.ll.ho$MEDIAN, tpm.gene.log2.ll.ho$TSS_RFD, c("Left-leading", "Head-on"), file.name, c(-1, 0), ylim, col=c("steelblue1", "red"), "topright")











# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
## > length(genes.rb1.q0.1)
## [1] 639
## > length(genes.rb1.q0.05)
## [1] 145

## RB1 status on D.E genes
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene[genes.rb1.q0.1, rownames(samples)]
pca.de <- getPCA(t(test))

trait <- as.numeric(samples[,"RB1_MUT"])
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

##
file.main <- "LCNEC RB1 status on 145 D.E. (Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_RB1_Q0.05_145DE", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)

file.main <- "LCNEC RB1 status on 639 D.E. (Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_RB1_Q0.1_639DE", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)
