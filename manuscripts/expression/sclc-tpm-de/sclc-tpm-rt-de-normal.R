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
handbooks  <- c("Commons.R", "Transcription.R", "TranscriptionReplicationConflict.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

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

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples.rna <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")
samples.wgs <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
overlaps <- intersect(rownames(samples.rna), rownames(samples.wgs))

samples  <- samples.wgs[overlaps,]
samples$M2 <- as.factor(samples$M2)
# > length(which(samples$M2 == 1))
# [1] 35
# > length(which(samples$M2 == 0))
# [1] 35

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 0.01)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 0.01)
nrow(tpm.gene.log2.m)
# [1] 23572
# [1] 19131

# -----------------------------------------------------------------------------
# Prepare Normal and SCLC samples
# Last Modified: 26/04/22
# -----------------------------------------------------------------------------
samples.sclc.normal <- toTable("", 2, length(samples), c("SAMPLE_ID", "NORMAL"))
samples.sclc.normal$SAMPLE_ID <- samples
rownames(samples.sclc.normal) <- samples
samples.sclc.normal[samples.normal,]$NORMAL <- 0
samples.sclc.normal[samples.sclc,  ]$NORMAL <- 1
samples.sclc.normal$NORMAL <- as.factor(samples.sclc.lung$NORMAL)

load(file=file.path(wd.de.data, paste0("sclc+normal", "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples.sclc.normal)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2 <- log2(tpm.gene + 1)

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric)
# Last Modified: 24/04/22; 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TR (1) vs UN (0) as factor
argv      <- data.frame(predictor="NORMAL", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47_NORMAL_wilcox_q_n122")
file.main <- paste0("NORMAL (n=41) vs SCLC (n=81) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples.sclc.normal, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
xlab.text <- expression("-Log" * ""[2] * " SCLC to normal fold change")
ylab.text <- expression("-Log" * ""[10] * "(" * italic("P") * ") significance")

plot.de <- file.path(wd.de.plots, "volcanoplot_normal+sclc_r5p47")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("Differentially expressed genes", "")
legends <- c("SCLC", "Normal lung")
plotVolcano(de.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text, ylab.text, legends, fold=1)

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
genes <- c("TERT", "PIF1", "BLM", "DHX36", "DNAJC2")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   tpm.1 <- as.numeric(tpm.gene.log2[id, samples.normal])
   tpm.2 <- as.numeric(tpm.gene.log2[id, samples.sclc])
   
   file.name <- file.path(wd.de.plots, paste0("boxplot_TPM_", genes[g]))
   plotBox(file.name, tpm.1, tpm.2, genes[g], names=c("Normal", "SCLC"), cols=c(green, red))
}



# -----------------------------------------------------------------------------
# D.E. LUNG vs. SCLC (ALL)
# Last Modified: 25/04/22; 24/02/21; 25/10/20; 01/09/20; 08/01/20
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(tpm.gene.log2), rownames(tpm.gene.lung.log2))
length(overlaps)
# [1] 14013
tpm.gene.log2.o      <- tpm.gene.log2[overlaps,]
tpm.gene.lung.log2.o <- tpm.gene.lung.log2[overlaps,]

colnames <- c("P", "FDR", "N", "T", "LOG2_FC")
de3 <- toTable(0, length(colnames), nrow(tpm.gene.log2.o), colnames)
rownames(de3) <- rownames(tpm.gene.log2.o)

## SRC
de3$P <- mapply(x = 1:nrow(tpm.gene.log2.o), function(x) testU(as.numeric(tpm.gene.log2.o[x,]), as.numeric(tpm.gene.lung.log2.o[x,])))

## Log2 fold change
de3$N <- mapply(x = 1:nrow(tpm.gene.log2.o), function(x) median(as.numeric(tpm.gene.lung.log2.o[x,])))
de3$T <- mapply(x = 1:nrow(tpm.gene.log2.o), function(x) median(as.numeric(tpm.gene.log2.o[x,])))
de3$LOG2_FC <- de3$T - de3$N

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de3$FDR <- p.adjust(de3$P, method="BH", n=length(de3$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de3 <- de3[order(de3$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de3.tpm.gene <- cbind(annot[rownames(de3),], de3)   ## BE EXTRA CAREFUL!!

save(de3.tpm.gene, samples, file=file.path(wd.de.data, "de3_sclc_tpm-gene-median1-lung_U_BH_n70+41.RData"))
writeTable(de3.tpm.gene, file.path(wd.de.data, "de3_sclc_tpm-gene-medain1-lung_U_BH_n70+41.txt"), colnames=T, rownames=F, sep="\t")
nrow(de3.tpm.gene)

###
##
plotNormal <- function(gene, id, tpm.gene.log2, tpm.gene.lung.log2) {
   file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_", gene)
   plotBox2(wd.de.plots, file.name, as.numeric(tpm.gene.lung.log2[id,]), as.numeric(tpm.gene.log2[id,]), main=gene, names=c("Lung", "SCLC"))
}

genes <- c("RSAD2", "IRF2", "IFI35", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFIT5")
genes <- c("BRD9", "SMARCD1")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotNormal(genes[g], id, tpm.gene.log2, tpm.gene.lung.log2)
}

# -----------------------------------------------------------------------------
# D.E. LUNG vs. SCLC
# Last Modified: 25/10/20; 01/09/20; 08/01/20
# -----------------------------------------------------------------------------
plotVolcano0 <- function(de, pvalue, genes, file.de, file.main, xlab.text, ymax=0, cols, fc) {
 #pvalue <- fdrToP(fdr, de)
 #fdr <- pvalueToFDR(pvalue, de)
 de.sig <- subset(de, P <= pvalue)
 de.sig$log10P <- -log10(de.sig$P)
 
 de$log10P <- -log10(de$P)
 xlim <- c(min(de$LOG2_FC), max(de$LOG2_FC))
 if (ymax ==0) ymax <- max(de$log10P)
 
 pdf(file.de, height=6, width=6)
 plot(de$LOG2_FC, de$log10P, pch=16, xlim=xlim, ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col="gray88", main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
 #text(xmax*-1 + 2*xmax/15, -log10(pvalue) - ymax/30, paste0("FDR=", fdr, "%"), cex=1.15)    ## SCLC (IZ)
 #text(xmax*-1 + 2*xmax/13, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)   ## SCLC (AA)
 #text(xmax*-1 + 2*xmax/9.5, -log10(pvalue) - ymax/30, paste0("BH=1.00E-16"), cex=1.1)
 #text(xmax*-1 + 2*xmax/30, -log10(pvalue) + ymax/30, paste0("***"), cex=2)
 
 de.up   <- subset(de.sig, LOG2_FC > log2(fc))
 points(de.up$LOG2_FC, de.up$log10P, pch=16, col=cols[1], cex=1.4)
 de.down <- subset(de.sig, LOG2_FC < -log2(fc))
 points(de.down$LOG2_FC, de.down$log10P, pch=16, col=cols[2], cex=1.4)
 
 #abline(v=log2(fc), lty=5, lwd=1.8, col=cols[1]) 
 #abline(v=-log2(fc), lty=5, lwd=1.8, col=cols[2])
 abline(h=c(-log10(pvalue)), lty=5, lwd=1.5)
 
 if (nrow(genes) != 0) {
  for (g in 1:nrow(genes)) {
   gene <- subset(de, external_gene_name == genes[g,]$GENE)
   gene <- cbind(gene, genes[g,])
   
   if (nrow(gene) > 0) {
    points(gene$LOG2_FC, gene$log10P, pch=1, col="black", cex=1.4)
    
    if (!is.na(gene$ADJ_1))
     if (is.na(gene$ADJ_2))
      text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=1.25)
    else
     text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.25)
    else
     if (gene$LOG2_FC > 0)
      text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.6), cex=1.25)
    else
     text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.6), cex=1.25)
   } else
    print(genes[g])
  }
 }
 
 axis(side=1, at=seq(-10, 10, by=4), labels=c(-10, -6, -2, 2, 6, 10), cex.axis=1.2)
 axis(side=1, at=seq(-8, 8, by=4), labels=c(-8, -4, 0, 4, 8), cex.axis=1.2)
 mtext(file.main[2], cex=1.25, line=0.3)
 legend("topleft", legend=c("Upregulated (SCLC > Lung)", "Downregulated (SCLC < Lung)"), col=cols, pch=19, cex=1.25)
 dev.off()
}

## Volcano
xlab.text <- "SCLC to lung fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_iz_e_cd_p1e-15_lung")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("Differential expression between SCLC and lung", paste0("Early IZ, CD genes (n=1,159)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano0(de2.tpm.gene, 1E-15, genes, file.de, file.main, xlab.text, ymax=23, c(yellow, lightblue), 1)









# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
load(file=file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_n70.RData"))
dim(de.tpm.gene)
# [1] 19131    13

nrds.RT.NRFD <- nrds.RT.NRFD.sclc.nl   ## 17/05/20 WHY nrds.RT.NRFD.sclc.nl? ## MUY MUY IMPORTANTE!!
tpm.gene.log2.m <- tpm.gene.log2.m[rownames(de.tpm.gene),]
nrow(tpm.gene.log2.m)
# [1] 23572
# [1] 19131

tpm.gene.log2.m.rfd.nl <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD)
nrow(tpm.gene.log2.m.rfd.nl)
# [1] 22041
# [1] 18014
nrow(subset(tpm.gene.log2.m.rfd.nl, TRC == 0))
# [1] 1
# [1] 1

###
## TTR and CTR (IZ +TZ)
file.name <- file.path(wd.de.data, "tpm_gene_median0_log2_m_rfd_nl.RData")
rfd <- 0.9
#setTRC(tpm.gene.log2.m.rfd, rfd=0.9, file.name)

## CTR 
tpm.gene.log2.m.rfd.ctr.nl <- subset(subset(tpm.gene.log2.m.rfd.nl, TSS_RFD < rfd), TSS_RFD > -rfd)

tpm.gene.log2.m.rfd.ctr.iz.nl <- subset(tpm.gene.log2.m.rfd.ctr.nl, TSS_NRFD > 0)
tpm.gene.log2.m.rfd.ctr.tz.nl <- subset(tpm.gene.log2.m.rfd.ctr.nl, TSS_NRFD < 0)

## 31/08/20




tpm.gene.log2.m.rfd.ctr.iz.cd.nl <- subset(tpm.gene.log2.m.rfd.ctr.iz.nl, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho.nl <- subset(tpm.gene.log2.m.rfd.ctr.iz.nl, TRC < 0)
tpm.gene.log2.m.rfd.ctr.tz.cd.nl <- subset(tpm.gene.log2.m.rfd.ctr.tz.nl, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho.nl <- subset(tpm.gene.log2.m.rfd.ctr.tz.nl, TRC < 0)

## TTR
tpm.gene.log2.m.rfd.ttr.nl <- rbind(subset(tpm.gene.log2.m.rfd.nl, TSS_RFD >= rfd), subset(tpm.gene.log2.m.rfd.nl, TSS_RFD <= -rfd))
save(tpm.gene.log2.m.rfd.nl, tpm.gene.log2.m.rfd.ctr.nl, tpm.gene.log2.m.rfd.ctr.iz.nl, tpm.gene.log2.m.rfd.ctr.iz.cd.nl, tpm.gene.log2.m.rfd.ctr.iz.ho.nl, tpm.gene.log2.m.rfd.ctr.tz.nl, tpm.gene.log2.m.rfd.ctr.tz.cd.nl, tpm.gene.log2.m.rfd.ctr.tz.ho.nl, tpm.gene.log2.m.rfd.ttr.nl, file=file.name)

## CTR
length(which(tpm.gene.log2.m.rfd.ctr.iz.nl$GENE_NRFD < 0))
# [1] 43
length(which(tpm.gene.log2.m.rfd.ctr.tz.nl$GENE_NRFD > 0))
# [1] 52











# -----------------------------------------------------------------------------
# RFD vs. TPM (Tables)
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd))
de.tpm.gene.rfd <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd$Q <- qvalue(de.tpm.gene.rfd$P)$qvalue
#writeTable(de.tpm.gene.rfd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz))
de.tpm.gene.rfd.iz <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd.iz$Q <- qvalue(de.tpm.gene.rfd.iz$P)$qvalue
#writeTable(de.tpm.gene.rfd.iz, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_iz_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.cd))
de.tpm.gene.rfd.iz.cd <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd.iz.cd$Q <- qvalue(de.tpm.gene.rfd.iz.cd$P)$qvalue
#writeTable(de.tpm.gene.rfd.iz.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_iz_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.ho))
de.tpm.gene.rfd.iz.ho <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd.iz.ho$Q <- qvalue(de.tpm.gene.rfd.iz.ho$P)$qvalue
#writeTable(de.tpm.gene.rfd.iz.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_iz_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

## TZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz))
de.tpm.gene.rfd.tz <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd.tz$Q <- qvalue(de.tpm.gene.rfd.tz$P)$qvalue
#writeTable(de.tpm.gene.rfd.tz, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_tz_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.cd))
de.tpm.gene.rfd.tz.cd <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd.tz.cd$Q <- qvalue(de.tpm.gene.rfd.tz.cd$P)$qvalue
#writeTable(de.tpm.gene.rfd.tz.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_tz_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.ho))
de.tpm.gene.rfd.tz.ho <- de.tpm.gene[overlaps,]
#de.tpm.gene.rfd.tz.ho$Q <- qvalue(de.tpm.gene.rfd.tz.ho$P)$qvalue
#writeTable(de.tpm.gene.rfd.tz.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_src_q_rfd_tz_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# RFD vs. TPM (median0)
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_TTR+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC expression", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)

file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_TTR+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC expression", names=c("TTR", "TZ"), cols=c("black", "blue"), ylim)

file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_IZ+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC expression", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

##
file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_HO+CD_IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.iz.cd, main="SCLC IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_HO+CD_TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.tz.cd, main="SCLC TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)

# -----------------------------------------------------------------------------
# RFD vs. TPM (r5p47)
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TTR+IZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="SCLC expressed genes", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)

file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_TTR+TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC expressed genes", names=c("TTR", "TZ"), cols=c("black", "blue"), ylim)

file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_IZ+TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="SCLC expressed genes", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

##
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_HO+CD_IZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.iz.cd, main="SCLC expressed IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_rfd_HO+CD_TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.tz.cd, main="SCLC expressed TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
genes <- c("PIF1", "MARS", "KIF18B", "BRCA2", "RAD9A", "TOR1AIP1", "TOR1A", "TOR1B", "TP53", "RB1")
genes <- c("RP3-407E4.3", "BRD9")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}

# -----------------------------------------------------------------------------
# Gene length vs. RFD slopes
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
tpm.gene.log2.m.rfd.ctr.iz <- cbind(de.tpm.gene.rfd.iz, tpm.gene.log2.m.rfd.ctr.iz[rownames(de.tpm.gene.rfd.iz),])
tpm.gene.log2.m.rfd.ctr.iz$length <- abs(tpm.gene.log2.m.rfd.ctr.iz$end_position - tpm.gene.log2.m.rfd.ctr.iz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.iz$length)))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD) * -1)
ylab.text <- "IZ efficiency"

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz$length), c("SCLC IZ genes", "n=2,565"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.cd$length), c("SCLC IZ, CD genes", "n=1,275"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.ho$length), c("SCLC IZ, HO genes", "n=1,289"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

##
tpm.gene.log2.m.rfd.ctr.tz <- cbind(de.tpm.gene.rfd.tz, tpm.gene.log2.m.rfd.ctr.tz[rownames(de.tpm.gene.rfd.tz),])
tpm.gene.log2.m.rfd.ctr.tz$length <- abs(tpm.gene.log2.m.rfd.ctr.tz$end_position - tpm.gene.log2.m.rfd.ctr.tz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.tz$length)))
ylim <- c(max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD)) * -1
ylab.text <- "TZ efficiency"

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz$length), c("SCLC, TZ-S genes", "n=1,577"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.cd$length), c("SCLC TZ, CD genes", "n=802"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.ho$length), c("SCLC TZ, HO genes", "n=775"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)









# -----------------------------------------------------------------------------
# 
# Last Modified: 01/21/20
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz), rownames(tpm.gene.log2.m.rfd.ctr.iz.nl))
tpm.gene.log2.m.rfd.ctr.iz.sp <- tpm.gene.log2.m.rfd.ctr.iz[setdiff(rownames(tpm.gene.log2.m.rfd.ctr.iz), overlaps),]

tpm.gene.log2.m.rfd.ctr.iz.sp <- tpm.gene.log2.m.rfd.ctr.iz.sp[intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.sp), rownames(de.tpm.gene)),]
tpm.gene.log2.m.rfd.ctr.iz.nsp <- tpm.gene.log2.m.rfd.ctr.iz[overlaps,]
tpm.gene.log2.m.rfd.ctr.iz.nsp <- tpm.gene.log2.m.rfd.ctr.iz.nsp[intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.nsp), rownames(de.tpm.gene)),]

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.sp))
de.tpm.gene.rfd.iz.sp <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.sp$Q <- qvalue(de.tpm.gene.rfd.iz.sp$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.sp, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_iz_sp_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ-SP (HO vs CD)  ## 2020/05/17
de.tpm.gene.rfd.iz.sp <- readTable(file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_iz_sp_n70.txt"), header=T, rownames=T, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene.rfd.iz.sp), rownames(subset(tpm.gene.log2.m.rfd.ctr.iz, TRC == 1)))
de.tpm.gene.rfd.iz.sp.cd <- cbind(de.tpm.gene.rfd.iz.sp[overlaps,], tpm.gene.log2.m.rfd.ctr.iz[overlaps,])
de.tpm.gene.rfd.iz.sp.cd$Q <- qvalue(de.tpm.gene.rfd.iz.sp.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.sp.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_iz_sp_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene.rfd.iz.sp), rownames(subset(tpm.gene.log2.m.rfd.ctr.iz, TRC == -1)))
de.tpm.gene.rfd.iz.sp.ho <- cbind(de.tpm.gene.rfd.iz.sp[overlaps,], tpm.gene.log2.m.rfd.ctr.iz[overlaps,])
de.tpm.gene.rfd.iz.sp.ho$Q <- qvalue(de.tpm.gene.rfd.iz.sp.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.sp.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_iz_sp_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ (CD vs. HO)
# > median(de.tpm.gene.rfd.iz.sp.cd$MEDIAN)
# [1] 3.47189
# > median(de.tpm.gene.rfd.iz.sp.ho$MEDIAN)
# [1] 3.397869
# > testU(de.tpm.gene.rfd.iz.sp.cd$MEDIAN, de.tpm.gene.rfd.iz.sp.ho$MEDIAN)
# [1] 0.5321021
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_IZ_SP_HO-vs-CD")
plotBox(wd.de.plots, file.name, de.tpm.gene.rfd.iz.sp.ho, de.tpm.gene.rfd.iz.sp.cd, main="IZ-S", names=c("HO", "CD"), cols=c("red", "red"), ylim)

## TZ-SP (HO vs CD)  ## 2020/05/17
de.tpm.gene.rfd.tz.sp <- readTable(file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_n70.txt"), header=T, rownames=T, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene.rfd.tz.sp), rownames(subset(tpm.gene.log2.m.rfd.ctr.tz, TRC == 1)))
de.tpm.gene.rfd.tz.sp.cd <- cbind(de.tpm.gene.rfd.tz.sp[overlaps,], tpm.gene.log2.m.rfd.ctr.tz[overlaps,])
de.tpm.gene.rfd.tz.sp.cd$Q <- qvalue(de.tpm.gene.rfd.tz.sp.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene.rfd.tz.sp), rownames(subset(tpm.gene.log2.m.rfd.ctr.tz, TRC == -1)))
de.tpm.gene.rfd.tz.sp.ho <- cbind(de.tpm.gene.rfd.tz.sp[overlaps,], tpm.gene.log2.m.rfd.ctr.tz[overlaps,])
de.tpm.gene.rfd.tz.sp.ho$Q <- qvalue(de.tpm.gene.rfd.tz.sp.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

## TZ (CD vs. HO)
# > median(de.tpm.gene.rfd.tz.sp.cd$MEDIAN)
# [1] 3.125119
# > median(de.tpm.gene.rfd.tz.sp.ho$MEDIAN)
# [1] 3.504559
# > testU(de.tpm.gene.rfd.tz.sp.cd$MEDIAN, de.tpm.gene.rfd.tz.sp.ho$MEDIAN)
# [1] 0.0518692
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_TZ_SP_HO-vs-CD")
plotBox(wd.de.plots, file.name, de.tpm.gene.rfd.tz.sp.ho, de.tpm.gene.rfd.tz.sp.cd, main="TZ-S", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)





###
##
overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.tz), rownames(tpm.gene.log2.m.rfd.ctr.tz.nl))
tpm.gene.log2.m.rfd.ctr.tz.sp <- tpm.gene.log2.m.rfd.ctr.tz[setdiff(rownames(tpm.gene.log2.m.rfd.ctr.tz), overlaps),]

tpm.gene.log2.m.rfd.ctr.tz.sp <- tpm.gene.log2.m.rfd.ctr.tz.sp[intersect(rownames(tpm.gene.log2.m.rfd.ctr.tz.sp), rownames(de.tpm.gene)),]

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.sp))
de.tpm.gene.rfd.tz.sp <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.sp$Q <- qvalue(de.tpm.gene.rfd.tz.sp$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_n70.txt"), colnames=T, rownames=F, sep="\t")

## IZ (Overlapping SCLC and NBL)
overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.sp), rownames(de.tpm.gene.rfd.iz.s))


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
   pvalue <- fdrToP(fdr, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #ymax <- max(de$log10P)
   ymax <- 4.5
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=xlab.text, ylab="-log10(p-value)", col="lightgray", main=file.main[1], cex=1.2, cex.axis=1.2, cex.lab=1.2, cex.main=1.25)

   abline(h=c(-log10(pvalue)), lty=5)
   text(xmax*-1 + 2*xmax/15, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)    ## SCLC (IZ)
   #text(xmax*-1 + 2*xmax/13, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)   ## SCLC (AA)
   #text(xmax*-1 + 2*xmax/9.5, -log10(pvalue) - ymax/30, paste0("BH=1.00E-16"), cex=1.1)
   
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
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.1")

## All genes
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC all genes", paste0("(n=18,007)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd, 0.14, genes, file.de, file.main, xlab.text)

## IZ
genes <- readTable(paste0(plot.de, "_iz.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", paste0("(n=2,565)"))
file.de <- paste0(plot.de, "_iz.pdf")
plotVolcano(de.tpm.gene.rfd.iz, 0.08, genes, file.de, file.main, xlab.text)

## IZ (SP)
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_rfd_fdr0.1")
genes <- readTable(paste0(plot.de, "_iz_sp.tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC-specific initiated genes (IZ-S)", paste0("(n=2,045)"))
file.de <- paste0(plot.de, "_iz_sp.pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp, 0.1, genes, file.de, file.main, xlab.text)

## IZ (Expressed, SP)
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08")
genes <- readTable(paste0(plot.de, "_iz_sp2.tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC-specificly initiated (IZ-S), expressed genes", paste0("(n=1,299)"))
file.de <- paste0(plot.de, "_iz_sp2.pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp, 0.08, genes, file.de, file.main, xlab.text)

## IZ (Expressed, SP, CD)
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08_iz_sp_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC, IZ-S and CD genes", paste0("n=642"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp.cd, 0.08, genes, file.de, file.main, xlab.text)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08_iz_sp_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC, IZ-S and HO genes", paste0("n=657"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp.ho, 0.08, genes, file.de, file.main, xlab.text)






















###
##
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47.rfd_IZ-NSP-vs-IZ-SP_3")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.nsp, tpm.gene.log2.m.rfd.ctr.iz.sp, main="SCLC-specific IZ (IZ-S)", names=c("IZ-NS", "IZ-S"), cols=c("red", "red"), ylim)

###
## IZ vs TZ
# > median(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN)
# [1] 3.455006
# > median(tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 3.236536
# > testU(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 0.0180636
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_IZ-vs-TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="CTR", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

## IZ (CD vs. HO)
# > median(tpm.gene.log2.m.rfd.ctr.iz.cd$MEDIAN)
# [1] 3.44142
# > median(tpm.gene.log2.m.rfd.ctr.iz.ho$MEDIAN)
# [1] 3.470568
# > testU(tpm.gene.log2.m.rfd.ctr.iz.cd$MEDIAN, tpm.gene.log2.m.rfd.ctr.iz.ho$MEDIAN)
# [1] 0.6792462
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_IZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.iz.cd, main="IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

## TZ (CD vs. HO)
# > median(tpm.gene.log2.m.rfd.ctr.tz.cd$MEDIAN)
# [1] 0.6792462
# > median(tpm.gene.log2.m.rfd.ctr.tz.ho$MEDIAN)
# [1] 2.994693
# > testU(tpm.gene.log2.m.rfd.ctr.tz.cd$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz.ho$MEDIAN)
# [1] 0.008081627
file.name <- paste0("boxplot_sclc_tpm.gene.rfd_TZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.tz.cd, main="TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)











# -----------------------------------------------------------------------------
# Gene length vs. RFD slopes
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
tpm.gene.log2.m.rfd.ctr.iz <- cbind(de.tpm.gene.rfd.iz.sp, tpm.gene.log2.m.rfd.ctr.iz[rownames(de.tpm.gene.rfd.iz.sp),])
tpm.gene.log2.m.rfd.ctr.iz$length <- abs(tpm.gene.log2.m.rfd.ctr.iz$end_position - tpm.gene.log2.m.rfd.ctr.iz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.iz$length)))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD) * -1)
ylab.text <- "IZ efficiency"

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz$length), c("SCLC, IZ-S genes", "n=1,299"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.cd$length), c("SCLC, IZ-S and CD genes", "n=642"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.ho$length), c("SCLC, IZ-S and HO genes", "n=657"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

##
tpm.gene.log2.m.rfd.ctr.tz <- cbind(de.tpm.gene.rfd.tz.sp, tpm.gene.log2.m.rfd.ctr.tz[rownames(de.tpm.gene.rfd.tz.sp),])
tpm.gene.log2.m.rfd.ctr.tz$length <- abs(tpm.gene.log2.m.rfd.ctr.tz$end_position - tpm.gene.log2.m.rfd.ctr.tz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.tz$length)))
ylim <- c(max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD)) * -1
ylab.text <- "TZ efficiency"

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz$length), c("SCLC, TZ-S genes", "n=874"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.cd$length), c("SCLC, TZ-S and CD genes", "n=446"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "SCLC-SP_TRC_LENGTH-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.ho$length), c("SCLC, TZ-S and HO genes", "n=428"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)




###
##
tpm.gene.log2.m.rfd.ctr.iz <- cbind(tpm.gene.log2.m.rfd.ctr.iz, de.tpm.gene[rownames(tpm.gene.log2.m.rfd.ctr.iz),])
tpm.gene.log2.m.rfd.ctr.iz$length <- abs(tpm.gene.log2.m.rfd.ctr.iz$end_position - tpm.gene.log2.m.rfd.ctr.iz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.iz$length)))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD) * -1)
ylab.text <- "Initiation efficiency"

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz$length), c("SCLC initiation zone (IZ)", "(n=2,565)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.cd$length), c("SCLC IZ (CD) genes", "(n=1,275)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.ho$length), c("SCLC IZ (HO) genes", "(n=1,289)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

##
tpm.gene.log2.m.rfd.ctr.tz <- cbind(tpm.gene.log2.m.rfd.ctr.tz, de.tpm.gene[rownames(tpm.gene.log2.m.rfd.ctr.tz),])
tpm.gene.log2.m.rfd.ctr.tz$length <- abs(tpm.gene.log2.m.rfd.ctr.tz$end_position - tpm.gene.log2.m.rfd.ctr.tz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.tz$length)))
ylim <- c(max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD)) * -1
ylab.text <- "Termination efficiency"

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz$length), c("SCLC termination zone (TZ) genes", "(n=1,577)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.cd$length), c("SCLC TZ (CD) genes", "(n=802)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "TRC_LENGTH-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.ho$length), c("SCLC TZ (HO) genes", "(n=775)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)
















###
## RHO
xlim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$RHO), max(tpm.gene.log2.m.rfd.ctr.iz$RHO))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), max(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD))

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.iz$RHO, c("Initiation zone (IZ)", "(n=2,565)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.iz.cd$RHO, c("IZ (CD)", "(n=1,275)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.iz.ho$RHO, c("IZ (HO)", "(n=1,289)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

##
xlim <- c(min(tpm.gene.log2.m.rfd.ctr.tz$RHO), max(tpm.gene.log2.m.rfd.ctr.tz$RHO))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD))

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.tz$RHO, c("Termination zone (TZ)", "(n=1,577)"), file.name, xlim, ylim, col=c("blue", "black"), "bottomright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.tz.cd$RHO, c("TZ (CD)", "(n=802)"), file.name, xlim, ylim, col=c("blue", "black"), "bottomright")

file.name <- file.path(wd.de.plots, "TRC_RHO-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD, tpm.gene.log2.m.rfd.ctr.tz.ho$RHO, c("TZ (HO)", "(n=775)"), file.name, xlim, ylim, col=c("blue", "black"), "bottomright")











ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), max(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD))
xlim <- c(min(-log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(-log10(tpm.gene.log2.m.rfd.ctr.iz$)))

## 
file.name <- file.path(wd.de.plots, "TRC_P-vs-G-NRFD_IZ")
plotTRC(-log10(de.tpm.gene.iz$P), tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, c("Initiation (IZ)", "(n=2,381)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_P-vs-G-NRFD_IZ_HO")
plotTRC(-log10(de.tpm.gene.iz.ho$P), tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, c("IZ (HO)", "(n=1,198)"), file.name, xlim, ylim, col=c("red", "black"), "topright")

file.name <- file.path(wd.de.plots, "TRC_P-vs-G-NRFD_IZ_CD")
plotTRC(-log10(de.tpm.gene.iz.cd$P), tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, c("IZ (CD)", "(n=1,182)"), file.name, xlim, ylim, col=c("red", "black"), "topright")



# -----------------------------------------------------------------------------
# D.E. LUNG vs. SCLC (ALL)
# Last Modified: 24/02/21; 25/10/20; 01/09/20; 08/01/20
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(tpm.gene.log2), rownames(tpm.gene.lung.log2))
length(overlaps)
# [1] 18265
tpm.gene.log2.o      <- tpm.gene.log2[overlaps,]
tpm.gene.lung.log2.o <- tpm.gene.lung.log2[overlaps,]

colnames <- c("P", "FDR", "N", "T", "LOG2_FC")
de3 <- toTable(0, length(colnames), nrow(tpm.gene.log2.o), colnames)
rownames(de3) <- rownames(tpm.gene.log2.o)

## SRC
de3$P <- mapply(x = 1:nrow(tpm.gene.log2.o), function(x) testU(as.numeric(tpm.gene.log2.o[x,]), as.numeric(tpm.gene.lung.log2.o[x,])))

## Log2 fold change
de3$N <- mapply(x = 1:nrow(tpm.gene.log2.o), function(x) median(as.numeric(tpm.gene.lung.log2.o[x,])))
de3$T <- mapply(x = 1:nrow(tpm.gene.log2.o), function(x) median(as.numeric(tpm.gene.log2.o[x,])))
de3$LOG2_FC <- de3$T - de3$N

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de3$FDR <- p.adjust(de3$P, method="BH", n=length(de3$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de3 <- de3[order(de3$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de3.tpm.gene <- cbind(annot[rownames(de3),], de3)   ## BE EXTRA CAREFUL!!

save(de3.tpm.gene, samples, file=file.path(wd.de.data, "de3_sclc_tpm-gene-r5p47-lung_U_BH_n70+41.RData"))
writeTable(de3.tpm.gene, file.path(wd.de.data, "de3_sclc_tpm-gene-r5p47-lung_U_BH_n70+41.txt"), colnames=T, rownames=F, sep="\t")
nrow(de3.tpm.gene)

###
##
plotNormal <- function(gene, id, tpm.gene.log2, tpm.gene.lung.log2) {
   file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_", gene)
   plotBox2(wd.de.plots, file.name, as.numeric(tpm.gene.lung.log2[id,]), as.numeric(tpm.gene.log2[id,]), main=gene, names=c("Lung", "SCLC"))
}

genes <- c("RSAD2", "IRF2", "IFI35", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFIT5")
genes <- c("BRD9", "SMARCD1")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotNormal(genes[g], id, tpm.gene.log2, tpm.gene.lung.log2)
}

# -----------------------------------------------------------------------------
# D.E. LUNG vs. SCLC
# Last Modified: 25/10/20; 01/09/20; 08/01/20
# -----------------------------------------------------------------------------
plotVolcano0 <- function(de, pvalue, genes, file.de, file.main, xlab.text, ymax=0, cols, fc) {
   #pvalue <- fdrToP(fdr, de)
   #fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xlim <- c(min(de$LOG2_FC), max(de$LOG2_FC))
   if (ymax ==0) ymax <- max(de$log10P)
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=xlim, ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col="gray88", main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   #text(xmax*-1 + 2*xmax/15, -log10(pvalue) - ymax/30, paste0("FDR=", fdr, "%"), cex=1.15)    ## SCLC (IZ)
   #text(xmax*-1 + 2*xmax/13, -log10(pvalue) - ymax/30, paste0("FDR=", fdr*100, "%"), cex=1.1)   ## SCLC (AA)
   #text(xmax*-1 + 2*xmax/9.5, -log10(pvalue) - ymax/30, paste0("BH=1.00E-16"), cex=1.1)
   #text(xmax*-1 + 2*xmax/30, -log10(pvalue) + ymax/30, paste0("***"), cex=2)
   
   de.up   <- subset(de.sig, LOG2_FC > log2(fc))
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=cols[1], cex=1.4)
   de.down <- subset(de.sig, LOG2_FC < -log2(fc))
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=cols[2], cex=1.4)
 
   #abline(v=log2(fc), lty=5, lwd=1.8, col=cols[1]) 
   #abline(v=-log2(fc), lty=5, lwd=1.8, col=cols[2])
   abline(h=c(-log10(pvalue)), lty=5, lwd=1.5)
   
   if (nrow(genes) != 0) {
      for (g in 1:nrow(genes)) {
         gene <- subset(de, external_gene_name == genes[g,]$GENE)
         gene <- cbind(gene, genes[g,])
   
         if (nrow(gene) > 0) {
            points(gene$LOG2_FC, gene$log10P, pch=1, col="black", cex=1.4)
    
            if (!is.na(gene$ADJ_1))
               if (is.na(gene$ADJ_2))
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=1.25)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.25)
            else
               if (gene$LOG2_FC > 0)
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.6), cex=1.25)
               else
                  text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.6), cex=1.25)
         } else
            print(genes[g])
      }
   }
 
   axis(side=1, at=seq(-10, 10, by=4), labels=c(-10, -6, -2, 2, 6, 10), cex.axis=1.2)
   axis(side=1, at=seq(-8, 8, by=4), labels=c(-8, -4, 0, 4, 8), cex.axis=1.2)
   mtext(file.main[2], cex=1.25, line=0.3)
   legend("topleft", legend=c("Upregulated (SCLC > Lung)", "Downregulated (SCLC < Lung)"), col=cols, pch=19, cex=1.25)
   dev.off()
}

## Volcano
xlab.text <- "SCLC to lung fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_iz_e_cd_p1e-15_lung")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("Differential expression between SCLC and lung", paste0("Early IZ, CD genes (n=1,159)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano0(de2.tpm.gene, 1E-15, genes, file.de, file.main, xlab.text, ymax=23, c(yellow, lightblue), 1)

## Volcano
xlab.text <- "SCLC to lung fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_p1e-15_lung_BRD9")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("Differential expression between SCLC and lung", paste0("Expressed genes (n=18,265)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano0(de3.tpm.gene, 1E-15, genes, file.de, file.main, xlab.text, ymax=23, c(yellow, lightblue), 8)

## Volcano
xlab.text <- "SCLC to lung fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_p1e-15_lung_n7_up_TONSL")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("Differential expression between SCLC and lung", paste0("Expressed genes (n=18,265)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano0(de3.tpm.gene, 1E-15, genes, file.de, file.main, xlab.text, ymax=23, c(yellow, lightblue), 8)

## Volcano
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_p1e-15_lung_n41_up")

ensembl <- readTable(paste0(plot.de, ".txt"), header=F, rownames=F, sep="\t")
ensembl2 <- intersect(ensembl, overlaps)

t <- de3.tpm.gene[ensembl2,]
t <- t[order(t$P), ]
subset(subset(t, P < 1E-15), LOG2_FC > 3)












## Volcano (S)
xlab.text <- "SCLC/Lung [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_median0_iz_e_cd_p1e-15_lung_UBE2I")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("Differential expression between SCLC and lung", paste0("SCLC early IZ, CD genes (n=1,178)"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano0(de2.tpm.gene, 1E-15, genes, file.de, file.main, xlab.text, ymax=21)

###
##
#tpm.gene.log2.0      <- tpm.gene.log2[rownames(de.tpm.gene.iz),]
#tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[rownames(de.tpm.gene.iz),]
load("/Users/tpyang/Work/uni-koeln/tyang2/LUSC/analysis/expression/kallisto/normal-tpm-de/data/lusc_kallisto_0.43.1_tpm.gene.RData")
tpm.gene.lung.log2   <- log2(tpm.gene + 0.01)
# > dim(tpm.gene.lung.log2)
# [1] 34908    41

overlaps <- intersect(intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.e.cd), rownames(tpm.gene.log2)), rownames(tpm.gene.lung.log2))
length(overlaps)
# [1] 1377
tpm.gene.log2.0      <- tpm.gene.log2[overlaps,]
tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[overlaps,]

colnames <- c("P", "FDR", "N", "T", "LOG2_FC")
de2 <- toTable(0, length(colnames), nrow(tpm.gene.log2.0), colnames)
rownames(de2) <- rownames(tpm.gene.log2.0)

## SRC
de2$P <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) testU(as.numeric(tpm.gene.log2.0[x,]), as.numeric(tpm.gene.lung.log2.0[x,])))

## Log2 fold change
de2$N <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) median(as.numeric(tpm.gene.lung.log2.0[x,])))
de2$T <- mapply(x = 1:nrow(tpm.gene.log2.0), function(x) median(as.numeric(tpm.gene.log2.0[x,])))
de2$LOG2_FC <- de2$T - de2$N

#de2.iz.sp <- de2[intersect(rownames(de2), rownames(tpm.gene.log2.m.rfd.ctr.iz.cd)),]
#de2.iz.sp$FDR <- p.adjust(de2.iz.sp$P, method="BH", n=length(de2.iz.sp$P))
#de2.iz.sp$Q <- qvalue(de2.iz.sp$P)$qvalue
#de2.iz.sp <- de2.iz.sp[order(de2.iz.sp$P),]

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de2$FDR <- p.adjust(de2$P, method="BH", n=length(de2$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de2 <- de2[order(de2$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de2.tpm.gene <- cbind(annot[rownames(de2),], de2)   ## BE EXTRA CAREFUL!!

save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-median0-lung_U_BH_n70+41_iz_e_cd.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-median0-lung_U_BH_n70+41_iz_e_cd.txt"), colnames=T, rownames=F, sep="\t")
nrow(de2.tpm.gene)
# [1] 1377

###
##
plotNormal <- function(gene, id, tpm.gene.log2, tpm.gene.lung.log2) {
   file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_", gene)
   plotBox2(wd.de.plots, file.name, as.numeric(tpm.gene.lung.log2[id,]), as.numeric(tpm.gene.log2[id,]), main=gene, names=c("Lung", "SCLC"))
}

genes <- c("PIF1", "KIF18B", "MARS", "AL049840.1", "GTPBP3", "EIF3B")
genes <- c("BRCA2")
genes <- c("IRF2", "PIAS3", "UBE2I", "AAAS")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
}

# -----------------------------------------------------------------------------
# D.E. LUNG vs. SCLC
# Last Modified: 25/10/20; 01/09/20; 08/01/20
# -----------------------------------------------------------------------------
   plotNormal(genes[g], id, tpm.gene.log2, tpm.gene.lung.log2)
load("/Users/tpyang/Work/uni-koeln/tyang2/LUSC/analysis/expression/kallisto/normal-tpm-de/data/lusc_kallisto_0.43.1_tpm.gene.RData")
tpm.gene.lung.log2   <- log2(tpm.gene + 0.01)
# > dim(tpm.gene.lung.log2)
# [1] 34908    41

de.tpm.gene.down <- subset(subset(de.tpm.gene, P < 0.01), LOG2_FC < 0)
overlaps <- intersect(rownames(de.tpm.gene.down), rownames(tpm.gene.lung.log2))
length(overlaps)
# [1] 267
tpm.gene.log2.0      <- tpm.gene.log2[overlaps,]
tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[overlaps,]

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

save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-median0-lung_U_BH_n70+41_down_n267.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-median0-lung_U_BH_n70+41_iz_e_cd.txt"), colnames=T, rownames=F, sep="\t")
nrow(de2.tpm.gene)
# [1] 267









overlaps <- intersect(rownames(de2.tpm.gene), rownames(tpm.gene))
# length(overlaps)
# [1] 2199
de2.tpm.gene <- de2.tpm.gene[overlaps,]
save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_expressed-in-lung.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_expressed-in-lung.txt"), colnames=T, rownames=F, sep="\t")

##
genes <- readTable(paste0(plot.de, "_lung.tab"), header=T, rownames=F, sep="\t")
file.main <- c("TTR + CTR", paste0("(n=17,311)"))
file.de <- paste0(plot.de, "_lung.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

###
## Volcano (All genes)
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_iz_sp_bh1e-16")

genes <- readTable(paste0(plot.de, "_lung2.tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC-specificly initiated (IZ-S), expressed genes", paste0("(n=1,245)"))
file.de <- paste0(plot.de, "_lung2.pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)

genes <- readTable(paste0(plot.de, "_lung_FOXH1.tab"), header=T, rownames=F, sep="\t")
file.main <- c("Initiation (IZ)", paste0("(n=1,299)"))
file.de <- paste0(plot.de, "_lung_FOXH1.pdf")
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
## B2M
median(as.numeric(tpm.gene.lung.log2["ENSG00000166710",]))
median(as.numeric(tpm.gene.log2["ENSG00000166710",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000166710",]), as.numeric(tpm.gene.log2["ENSG00000166710",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_B2M_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000166710",], tpm.gene.log2["ENSG00000166710",], main="B2M", names=c("Lung", "SCLC"))

## MARS
median(as.numeric(tpm.gene.lung.log2["ENSG00000166986",]))
median(as.numeric(tpm.gene.log2["ENSG00000166986",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000166986",]), as.numeric(tpm.gene.log2["ENSG00000166986",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_MARS_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000166986",], tpm.gene.log2["ENSG00000166986",], main="MARS", names=c("Lung", "SCLC"))

## KIF18B
median(as.numeric(tpm.gene.lung.log2["ENSG00000186185",]))
median(as.numeric(tpm.gene.log2["ENSG00000186185",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000186185",]), as.numeric(tpm.gene.log2["ENSG00000186185",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_KIF18B_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000186185",], tpm.gene.log2["ENSG00000186185",], main="KIF18B", names=c("Lung", "SCLC"))

## BRCA2
median(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]))
median(as.numeric(tpm.gene.log2["ENSG00000139618",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]), as.numeric(tpm.gene.log2["ENSG00000139618",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRCA2_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139618",], tpm.gene.log2["ENSG00000139618",], main="BRCA2", names=c("Lung", "SCLC"))

## TP53
median(as.numeric(tpm.gene.lung.log2["ENSG00000141510",]))
median(as.numeric(tpm.gene.log2["ENSG00000141510",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000141510",]), as.numeric(tpm.gene.log2["ENSG00000141510",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TP53_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000141510",], tpm.gene.log2["ENSG00000141510",], main="TP53", names=c("Lung", "SCLC"))

## RB1
median(as.numeric(tpm.gene.lung.log2["ENSG00000139687",]))
median(as.numeric(tpm.gene.log2["ENSG00000139687",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139687",]), as.numeric(tpm.gene.log2["ENSG00000139687",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_RB1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139687",], tpm.gene.log2["ENSG00000139687",], main="RB1", names=c("Lung", "SCLC"))



## TONSL
median(as.numeric(tpm.gene.lung.log2["ENSG00000160949",]))
median(as.numeric(tpm.gene.log2["ENSG00000160949",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000160949",]), as.numeric(tpm.gene.log2["ENSG00000160949",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TONSL_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000160949",], tpm.gene.log2["ENSG00000160949",], main="TONSL", names=c("Lung", "SCLC"))

## RECQL4
median(as.numeric(tpm.gene.lung.log2["ENSG00000160957",]))
median(as.numeric(tpm.gene.log2["ENSG00000160957",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000160957",]), as.numeric(tpm.gene.log2["ENSG00000160957",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_RECQL4_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000160957",], tpm.gene.log2["ENSG00000160957",], main="RECQL4", names=c("Lung", "SCLC"))

## CPSF1
median(as.numeric(tpm.gene.lung.log2["ENSG00000071894",]))
median(as.numeric(tpm.gene.log2["ENSG00000071894",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000071894",]), as.numeric(tpm.gene.log2["ENSG00000071894",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_CPSF1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000071894",], tpm.gene.log2["ENSG00000071894",], main="CPSF1", names=c("Lung", "SCLC"))

## FOXH1
median(as.numeric(tpm.gene.lung.log2["ENSG00000160973",]))
median(as.numeric(tpm.gene.log2["ENSG00000160973",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000160973",]), as.numeric(tpm.gene.log2["ENSG00000160973",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_FOXH1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000160973",], tpm.gene.log2["ENSG00000160973",], main="FOXH1", names=c("Lung", "SCLC"))

## PIF1
median(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]))
median(as.numeric(tpm.gene.log2["ENSG00000140451",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000140451",]), as.numeric(tpm.gene.log2["ENSG00000140451",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_PIF1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000140451",], tpm.gene.log2["ENSG00000140451",], main="PIF1", names=c("Lung", "SCLC"))

## TOR1AIP1
median(as.numeric(tpm.gene.lung.log2["ENSG00000143337",]))
median(as.numeric(tpm.gene.log2["ENSG00000143337",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000143337",]), as.numeric(tpm.gene.log2["ENSG00000143337",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_TOR1AIP1_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000143337",], tpm.gene.log2["ENSG00000143337",], main="TOR1AIP1", names=c("Lung", "SCLC"))

## BRCA2
median(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]))
median(as.numeric(tpm.gene.log2["ENSG00000139618",]))
testU(as.numeric(tpm.gene.lung.log2["ENSG00000139618",]), as.numeric(tpm.gene.log2["ENSG00000139618",]))
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_BRCA2_")
plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2["ENSG00000139618",], tpm.gene.log2["ENSG00000139618",], main="BRCA2", names=c("Lung", "SCLC"))

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
