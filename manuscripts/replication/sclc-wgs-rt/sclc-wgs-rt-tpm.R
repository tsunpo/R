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
#overlaps <- intersect(rownames(tpm.gene), rownames(wgs.gene))   ## After sclc-wgs-de-cm2.R
#tpm.gene <- tpm.gene[overlaps,]
#tpm.gene.log2 <- log2(tpm.gene + 0.01)
tpm.gene.log2 <- getLog2andMedian(tpm.gene, 0.01)

# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
tpm.gene.log2 <- cbind(annot[rownames(tpm.gene.log2),], tpm.gene.log2[, "MEDIAN"], ensGene.bed[rownames(tpm.gene.log2),]) 
colnames(tpm.gene.log2)[8] <- "MEDIAN"

tpm.gene.log2$TSS_RFD <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$RFD
#tpm.gene.log2$TTS_RFD <- nrds.RT.NRFD[tpm.gene.log2$TTS,]$RFD
tpm.gene.log2$TSS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$SPLINE
#tpm.gene.log2$TTS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2$TTS,]$SPLINE

tpm.gene.log2 <- tpm.gene.log2[!is.na(tpm.gene.log2$TSS_RFD),]

##
tpm.gene.log2$TRC <- tpm.gene.log2$strand * tpm.gene.log2$TSS_RFD
tpm.gene.log2$TRC[which(tpm.gene.log2$TRC > 0)] <- 1
tpm.gene.log2$TRC[which(tpm.gene.log2$TRC < 0)] <- -1
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
de$Q   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_n70.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.cd))
de.tpm.gene.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.cd$Q <- qvalue(de.tpm.gene.cd$P)$qvalue
writeTable(de.tpm.gene.cd, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_cd_n70.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ho))
de.tpm.gene.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.ho$Q <- qvalue(de.tpm.gene.ho$P)$qvalue
writeTable(de.tpm.gene.ho, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47_cor_src_q_ho_n70.txt"), colnames=T, rownames=F, sep="\t")

testU(tpm.gene.log2.cd$MEDIAN, tpm.gene.log2.ho$MEDIAN)

# -----------------------------------------------------------------------------
# 
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
plotCYS <- function(gene, cn, snr, pch, col, pos) {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   xlab.text <- "Cell cycle status"
   ylab.text <- "log2(TPM + 0.01)"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.rt.plots, paste0("TPM-vs-CYCLE_", genes[g]))
   
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
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}




file.name <- file.path(wd.rt.plots, "TPM-vs-CYCLE_GTF3C2")
plotCN(as.numeric(tpm.gene.log2["ENSG00000115207",]), samples$COR, 1, file.name, "black", "bottomright")











# -----------------------------------------------------------------------------
# D.E. using non-parametric test
# Last Modified: 04/04/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: M2 vs M1 (as factor)
argv      <- data.frame(predictor="M2", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47-wgs_m2_wilcox_q_n70")
file.main <- paste0("M2 (n=35) vs M1 (n=35) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_sclc_tpm-gene-r5p47-wgs_m2_wilcox_q_n70_weight")
writeRNKformat(de.tpm.gene, wd.de.gsea, file.name)

## Tirosh et al 2016 
load(file.path(wd.src.ref, "cycle.RData"))                           ## See guide-to-the/cycle.R
core.G1S <- intersect(core.G1S, rownames(tpm.gene.log2))             ## 41/43/43 are expressed in the dataset
core.G2M <- intersect(core.G2M, rownames(tpm.gene.log2))             ## 54/54/55
core.Stemness <- intersect(core.Stemness, rownames(tpm.gene.log2))   ## 54/58/63

writeGRPformat(core.G1S, wd.de.gsea, "core.G1-S")
writeGRPformat(core.G2M, wd.de.gsea, "core.G2-M")
writeGRPformat(core.Stemness, wd.de.gsea, "core.Stemness")

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 262/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 792/876

writeGRPformat(periodic.G1S, wd.de.gsea, "G1-S")
writeGRPformat(periodic.G2M, wd.de.gsea, "G2-M")

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

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 15/08/18
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, pvalue, genes, file.de, file.main, xlab.text) {
   #pvalue <- fdrToP(fdr, de)
   fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #xmax <- 1
   ymax <- max(de$log10P)
   
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab=xlab.text, ylab="-log10(p-value)", col="darkgray", main=file.main[1])

   abline(h=c(-log10(pvalue)), lty=5)
   text(xmax*-1 + 2*xmax/25, -log10(pvalue) + ymax/45, paste0("FDR=", fdr, "%"), cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="red")
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:nrow(genes)) {
      gene <- subset(de, external_gene_name == genes[g,]$GENE)
      gene <- cbind(gene, genes[g,])
      
      if (nrow(gene) > 0) {
         points(gene$LOG2_FC, gene$log10P, pch=1, col="black")
         
         if (!is.na(gene$ADJ_1))
            if (is.na(gene$ADJ_2))
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=0.75)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=0.75)
         else
            if (gene$LOG2_FC > 0)
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
   
   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topleft", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

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
