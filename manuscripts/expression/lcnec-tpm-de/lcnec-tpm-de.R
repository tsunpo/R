# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/lcnec-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "LCNEC"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "MAX_INSERT_SIZE", "AVG_FRAGMENT_LENGTH", "RB1_MUT")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# D.E. using non-parametric test (20 RB1 vs 34 WT; n=69-15NA)
# Last Modified: 22/05/18
# -----------------------------------------------------------------------------
## Parameters for this test
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RB1_MUT (1) vs RB1_WT(0) as factor
argv      <- data.frame(predictor="RB1_MUT", predictor.wt=0, test="Wilcox", test.fdr="Q", fdr=0.1, effect=0, stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm_gene_r5_p47_rb1_wilcox_q_n54")
file.main <- paste0("RB1 MUT (n=20) vs WT (n=34) in ", BASE)

de.tpm.gene <- pipeDE(tpm.gene.log2, samples, argv, ensGene)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
genes <- rownames(subset(de.tpm.gene, FDR <= 0.1))
pca.de <- getPCA(t(tpm.gene[genes,]))   ## BUG FIX 13/02/17: Perform PCA using normalised data

###
## RB1 status on D.E genes
pheno.expr <- pheno[samples.expr,]

#trait <- pheno.expr[,"RB1NEW"]
trait <- as.numeric(samples[,"RB1_MUT"])
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

file.main <- "LCNEC RB1 status on 510 D.E. (FDR=10%) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "test2_RB1_510DE_TEST", file.main, NA, NA, c("red", "dodgerblue"))

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 01/10/17
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, fdr, genes, file.de, file.main, isIg) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="Median fold change (log2 RB1/WT)", ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(fdrToP(0.05, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(fdrToP(0.05, de)) + ymax/42, "FDR=5%", col="darkgray")
   
   if (isIg) {   ## ADD 02/11/17: Immunoglobulin (Ig) variable chain and T-cell receptor (TcR) genes
      de.ig <- de[grep("IG", de$gene_biotype),]
      points(de.ig$Effect, de.ig$log10P, pch=16, col="gold")
      de.tr <- de[grep("TR", de$gene_biotype),]
      points(de.tr$Effect, de.tr$log10P, pch=16, col="forestgreen")
   }
   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")

         if (genes[g] == "CD274")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.04, 1.3), cex=0.75)   
         else if (genes[g] == "NOTCH4" || genes[g] == "BRD4" || genes[g] == "TOP1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.3), cex=0.75)
         else if (genes[g] == "REST" || genes[g] == "CDK9")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1.05, 1.3), cex=0.75)
         else if (genes[g] == "NOTCH3" || genes[g] == "IL6R" || genes[g] == "IL6ST")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
         else if (genes[g] == "TUBA1A")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.08, cex=0.75)
         else if (genes[g] == "CHEK2" || genes[g] == "NOTCH2" || genes[g] == "CYP1B1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.1, cex=0.75)
         else if (genes[g] == "TNFRSF1A")
             text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.08, cex=0.75)
         else if (genes[g] == "TOP2A" || genes[g] == "STMN1" || genes[g] == "BRCA2" || genes[g] == "BRIP1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.1, cex=0.75)
         else if (genes[g] == "ASCL1" || genes[g] == "MYCN")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.13, cex=0.75)
         else if (genes[g] == "E2F2" || genes[g] == "RAD51AP1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.15, cex=0.75)  
         else if (genes[g] == "BLM" || genes[g] == "ATM")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.17, cex=0.75)
         else if (genes[g] == "RNASEH2A" || genes[g] == "BRCA1" || genes[g] == "HMGB2" || genes[g] == "XRCC1" || genes[g] == "XRCC4" || genes[g] == "NOTCH1" || genes[g] == "XRCC4" || genes[g] == "HES1" || genes[g] == "MYD88")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
         else if (genes[g] == "PARP1" || genes[g] == "PCNA" || genes[g] == "RBL1" || genes[g] == "RBL2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.55), cex=0.75)
         else
            if (gene$Effect > 0)
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
   
   if (isIg)
      legend("topleft", legend=c("Up-regulation", "Down-regulation", "Ig VJC genes (no D)", "TCR VC genes (no DJ)"), col=c("red", "dodgerblue", "gold", "forestgreen"), pch=19)
   else
      legend("topleft", legend=c("Up-regulation", "Down-regulation"), col=c("red", "dodgerblue"), pch=19)
   
   abline(h=c(-log10(p)), lty=5)
   text(xmax*-1 + 2*xmax/24, -log10(p) + ymax/42, "FDR=10%")
   dev.off()
}

## Volcano plot
de.lcnec <- de.lcnec.tpm.gene
de.lcnec <- de.lcnec[setdiff(rownames(de.lcnec), "ENSG00000269846"),]

plot.main <- "RB1mut-associated differential expression in LCNEC"
plot.de <- paste0(wd.lcnec.de, "volcanoplot_lcnec_tpm-gene_rb1_fdr10")

##
genes <- c("MYCL", "MYCN", "NEUROD1", "ASCL1", "SYP", "NCAM1", "ELAVL4", "CD44", "BMP4", "ATXN1")
file.main <- c(plot.main, "", "Non-/Neuroendocrine marker", "")
file.de   <- paste0(plot.de, "_Neuroendocrine.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, F)

## Cycle
#genes <- c("CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN2C", "PARP1", "PARP2", "BARD1", "BRIP1", "BRCA1", "BRCA2", "CDC7", "CCNE2", "CCND1", "E2F2", "E2F7", "E2F1", "ATM", "ATR", "CHEK1", "CHEK2", "CDK2", "POLD4", "TOPBP1", "CLSPN", "TOP2A", "TOP1", "RB1", "RBL1", "RBL2", "TP73", "TP53", "TP53BP1", "PCNA")
#genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
genes <- c("CHAF1B", "ASF1B", "PPM1D", "TMSB15A", "TMSB15B", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN2C", "PARP1", "PARP2", "BARD1", "BRIP1", "BRCA1", "BRCA2", "CDC7", "CCNE2", "CCND1", "E2F2", "E2F7", "E2F1", "ATM", "ATR", "CHEK1", "CHEK2", "CDK2", "POLD4", "TOPBP1", "CLSPN", "TOP2A", "TOP1", "RB1", "RBL1", "RBL2", "TP73", "TP53", "TP53BP1", "PCNA", "STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2")
file.main <- c(plot.main, "", "Cell cycle pathway", "")
file.de <- paste0(plot.de, "_Cycle.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, F)

## 53BP1
#genes <- c("SIRT1", "MSH2", "MSH6", "DNMT1", "E2F2", "E2F7", "E2F1", "SMARCA4", "MAX", "SMARCB1", "HELLS", "ARID1A", "NHEJ1", "LIG4", "POLQ", "DNA2", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
#genes <- c("CHAF1B", "ASF1B", "MSH2", "MSH6", "DNMT1", "E2F2", "E2F7", "E2F1", "SMARCA4", "SMARCB1", "ARID1A", "NHEJ1", "LIG1", "CHEK2", "LIG4", "POLQ", "DNA2", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
genes <- c("RAD50", "NBN", "MRE11A", "MSH2", "MSH6", "TOPBP1", "SMARCAD1", "SMARCA4", "SMARCB1",  "NHEJ1", "LIG1", "LIG4", "E2F4", "POLQ", "DNA2", "SETD2", "RNF8", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "RIF1", "MUM1", "MDC1", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
file.main <- c(plot.main, "", "DNA repair pathway", "")
file.de <- paste0(plot.de, "_53BP1.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, F)

## NOTCH
#genes <- c("CDK9", "KRAS", "MYC", "MYCL", "MYCN", "NEUROD1", "ASCL1", "PSIP1", "HMGB2", "BRD4", "EZH2", "UCHL1", "HES1", "REST", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "EGFR", "ERBB2", "ERBB3", "DPYSL5", "CRMP1", "DPYSL3")
genes <- c("CDK9", "KRAS", "MYC", "MAX", "MYCL", "MYCN", "NEUROD1", "ASCL1", "DLL3", "PSIP1", "HMGB2", "BRD4", "EZH2", "UCHL1", "HES1", "REST", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "EGFR", "ERBB2", "ERBB3", "DPYSL5", "CRMP1", "DPYSL3")
file.main <- c(plot.main, "", "Neuroendocrine marker", "")
file.de <- paste0(plot.de, "_Neuroendocrine.pdf")
plotVolcano(de.lcnec, 0.05, genes, file.de, file.main, F)

# -----------------------------------------------------------------------------
# Volcano plots Lite
# Last Modified: 25/10/17
# -----------------------------------------------------------------------------
plotVolcanoLite <- function(de, fdr, genes, file.de, file.main) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="Median fold change (log2 RB1/WT)", ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(fdrToP(0.05, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(fdrToP(0.05, de)) + ymax/42, "FDR=5%", col="darkgray")
 
   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
 
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")
    
         if (genes[g] == "XRCC2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.10, cex=0.75)
         else if (genes[g] == "RAD52")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.11, cex=0.75)
         else if (genes[g] == "CDK6" || genes[g] == "E2F4" || genes[g] == "MDM2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.13, cex=0.75)
         else if (genes[g] == "RNASEH2A")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.05, cex=0.75)
         else if (genes[g] == "MAD2L2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         #else if (genes[g] == "RBL2")
         #   text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         else if (genes[g] == "PARP1" || genes[g] == "XRCC4" || genes[g] == "RNF138" || genes[g] == "BRCA2" || genes[g] == "ASF1A" || genes[g] == "HMGB1" || genes[g] == "CHAF1A")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.1, cex=0.75)
         else if (genes[g] == "PSIP1" || genes[g] == "RBBP4")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.12, cex=0.75)
         else if (genes[g] == "E2F8")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.14, cex=0.75)
         else if (genes[g] == "LIG1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.16, cex=0.75)
         else if (genes[g] == "HMGB2" )
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.3), cex=0.75)
         else if (genes[g] == "CREBBP")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.05, 1.3), cex=0.75)
         else
            if (gene$Effect > 0)
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   abline(h=c(-log10(p)), lty=5)
   text(xmax*-1 + 2*xmax/24, -log10(p) + ymax/42, "FDR=10%")
   legend("topleft", legend=c("Up-regulation", "Down-regulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

plot.main <- "Differential effect of RB1 loss in LCNEC"
plot.de <- paste0(wd.lcnec.de, "volcanoplot_lcnec_tpm-gene_rb1_fdr10")

###
## Figure 1 (Cell cycle)
genes <- c("TP53", "MDM2", "RB1", "E2F1", "E2F7", "E2F8", "E2F4", "E2F5", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN1B", "CDKN2C", "BARD1", "CDC7", "CDK2", "CCNE2", "CCND1", "POLD4", "TOPBP1", "CHEK2", "ATR", "ATM")
genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
file.main <- c(plot.main, "", "Cell cycle control", "")
file.de <- paste0(plot.de, "_Figure_1_Cycle2.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)

## Figure 1 (DSB repair pathway choice)
genes <- c("XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "FEN1", "LIG1", "LIG3", "LIG4", "TP53BP1", "PARP1", "RB1",  "RBL1", "RBL2", "BRCA1", "BRCA2", "BARD1", "TOPBP1", "PSIP1", "RBBP8", "RNF8", "RNF168", "RNF138", "RAD51AP1", "RAD52")   ##"CHAF1A", "ASF1A", "SMARCAD1", "RAD51", "EXO1", "DNA2", "MDC1", "RIF1", "MAD2L2", "CDH1", "RNF8")
file.main <- c(plot.main, "", "DSB repair pathway choice", "")
file.de <- paste0(plot.de, "_Figure_1_53BP1.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)

## Figure 1 (Histone chaperone and DSB repair pathway choice)
genes <- c("TP53BP1", "RB1",  "RBL1", "RBL2", "RNF8", "RNF168", "RNF138", "ASF1B", "CHAF1B", "CHAF1A", "ASF1A", "TOPBP1", "PSIP1", "RBBP8", "RBBP4", "HMGB2", "HMGB1", "SMARCAD1", "SMARCB1", "EP300", "CREBBP", "DNMT1")   ## "RAD51", "EXO1", "DNA2", "MDC1", "RIF1", "MAD2L2", "CDH1", "RNF8")
file.main <- c(plot.main, "", "Histone chaperone and chromatin remodeler", "")
file.de <- paste0(plot.de, "_Figure_1_53BP1_ASF1A.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)

## Figure 1
genes <- c("XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "FEN1", "LIG1", "LIG4", "TP53", "MDM2", "TP53BP1", "PARP1", "RB1", "RBL1", "RBL2", "E2F1", "E2F7", "E2F8", "E2F4", "E2F5", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN1A", "CDKN1B", "BARD1", "CDC7", "CDK2", "CCNE2", "CCND1", "POLD4", "TOPBP1", "CHEK2", "BRCA1", "BRCA2", "PSIP1", "RBBP8", "RNF168", "BARD1")
file.main <- c(plot.main, "", "DSB repair pathway choice & Cell cycle control", "")
file.de <- paste0(plot.de, "_Figure_1.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)
