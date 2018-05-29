# =============================================================================
# Manuscript   : Loss of RB1 escapes cell cycle arrest
# Chapter      : RB-loss differential effect in LCNEC
# Name         : manuscripts/expression/lcnec-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 02/11/17
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/re/home/tyang2/dev/R"                ## tyang2@gauss
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
setwd(wd.de)

samples <- readTable(file.path(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "MAX_INSERT_SIZE", "RB1_MUT")

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
load(file.path(wd, "ALL", "analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/data/all_kallisto_0.43.1_tpm.gene_r5_p47.RData"))
tpm.gene      <- tpm.gene[,rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# D.E. using non-parametric test (20 RB1 vs 34 WT; n=69-15NA)
# Last Modified: 22/05/18
# -----------------------------------------------------------------------------
## Parameters for this test
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RB1_MUT (1) vs RB1_WT(0) as factor
argv      <- data.frame(predictor="RB1_MUT", predictor.wt=0, test="Wilcox", test.fdr="Q", fdr=0.05, effect=0, stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm_gene_rb1_wilcox_q_n54")
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
# PCA (with SCLC and HeLa)
# Last Modified: 14/12/17
# -----------------------------------------------------------------------------
signature.rb1 <- rownames(subset(de.lcnec.tpm.gene, FDR <= 0.1))
signature.rb1.overlaps <- intersect(signature.rb1, rownames(hela.tpm.gene.log2))
## > length(signature.rb1.overlaps)
## [1] 464

test <- cbind(all.tpm.gene.log2[signature.rb1.overlaps, c(samples.lcnec, samples.sclc)], hela.tpm.gene.log2[signature.rb1.overlaps, hela.samples])
pca.de <- getPCA(t(test))

trait <- pheno.expr[,"RB1NEW"]
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

trait <- c(trait, rep("SCLC", length(samples.sclc)))
trait <- c(trait, rep("HeLa", length(samples.hela)))

file.main <- "SCLC and HeLa (n=14) on LCNEC 464 D.E. (RB1 v WT; FDR=10%) genes"
plotPCA(1, 2, pca.de, trait, wd.de, "pca_RB1_464DE_PC1-PC2_SCLC+HeLa", file.main, NA, -15, c("purple", "lightgray", "red", "orange", "dodgerblue"))

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

## PPM1D   ## MAPK14 (p38), HSPB1 (Hsp27), RBM38, CBX5 (HP1)
genes <- c("TP53", "RB1", "RBL1", "RBL2", "CDKN2A", "CCND1", "PPM1D", "MAPK14", "HSPB1", "RBM38", "CBX5", "EZH2", "SLFN11", "DNMT1", "DNMT3A", "DNMT3B", "BRCA1", "TP53BP1", "ATM")
file.main <- c(plot.main, "", "DNA methylation-mediated silencing", "")
file.de   <- paste0(plot.de, "_Cycle_PPM1D.pdf")
plotVolcano(de.lcnec, 0.05, genes, file.de, file.main, F)

## 53BP1
#genes <- c("SIRT1", "MSH2", "MSH6", "DNMT1", "E2F2", "E2F7", "E2F1", "SMARCA4", "MAX", "SMARCB1", "HELLS", "ARID1A", "NHEJ1", "LIG4", "POLQ", "DNA2", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
#genes <- c("CHAF1B", "ASF1B", "MSH2", "MSH6", "DNMT1", "E2F2", "E2F7", "E2F1", "SMARCA4", "SMARCB1", "ARID1A", "NHEJ1", "LIG1", "CHEK2", "LIG4", "POLQ", "DNA2", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
genes <- c("RAD50", "NBN", "MRE11A", "MSH2", "MSH6", "TOPBP1", "SMARCAD1", "SMARCA4", "SMARCB1",  "NHEJ1", "LIG1", "LIG4", "E2F4", "POLQ", "DNA2", "SETD2", "RNF8", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "RIF1", "MUM1", "MDC1", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
file.main <- c(plot.main, "", "DNA repair pathway", "")
file.de <- paste0(plot.de, "_53BP1.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, F)

## Immune response
genes <- c("CD274", "PDCD1LG2", "MYD88", "NFKB1", "NFKB2", "TLR2", "BCL3", "IL6", "IL6ST", "IL6R")
file.main <- c(plot.main, "", "Immune response", "")
file.de <- paste0(plot.de, "_Immune.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, T)

## MHC I
genes <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F")
file.main <- c(plot.main, "", "MHC class I", "")
file.de <- paste0(plot.de, "_MHCI.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, T)

## MHC II
genes <- c("CIITA", "HLA-DRA", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DOA", "HLA-DOB", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1")
file.main <- c(plot.main, "", "MHC class II", "")
file.de <- paste0(plot.de, "_MHCII.pdf")
plotVolcano(de.lcnec, 0.1, genes, file.de, file.main, T)

###
## Immunoglobulin class switching
genes <- c("IL4R", "IFNGR2", "IFNGR1", "TGFBR2", "TGFB1", "TGFB3", "TGFB2", "TGFBR1", "TGFBR3")
file.main <- c(plot.main, "", "Immunoglobulin class switching", "")
file.de <- paste0(plot.de, "_Ig.pdf")
plotVolcano(de.lcnec, 0.05, genes, file.de, file.main, T)

## Innate immunity
genes <- c("IL6R", "GSK3A", "GSK3B")
file.main <- c(plot.main, "", "Innate immunity", "")
file.de <- paste0(plot.de, "_Innate.pdf")
plotVolcano(de.lcnec, 0.05, genes, file.de, file.main, T)
## XRCC   ##"UBE2N", "RNF8" 
genes <- c("SMARCA4", "SMARCB1", "NHEJ1", "LIG4", "POLQ", "SETD2", "RNF138", "RNF168", "PSIP1", "HMGB1", "HMGB2", "PARP1", "PARP2", "PARP3", "BARD1", "BRIP1", "BRCA1", "BRCA2", "RB1", "RBL1", "RBL2", "TP53", "TP53BP1", "XRCC1", "EZH2", "XRCC2", "XRCC3", "XRCC4", "RAD51AP1", "RBBP8", "FEN1", "EXO1", "BLM", "WRN")
file.main <- c(plot.main, "", "DNA repair pathway", "")
file.de   <- paste0(plot.de, "_XRCC.pdf")
plotVolcano(de.lcnec, 0.05, genes, file.de, file.main, F)

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


## PPM1D   ## MAPK14 (p38), HSPB1 (Hsp27), CBX1/CBX5 (HP1), MAP2K1 (MEK1), MAP2K2 (MEK2), MAPK3 (ERK1), MAPK1 (ERK2)
# genes <- c("TP53", "RB1", "RBL1", "RBL2", "CDKN2A", "CCND1", "CBX1", "PPM1D", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MAPK14", "HSPB1", "CBX5", "EZH2", "SLFN11", "DNMT1", "DNMT3A", "DNMT3B", "BRCA1", "TP53BP1", "ATM")
genes <- c("TP53", "PPM1D", "DNMT1", "MAPK14", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "ATM", "CHEK2")
file.main <- c(plot.main, "", "DNA methylation-induced silencing", "")
file.de   <- paste0(plot.de, "_Figure_1_Cycle_PPM1D.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)









## Volcano plot
file.de <- paste0(wd.de, "LCNEC_DE_kallisto_tpm-gene-pcg-exp_fdr5.pdf")
file.main <- c("Differential Gene Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")
genes <- c("MCM2", "RNASEH2A", "CDKN2A", "CDKN2C", "PARP1", "PARP2", "BARD1", "BRIP1", "BRCA1", "BRCA2", "CDC7", "CCNE2", "CCND1", "E2F7", "ATM", "ATR", "CHEK1", "CHEK2", "CDK2", "POLD4", "POLH", "TOPBP1", "CLSPN", "TOP2A", "TOP1", "RB1", "TP53", "TP73", "TP53BP1", "XRCC1", "XRCC4", "PCNA", "RAD51AP1", "RBBP8", "FEN1", "BLM")
genes <- c("STMN1", "FBXO5", "TUBA1A", "PSIP1", "TMPO", "HMGB2", genes)

#colnames(de.tpm.gene.pcg) <- c("ensembl_gene_id", "external_gene_name", "SRC_RHO", "SRC_P", "SRC_FDR", "ANOVA_P", "ANOVA_FDR", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")
plotVolcano(de.tpm.gene.pcg.exp2, 0.05, genes, file.de, file.main)

## NOTCH
file.de <- paste0(wd.de, "LCNEC_DE_kallisto_tpm-gene-pcg-exp_fdr5_NOTCH.pdf")
genes <- c("CDK9", "MYC", "MYCL", "MYCN", "NEUROD1", "ASCL1", "PSIP1", "HMGB2", "BRD4", "EZH2", "UCHL1", "HES1", "REST", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
plotVolcano(de.tpm.gene.pcg.exp2, 0.05, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 19/09/17
# -----------------------------------------------------------------------------
pca.de <- getPCA(t(tpm.gene.pcg.exp2.log2[rownames(subset(de.tpm.gene.pcg.exp2, FDR <= 0.05)),]))   ## BUG FIX 13/02/17: Perform PCA using normalised data

###
## RB1 status on D.E genes
pheno.expr <- pheno[samples.expr,]

trait <- pheno.expr[,"RB1NEW"]
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

file.main <- "LCNEC RB1 Status on 129 D.E. Genes"
plotPCA(1, 2, pca.de, trait, wd.de, "RB1_129DE", file.main, NA, NA, c("purple", "red", "dodgerblue", "blue"))












# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 03/07/17
# -----------------------------------------------------------------------------
plotVolcanoLCNEC <- function(de, de.sig, genes, file.de, file.main) {
   de$log10P <- -log10(de$P)
   de.sig$log10P <- -log10(de.sig$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="Effect size (log2 FC)", ylab="Significance (-log10 P)", col="darkgray", main=file.main)
 
   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
   for (g in 1:length(genes)) {
      gene <- subset(de.up, external_gene_name == genes[g])
      if (nrow(gene) != 0) {
         #if (genes[g] == "MCM6" || genes[g] == "MCM3" || genes[g] == "BRIP1" || genes[g] == "CLSPN" || genes[g] == "RFC5" || genes[g] == "RFC4" || genes[g] == "RNASEH2A")
         #points(gene$Effect, gene$log10P, pch=1, col="dimgray")
         points(gene$Effect, gene$log10P, pch=1, col="black")
         
         if (genes[g] == "CDK2" || genes[g] == "RFC4")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.1, 1.1), cex=0.75)   ##adj=c(1, 1.5)
         else if (genes[g] == "")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, 1.5), cex=0.75)
         else if (genes[g] == "MCM6" || genes[g] == "RFC5")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
         else if (genes[g] == "RNASEH2A")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.05, cex=0.75)
         else
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
      }
   }
   
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="blue")
   for (g in 1:length(genes)) {
      gene <- subset(de.down, external_gene_name == genes[g])
      if (nrow(gene) != 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")
         text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      }
   }
   
   genes <- c("TOP1", "POLA1", "POLE2", "POLD4", "CCND1", "CDK4", "CDK6", "CHEK1", "CHEK2", "BRCA1", "HELLS", "ARID1A", "SMARCA4", "POLH", "RB1", "TP53", "RAD51AP1")
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      points(gene$Effect, gene$log10P, pch=1, col="black")
      if (genes[g] == "TOP1" || genes[g] == "POLD4" || genes[g] == "CCND1" || genes[g] == "CDK6" || genes[g] == "CHGA" || genes[g] == "CHEK2")
         text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      else if (genes[g] == "POLA1")
         text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.1, 1.1), cex=0.75)   ##adj=c(0, 1.5)
      else if (genes[g] == "CDK1")
         text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, 1.5), cex=0.75)
      else
         text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
   }
   
   abline(h=c(-log10(p)), lty=5)
   legend("topleft", legend=c("Upregulation", "Downregulation"), col=c("red", "blue"), pch=19)
   dev.off()
}
#plotVolcanoLCNEC(de.tpm.gene.pcg.exp2, de.tpm.gene.pcg.exp2.fdr5, genes, file.de, file.main)

## Volcano plot
file.de <- paste0(wd.de, "LCNEC_DE_kallisto_tpm-gene-pcg-exp-fdr5.pdf")
file.main <- c("Differential Gene Expression in LCNEC", "Mutant RB1 (n=20) vs WT (n=34)")

de.tpm.gene.pcg.exp2.fdr5 <- subset(de.tpm.gene.pcg.exp2, FDR <= 0.05)
genes <- c("MCM2", "MCM3", "MCM6", "CDKN2A", "BARD1", "CLSPN", "TOPBP1", "TOP2A", "CDK2", "CDC7", "CCNE2", "E2F7", "RFC5", "RFC4", "RNASEH2A", "CHEK1", "CHEK2", "BRCA1", "HELLS", "ARID1A", "SMARCA4", "POLH", "RB1", "TP53", "TP73", "RAD51AP1")   ##"TMPO"
plotVolcanoLCNEC(de.tpm.gene.pcg.exp2, de.tpm.gene.pcg.exp2.fdr5, genes, file.de, file.main)





# -----------------------------------------------------------------------------
# Arturo
# Last Modified: 03/07/17
# -----------------------------------------------------------------------------
arturo <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg-exp/Arturo.txt", header=F, rownames=F, sep="\t")
out <- de.tpm.gene.pcg.exp2[match(arturo, de.tpm.gene.pcg.exp2$external_gene_name),]
out <- out[order(out$P),]
writeTable(out, "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de-tpm-gene-pcg-exp/Arturo.tsv", colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 10/02/17
# -----------------------------------------------------------------------------
pca.de <- getPCA(t(tpm.gene.pcg.exp2.log2[rownames(significantAndVariable(de.tpm.gene.pcg.exp2,      dao$effect, dao$fdr)),]))   ## BUG FIX 13/02/17: Perform PCA using normalised data

###
## RB1 status on D.E genes
pheno.expr <- pheno[samples.expr,]

trait <- pheno.expr[,"RB1NEW"]
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

file.main <- "LCNEC RB1 Status on 54 D.E. Genes"
plotPCA(1, 2, pca.de, trait, wd.de, "RB1_54DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))






# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results
# -----------------------------------------------------------------------------
wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de.tpm.gene.pcg.exp2.res2/pathway/"
#wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/kallisto/lcnec-rnaseq-de-kallisto-ensembl/de.tpm.gene.exp/pathway/"
list <- readTable(paste0(wd.reactome, "untitled.txt"), header=T, rownames=T, sep="\t")
colnames(list) <- c("ensembl_gene_id",	"external_gene_name")

reactome <- read.csv(paste0(wd.reactome, "result.csv"))
colnames(reactome) <- gsub("X.", "", colnames(reactome))
reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
for (r in 1:nrow(reactome)) {
   ids <- as.vector(reactome$Submitted.entities.found[r])
   ids <- unlist(strsplit(ids, ";"))
 
   for (i in 1:length(ids))
      if (nrow(list[ids[i],]) != 0)
         ids[i] <- list[ids[i],]$external_gene_name
 
   reactome$Submitted.entities.found[r] <- paste(ids, collapse=";")
}
writeTable(reactome, paste0(wd.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")










## Find outliers
scores <- pcaScores(pca.de)

outliers <- subset(scores[,c("PC1", "PC2")], PC1 < -2.5)
subset(cbind(outliers, pheno.expr[rownames(outliers),]), RB1NEW == 0)

outliers <- subset(scores[,c("PC1", "PC2")], PC1 > -2.5)
subset(cbind(outliers, pheno.expr[rownames(outliers),]), RB1NEW == 1)

# -----------------------------------------------------------------------------
# Transcription arrest: 20/04/17
# ----------------------------------------------------------------------------- 
de.gene.down <- subset(de.tpm.gene.pcg.exp2, Effect <= -1)
de.gene.up <- subset(de.tpm.gene.pcg.exp2, Effect >= 1)

##
tpm.gene.log2.median <- median0(tpm.gene.pcg.exp2.log2)
tpm.gene.log2.median.down <- median0(tpm.gene.pcg.exp2.log2[rownames(de.gene.down),])
tpm.gene.log2.median.up   <- median0(tpm.gene.pcg.exp2.log2[rownames(de.gene.up),])

xr <- range(density(tpm.gene.log2.median)$x)
yr <- range(density(tpm.gene.log2.median)$y)

pdf(paste0(wd.de, "density_kallisto_tpm.gene.pcg.exp2.pdf"))
plot(density(tpm.gene.log2.median), xlim=xr, ylim=yr, main="")
lines(density(tpm.gene.log2.median.up), xlim=xr, ylim=yr, col="red")
lines(density(tpm.gene.log2.median.down), xlim=xr, ylim=yr, col="blue")
#legend(x=legend.x, y=ymax, legend=c("Overexpressed", "Underexpressed"), col=c("red", "blue"), pch=21)
dev.off()






# -----------------------------------------------------------------------------
# Transcription arrest: 04/04/17
# ----------------------------------------------------------------------------- 
hist(median0(tpm.gene.pcg.exp2.log2[,rownames(subset(pheno.expr.final, RB1NEW == 0))]))
hist(median0(expr.res2[,rownames(subset(pheno.expr.final, RB1NEW == 1))]))

##
pdf(paste0(wd.de, "hist_kallisto_effect.pdf"))
hist(de.tpm.gene.pcg.exp2$Effect, xlab="Effect size (log2 FC)", ylab="Frequency")
dev.off()

##
de.down <- subset(de.tpm.gene.pcg.exp2, Effect <= 0)
de.down.chr1 <- subset(subset(subset(ensGene.gene[rownames(de.down),], chromosome_name == "chr1"), start_position > 35000000), start_position < 47500000)
#de.down.chr1 <- subset(ensGene.gene[rownames(de.down),], chromosome_name == "chr1")
de.down.chr1 <- cbind(de.down.chr1[,-6], de.tpm.gene.pcg.exp2[rownames(de.down.chr1),4:8])
de.down.chr1 <- de.down.chr1[order(de.down.chr1$start_position),]

de.up <- subset(de.tpm.gene.pcg.exp2, Effect >= 0)
de.up.chr1 <- subset(subset(subset(ensGene.gene[rownames(de.up),], chromosome_name == 1), start_position > 35000000), start_position < 47500000)
de.up.chr1 <- cbind(de.up.chr1[,-6], de.tpm.gene.pcg.exp2[rownames(de.up.chr1),4:8])
de.up.chr1 <- de.up.chr1[order(de.up.chr1$start_position),]

##
pdf(paste0(wd.de, "effect_chr1_MYCL_pcg_exp2.pdf"))
plot(de.up.chr1[,c(4,11)], type="o", col="red", xlim=c(35000000,47500000), ylim=c(-1.25,2))
lines(de.down.chr1[,c(4,11)], type="o", col="blue")
dev.off()

# -----------------------------------------------------------------------------
# Transcription arrest: 27/04/17
# ----------------------------------------------------------------------------- 
tpm.gene.pcg.log2 <- log2(tpm.gene.pcg + 0.01) 
de.tpm.gene.pcg <- pipeDE(tpm.gene.pcg.log2, pheno.expr, dao, file.de, file.main, annot.gene)

start <- 34600000
end   <- 46800000
bands <- c(start, 40100000, 44100000, end)
start <- start - 1000000
end   <- end   + 1000000

de.down <- subset(de.tpm.gene.pcg, Effect <= 0)
de.down.chr1 <- subset(subset(subset(ensGene.gene[rownames(de.down),], chromosome_name == 1), start_position >= start), start_position <= end)
#de.down.chr1 <- subset(ensGene.gene[rownames(de.down),], chromosome_name == 1)
de.down.chr1 <- cbind(de.down.chr1[,-6], de.tpm.gene.pcg[rownames(de.down.chr1),4:8])
de.down.chr1 <- de.down.chr1[order(de.down.chr1$start_position),]

de.up <- subset(de.tpm.gene.pcg, Effect >= 0)
de.up.chr1 <- subset(subset(subset(ensGene.gene[rownames(de.up),], chromosome_name == 1), start_position >= start), start_position <= end)
#de.up.chr1 <- subset(ensGene.gene[rownames(de.up),], chromosome_name == 1)
de.up.chr1 <- cbind(de.up.chr1[,-6], de.tpm.gene.pcg[rownames(de.up.chr1),4:8])
de.up.chr1 <- de.up.chr1[order(de.up.chr1$start_position),]

###
##
start <- 34600000
end   <- 46800000

pdf(paste0(wd.de, "effect_tpm.gene.pcg_1p34.3-2-1.pdf"), height=4, width=10)
plot(de.up.chr1[,c(4,11)], type="o", col="red", xlim=c(start, end), ylim=c(-1.25,2), main="Effect Size", xlab="Chromosome 1", ylab="log2 FC")
#plot(de.up.chr1[,c(4,11)], type="o", col="red", ylim=c(-2.5,3))
lines(de.down.chr1[,c(4,11)], type="o", col="blue")
legend(45900000, 2, c("Up","Down"), cex=0.8, col=c("red","blue"), pch=21:21, lty=1:2)

abline(v=c(bands), lty=5)
abline(h=1, col="red", lty=5)
abline(h=-1, col="blue", lty=5)
dev.off()






##
pdf(paste0(wd.de, "arrest_expression_chr1_pcg.pdf"), height=4, width=10)
plot(de.up.chr1[,c(4,10)], type="o", pch=22, col="red", xlim=c(35000000,47500000))
#plot(de.up.chr1[,c(4,11)], type="o", col="red", ylim=c(-2.5,3))
lines(de.up.chr1[,c(4,9)], type="o", pch=22, col="blue")
dev.off()

pdf(paste0(wd.de, "arrest_expression_chr1_pcg_down.pdf"), height=4, width=10)
plot(de.down.chr1[,c(4,10)], type="o", pch=22, col="red", xlim=c(35000000,47500000))
#plot(de.up.chr1[,c(4,11)], type="o", col="red", ylim=c(-2.5,3))
lines(de.down.chr1[,c(4,9)], type="o", pch=22, col="blue")
dev.off()

###
##
de.all.chr1 <- rbind(de.down.chr1, de.up.chr1)
de.all.chr1 <- de.all.chr1[order(de.all.chr1$start_position),]

pdf(paste0(wd.de, "arrest_expression_chr1_pcg_all.pdf"), height=4, width=10)
plot(de.all.chr1[,c(4,10)], type="o", pch=21, col="red", xlim=c(35000000,47500000), main="Gene Expression Level", xlab="Chromosome 1", ylab="log2(TPM + 0.01)")
#plot(de.up.chr1[,c(4,11)], type="o", col="red", ylim=c(-2.5,3))
lines(de.all.chr1[,c(4,9)], type="o", pch=22, col="blue")

#legend(x=legend.x, y=ymax, legend=c("Overexpressed", "Underexpressed"), col=c("red", "blue"), pch=21)
legend(34650000, 10.1, c("RB1","WT"), cex=0.8, col=c("red","blue"), pch=21:22, lty=1:2)
dev.off()





# -----------------------------------------------------------------------------
# Last Modified: 04/04/17
# -----------------------------------------------------------------------------
#wd.src <- "/Users/tpyang/Work/local/R"
wd.src <- "/re/home/tyang2/local/R"
sources <- c("DailyMeal.R", "DifferentialAnalysis.R")
invisible(sapply(sources, function(s) source(file.path(wd.src, "handbook-of", s))))
load(file.path(wd.src, "guide-to-the", "hg19.transcript.RData"))

##
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/rearrangements/"
wd <- "/re/home/tyang2/LCNEC/analysis/expression/"
#wd.ngs <- "/ngs/cangen/tyang2/ngs/LCNEC/WGS/"
wd.ngs <- "/ngs/ngs2/DNAseq/LCNEC/LCNEC/"
setwd(wd)

#samples.ngs <- rownames(pheno.expr.final)
#save(samples.ngs, file="samples.ngs.RData")
load(file="samples.ngs.RData")
genes <- c("TP53", "RB1", "XRN1", "XRN2", "ATRX", "BRCA1", "BRCA2", "BARD1", "SMARCA4", "NEDD4")

results <- toTable(0, length(genes), length(samples.ngs), genes)
rownames(results) <- samples.ngs
for (s in 1:length(samples.ngs)) {
   sample <- samples.ngs[s]
 
   ## INPUT: INV
   #snv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_muts.txt")
   snv.file <- paste0(wd.ngs, sample, "_muts.txt")
   snv <- readTable(snv.file, header=T, rownames=F, sep="\t")
   snv <- subset(snv, Type_1 != "silent")

   for (g in 1:length(genes)) {
      gene <- genes[g]
      mutations <- subset(snv, Gene_Hugo == gene)
      results[s, g] <- nrow(mutations)
   }
}
results[order(results$RB1),]











# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# http://www.bioconductor.org/help/workflows/rnaseqGene/#eda
# https://onlinecourses.science.psu.edu/stat505/node/55
# Last Modified: 29/03/17
# -----------------------------------------------------------------------------
pca <- getPCA(t(expr.magic))
#summary(pca)

##
wd.pca <- paste0(wd, "analysis/expression/data/")
#samples.expr <- names(expr.d.e)
#pheno.expr <- pheno[samples.expr,]

traits <- pheno.expr[,c("Sex", "Smoking_Hx", "Stage_UICC")]
colnames(traits) <- c("Sex", "Smoking History", "UICC Stage")
plotPCAs(1, 2, pca, traits, wd.pca, NA, 75, -50)   ## TO-DO

traits <- pheno.expr[,c("Classification", "Survival_censor")]
colnames(traits) <- c("Histological Subtype", "Survival")
plotPCAs(1, 2, pca, traits, wd.pca, NA, 75, -50, c("deepskyblue","purple3","red","blue"))   ## TO-DO

traits <- pheno.expr[,c("Age", "Classification2")]
colnames(traits) <- c("Age", "Histological Subtype")
plotPCAs(1, 2, pca, traits, wd.pca, NA, 60, -50)   ## TO-DO











## Find outliers
outliers <- subset(pca.de[,c("PC1", "PC2")], PC1 > 2.5)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 0)

outliers <- subset(subset(pca.de[,c("PC1", "PC2")], PC2 > 2), PC1 > 0.5)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW != "NA")

outliers <- subset(pca.de[,c("PC1", "PC2")], PC1 < 0.5)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 1)

outliers <- pca.de[c("S00714", "S02229"), c("PC1", "PC2")]
subset(cbind(outliers, pheno[rownames(outliers),]))

## Find outliers
outliers <- subset(pca.de.resid[,c("PC1", "PC2")], PC1 < 1)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 1)

outliers <- subset(subset(pca.de.resid[,c("PC1", "PC2")], PC1 > 0), PC1 < 1)
subset(cbind(outliers, pheno[rownames(outliers),]), RB1NEW == 0)

outliers <- pca.de.resid[c("S00714", "S02229"), c("PC1", "PC2")]
subset(cbind(outliers, pheno[rownames(outliers),]))

###
## LCNEC subtypes on D.E genes
trait <- pheno.expr[,"Classification"]

file.main <- c("LCNEC subtypes on 51 D.E. Genes", "FPKM Expression")
plotPCA(1, 2, pca.de, trait, wd.de, "Classification_51DE", file.main, 5, -2.2, c("deepskyblue","gold","blue","red"))

file.main <- c("LCNEC subtypes on 18 D.E. Genes", "2-Subtype Residuals")
plotPCA(1, 2, pca.de.resid2, trait, wd.de, "Classification_18DE", file.main, 2.97, -2.87, c("deepskyblue","gold","blue","red"))

file.main <- c("LCNEC subtypes on 4 D.E. Genes", "4-Subtype Residuals")
plotPCA(1, 2, pca.de.resid, trait, wd.de, "Classification_4DE", file.main, 1.427, -1.365, c("deepskyblue","gold","blue","red"))

## Find outliers
outliers <- subset(subset(pca.de[,c("PC1", "PC2")], PC2 > 2.5), PC1 > 0)
subset(cbind(outliers, pheno[rownames(outliers),]), Classification != "LCNEC")

# -----------------------------------------------------------------------------
# Fisher's exact (RB1NEW vs WT)
# Last Modified: 08/02/17
# -----------------------------------------------------------------------------
trait1 <- "RB1NEW"
trait2 <- "Classification2"
table <- getFishersTable(pheno.expr, trait1, trait2) 
table
#                           RB1NEW_0 RB1NEW_1
# Classification2_LCNEC           30       14
# Classification2_LCNEC+Mix        4        6
testFishers(table)
# [1] 0.1470132

nrow(subset(subset(pheno.expr, Classification == "LCNEC"), RB1NEW != 1))
# [1] 30
nrow(subset(subset(pheno.expr, Classification == "LCNEC"), RB1NEW == 1))
# [1] 14
nrow(subset(subset(pheno.expr, Classification != "LCNEC"), RB1NEW != 1))
# [1] 4
nrow(subset(subset(pheno.expr, Classification != "LCNEC"), RB1NEW == 1))
# [1] 6
fisher.test(rbind(c(30, 14), c(4, 6)))$p.value
# [1] 0.1470132

# -----------------------------------------------------------------------------
# Differential expression (LCNEC vs LCNEC-Mix)
# Last Modified: 28/02/17
# -----------------------------------------------------------------------------
samples <- sort(rownames(pheno))
samples.expr <- intersect(colnames(expr.d.e.log2), samples)
pheno.expr <- pheno[samples.expr,]

wd.de <- paste0(wd, "analysis/expression/de/")
test <- "Wilcox"   ## Wilcoxon/Wilcox/U/Students/ttest
test.fdr <- "Q"    ## Q/BH
effect <- 1
fdr <- 0.05

## DE (LCNEC vs Mix)
predictor <- "Classification2"
predictor.wt <- "LCNEC"
file.de.sub <- paste0(wd.de, "LCNEC_DE_SUB_RPKM_", predictor, "_", test, "_", test.fdr)

results.de.sub <- differentialAnalysis(expr.d.e.log2, pheno.expr, predictor, predictor.wt, test, test.fdr)
writeTable(results.de.sub, paste0(file.de.sub, ".txt"), colnames=T, rownames=T, sep="\t")

## Only varible genes
writeTable(onlyVariable(results.de.sub, effect), paste0(file.de, "_Effect>", effect, ".txt"), colnames=T, rownames=T, sep="\t")





# -----------------------------------------------------------------------------
# Differential expression on the residuals of expression (LCNEC vs Mix; n=69-15NA)
# Last Modified: 01/03/17
# -----------------------------------------------------------------------------
covariate <- "RB1NEW"

## DE on the residuals (LCNEC vs Mix)
expr.resid <- residualsOf(expr.d.e.log2, pheno.expr, covariate)
results.de.resid <- differentialAnalysis(expr.resid, pheno.expr, predictor, predictor.wt, test, test.fdr)

## Use the original FPKM effect size
results.de.resid[,8:10] <- results.de.sub[rownames(results.de.resid), 8:10]
writeTable(results.de.resid, paste0(file.de, "_Resid.txt2"), colnames=T, rownames=T, sep="\t")

## Only varible genes
writeTable(onlyVariable(results.de.resid, effect), paste0(file.de, "_Resid_Effect>", effect, ".txt2"), colnames=T, rownames=T, sep="\t")





# -----------------------------------------------------------------------------
# Sandbox: Plot un-detected, un-expressed genes
# Last Modified: 07/03/17
# -----------------------------------------------------------------------------
expr.log2     <- log2(1 + expr)
expr.d.log2   <- log2(1 + expr.d)
expr.d.e.log2

undetected <- setdiff(rownames(expr.log2), rownames(expr.d.log2))
unexpressed <- setdiff(rownames(expr.d.log2), rownames(expr.d.e.log2))

expr.unexpressed.log2 <- expr.d.log2[unexpressed,]
medians.unexpressed.log2 <- median0(expr.unexpressed.log2)
unexpressed.1 <- unexpressed[which(medians.unexpressed.log2 > 1)]
unexpressed.2 <- unexpressed[which(medians.unexpressed.log2 < 1)]

expr.removed.log2 <- expr.log2[c(undetected, unexpressed.1, unexpressed.2),]










# -----------------------------------------------------------------------------
# Plot differential analysis (RB1NEW vs WT)
# Last Modified: 09/02/17
# -----------------------------------------------------------------------------
plotDE <- function(expr, samples.expr.mut, name, name2, legend, legend.y) {
   red  <- rgb(1, 0, 0, 0.5)
   blue <- rgb(0, 0, 1, 0.5)
   samples.expr.wt <- setdiff(colnames(expr), samples.expr.mut)
   xmin <- min(expr[name,])
   xmax <- max(as.numeric(expr[name,]))
   
   pdf(paste0("fpkm-d-e-log2_", name2, "_", name, ".pdf"))
   hist(as.numeric(expr[name, samples.expr.mut]), col=red, xlim=c(xmin, xmax+0.1), ylim=c(0, legend.y), main=paste0(name2, " (", name, ")"), xlab="log2(1 + FPKM)")
   hist(as.numeric(expr[name, samples.expr.wt]), col=blue, add=T)
 
   box()
   legend(x=xmin, y=legend.y, legend=legend, fill=c(red, blue), border=c(red, blue))
   dev.off()
}

legend <- c("RB1 (n=20)", "WT (n=49)")
plotDE(expr.d.e.log2, samples.expr.mut, "NM_004526", "MCM2", legend, 12)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_053056", "CCND1", legend, 14)

plotDE(expr.d.e.log2, samples.expr.mut, "NM_000465", "BARD1", legend, 11)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_203394", "E2F7", legend, 13)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_001204192", "TP73", legend, 30)
plotDE(expr.d.e.log2, samples.expr.mut, "NM_058195", "CDKN2A", legend, 11)

# -----------------------------------------------------------------------------
# Manual validation
# Last Modified: 09/02/17
# -----------------------------------------------------------------------------
samples.expr.mut
#  [1] "S00550" "S00567" "S00587" "S00609" "S00710" "S00718" "S00732" "S00739" "S00794" "S00795" "S01407" "S01496" "S01499" "S01505" "S01508" "S01574" "S01581" "S02236"
# [19] "S02259" "S02266"

getRefGene("name2", "MCM2")
# name chrom strand   txStart     txEnd  cdsStart    cdsEnd exonCount
# NR_073375 NR_073375  chr3      + 127317199 127341278 127341278 127341278        16
# NM_004526 NM_004526  chr3      + 127317199 127341278 127317309 127340616        16

wilcox.test(as.numeric(expr.d.e.log2["NM_004526",]) ~ pheno.expr[,"RB1NEW"], exact=F)$p.value
# [1] 2.583733e-07
median(as.numeric(expr.d.e.log2["NM_004526", samples.expr.wt]))
# [1] 4.515555
median(as.numeric(expr.d.e.log2["NM_004526", samples.expr.mut]))
# [1] 5.656553

# > wilcox.test(as.numeric(expr["NM_004526",]) ~ pheno.expr[,trait], exact=F)
# Wilcoxon rank sum test with continuity correction
# data:  as.numeric(expr["NM_004526", ]) by pheno.expr[, trait]
# W = 100, p-value = 2.584e-07
# alternative hypothesis: true location shift is not equal to 0
# 
# > wilcox.test(as.numeric(expr["NM_004526",]) ~ pheno.expr[,trait])
# Wilcoxon rank sum test
# data:  as.numeric(expr["NM_004526", ]) by pheno.expr[, trait]
# W = 100, p-value = 1.593e-08
# alternative hypothesis: true location shift is not equal to 0

# -----------------------------------------------------------------------------
# Differential expression (LCNEC vs LCNEC-Mix; n=69)
# Last Modified: 06/02/17
# -----------------------------------------------------------------------------
wd.de <- paste0(wd, "analysis/expression/de/")
test <- "Wilcox"   ## Wilcoxon/Wilcox/U/Students/ttest
test.fdr <- "Q"    ## Q/BH
effect <- 1
fdr <- 0.05

###
## DE (LCNEC vs LCNEC-Mix)
pheno.expr <- pheno[samples.expr,]   ## NOTE
predictor <- "Classification2"
predictor.wt <- "LCNEC"
file.de <- paste0(wd.de, "LCNEC_DE_RPKM_", predictor, "_", test, "_", test.fdr)

results.de <- differentialAnalysis(expr.d.e.log2, pheno.expr, predictor, predictor.wt, test, test.fdr)
writeTable(results.de, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")

## Only varible genes
pdf(paste0(wd.de, "hist_expr-d-e_effect.pdf"))
hist(abs(results.de$Effect))
dev.off()

writeTable(onlyVariable(results.de, effect), paste0(file.de, "_Effect>", effect, ".txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots (LCNEC vs LCNEC-Mix)
# Last Modified: 07/02/17
# -----------------------------------------------------------------------------
file.vol <- paste0(file.de, "_Effect>", effect, ".pdf")
file.main <- "55 LCNEC vs 14 LCNEC-Mix"
plotVolcano(results.de, fdr, effect, file.vol, file.main, 1.8)
