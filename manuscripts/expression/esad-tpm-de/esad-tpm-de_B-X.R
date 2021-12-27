# =============================================================================
# Manuscript   : 
# Name         : manuscripts/expression/esad-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 21/08/20
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "ESAD"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt2"), header=T, rownames=T, sep="\t")
samples <- samples[!is.na(samples[, "tumorcontent_per_sample"]),]
samples$tumorcontent_per_sample <- as.numeric(sub("%", "", samples$tumorcontent_per_sample))/100

b <- rownames(subset(samples, GROUP_ID == 0))
x <- rownames(subset(samples, GROUP_ID == 1))
samples$GROUP_ID2 <- 0
samples[b,]$GROUP_ID2 <- 0
samples[x,]$GROUP_ID2 <- 1
samples$GROUP_ID2 <- as.numeric(samples$GROUP_ID2)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=0.01

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)
# [1] 15658

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : X (1) vs B (0) as factor
argv      <- data.frame(predictor="GROUP_ID", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("DE_EAC_tpm-gene-median0_B-vs-X_wilcox_q_n46")
file.main <- paste0("X (n=22) vs B (n=23) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

de.tpm.gene.BX.up   <- subset(de.tpm.gene, LOG2_FC >= 0)
de.tpm.gene.BX.down <- subset(de.tpm.gene, LOG2_FC < 0)
#genes.BX      <- subset(de.tpm.gene, P < 1E-6)
#genes.BX.up   <- subset(genes.BX, LOG2_FC >= 0)
#genes.BX.down <- subset(genes.BX, LOG2_FC < 0)

# -----------------------------------------------------------------------------
# N, B, X
# Last Modified: 20/12/21
# -----------------------------------------------------------------------------
de.tpm.gene.BX.NBX <- de.tpm.gene.BX[c(genes.NBX.up, genes.NBX.down),]

genes.BX.NBX.up   <- rownames(subset(subset(de.tpm.gene.BX.NBX, P < 1E-6), LOG2_FC >= 0))
genes.BX.NBX.down <- rownames(subset(subset(de.tpm.gene.BX.NBX, P < 1E-6), LOG2_FC < 0))
genes.BX.NBX <- c(genes.BX.NBX.up, genes.BX.NBX.down)

##
genes <- genes.BX.NBX
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("B vs. X genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_B-vs-X", size=6, file.main, "topleft", trait.v, cols, NULL, flip.x=1, flip.y=1, legend.title=NA)

# -----------------------------------------------------------------------------
# N, X, B
# Last Modified: 21/12/21
# -----------------------------------------------------------------------------
de.tpm.gene.BX.NXB <- de.tpm.gene.BX[c(genes.NXB.up, genes.NXB.down),]

genes.XB.NXB.up   <- rownames(subset(subset(de.tpm.gene.BX.NXB, P < 1E-6), LOG2_FC < 0))
genes.XB.NXB.down <- rownames(subset(subset(de.tpm.gene.BX.NXB, P < 1E-6), LOG2_FC >= 0))
genes.XB.NXB <- c(genes.XB.NXB.up, genes.XB.NXB.down)

##
genes <- genes.XB.NXB
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("X vs. B genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-X-B_n68_X-vs-B", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=-1, flip.y=1, legend.title=NA)










###
##
plot.main <- "387 differentially expressed genes in EAC"
plot.de <- file.path(wd.de.plots, "volcanoplot_DE_EAC_median0_B-vs-X_p1e-6_CNA")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Before (B) vs. After treatment (X)")
file.de <- paste0(plot.de, ".pdf")
de.tpm.gene$LOG2_FC <- de.tpm.gene$GROUP_ID - de.tpm.gene$GROUP_ID_WT
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main, ymax=8.5)




# -----------------------------------------------------------------------------
# PCA vs. Purites
# Last Modified: 02/06/21
# -----------------------------------------------------------------------------
#rownames(samples) <- samples$PATIENT_ID2
#overlaps <- intersect(samples$PATIENT_ID2, purities$sample_name)
#samples <- samples[overlaps,]
#samples$purity <- purities[overlaps,]$purity
#rownames(samples) <- samples$SAMPLE_ID

#load("/Users/tpyang/Work/uni-koeln/tyang2/ESAD/analysis/expression/kallisto/esad-tpm-de/data/PCA_EAC_N-B-X.RData")
scores <- pcaScores(pca.de)
scores <- scores[rownames(samples),]

# -----------------------------------------------------------------------------
# 
# Last Modified: 03/06/21
# -----------------------------------------------------------------------------
plotFACS3 <- function(p, pc, file.name, main.text, xlab.text, ylab.text, col, col2) {
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(pc ~ p, ylab="", xlab=xlab.text, main=main.text[1], col="white", pch=19, cex=2, lwd=0, cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
   points(p[ 1:23], pc[ 1:23], col=col2[1], pch=19, cex=1.7)
   points(p[24:44], pc[24:44], col=col2[2], pch=19, cex=1.7)
   
   lm.fit <- lm(pc ~ p)
   abline(lm.fit, col=col, lwd=4)
 
   cor3 <- cor.test(p, pc, method="spearman", exact=F)
   legend("bottomright", c(paste0("rho = ", round0(cor3[[4]], digits=2)), paste0("p-value = ", scientific(cor3[[3]]))), text.col=cols, text.font=2, bty="n", cex=1.5)
 
   legend("topleft", c("B", "X"), col=col2, pch=19, pt.cex=2.5, cex=1.6)   ##bty="n")
   #axis(side=2, at=seq(-0.2, 0.2, by=0.2), labels=c(-0.2, 0, 0.2), cex.axis=1.5)
   #axis(side=2, at=seq(-0.3, 0.1, by=0.2), labels=c(-0.3, -0.1, 0.1), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.75, cex=1.6)
   #mtext(main.text[2], line=0.3, cex=1.6)
   dev.off()
}

file.name <- file.path(wd.de.plots, "ESAD_Content-vs-PC1_n44_NEW")
main.text <- c(paste("Tumour content vs. PC1 (n=44)"), "")
ylab.text <- paste0("PC", 1, " (", pcaProportionofVariance(pca.de, 1), "%)")
xlab.text <- "Tumour content [%]"                                                                         ## "#619CFF", "#F8766D", "#00BA38"      "skyblue3", "lightcoral", "#59a523"
cols <- "black"
cols2 <- c(red, purple)
plotFACS3(samples$tumorcontent_per_sample, scores$PC1*-1, file.name, main.text, xlab.text, ylab.text, cols, cols2)

# -----------------------------------------------------------------------------
# Spearman's rank correlation
# Last Modified: 03/-6/21; 22/08/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "FDR", "X", "B", "FC_B_X")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$tumorcontent_per_sample, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$tumorcontent_per_sample, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 0)))
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 1)))
de$FC_B_X <- de$B - de$X

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "SRC_EAC_tpm-gene-median0_B-X_vs_content_q_n44.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "SRC_EAC_tpm-gene-median0_B-X_vs_content_q_n44.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
plotSRC <- function(gene, cn, snr, pch, col, pos, xlab.text="") {
   unit <- (max(snr) - min(snr))/10
   xlim <- c(min(snr) - unit, max(snr) + unit)
   unit <- (max(cn) - min(cn))/10
   ylim <- c(min(cn) - unit, max(cn) + unit)
 
   if (xlab.text == "")
      xlab.text <- expression(italic('In silico')~'sorting [rho]')
   ylab.text <- "log2(TPM + 1)"
   id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id
   file.name <- file.path(wd.de.plots, paste0("TPM-vs-SORTING_", genes[g], ""))
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   #plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xaxt="n", xlab=xlab.text, main=paste0(gene, " (", id, ")"), col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
   plot(cn ~ snr, ylim=ylim, xlim=xlim, ylab="", xlab="", main=gene, col="black", pch=pch, cex=2, cex.axis=1.7, cex.lab=1.9, cex.main=2)
 
   lm.fit <- lm(cn ~ snr)
   abline(lm.fit, col=col, lwd=7)
 
   cor <- cor.test(cn, snr, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col, text.font=2, bty="n", cex=1.9)
 
   #axis(side=1, at=seq(-0.5, 0.5, by=0.5), labels=c(-0.5, 0, 0.5), cex.axis=1.7)
   #axis(side=2, at=seq(6, 8, by=1), labels=c(6, 7, 8), cex.axis=1.7)   ## MARS
   #axis(side=2, at=seq(4, 6, by=1), labels=c(4, 5, 6), cex.axis=1.7)   ## GTPBP3
   mtext(ylab.text, side=2, line=2.74, cex=1.85)
   mtext(xlab.text, side=1, line=3.5, cex=1.9)
   #mtext(main.text[2], cex=1.2, line=0.3)
   dev.off()
}

genes <- c("CLSPN", "ECHS1", "ERBB3", "EZH2")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples$tumorcontent_per_sample, 1, red, "bottomright", xlab.text="Tumour content [%]")
}

genes <- c("DDR2", "SMARCA2", "MMP2", "IGF1")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples$tumorcontent_per_sample, 1, purple, "bottomright", xlab.text="Tumour content [%]")
}











# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt2"), header=T, rownames=T, sep="\t")

n <- rownames(subset(samples, GROUP_ID == 2))
x <- rownames(subset(samples, GROUP_ID == 1))
b <- rownames(subset(samples, GROUP_ID == 0))
samples$GROUP_ID2     <- 0
samples[b,]$GROUP_ID2 <- 1
samples[x,]$GROUP_ID2 <- 2
samples$GROUP_ID2 <- as.numeric(samples$GROUP_ID2)

trait <- as.numeric(samples[, "GROUP_ID2"])
trait[which(trait == 0)] <- "N"
trait[which(trait == 1)] <- "B"
trait[which(trait == 2)] <- "X"
trait.v <- c("B", "X", "N")
cols    <- c(red, purple, blue)

###
##
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("Expressed genes (n=15,658)", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_ALL", size=6, file.main, "topright", trait.v, cols, NULL, flip.x=-1, flip.y=1, legend.title=NA)

###
##
genes <- intersect(genes.G1S, rownames(tpm.gene))
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("EAC samples (n=68)", "G1/S genes (n=304)")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_G1-S", size=6, file.main, "topright", trait.v, cols, NULL, flip.x=-1, flip.y=1, legend.title=NA)

genes <- intersect(genes.G2M, rownames(tpm.gene))
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("EAC samples (n=68)", "G2/M genes (n=876)")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_G2-M", size=6, file.main, "topright", trait.v, cols, NULL, flip.x=-1, flip.y=1, legend.title=NA)

###
##
genes <- genes.Content
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("Tumour content genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_Content", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=-1, flip.y=1, legend.title=NA)

genes <- genes.BX
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("B vs. X genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_B-vs-X", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=-1, flip.y=-1, legend.title=NA)

genes <- genes.NB
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("N vs. B genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_N-vs-B", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=-1, flip.y=-1, legend.title=NA)

genes <- genes.NX
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("N vs. X genes (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_N-vs-X", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=1, flip.y=-1, legend.title=NA)

genes <- genes.BandXcontent
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c(paste0("B vs. X, corrected for content (n=", length(genes), ")"), "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_n68_B-vs-X-content", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=1, flip.y=1, legend.title=NA)

###
##
genes <- c("CXCL14", "AMIGO2", "GPX2", "D4S234E", "PTPRZ1", "AKR1C2", "CSTA", "KRT16", "DSG1", "IL1A", "GJB5", "DSC3", "PKP1", "CLCA2", "IVL", "TP63", "KRT6A", "SPRR1A", "SPRR1B", "KRT75", "RHCG", "S100A7", "PTHLH", "MMP10", "S100A2", "TRIM29")
genes <- rownames(getGenes(genes))
genes <- intersect(genes, rownames(tpm.gene))
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("EAC samples (n=68)", paste0("HUPER_BREAST_BASAL_VS_LUMINAL (n=", length(genes), ")"))
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_Guo_Huper", size=6, file.main, "bottomleft", trait.v, cols, NULL, flip.x=1, flip.y=-1, legend.title=NA)

genes <- c("TRIM29", "CDA", "S100A2", "CSTA", "ARL4D", "LY6G6C", "KRT6A", "MAL", "KRT4", "SLURP1", "SERPINB3", "PAX9", "SULT2B1", "NFRKB", "VAT1")
genes <- rownames(getGenes(genes))
genes <- intersect(genes, rownames(tpm.gene))
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("EAC samples (n=68)", paste0("WANG_BE_AND_EAC_DN (n=", length(genes), ")"))
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_Guo_Wang", size=6, file.main, "bottomleft", trait.v, cols, NULL, flip.x=1, flip.y=1, legend.title=NA)

###
##
genes <- intersect(genes.NBX.25, rownames(tpm.gene))
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("N, B, X genes (n=25)", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_NBX_25", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=1, flip.y=1, legend.title=NA)

genes <- intersect(genes.NBX.199, rownames(tpm.gene))
test <- tpm.gene[genes, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
file.main <- c("N, B, X genes (n=199)", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_NBX_199", size=6, file.main, "bottomright", trait.v, cols, NULL, flip.x=1, flip.y=1, legend.title=NA)

# -----------------------------------------------------------------------------
# Outputs
# Last Modified: 21/12/21; 10/12/21
# -----------------------------------------------------------------------------
colnames <- c("N", "B", "X")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

de$N <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 0)))
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 1)))
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 2)))

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

de.tpm.gene.1 <- de.tpm.gene[genes.NBX.down.pe3,]
de.tpm.gene.1 <- de.tpm.gene.1[order(de.tpm.gene.1$X, decreasing=F),]
writeTable(de.tpm.gene.1, file.path(wd.de.data, "PCA_EAC_N-B-X_n68_down_n14.txt"), colnames=T, rownames=F, sep="\t")




de.tpm.gene.1 <- de.tpm.gene[genes.1,]
de.tpm.gene.1 <- de.tpm.gene.1[order(de.tpm.gene.1$X, decreasing=F),]
writeTable(de.tpm.gene.1, file.path(wd.de.data, "PCA_EAC_N-B-X_n68_Guo_Huper_n22.txt"), colnames=T, rownames=F, sep="\t")

de.tpm.gene.2 <- de.tpm.gene[genes.2,]
de.tpm.gene.2 <- de.tpm.gene.2[order(de.tpm.gene.2$X, decreasing=F),]
writeTable(de.tpm.gene.2, file.path(wd.de.data, "PCA_EAC_N-B-X_n68_Guo_Wang_n15.txt"), colnames=T, rownames=F, sep="\t")












# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure S1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
plotVolcano <- function(de, pvalue, genes, file.de, file.main, ymax) {
   #pvalue <- fdrToP(fdr, de)
   #fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- min(de$LOG2_FC) * -1
   #ymax <- max(de$log10P)
   #ymax <- 8.5
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="X to B fold change [log2]", ylab="P-value significance [-log10]", col="lightgray", main=file.main[1])

   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=purple)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=red)
 
   abline(v=c(-log2(2), log2(2)), lty=5, col="darkgray")
   abline(h=c(-log10(pvalue)), lty=5, col="black")
   
   if(nrow(genes) != 0) {
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
   }
   
   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topleft", legend=c("Up-regulated (B < X)", "Down-regulated (B > X)"), col=c(purple, red), pch=19)
   dev.off()
}

##
plot.main <- "205 genes significantly correlated with tumour content"
plot.de <- file.path(wd.de.plots, "volcanoplot_SRC_EAC_median0_B-X_vs_content_p1e-6")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Before (B) vs. After treatment (X)")
file.de <- paste0(plot.de, ".pdf")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_X_B
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main, ymax=10.5)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-X_content_p1e-6_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-X_content_p1e-6_down")

list <- ensGene[,c("ensembl_gene_id",	"external_gene_name")]

reactome <- read.csv(file.path(wd.de.reactome, "result.csv"))
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
writeTable(reactome, file.path(wd.de.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-X_content_p1e-6_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[5, 2] <- "FOXO-mediated oxidative stress, metabolic and neuronal genes"

reactome.up <- subset(reactome, Entities.pValue <= 1e-3)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated after treatment", "135 genes")
file.de <- file.path(wd.de.reactome, "genes_B-X_content_p1e-6_n135_up.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,25,4,2))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 9), xaxt="n", names.arg=reactome.up$Pathway.name, col=purple, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=3, lty=5)

axis(side=1, at=seq(0, 8, by=2))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_B-X_content_p1e-6_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[2, 2] <- "Transcription of E2F targets under negative control by DREAM"
reactome[3, 2] <- "HDR through Homologous Recombination or Single Strand Annealing"

reactome.down <- subset(reactome, Entities.pValue <= 1e-3)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated after treatment", "70 genes")
file.de <- file.path(wd.de.reactome, "genes_B-X_content_p1e-6_n70_down.pdf")

pdf(file.de, height=2.8, width=7.5)
par(mar=c(4,2,4,25))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-10, 0), xaxt="n", names.arg="", col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-3, lty=5)

axis(side=1, at=seq(-8, 0, by=2), labels=c(8, 6, 4, 2, 0))
mtext(main.text[2], line=0.3)
dev.off()












# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "ESAD"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/3RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "esad_3rna_n68.txt2"), header=T, rownames=T, sep="\t")
samples <- samples[rownames(subset(samples, GROUP_ID != 2)),]

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=0.01

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)
# [1] 15658








# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
test <- tpm.gene[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.de.data, paste0("PCA_EAC_N-B-X.RData")))

trait <- as.numeric(samples[, "GROUP_ID2"])
trait[which(trait == 0)] <- "N"
trait[which(trait == 1)] <- "B"
trait[which(trait == 2)] <- "X"

file.main <- c("EAC samples (n=68)", "")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "PCA_EAC_N-B-X_n68_EAC", size=6.5, file.main, "topright", c(green, red, purple), NULL, flip.x=-1, flip.y=1, legend.title=NA)

# -----------------------------------------------------------------------------
# PCA (Controlled for PC1)
# Last Modified: 26/03/21; 23/01/21
# -----------------------------------------------------------------------------
test <- tpm.gene.res[, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de.res <- getPCA(t(test))
save(pca.de.res, file=file.path(wd.de.data, paste0("PCA_EAC_N-B-X_RES.RData")))

trait <- as.numeric(samples[, "GROUP_ID2"])
trait[which(trait == 0)] <- "N"
trait[which(trait == 1)] <- "B"
trait[which(trait == 2)] <- "X"

file.main <- c("EAC samples (Controlled for PC1)", "")
plotPCA(1, 2, pca.de.res, trait, wd.de.plots, "PCA_EAC_B-X-N_n68_RES", size=6.5, file.main, "topright", c(green, red, purple), NULL, flip.x=-1, flip.y=1, legend.title=NA)
