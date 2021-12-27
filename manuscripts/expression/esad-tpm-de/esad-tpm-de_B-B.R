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
samples$tumorcontent_per_sample <- as.numeric(sub("%", "", samples$tumorcontent_per_sample))/100
samples <- subset(samples, GROUP_ID == 0)

##
type1 <- c("P-0014-B", "P-0023-B", "P-0002-B", "P-0009-B", "P-0018-B", "P-0019-B", "P-0003-B", "P-0007-B", "P-0015-B", "P-0021-B", "P-0010-B", "P-0013-B")
samples$TYPE <- 1
rownames(samples) <- samples$PATIENT_ID2
samples[type1, ]$TYPE <- 0
rownames(samples) <- samples$SAMPLE_ID
samples$TYPE <- as.factor(samples$TYPE)

##
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 1)   ## Use pseudocount=0.01

tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
nrow(tpm.gene.log2)

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 22/08/20
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : X (1) vs B (0) as factor
argv      <- data.frame(predictor="TYPE", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("DE_EAC_tpm-gene-median0_B-vs-B_wilcox_q_n23")
file.main <- paste0("B (n=12) vs B (n=11) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

###
genes.BB.up.pe3   <- rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC >= 0))
genes.BB.down.pe3 <- rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC < 0))

file.name <- paste0("genes_ANOVA_N-B-X_p1e-9_n37_up")
writeTable(genes.anova.NBX.up.pe9, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=F, rownames=F, sep="\t")
file.name <- paste0("genes_ANOVA_N-B-X_p1e-9_n84_down")
writeTable(genes.anova.NBX.down.pe9, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=F, rownames=F, sep="\t")






# -----------------------------------------------------------------------------
# ANNOVA
# Last Modified: 22/12/21
# -----------------------------------------------------------------------------
testANOVA <- function(expr, group) {
   fit1 <- lm(expr ~ group + 1)
   fit2 <- lm(expr ~ 1)
   a1 <- anova(fit1, fit2)
 
   return(a1$Pr[2])
}

colnames <- c("P", "FDR", "N", "B", "X", "FC_B_N", "FC_X_N")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## ANNOVA
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) testANOVA(as.numeric(tpm.gene.log2[x,]), samples$GROUP_ID2))

## Log2 fold change
de$N <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 0)))
de$B <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 1)))
de$X <- median00(tpm.gene.log2, rownames(subset(samples, GROUP_ID2 == 2)))
de$FC_B_N <- de$B - de$N
de$FC_X_N <- de$X - de$N

## FDR
library(qvalue)
de$FDR   <- qvalue(de$P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "ANOVA_EAC_tpm-gene-median0_q_n44.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "ANOVA_EAC_tpm-gene-median0_q_n44.txt"), colnames=T, rownames=F, sep="\t")

genes.anova.up.pe9   <- rownames(subset(subset(de.tpm.gene, P < 1E-9), FC_X_N >= 0))
genes.anova.down.pe9 <- rownames(subset(subset(de.tpm.gene, P < 1E-9), FC_X_N < 0))

file.name <- paste0("genes_ANOVA_p1e-9_n463_up")
writeTable(genes.NXB.up.pe4, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=F, rownames=F, sep="\t")
file.name <- paste0("genes_ANOVA_p1e-9_n414_down")
writeTable(genes.NXB.down.pe4, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=F, rownames=F, sep="\t")

###
## FDR (N, B, X)
de.NBX <- de[c(genes.NBX.up, genes.NBX.down),]

library(qvalue)
de.NBX$FDR   <- qvalue(de.NBX$P)$qvalue
de.NBX <- de.NBX[order(de.NBX$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de.NBX),], de.NBX)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "ANOVA_EAC_tpm-gene-median0_N-B-X_q_n44.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "ANOVA_EAC_tpm-gene-median0_N-B-X_q_n44.txt"), colnames=T, rownames=F, sep="\t")

genes.anova.NBX.up.pe9   <- rownames(subset(subset(de.tpm.gene, P < 1E-9), FC_X_N >= 0))
genes.anova.NBX.down.pe9 <- rownames(subset(subset(de.tpm.gene, P < 1E-9), FC_X_N < 0))

file.name <- paste0("genes_ANOVA_N-B-X_p1e-9_n37_up")
writeTable(genes.anova.NBX.up.pe9, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=F, rownames=F, sep="\t")
file.name <- paste0("genes_ANOVA_N-B-X_p1e-9_n84_down")
writeTable(genes.anova.NBX.down.pe9, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=F, rownames=F, sep="\t")

###
##
plot.main <- "28 differentially expressed genes between N, B and X"
plot.de <- file.path(wd.de.plots, "volcanoplot_ANOVA_EAC_median0_N-B-X_p1e-12_n=14")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
overlaps <- intersect(genes$GENE, de.tpm.gene$external_gene_name)
genes <- genes[overlaps,]
file.main <- c(plot.main, "")
file.de <- paste0(plot.de, ".pdf")
de.tpm.gene$LOG2_FC <- de.tpm.gene$FC_X_N
plotVolcano(de.tpm.gene, 1.00E-12, genes, file.de, file.main, ymax=19.5)

# -----------------------------------------------------------------------------
# Heatmap
# Last Modified: 22/12/21; 15/04/21; 11/01/20
# -----------------------------------------------------------------------------
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)

## B subgroup
samples.b <- samples[rownames(subset(samples, GROUP_ID2 == 1)),]
tpm.gene.b <- tpm.gene[, rownames(samples.b)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.b.log2   <- log2(tpm.gene.b + 1)
tpm.gene.b.log2.m <- getLog2andMedian(tpm.gene.b, 1)

## X subgroup
samples.x <- samples[rownames(subset(samples, GROUP_ID2 == 2)),]
tpm.gene.x <- tpm.gene[, rownames(samples.x)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.x.log2   <- log2(tpm.gene.x + 1)
tpm.gene.x.log2.m <- getLog2andMedian(tpm.gene.x, 1)

# -----------------------------------------------------------------------------
# Heatmap
# Links: https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
# Last Modified: 13/01/20
# -----------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#install.packages("gplots")

library("DESeq")
library(gplots)

#genes.rownames <- c("TRIM29", "CDA", "S100A2", "CSTA", "ARL4D", "ALOX12", "LY6G6C", "KRT1", "LGALS7", "KRT6A", "MAL", "KRT4", "SLURP1", "SERPINB3", "SPRR2C", "PAX9", "SULT2B1", "NFRKB", "VAT1")
#genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L", "HDGF", "SPDL1", "ITPR3", "IQCE", "KIF11", "AURKA", "DSG2", "NCEH1", "FRK", "HMMR", "HNRNPC", "TIMELESS")
#genes.rownames <- c("CLSPN", "CENPF", "PDIA6", "TPX2", "MLEC", "EIF2AK1", "KIF20B", "MKI67", "HELLS", "HNRNPH1", "C6orf62", "CBX3", "NT5C3A", "ASPM", "DIAPH3", "KIF14", "PRC1", "TOPBP1", "TOP2A", "ARHGAP11A", "DTX3L", "HDGF", "SPDL1", "ITPR3", "IQCE", "KIF11", "AURKA", "DSG2", "NCEH1", "FRK", "HMMR", "HNRNPC", "TIMELESS", "KIF4A", "C10orf12", "KNTC1", "SLC30A4", "TPD52", "ASTN2", "BUB1", "ARHGAP11B", "NUSAP1", "RP11-366L20.2", "NOP56", "TBL1XR1", "ACBD5", "MCM8", "ESCO2", "SSRP1", "CDK1", "KIAA0319L", "LAMC2", "KCNE3", "SPINK9", "ABCC3", "RAB3IP", "CGN", "CDCA2")

#genes.rownames <- c("PDIA6", "CLSPN", "MLEC", "C6orf62", "HNRNPH1", "TPX2", "CBX3", "CENPF", "TOPBP1", "EIF2AK1", "ITPR3", "KIF14", "DTX3L", "HELLS", "NCEH1", "SPDL1", "C10orf12", "KIF18A", "FRK", "HNRNPC", "IQCE", "HMMR", "MKI67", "ASTN2", "SPINK9", "NT5C3A", "KIF11", "KIF4A", "KNTC1", "PRC1", "RP11-366L20.2", "CDCA2", "MCM8", "ACBD5", "ARHGAP11A", "DIAPH3", "KIF20B", "TBL1XR1", "TOP2A", "TIMELESS")
#genes.rownames <- c("PDIA6", "CLSPN", "MLEC", "C6orf62", "HNRNPH1", "TPX2", "CBX3", "CENPF", "TOPBP1", "EIF2AK1", "ITPR3", "KIF14", "DTX3L", "HELLS", "NCEH1", "SPDL1", "C10orf12", "KIF18A", "FRK", "HNRNPC", "IQCE", "HMMR", "MKI67", "ASTN2", "SPINK9", "NT5C3A", "KIF11", "KIF4A", "KNTC1", "PRC1", "RP11-366L20.2", "CDCA2", "MCM8", "ACBD5", "ARHGAP11A", "DIAPH3", "KIF20B", "TBL1XR1", "TOP2A", "TIMELESS", "CDCA8", "BUD31", "TMSB10", "KCNE3", "KDELR1", "SSRP1", "PARP9", "CGN", "CDC6", "KIAA0319L", "EML4", "GMNN", "TNS3", "NUSAP1", "POLA1", "PARP4", "WHSC1", "TMPO", "ESCO2", "DSG2", "ARHGAP11B", "EDEM3", "NCAPG2", "AURKA", "FBXO5", "SSFA2", "OTUD6B", "TMEM63A", "CDK1", "TRIM11", "C9orf64", "HOXB6", "IFIH1", "RASEF", "RAD51AP1", "FAM83A", "HOXC8", "RAB3IP")

genes <- c(genes.NBX.up.pe3, genes.NBX.down.pe3)
genes <- c("ENSG00000133110", "ENSG00000145685", "ENSG00000101825", "ENSG00000088882", "ENSG00000168542", "ENSG00000154096", "ENSG00000183508", "ENSG00000211896", "ENSG00000204262", "ENSG00000108821", "ENSG00000099958", "ENSG00000099875", "ENSG00000170477", "ENSG00000186806", "ENSG00000136810", "ENSG00000198807", "ENSG00000183347", "ENSG00000177508", "ENSG00000196344", "ENSG00000133710", "ENSG00000187642", "ENSG00000163209", "ENSG00000124466", "ENSG00000170423", "ENSG00000226894")
#genes <- c("ENSG00000099875", "ENSG00000170477", "ENSG00000186806", "ENSG00000136810", "ENSG00000198807", "ENSG00000183347", "ENSG00000177508", "ENSG00000196344", "ENSG00000133710", "ENSG00000187642", "ENSG00000163209", "ENSG00000124466", "ENSG00000170423", "ENSG00000226894")

#wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-9_up")
#genes.up <- readTable(file.path(wd.de.reactome, "genes_ANOVA_N-B-X_p1e-9_n37_up.txt"), header=F, rownames=F, sep="\t")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-9_down")
genes <- readTable(file.path(wd.de.reactome, "genes_ANOVA_N-B-X_p1e-9_n84_down.txt"), header=F, rownames=F, sep="\t")
#genes <- c(genes.up, genes.down)

wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-12_up")
genes.up <- readTable(file.path(wd.de.reactome, "genes_ANOVA_N-B-X_p1e-12_n15_up.txt"), header=F, rownames=F, sep="\t")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-12_down")
genes.down <- readTable(file.path(wd.de.reactome, "genes_ANOVA_N-B-X_p1e-12_n13_down.txt"), header=F, rownames=F, sep="\t")
genes <- c(genes.up, genes.down)

###
##
b <- tpm.gene.b.log2[genes,]
colnames(b) <- gsub("P-00", "P", samples.b$PATIENT_ID2)
rownames(b) <- ensGene[genes, ]$external_gene_name
b <- data.matrix(b)

x <- tpm.gene.x.log2[genes,]
colnames(x) <- gsub("P-00", "P", samples.x$PATIENT_ID2)
rownames(x) <- ensGene[genes, ]$external_gene_name
x <- data.matrix(x)

xb <- cbind(x, b)

###
##
heatmap(b)

## Returning the values used for the heatmap
test <- heatmap.2(b, scale="row")
b[rev(test$rowInd), test$colInd]

hr <- hclust(as.dist(1-cor(t(b), method="pearson")),  method="complete")
hc <- hclust(as.dist(1-cor(b,    method="spearman")), method="complete")

colfunc <- colorRampPalette(c("red", "black", "green"))
heatmap.2(b, col=colfunc(7), scale="row", trace="none")

###
##
heatmap(x)

## Returning the values used for the heatmap
test <- heatmap.2(x, scale="row")
x[rev(test$rowInd), test$colInd]

hr <- hclust(as.dist(1-cor(t(x), method="pearson")),  method="complete")
hc <- hclust(as.dist(1-cor(x,    method="spearman")), method="complete")

colfunc <- colorRampPalette(c("yellow", "black", "green"))
heatmap.2(x, col=colfunc(7), scale="row", trace="none")

###
##
heatmap(xb)

## Returning the values used for the heatmap
test <- heatmap.2(xb, scale="row")
xb[rev(test$rowInd), test$colInd]

hr <- hclust(as.dist(1-cor(t(xb), method="pearson")),  method="complete")
hc <- hclust(as.dist(1-cor(xb,    method="spearman")), method="complete")

colfunc <- colorRampPalette(c("red", "black", "green"))
heatmap.2(xb, col=colfunc(7), scale="row", trace="none")

# -----------------------------------------------------------------------------
# Tumour contents
# Last Modified: 23/12/21
# -----------------------------------------------------------------------------
type1 <- c("P-0014-B", "P-0023-B", "P-0002-B", "P-0009-B", "P-0018-B", "P-0019-B", "P-0003-B", "P-0007-B", "P-0015-B", "P-0021-B", "P-0010-B", "P-0013-B")

samples.b$tumorcontent_per_sample <- as.numeric(sub("%", "", samples.b$tumorcontent_per_sample))/100
samples.b$TYPE <- 1
rownames(samples.b) <- samples.b$PATIENT_ID2
samples.b[type1, ]$TYPE <- 0
rownames(samples.b) <- samples.b$SAMPLE_ID

trait <- as.factor(samples.b$TYPE)
wilcox.test(samples.b$tumorcontent_per_sample ~ trait, exact=F)$p.value











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
trait.v <- c("X", "B", "N")
cols    <- c(purple, red, green)

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
   #xmax <- min(de$LOG2_FC) * -1
   #ymax <- max(de$log10P)
   #ymax <- 8.5
 
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-12, max(de$LOG2_FC)), ylim=c(0, ymax), xlab="X to N fold change [log2]", ylab="P-value significance [-log10]", col="lightgray", main=file.main[1])

   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, paste0("FDR=", fdr, "%"), cex=0.85)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=purple)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=blue)
 
   #abline(v=c(-log2(2), log2(2)), lty=5, col="darkgray")
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
   legend("topleft", legend=c("Up-regulated (N < B < X)", "Down-regulated (N > B > X)"), col=c(purple, blue), pch=19)
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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-12_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-12_down")

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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-12_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[5, 2] <- "FOXO-mediated oxidative stress, metabolic and neuronal genes"

reactome.up <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("N < B < X genes", "37 genes")
file.de <- file.path(wd.de.reactome, "genes_N-B-X_anova_p1e-9_n37_up.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,25,4,2))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 13), xaxt="n", names.arg=reactome.up$Pathway.name, col=purple, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=4, lty=5)

axis(side=1, at=seq(0, 12, by=4))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-12_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[2, 2] <- "Transcription of E2F targets under negative control by DREAM"
#reactome[3, 2] <- "HDR through Homologous Recombination or Single Strand Annealing"

reactome.down <- subset(reactome, Entities.pValue <= 5e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("N > B > X genes", "13 genes")
file.de <- file.path(wd.de.reactome, "genes_B-X_content_p1e-12_n13_down.pdf")

pdf(file.de, height=3, width=7.5)
par(mar=c(4,2,4,25))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-2, lty=5)

axis(side=1, at=seq(-4, 0, by=1), labels=c(4, 3, 2, 1, 0))
mtext(main.text[2], line=0.3)
dev.off()








# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-9_up")
#wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-9_down")

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
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-9_up")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[5, 2] <- "FOXO-mediated oxidative stress, metabolic and neuronal genes"

reactome.up <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("N < B < X genes", "37 genes")
file.de <- file.path(wd.de.reactome, "genes_N-B-X_anova_p1e-9_n37_up.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,25,4,2))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 13), xaxt="n", names.arg=reactome.up$Pathway.name, col=purple, xlab="-log10(p-value)", width=1)   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
abline(v=4, lty=5)

axis(side=1, at=seq(0, 12, by=4))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_N-B-X_anova_p1e-9_down")
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[2, 2] <- "Transcription of E2F targets under negative control by DREAM"
#reactome[3, 2] <- "HDR through Homologous Recombination or Single Strand Annealing"

reactome.down <- subset(reactome, Entities.pValue <= 5e-3)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("N > B > X genes", "84 genes")
file.de <- file.path(wd.de.reactome, "genes_B-X_content_p1e-6_n84_down.pdf")

pdf(file.de, height=2.8, width=7.5)
par(mar=c(4,2,4,25))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-13, 0), xaxt="n", names.arg="", col=blue, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-4, lty=5)

axis(side=1, at=seq(-12, 0, by=4), labels=c(12, 8, 4, 0))
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
