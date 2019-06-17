# =============================================================================
# Manuscript   : 
# Chapter I    : 
# Figure(s)    : 
# Name         : manuscripts/expression/sclc-wgs-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/05/19
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "NBL"
base <- tolower(BASE)

wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-wgs-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
wd.de.gsea  <- file.path(wd.de, "gsea")

wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")

samples <- readTable(file.path(wd.wgs, "nbl_wgs_n57-1.cm2"), header=T, rownames=T, sep="")[,-1]

##
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
samples.rna <- readTable(file.path(wd.rna, "nbl_rna_n54.list"), header=F, rownames=T, sep="")
tpm.gene <- tpm.gene[, rownames(samples.rna)]
colnames(tpm.gene) <- samples.rna$V2

rownames(samples.rna) <- samples.rna$V2
overlaps <- intersect(rownames(samples), rownames(samples.rna))
samples <- samples[overlaps,]
tpm.gene <- tpm.gene[, overlaps]

samples <- readTable(file.path(wd.wgs, "nbl_wgs_n53.cm2"), header=T, rownames=T, sep="")[,-1]

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; CM2 vs CM1)
# Last Modified: 22/05/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : RB1_MUT (1) vs RB1_WT (0) as factor
file.name <- paste0("de_", base, "_wgs-gene_cm2_wilcox_q_n53")

wgs.gene <- c()
de <- c()
for (c in 1:22) {
   chr <- chrs[c]
   samples[, chr] <- as.factor(samples[, chr])
   
   wgs.gene.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.ensGene_", chr, ".txt.gz")), header=T, rownames=T, sep="\t")[, rownames(samples)]
   wgs.gene.chr <- wgs.gene.chr[intersect(rownames(wgs.gene.chr), expressed),]
   wgs.gene.chr.log2 <- log2(wgs.gene.chr + 0.01)
   wgs.gene <- rbind(wgs.gene, wgs.gene.chr)
   
   argv <- data.frame(predictor=chr, predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
   de.chr <- differentialAnalysis(wgs.gene.chr.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
   colnames(de.chr)[3:4] <- c("CM1", "CM2")
   de <- rbind(de, de.chr)
}
de$FDR <- testFDR(de$P, argv$test.fdr)
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.wgs.gene <- cbind(annot[rownames(de),], de)

save(de.wgs.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.wgs.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(tpm.gene), rownames(wgs.gene))
tpm.gene <- tpm.gene[overlaps,]
wgs.gene <- wgs.gene[overlaps, colnames(tpm.gene)]   ## ADD 11/06/19
wgs.gene.log2 <- log2(wgs.gene + 0.01)
tpm.gene.log2 <- log2(tpm.gene + 0.01)
save(wgs.gene, file=file.path(wd.de.data, paste0(base, "_wgs.gene_n53.RData")))
save(tpm.gene, file=file.path(wd.de.data, paste0(base, "_tpm.gene_n53.RData")))
# > dim(tpm.gene)
# [1] 18015    53
# > dim(wgs.gene)
# [1] 18015    53

# -----------------------------------------------------------------------------
# SRC between WGS and TPM
# Last Modified: 23/05/19
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "TPM_P", "TPM_LOG2_FC", "TPM_WEIGHT", "WGS_P", "WGS_LOG2_FC", "WGS_WEIGHT", "KEEP")
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), as.numeric(wgs.gene.log2[x,]), method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), as.numeric(wgs.gene.log2[x,]), method="spearman", exact=F)[[3]])

## Pearson
#de$PCC_R <- mapply(x = 1:nrow(tpm.gene.o.log2), function(x) cor.test(as.numeric(tpm.gene.o.log2[x,]), as.numeric(wgs.gene.o.log2[x,]), method="pearson")$estimate)
#de$PCC_P <- mapply(x = 1:nrow(tpm.gene.o.log2), function(x) cor.test(as.numeric(tpm.gene.o.log2[x,]), as.numeric(wgs.gene.o.log2[x,]), method="pearson")$p.value)

## Log2 fold change
de$TPM_P <- de.tpm.gene[rownames(de),]$P
de$WGS_P <- de.wgs.gene[rownames(de),]$P
de$TPM_LOG2_FC <- de.tpm.gene[rownames(de),]$LOG2_FC
de$WGS_LOG2_FC <- de.wgs.gene[rownames(de),]$LOG2_FC
de$TPM_WEIGHT <- -log10(de$TPM_P) * de$TPM_LOG2_FC
de$WGS_WEIGHT <- -log10(de$WGS_P) * de$WGS_LOG2_FC
de$KEEP <- de$TPM_LOG2_FC * de$WGS_LOG2_FC

## FDR
library(qvalue)
de$Q <- qvalue(de$P)$qvalue
#de$PCC_Q <- qvalue(de$PCC_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.wgs.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.wgs.gene, samples, file=file.path(wd.de.data, "de_nbl_tpm-wgs-gene_cm2_src_q_n53.RData"))
writeTable(de.tpm.wgs.gene, file.path(wd.de.data, "de_nbl_tpm-wgs-gene_cm2_src_q_n53.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# 
# Last Modified: 11/06/19
# -----------------------------------------------------------------------------
ylab.text <- "Differential replication"
xlab.text <- "Differential expression"

de.all <- de.tpm.wgs.gene
de.all <- subset(de.all, P < 1E-03)
de.all <- subset(de.all, TPM_P < 5E-03)
#de.all$TPM_WEIGHT <- de.all$TPM_LOG2_FC
#de.all$WGS_WEIGHT <- de.all$WGS_LOG2_FC
de.tpm.wgs.gene.con <- subset(de.all, KEEP > 0)
#consistents <- rownames(de.tpm.wgs.gene.con)
de.tpm.wgs.gene.not <- subset(de.all, KEEP < 0)
#inconsistents <- rownames(de.tpm.wgs.gene.not)

###
## Consistents
de <- subset(de.tpm.wgs.gene.con, TPM_LOG2_FC > 0)
file.name <- file.path(wd.de.plots, "plot_TPM-WGS_SCLC_SRC-P1E03-P5E03_WEIGHT_UP.pdf")
main.text <- c("Differential genes in SCLC (n=70)", paste0(nrow(de), " up-expressed and up-replicated"))
pdf(paste0(file.name, ".pdf"), height=5, width=5)
plot(de$WGS_WEIGHT ~ de$TPM_WEIGHT, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=adjustcolor("red", alpha.f=0.3), xlim=c(0, max(de.all$TPM_WEIGHT)), ylim=c(0, -min(de.all$WGS_WEIGHT)))
mtext(main.text[2], cex=1.2, line=0.3)

cor <- cor.test(de$WGS_WEIGHT, de$TPM_WEIGHT, method="spearman", exact=F)
lm.fit <- lm(de$WGS_WEIGHT ~ de$TPM_WEIGHT)
abline(lm.fit, col="red", lwd=3)
legend("topright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col="black", bty="n", cex=1.1)
#legend("topright", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col="black", bty="n", cex=1.1)
dev.off()

##
de <- subset(de.tpm.wgs.gene.con, TPM_LOG2_FC < 0)
file.name <- file.path(wd.de.plots, "plot_TPM-WGS_SCLC_SRC-P1E03-P5E03_WEIGHT_DOWN.pdf")
main.text <- c("Differential genes in SCLC (n=70)", paste0(nrow(de), " down-expressed and down-replicated"))
pdf(paste0(file.name, ".pdf"), height=5, width=5)
plot(de$WGS_WEIGHT ~ de$TPM_WEIGHT, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=adjustcolor("blue", alpha.f=0.3), xlim=c(-max(de.all$TPM_WEIGHT), 0), ylim=c(min(de.all$WGS_WEIGHT), 0))
mtext(main.text[2], cex=1.2, line=0.3)

cor <- cor.test(de$WGS_WEIGHT, de$TPM_WEIGHT, method="spearman", exact=F)
lm.fit <- lm(de$WGS_WEIGHT ~ de$TPM_WEIGHT)
abline(lm.fit, col="blue", lwd=3)
legend("bottomleft", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col="black", bty="n", cex=1.1)
#legend("bottomleft", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col="black", bty="n", cex=1.1)
dev.off()

###
## Inconsistents
de <- subset(de.tpm.wgs.gene.not, TPM_LOG2_FC > 0)
file.name <- file.path(wd.de.plots, "plot_TPM-WGS_SCLC_SRC-P1E03-P5E03_WEIGHT_UP-DOWN.pdf")
main.text <- c("Differential genes in SCLC (n=70)", paste0(nrow(de), " up-expressed and down-replicated"))
pdf(paste0(file.name, ".pdf"), height=5, width=5)
plot(de$WGS_WEIGHT ~ de$TPM_WEIGHT, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=adjustcolor("purple", alpha.f=0.3), xlim=c(0, max(de.all$TPM_WEIGHT)), ylim=c(min(de.all$WGS_WEIGHT), 0))
mtext(main.text[2], cex=1.2, line=0.3)

cor <- cor.test(de$WGS_WEIGHT, de$TPM_WEIGHT, method="spearman", exact=F)
lm.fit <- lm(de$WGS_WEIGHT ~ de$TPM_WEIGHT)
abline(lm.fit, col="purple", lwd=3)
legend("bottomright", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col="black", bty="n", cex=1.1)
#legend("bottomright", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col="black", bty="n", cex=1.1)
dev.off()

##
de <- subset(de.tpm.wgs.gene.not, TPM_LOG2_FC < 0)
file.name <- file.path(wd.de.plots, "plot_TPM-WGS_SCLC_SRC-P1E03-P5E03_WEIGHT_DOWN-UP.pdf")
main.text <- c("Differential genes in SCLC (n=70)", paste0(nrow(de), " down-expressed and up-replicated"))
pdf(paste0(file.name, ".pdf"), height=5, width=5)
plot(de$WGS_WEIGHT ~ de$TPM_WEIGHT, ylab=ylab.text, xlab=xlab.text, main=main.text[1], col=adjustcolor("purple", alpha.f=0.3), xlim=c(-max(de.all$TPM_WEIGHT), 0), ylim=c(0, -min(de.all$WGS_WEIGHT)))
mtext(main.text[2], cex=1.2, line=0.3)

cor <- cor.test(de$WGS_WEIGHT, de$TPM_WEIGHT, method="spearman", exact=F)
lm.fit <- lm(de$WGS_WEIGHT ~ de$TPM_WEIGHT)
abline(lm.fit, col="purple", lwd=3)
legend("topleft", c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]]))), text.col="black", bty="n", cex=1.1)
#legend("topleft", c(paste0("R^2 = ", round0(summary(lm.fit)$r.squared, digits=2)), paste0("p-value = ", scientific(summary(lm.fit)$coefficients[2, 4]))), text.col="black", bty="n", cex=1.1)
dev.off()










#de.wgs.tpm.gene.out <- subset(de.wgs.tpm.gene.out, P <= 1E-03)

consistents <- intersect(rownames(de.tpm.wgs.gene.out), rownames(subset(de.tpm.gene, P <= 1E-03)))
overlaps <- intersect(overlaps, rownames(subset(de.wgs.gene, P <= 1E-06)))
length(overlaps)
# [1] 20

de.wgs.tpm.gene.out <- de.wgs.gene[overlaps,]
colnames(de.wgs.tpm.gene.out)[8:12] <- paste0("WGS_", colnames(de.wgs.tpm.gene.out)[8:12])
de.wgs.tpm.gene.out <- cbind(de.wgs.tpm.gene.out, de.tpm.gene[overlaps, 8:12])
colnames(de.wgs.tpm.gene.out)[13:17] <- paste0("RNA_", colnames(de.wgs.tpm.gene.out)[13:17])

save(de.wgs.tpm.gene.out, samples, file=file.path(wd.de.data, "de_sclc_wgs-tpm-gene-out_rt_p1e-06_p1e-02_n70.RData"))
writeTable(de.wgs.tpm.gene.out, file.path(wd.de.data, "de_sclc_wgs-tpm-gene-out_rt_p1e-06_p1e-02_n70.txt"), colnames=T, rownames=F, sep="\t")

###
## PCA of SCLC
#overlaps <- rownames(de.wgs.tpm.gene.out)
tpm.gene <- tpm.gene.sclc[overlaps, samples.sclc$SAMPLE_ID]

test <- tpm.gene
pca.de <- getPCA(t(test))

##
BASE <- "SCLC"
base <- "sclc"
traits <- samples.sclc$PCA
cols <- c("blue", "gray", "red")

file.main <- paste0(BASE, " on ", length(overlaps), " DRD and DE genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base, "_DRD+DE_20genes"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base, "_DRD+DE_20genes"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base, "_DRD+DE_20genes"), size=6.5, file.main, "topleft", cols, NA, 1, 1)

###
## PCA of combine HeLa and SCLC
#samples.hela$PCA <- samples.hela$CYCLE_PCA_2
#samples.sclc$PCA <- paste0("SCLC (", samples.sclc$G1S, ")")
samples <- rbind(samples.hela[,c("SAMPLE_ID", "PCA")], samples.sclc[,c("SAMPLE_ID", "PCA")])

overlaps <- intersect(overlaps, rownames(tpm.gene.hela))
# > setdiff(overlaps, rownames(tpm.gene))
# [1] "ENSG00000164114"
# > ensGene["ENSG00000164114",]
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name
# ENSG00000164114 ENSG00000164114            chr4     -1      156263810    156298122 protein_coding               MAP9

#overlaps <- intersect(rownames(de.wgs.tpm.gene.out), rownames(tpm.gene.hela))
tpm.gene <- cbind(tpm.gene.hela[overlaps,], tpm.gene.sclc[overlaps,])
tpm.gene <- tpm.gene[, samples$SAMPLE_ID]

test <- tpm.gene
pca.de <- getPCA(t(test))

##
BASE <- "SCLC and HeLa"
base <- "sclc+hela"
traits <- samples$PCA
cols <- c("purple3", "forestgreen", "gold", "blue", "gray", "red")

file.main <- paste0(BASE, " on ", nrow(de.wgs.tpm.gene.out), " DRD and DE genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base, "_DRD+DE_19genes"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base, "_DRD+DE_19genes"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base, "_DRD+DE_19genes"), size=6.5, file.main, "topleft", cols, NA, 1, 1)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
plotBox <- function(gene, wd.de.plots, tpm.gene.log2, de.tpm.gene, samples, type, ylim=NULL) {
   ids <- subset(ensGene, external_gene_name == gene)   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
   chr <- ids$chromosome_name[1]
   
   for (i in 1:nrow(ids)) {
      id <- ids$ensembl_gene_id[i]

      ##
      file.name <- paste0("boxplot_", tolower(type), ".gene.log2_", gene)
      if (nrow(ids) != 1)
         file.name <- paste0(file.name, "_", id)
      
      gene.tpm <- cbind(t(tpm.gene.log2[id, rownames(samples.o)]), samples)
      colnames(gene.tpm)[1] <- "MEDIAN"
  
      pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.5)
      boxplot(gene.tpm$MEDIAN ~ gene.tpm[, chr], outline=T, names=c("CM1", "CM2"), ylab=paste0("log2(", type, "+0.01)"), main=paste0(gene, " (", id, ")"))
      
      mtext(paste0("p-value = ", scientific(de.tpm.gene[id,]$P)), cex=1.2, line=0.3)
      dev.off()
   }
}

plotScatter <- function(gene, wd.de.plots, tpm.gene.log2, wgs.gene.log2, samples, ylim=NULL) {
   ids <- subset(ensGene, external_gene_name == gene)   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
   chr <- ids$chromosome_name[1]
 
   for (i in 1:nrow(ids)) {
      id <- ids$ensembl_gene_id[i]
  
      ##
      file.name <- paste0("plot_tpm.wgs.gene.log2_", gene)
      if (nrow(ids) != 1)
         file.name <- paste0(file.name, "_", id)
  
      gene.tpm <- cbind(t(tpm.gene.log2[id, rownames(samples)]), samples)
      gene.wgs <- cbind(t(wgs.gene.log2[id, rownames(samples)]), samples)
      colnames(gene.tpm)[1] <- "MEDIAN"
      colnames(gene.wgs)[1] <- "MEDIAN"
      
      pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=5, width=5)
      plot(gene.wgs$MEDIAN ~ gene.tpm$MEDIAN, xlab="log2(TPM+0.01)", ylab="log2(WGS+0.01)", main=paste0(gene, " (", id, ")"))

      cor <- cor.test(gene.wgs$MEDIAN, gene.tpm$MEDIAN, method="spearman", exact=F)
      lm.fit <- lm(gene.wgs$MEDIAN ~ gene.tpm$MEDIAN)
      abline(lm.fit, col="black", lwd=3)
      mtext(paste0("rho = ", round0(cor[[4]], digits=2), "; p-value = ", scientific(scientific(cor[[3]]))), cex=1.2, line=0.3)
      dev.off()
   }
}

##
genes <- c("RP11-141C7.3", "BRD9", "MAP9", "RAD9A", "POLE", "MCM10", "DNA2", "CCAR1", "SMARCD1", "E2F3", "ERCC8")
genes <- c("NRG1", "GULP1", "NIPSNAP3A", "CCAR1", "SMARCD1")
genes <- c("RP11-141C7.3", "AGBL4")

for (g in 1:length(genes)) {
   plotBox(genes[g], wd.de.plots, tpm.gene.log2, de.tpm.gene, samples, "TPM")
   plotBox(genes[g], wd.de.plots, wgs.gene.log2, de.wgs.gene, samples, "WGS")
   plotScatter(genes[g], file.path(wd.de.plots, "boxplots"), tpm.gene.log2, wgs.gene.log2, samples)
}






# -----------------------------------------------------------------------------
# Cannoli/Kilonova plot
# -----------------------------------------------------------------------------
plotCannoli <- function(de.1, de.2, de.sig, genes, file.de, file.main, xlim=NA, ylim=NA) {
   overlaps <- intersect(rownames(de.1), rownames(de.2))
   de.1 <- de.1[overlaps,]
   de.2 <- de.2[overlaps,]
   
   de <- de.1
   de$Effect1 <- de$LOG2_FC
   de$Effect2 <- de.2[overlaps,]$LOG2_FC
   de$Effect1 <- -log10(de.1$P) * de$Effect1
   de$Effect2 <- -log10(de.2$P) * de$Effect2
   de.positive <- subset(subset(subset(de, Effect1 > 0), Effect2 > 0), P < 1E-06)
   de.negative <- subset(subset(subset(de, Effect1 < 0), Effect2 < 0), P < 1E-06)
   de.sig <- de[rownames(de.sig),]
   
   xlab.text <- "Differential replication timing (CM2/CM1)"
   ylab.text <- "Differential gene expression (CM2/CM1)"
   pdf(file.de, height=7, width=7)
   if (is.na(xlim) && is.na(ylim)) {
      plot(de$Effect1, de$Effect2, pch=16, xlab=xlab.text, ylab=ylab.text, col=adjustcolor.gray, main=file.main[1])
   } else
      plot(de$Effect1, de$Effect2, pch=16, xlab=xlab.text, ylab=ylab.text, col=adjustcolor.gray, main=file.main[1], xlim=xlim, ylim=ylim)
   
   points(de.positive$Effect1, de.positive$Effect2, pch=16, col=adjustcolor.red)
   points(de.negative$Effect1, de.negative$Effect2, pch=16, col=adjustcolor.blue)
 
   for (g in 1:nrow(de.sig)) {
      de.sig.gene <- de.sig[g,]
      
      if (de.sig.gene$Effect1 > 0) {
         points(de.sig.gene$Effect1, de.sig.gene$Effect2, pch=16, col="red")
      } else
         points(de.sig.gene$Effect1, de.sig.gene$Effect2, pch=16, col="dodgerblue")
   }
   
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect1, gene$Effect2, pch=1, col="black")
       
         if (gene$Effect1 > 0)
            text(gene$Effect1, gene$Effect2, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
         else
            text(gene$Effect1, gene$Effect2, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   legend("topleft", legend=c("Early (P < 1E-06) & Up (P < 1.00E-03)", "Late  (P < 1E-06) & Down (P < 1.00E-03)", "Early (P < 1E-06) & Up", "Late  (P < 1E-06) & Down", "Weighted effect (-log10P * log2FC)"), col=c("red", "dodgerblue", adjustcolor.red, adjustcolor.blue, adjustcolor.gray), pch=19)
 
   abline(h=0, lty=5)
   abline(v=0, lty=5)
   dev.off()
}

## Cannoli plot
#overlaps <- intersect(rownames(de.wgs.gene), rownames(de.tpm.gene))
#de.cannoli <- de.wgs.gene[overlaps,]
#de.cannoli$Effect1 <- de.cannoli$LOG2_FC
#de.cannoli$Effect2 <- de.tpm.gene[overlaps,]$LOG2_FC

plot.main <- "Consensual differential replication and gene expression"
plot.de <- file.path(wd.de.plots, "cannoliplot_weight_sclc_")

###
## Immunoglobulin class switching
#genes <- c("DDX55", "KNTC1", "NSUN2", "AL359195.1", "EVPL", "PHKG2", "RCE1", "GTPBP4", "PAPD7", "IGHMBP2", "ZNF598", "GTF3C2", "RBM6", "AARSD1", "KAT2A")
genes <- c("BRD9", "SMARCD1", "MCM10", "RAD9A", "MAP9", "RCHY1")
file.main <- c(plot.main, "Cell cycle")

file.de <- paste0(plot.de, "Cycle.pdf")
plotCannoli(de.wgs.gene, de.tpm.gene, de.wgs.tpm.gene.out, genes, file.de, file.main)

file.de <- paste0(plot.de, "Cycle_xlim-ylim5.pdf")
plotCannoli(de.wgs.gene, de.tpm.gene, de.wgs.tpm.gene.out, genes, file.de, file.main, xlim=c(-2, 2), ylim=c(-5, 5))





## 
genes <- ""
file.main <- c("Consensual fold change (log2FC CM2/CM1)", "", "", "")
file.de <- paste0(wd.de, "cannoliplot_sclc_de-wgs-tpm_121genes.pdf")
plotKilonova(de1, de2, genes, file.de, file.main)






# -----------------------------------------------------------------------------
# Test 2
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-wgs-de/data/de_sclc_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101.RData")
de.tpm.gene.wgs <- de.tpm.gene
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-tpm-de/data/de_sclc_tpm-gene-r5p47_rt_wilcox_q_n70.RData")
de.tpm.gene <- de.tpm.gene[rownames(de.tpm.gene.wgs),]

idx <- which(de.tpm.gene.wgs$LOG2_FC * de.tpm.gene$LOG2_FC > 0)
de.tpm.gene.wgs.rna <- de.tpm.gene.wgs[idx,]

file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101-n70")
save(de.tpm.gene.wgs.rna, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.wgs.rna, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
writeRNKformat(de.tpm.gene.wgs.rna, wd.de.gsea, file.name)

de.tpm.gene.rna.wgs <- de.tpm.gene[rownames(de.tpm.gene.wgs.rna),]
de.tpm.gene.rna.wgs <- de.tpm.gene.rna.wgs[order(de.tpm.gene.rna.wgs$P),]

file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n70-n101")
writeTable(de.tpm.gene.rna.wgs, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
writeRNKformat(de.tpm.gene.rna.wgs, wd.de.gsea, file.name)

# -----------------------------------------------------------------------------
# Meta-analysis
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene.rna.wgs), rownames(de.tpm.gene.rna.wgs.nbl))
de.tpm.gene.rna.wgs.o      <- de.tpm.gene.rna.wgs[overlaps,]
de.tpm.gene.rna.wgs.nbl.o  <- de.tpm.gene.rna.wgs.nbl[overlaps,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 4882

idx <- which(de.tpm.gene.rna.wgs.o$LOG2_FC * de.tpm.gene.rna.wgs.nbl.o$LOG2_FC > 0)
de.tpm.gene.rna.wgs.o     <- de.tpm.gene.rna.wgs.o[idx,]
de.tpm.gene.rna.wgs.nbl.o <- de.tpm.gene.rna.wgs.nbl.o[idx,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 3659

de.tpm.gene.rna.wgs.o$FISHERS_P <- fishers(de.tpm.gene.rna.wgs.o$P, de.tpm.gene.rna.wgs.nbl.o$P)
library(qvalue)
de.tpm.gene.rna.wgs.o$FISHERS_FDR <- qvalue(de.tpm.gene.rna.wgs.o$FISHERS_P)$qvalue

de.tpm.gene.rna.wgs.o <- de.tpm.gene.rna.wgs.o[order(de.tpm.gene.rna.wgs.o$FISHERS_P),]
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101-n70_fishers_nbl")
save(de.tpm.gene.rna.wgs.o, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.rna.wgs.o, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Meta-analysis
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene.rna.wgs), rownames(de.tpm.gene.rna.wgs.nbl))
de.tpm.gene.rna.wgs.o      <- de.tpm.gene.rna.wgs[overlaps,]
de.tpm.gene.rna.wgs.nbl.o  <- de.tpm.gene.rna.wgs.nbl[overlaps,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 4882

idx <- which(de.tpm.gene.rna.wgs.o$LOG2_FC * de.tpm.gene.rna.wgs.nbl.o$LOG2_FC > 0)
de.tpm.gene.rna.wgs.o     <- de.tpm.gene.rna.wgs.o[idx,]
de.tpm.gene.rna.wgs.nbl.o <- de.tpm.gene.rna.wgs.nbl.o[idx,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 3659

de.tpm.gene.rna.wgs.o$FISHERS_P <- fishers(de.tpm.gene.rna.wgs.o$P, de.tpm.gene.rna.wgs.nbl.o$P)
library(qvalue)
de.tpm.gene.rna.wgs.o$FISHERS_FDR <- qvalue(de.tpm.gene.rna.wgs.o$FISHERS_P)$qvalue

de.tpm.gene.rna.wgs.o <- de.tpm.gene.rna.wgs.o[order(de.tpm.gene.rna.wgs.o$FISHERS_P),]
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101-n70_fishers_nbl")
save(de.tpm.gene.rna.wgs.o, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.rna.wgs.o, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Meta-analysis
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene.rna.wgs.o), rownames(de.tpm.gene.rna.wgs.cll))
de.tpm.gene.rna.wgs.o      <- de.tpm.gene.rna.wgs.o[overlaps,]
de.tpm.gene.rna.wgs.cll.o  <- de.tpm.gene.rna.wgs.cll[overlaps,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 4882

idx <- which(de.tpm.gene.rna.wgs.o$LOG2_FC * de.tpm.gene.rna.wgs.cll.o$LOG2_FC > 0)
de.tpm.gene.rna.wgs.o     <- de.tpm.gene.rna.wgs.o[idx,]
de.tpm.gene.rna.wgs.cll.o <- de.tpm.gene.rna.wgs.cll.o[idx,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 1101

de.tpm.gene.rna.wgs.o$FISHERS_P3 <- fishers(de.tpm.gene.rna.wgs.o$FISHERS_P, de.tpm.gene.rna.wgs.cll.o$P)
library(qvalue)
de.tpm.gene.rna.wgs.o$FISHERS_FDR3 <- qvalue(de.tpm.gene.rna.wgs.o$FISHERS_P3)$qvalue

de.tpm.gene.rna.wgs.o <- de.tpm.gene.rna.wgs.o[order(de.tpm.gene.rna.wgs.o$FISHERS_P3),]
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101-n70-n71_fishers_nbl+cll")
save(de.tpm.gene.rna.wgs.o, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.rna.wgs.o, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Meta-analysis on wgs.de first
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-wgs-de/data/de_sclc_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101.RData")
de.tpm.gene.wgs <- de.tpm.gene

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-wgs-de/data/de_nbl_tpm-gene-r30p47-r5p47_rt_wilcox_q_n56.RData")
de.tpm.gene.wgs.nbl <- de.tpm.gene

load("/Users/tpyang/Work/uni-koeln/tyang2/CLL/analysis/expression/kallisto/cll-wgs-de/data/de_cll_tpm-gene-r30p47-r5p47_rt_wilcox_q_n93.RData")
de.tpm.gene.wgs.cll <- de.tpm.gene

overlaps <- intersect(intersect(rownames(de.tpm.gene.wgs), rownames(de.tpm.gene.wgs.nbl)), rownames(de.tpm.gene.wgs.cll))
de.tpm.gene.wgs.o      <- de.tpm.gene.wgs[overlaps,]
de.tpm.gene.wgs.nbl.o  <- de.tpm.gene.wgs.nbl[overlaps,]
de.tpm.gene.wgs.cll.o  <- de.tpm.gene.wgs.cll[overlaps,]
# > length(overlaps)
# [1] 15939

idx1 <- which(de.tpm.gene.wgs.o$LOG2_FC * de.tpm.gene.wgs.nbl.o$LOG2_FC > 0)
idx2 <- which(de.tpm.gene.wgs.o$LOG2_FC * de.tpm.gene.wgs.cll.o$LOG2_FC > 0)
idx <- intersect(idx1, idx2)
de.tpm.gene.wgs.o     <- de.tpm.gene.wgs.o[idx,]
de.tpm.gene.wgs.nbl.o <- de.tpm.gene.wgs.nbl.o[idx,]
de.tpm.gene.wgs.cll.o <- de.tpm.gene.wgs.cll.o[idx,]
# > length(idx)
# [1] 9982

library(qvalue)
de.tpm.gene.wgs.o$FISHERS_P2 <- fishers(de.tpm.gene.wgs.o$P, de.tpm.gene.wgs.nbl.o$P)
de.tpm.gene.wgs.o$FISHERS_FDR2 <- qvalue(de.tpm.gene.wgs.o$FISHERS_P2)$qvalue
de.tpm.gene.wgs.o$FISHERS_P3 <- fishers(de.tpm.gene.wgs.o$FISHERS_P2, de.tpm.gene.wgs.cll.o$P)
de.tpm.gene.wgs.o$FISHERS_FDR3<- qvalue(de.tpm.gene.wgs.o$FISHERS_P3)$qvalue

de.tpm.gene.wgs.o <- de.tpm.gene.wgs.o[order(de.tpm.gene.wgs.o$FISHERS_P3),]
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_fishers_sclc+nbl+cll")
save(de.tpm.gene.wgs.o, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.wgs.o, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Meta-analysis
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene.rna.wgs.o), rownames(de.tpm.gene.rna.wgs.cll))
de.tpm.gene.rna.wgs.o      <- de.tpm.gene.rna.wgs.o[overlaps,]
de.tpm.gene.rna.wgs.cll.o  <- de.tpm.gene.rna.wgs.cll[overlaps,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 4882

idx <- which(de.tpm.gene.rna.wgs.o$LOG2_FC * de.tpm.gene.rna.wgs.cll.o$LOG2_FC > 0)
de.tpm.gene.rna.wgs.o     <- de.tpm.gene.rna.wgs.o[idx,]
de.tpm.gene.rna.wgs.cll.o <- de.tpm.gene.rna.wgs.cll.o[idx,]
# > nrow(de.tpm.gene.rna.wgs.o)
# [1] 1101

de.tpm.gene.rna.wgs.o$FISHERS_P3 <- fishers(de.tpm.gene.rna.wgs.o$FISHERS_P, de.tpm.gene.rna.wgs.cll.o$P)
library(qvalue)
de.tpm.gene.rna.wgs.o$FISHERS_FDR3 <- qvalue(de.tpm.gene.rna.wgs.o$FISHERS_P3)$qvalue

de.tpm.gene.rna.wgs.o <- de.tpm.gene.rna.wgs.o[order(de.tpm.gene.rna.wgs.o$FISHERS_P3),]
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101-n70-n71_fishers_nbl+cll")
save(de.tpm.gene.rna.wgs.o, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.rna.wgs.o, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")












# -----------------------------------------------------------------------------
# Principal component analysis (PCA) of LCNEC samples on RB1-loss DE genes
# Figure(s)    : Figure 1 (B)
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
#load(file.path("/Users/tpyang/Work/uni-koeln/tyang2", "LCNEC", "analysis/expression/kallisto", paste0("lcnec", "-tpm-de/data/", "de_", "lcnec", "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rt.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
# > length(genes.rb1.q0.05)
# [1] 145

##
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- as.numeric(samples[,"RB1_MUT"])
trait[which(trait == 0)] <- "LCNEC (WT)"
trait[which(trait == 1)] <- "LCNEC (RB1)"
trait[which(is.na(trait))] <- "LCNEC (NA)"

file.main <- c("LCNEC samples on top 145 DE genes", "RB1-loss differential effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcnec_rb1_q0.05_145de", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)


# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
nrow(ensGene[intersect(genes.rb1.q0.05, core.G1S),])   ## 10/41 = 23.4%
nrow(ensGene[intersect(genes.rb1.q0.05, core.G2M),])   ##  4/54 = 7.4%
nrow(ensGene[intersect(genes.rb1.q0.05, core.SC),])    ##  0/54 = 0%

nrow(ensGene[intersect(genes.rb1.q0.05, genes.G1S),])   ## 15
nrow(ensGene[intersect(genes.rb1.q0.05, genes.G2M),])   ##  9

nrow(ensGene[intersect(genes.rb1.q0.05, unique(c(core.G1S, genes.G1S))),])   ## 19
nrow(ensGene[intersect(genes.rb1.q0.05, unique(c(core.G2M, genes.G2M))),])   ## 11

## Output
genes.rb1.n.s   <- rownames(subset(de.tpm.gene, FDR > 0.5))
genes.rb1.q0.5  <- rownames(subset(de.tpm.gene, FDR <= 0.5))
file.name <- paste0("de_", base, "_tpm-gene-r5p47_rb1_wilcox_q_n54_n.s.0.5")
#genes.rb1 <- genes.rb1.n.s
#genes.rb1 <- genes.rb1.q0.5
genes.rb1 <- genes.rb1.q0.1

output <- de.tpm.gene[intersect(genes.rb1, unique(c(core.G1S, genes.G1S))), c(1,2,8,9,12)]
output$Tirosh_2016 <- ""
output$Dominguez_2016 <- ""
output[intersect(genes.rb1, core.G1S),]$Tirosh_2016 <- "yes"
output[intersect(genes.rb1, genes.G1S),]$Dominguez_2016 <- "yes"
writeTable(output, file.path(wd.de.data, paste0(file.name, "_G1-S.txt")), colnames=T, rownames=F, sep="\t")

output <- de.tpm.gene[intersect(genes.rb1, unique(c(core.G2M, genes.G2M))), c(1,2,8,9,12)]
output$Tirosh_2016 <- ""
output$Dominguez_2016 <- ""
output[intersect(genes.rb1, core.G2M),]$Tirosh_2016 <- "yes"
output[intersect(genes.rb1, genes.G2M),]$Dominguez_2016 <- "yes"
writeTable(output, file.path(wd.de.data, paste0(file.name, "_G2-M.txt")), colnames=T, rownames=F, sep="\t")
