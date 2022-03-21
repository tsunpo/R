# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/expression/icgc-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 12/03/22
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R", "TranscriptionReplicationInteraction.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
#load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/RNA")
#wd.rna.raw <- file.path(wd.rna, "kallisto_hg19.ensembl_quant-b100--bias")

wd.anlys <- file.path(wd, BASE, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

load(file.path(wd.de.data, paste0(base, "_kallisto_0.42.1_tpm.gene.median0.RData")))
dim(tpm.gene)
# [1] 20720  1359
# [1] 18502  1359

icgc.pos <- subset(icgc, survival_rho > 0)
samples.surv.pos <- subset(samples.surv, histology_abbreviation %in% rownames(icgc.pos))
nrow(samples.surv.pos)
# [1] 905

overlaps <- intersect(samples.surv.pos$icgc_specimen_id, colnames(tpm.gene))
length(overlaps)
# [1] 236
samples.tpm <- samples.pos[overlaps,]
tpm.gene <- tpm.gene[, overlaps]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2   <- log2(tpm.gene + 1)
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 1)
tpm.gene.log2 <- tpm.gene.log2[rownames(subset(tpm.gene.log2.m, MEDIAN != 0)),]
nrow(tpm.gene.log2)
# [1] 

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric)
# Last Modified: 12/03/22
# -----------------------------------------------------------------------------
## Test: Wilcoxon
## FDR : Q
argv      <- data.frame(predictor="SG1", predictor.wt="G1", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0_SG1_wilcox_q_RHO>0_n236")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples.tpm, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
# Samples with MUT SG1: 62
# Samples with  WT SG1: 174

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
nrow(de.tpm.gene)
# [1] 18432

pos.bad <- rownames(subset(subset(de.tpm.gene, P <= 1E-6), LOG2_FC < 0))
pos.good <- rownames(subset(subset(de.tpm.gene, P <= 1E-6), LOG2_FC > 0))

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=45, 22 TR vs 23 UN)
# Last Modified: 12/03/22
# -----------------------------------------------------------------------------
## Test: Wilcoxon
## FDR : Q
argv      <- data.frame(predictor="donor_sex", predictor.wt="male", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0_SEX_wilcox_q_RHO>0_n236")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples.tpm, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
# Samples with MUT donor_sex: 160
# Samples with  WT donor_sex: 76

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)
de.tpm.gene <- subset(de.tpm.gene, chromosome_name %in% chrs[-c(23, 24)])
de.tpm.gene$FDR <- qvalue(de.tpm.gene$P)$qvalue

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
nrow(de.tpm.gene)
# [1] 17751

pos.bad <- rownames(subset(subset(de.tpm.gene, P <= 1E-9), LOG2_FC > 0))
pos.good <- rownames(subset(subset(de.tpm.gene, P <= 1E-9), LOG2_FC < 0))




plotVolcano <- function(de, pvalue, genes, file.de, file.main, xlab.text, ymax=0) {
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   xmin <- min(de$LOG2_FC)
   if (ymax ==0) ymax <- max(de$log10P)
   #ymax <- 7
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(xmin, xmax), ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col="lightgray", main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   #abline(v=c(-log2(1.5), log2(1.5)), lty=5, col="darkgray")
 
   abline(h=c(-log10(pvalue)), lty=5)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red.lighter, cex=1.4)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=blue.lighter, cex=1.4)
 
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
            print(genes[g,])
      }
   }
 
   axis(side=1, at=seq(-12, 12, by=4), labels=c(-12, -8, -4, 0, 4, 8, 12), cex.axis=1.2)
   mtext(file.main[2], cex=1.25, line=0.3)
   legend("topleft", legend=c("Up-regulated (L < H)", "Down-regulated (L > H)"), col=c(red.lightest, blue.lightest), pch=19, cex=1.25)
   dev.off()
}

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE_p1e-3_TERT")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL risk group-associated genes", "Low (L) vs. High (H)")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# Last Modified: 12/08/21
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_median0")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_Risk_p1e-3_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_Risk_p1e-3_down")

setTable <- function(wd.de.reactome) {
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
}

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE_p1e-3_up")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 1e-7)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated (L < H)", "309 genes")
file.de <- file.path(wd.de.reactome, "genes_DE_p1e-3_n309_up.pdf")

pdf(file.de, height=3.7, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 12), xaxt="n", names.arg=reactome.up$Pathway.name, col=red.lightest, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 12, by=4))
mtext(main.text[2], line=0.3)
abline(v=7, lty=5)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE_p1e-3_down")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1.7e-3)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Down-regulated (L > H)", "727 genes")
file.de <- file.path(wd.de.reactome, "genes_DE_p1e-3_n727_down.pdf")

pdf(file.de, height=3, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col=blue.lightest, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-7, lty=5)

axis(side=1, at=seq(-12, 0, by=4))
mtext(main.text[2], line=0.3)
dev.off()








# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/03/22
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "G1", "S", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(src) <- rownames(tpm.gene.log2)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.tpm$COR, method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples.tpm$COR, method="spearman", exact=F)[[3]])
#src <- src[!is.na(src$P),]

## Log2 fold change
src$G1 <- median00(tpm.gene.log2, rownames(subset(samples.tpm, SG1 == "G1")))
src$S <- median00(tpm.gene.log2, rownames(subset(samples.tpm, SG1 == "S")))
src$LOG2_FC <- src$S - src$G1

## FDR
library(qvalue)
src$Q <- qvalue(src$P)$qvalue
src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, "icgc_tpm-gene-median0_src_q_RHO<0_n368.txt"), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples.tpm, file=file.path(wd.de.data, "icgc_tpm-gene-median0_src_q_RHO<0_n368.RData"))
nrow(src.tpm.gene)
# [1] 










# -----------------------------------------------------------------------------
# Correlation bwteen TPM and in-silico sorting
# Last Modified: 12/12/19; 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "M1", "M2", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
src <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(src) <- rownames(tpm.gene.log2)

## SRC
src$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[4]])
src$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[3]])
src <- src[!is.na(src$P),]

## Log2 fold change
src$M1 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 0)))
src$M2 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 1)))
src$LOG2_FC <- src$M2 - src$M1

## FDR
library(qvalue)
src$Q <- qvalue(src$P)$qvalue
src <- src[order(src$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
src.tpm.gene <- cbind(annot[rownames(src),], src)   ## BE EXTRA CAREFUL!!

writeTable(src.tpm.gene, file.path(wd.de.data, "src_nbl_tpm-gene-median0+1_src_q_n53.txt"), colnames=T, rownames=F, sep="\t")
save(src.tpm.gene, samples, file=file.path(wd.de.data, "src_nbl_tpm-gene-median0+1_src_q_n53.RData"))
nrow(src.tpm.gene)
# [1] 22899
# [1] 18764

# -----------------------------------------------------------------------------
# Last Modified: 15/08/21
# -----------------------------------------------------------------------------
plotVolcanoFDR <- function(de, fdr, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e.cd, legend, col) {
   pvalue <- fdrToP(fdr/100, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
   adjustcolor.gray  <- adjustcolor("lightgray", alpha.f=0.125)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   xmin <- min(de$LOG2_FC)
   ymax <- max(de$log10P)

   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(xmin, xmax), ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col=adjustcolor.gray, main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   #abline(v=c(-log2(1.5), log2(1.5)), lty=5, col="darkgray")
 
   abline(h=c(-log10(pvalue)), lty=5)
   text(xmin + (xmax-xmin)/13, -log10(pvalue) + ymax/26, paste0("FDR = ", fdr, "%"), cex=1.2)
   
   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red.lightest, cex=1.4)
 
   if (is.null(nrow(overlaps.up.pe3))) {
      de.up <- de[overlaps.up.pe3,]
      points(de.up$LOG2_FC, de.up$log10P, pch=16, col=flowjo.red, cex=1.4)
   }
   if (is.null(nrow(overlaps.up.pe3.iz.e.cd))) {
      de.up <- de[overlaps.up.pe3.iz.e.cd,]
      points(de.up$LOG2_FC, de.up$log10P, pch=16, col="red", cex=1.4)
   }
 
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
            print(genes[g,])
      }
   }
 
   axis(side=1, at=seq(-4, 4, by=2), labels=c(-4, -2, 0, 2, 4), cex.axis=1.2)
   mtext(file.main[2], cex=1.25, line=0.3)
   #if (!is.null(nrow(overlaps.up.pe3)) && !is.null(nrow(overlaps.up.pe3.iz.e.cd))) {
      #legend("topleft", , pch=19, cex=1.25) 
   #} else if (!is.null(nrow(overlaps.up.pe3)) && is.null(nrow(overlaps.up.pe3.iz.e.cd))) {
      #legend("topleft", legend=c("309 L < H", "158 S-like"), col=c(red.lightest, flowjo.red), pch=19, cex=1.25)  
   #} else {
      legend("topleft", legend=legend, col=col, pch=19, cex=1.25)
   #}
 
   dev.off()
}

### 17/09/21
##  DE
overlaps.up   <- intersect(rownames(subset(de.tpm.gene, LOG2_FC >= 0)), rownames(subset(src.tpm.gene, LOG2_FC >= 0)))
overlaps.down <- intersect(rownames(subset(de.tpm.gene, LOG2_FC < 0)),  rownames(subset(src.tpm.gene, LOG2_FC < 0)))
overlaps      <- c(overlaps.up, overlaps.down)

##  DE + SRC
overlaps.up.fdr1   <- intersect(rownames(subset(subset(de.tpm.gene, FDR <= 1E-2), LOG2_FC >= 0)), rownames(subset(subset(src.tpm.gene, Q <= 1E-2), LOG2_FC >= 0)))
overlaps.down.fdr1 <- intersect(rownames(subset(subset(de.tpm.gene, FDR <= 1E-2), LOG2_FC < 0)),  rownames(subset(subset(src.tpm.gene, Q <= 1E-2), LOG2_FC < 0)))
overlaps.fdr1      <- c(overlaps.up.fdr1, overlaps.down.fdr1)

de.tpm.gene.overlaps.fdr1 <- de.tpm.gene[overlaps.fdr1,]
de.tpm.gene.overlaps.fdr1 <- de.tpm.gene.overlaps.fdr1[order(de.tpm.gene.overlaps.fdr1$P),]

writeTable(overlaps.up.fdr1,   file.path(wd.de.data, "genes_DE+SRC_p1e-3_n125_up.txt"), colnames=F, rownames=F, sep="")
writeTable(overlaps.down.fdr1, file.path(wd.de.data, "genes_DE+SRC_p1e-3_n397_down.txt"), colnames=F, rownames=F, sep="")

## DE + SER + Early IZ genes
overlaps.up.fdr1.iz.e <- intersect(overlaps.up.fdr1, rownames(tpm.gene.log2.m.rfd.ctr.iz.e))
overlaps.down.fdr1.iz.e <- intersect(overlaps.down.fdr1, rownames(tpm.gene.log2.m.rfd.ctr.iz.e))

###
## 3rd tier
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE+SRC+IZ+E_fdr1_ICT1")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high-risk, proliferative, early IZ genes", "")
plotVolcanoFDR(de.tpm.gene, 1, genes, file.de, file.main, xlab.text, overlaps.up.fdr1, overlaps.up.fdr1.iz.e, c("215 Low < High", "125 Proliferative", "  20 Early IZ"), c(red.lightest, flowjo.red, "red"))

writeTable(de.tpm.gene[overlaps.up.fdr1.iz.e,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC+IZ+E_fdr1_n53_HIGH.txt"), colnames=T, rownames=F, sep="\t")
writeTable(de.tpm.gene[overlaps.down.fdr1.iz.e,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC+IZ+E_fdr1_n53_LOW.txt"), colnames=T, rownames=F, sep="\t")
save(overlaps.up.fdr1.iz.e, overlaps.down.fdr1.iz.e, file=file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC+IZ+E_fdr1_n53.RData"))

## 2nd tier
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE+SRC_fdr1_FKBP3")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high-risk and proliferative genes", "")
plotVolcanoFDR(de.tpm.gene, 1, genes, file.de, file.main, xlab.text, overlaps.up.fdr1, NA, c("215 Low < High", "125 Proliferative"), c(red.lightest, flowjo.red))

writeTable(de.tpm.gene[overlaps.up.fdr1,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC_fdr1_n53_HIGH.txt"), colnames=T, rownames=F, sep="\t")
writeTable(de.tpm.gene[overlaps.down.fdr1,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC_fdr1_n53_LOW.txt"), colnames=T, rownames=F, sep="\t")
save(overlaps.up.pe3, overlaps.down.fdr1, file=file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC_fdr1_n53.RData"))

## 1st tier
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE_fdr1_TERT")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high risk-associated genes", "")
plotVolcanoFDR(de.tpm.gene, 1, genes, file.de, file.main, xlab.text, NA, NA, "215 Low < High", red.lightest)






plotVolcano2 <- function(de, pvalue, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e.cd, legend) {
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
   adjustcolor.gray  <- adjustcolor("lightgray", alpha.f=0.125)
   
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #xmin <- -xmax
   xmin <- min(de$LOG2_FC)
   #if (ymax ==0) ymax <- max(de$log10P)
   ymax <- max(de$log10P) + 0.2
   #ymax <- 8.5
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(xmin, xmax), ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col=adjustcolor.gray, main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   #abline(v=c(-log2(1.5), log2(1.5)), lty=5, col="darkgray")
 
   abline(h=c(-log10(pvalue)), lty=5)
 
   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red.lightest, cex=1.4)
   
   if (is.null(nrow(overlaps.up.pe3))) {
      de.up <- de[overlaps.up.pe3,]
      points(de.up$LOG2_FC, de.up$log10P, pch=16, col=flowjo.red, cex=1.4)
   }
   if (is.null(nrow(overlaps.up.pe3.iz.e.cd))) {
      de.up <- de[overlaps.up.pe3.iz.e.cd,]
      points(de.up$LOG2_FC, de.up$log10P, pch=16, col="red", cex=1.4)
   }
 
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
            print(genes[g,])
      }
   }
 
   axis(side=1, at=seq(-4, 4, by=2), labels=c(-4, -2, 0, 2, 4), cex.axis=1.2)
   mtext(file.main[2], cex=1.25, line=0.3)
   #if (!is.null(nrow(overlaps.up.pe3)) && !is.null(nrow(overlaps.up.pe3.iz.e.cd))) {
   #legend("topleft", legend=c("309 L < H", "158 S-like", "  22 Early IZ"), col=c(red.lightest, flowjo.red, "red"), pch=19, cex=1.25) 
   #} else if (!is.null(nrow(overlaps.up.pe3)) && is.null(nrow(overlaps.up.pe3.iz.e.cd))) {
   #legend("topleft", legend=c("309 L < H", "158 S-like"), col=c(red.lightest, flowjo.red), pch=19, cex=1.25)  
      #} else {
         legend("topleft", legend=legend, col=c(red.lightest), pch=19, cex=1.25)
      #}
   
   #legend("topleft", legend=c("309 L < H", "158 S-like", "  16 Universal"), col=c(red.lightest, flowjo.red, "red"), pch=19, cex=1.25) 
   #legend("topleft", legend=c("426 LR < TERT"), col=c(red.lightest), pch=19, cex=1.25)
   #legend("topleft", legend=c("265 LR < HR"),  col=c(red.lightest), pch=19, cex=1.25)
   dev.off()
}

### 16/08/21
## NBL expressed genes 
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC+E+CD_p1e-3_TERT")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".google.red.pdf")
file.main <- c("11 NBL early IZ, CD genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e.cd)

### 30/08/21
## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE+SRC+IZ+E_p1e-3_ICT1")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high-risk, proliferative, early IZ genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e, "")

writeTable(de.tpm.gene[overlaps.up.pe3.iz.e,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC+IZ+E_pe1-3_n53_HIGH.txt"), colnames=T, rownames=F, sep="\t")
writeTable(de.tpm.gene[overlaps.down.pe3.iz.e,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC+IZ+E_pe1-3_n53_LOW.txt"), colnames=T, rownames=F, sep="\t")
save(overlaps.up.pe3.iz.e, overlaps.down.pe3.iz.e, file=file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC+IZ+E_pe1-3_n53.RData"))

### 29/08/21
## NBL expressed genes 
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE+SRC_p1e-3_FKBP3")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high-risk and proliferative genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, NA, "")

writeTable(de.tpm.gene[overlaps.up.pe3,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC_pe1-3_n53_HIGH.txt"), colnames=T, rownames=F, sep="\t")
writeTable(de.tpm.gene[overlaps.down.pe3,], file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC_pe1-3_n53_LOW.txt"), colnames=T, rownames=F, sep="\t")
save(overlaps.up.pe3, overlaps.down.pe3, file=file.path(wd.de.data, "NBL_tpm-gene-median0_DE+SRC_pe1-3_n53.RData"))

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE_p1e-3_TERT")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high risk-associated genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, NA, NA, "309 L < H")





## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC+E+CD_p1e-3_22q11.22")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL early IZ, CD genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e.cd)





###
##
overlaps.up <- intersect(rownames(subset(de.tpm.gene, LOG2_FC >= 0)), rownames(subset(src.tpm.gene, LOG2_FC >= 0)))
overlaps.down <- intersect(rownames(subset(de.tpm.gene, LOG2_FC < 0)), rownames(subset(src.tpm.gene, LOG2_FC < 0)))
overlaps <- c(overlaps.up, overlaps.down)

## NBL expressed genes
#xlab.text <- "High/Low risk fold change [log2]"
#plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC_p1e-3_TERT")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.de <- paste0(plot.de, ".pdf")
#file.main <- c("NBL high risk, proliferation-associated genes", "Low (L) vs. High (H)")
#plotVolcano(de.tpm.gene[overlaps,], 0.001, genes, file.de, file.main, xlab.text)

##
overlaps.up.pe3 <- intersect(rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC >= 0)), rownames(subset(subset(src.tpm.gene, P < 1E-3), LOG2_FC >= 0)))
overlaps.down.pe3 <- intersect(rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC < 0)), rownames(subset(subset(src.tpm.gene, P < 1E-3), LOG2_FC < 0)))
overlaps.pe3 <- c(overlaps.up.pe3, overlaps.down.pe3)

de.tpm.gene.overlaps.pe3 <- de.tpm.gene[overlaps.pe3,]
de.tpm.gene.overlaps.pe3 <- de.tpm.gene.overlaps.pe3[order(de.tpm.gene.overlaps.pe3$P),]

writeTable(overlaps.up.pe3, file.path(wd.de.data, "genes_DE+SRC_p1e-3_n158_up.txt"), colnames=F, rownames=F, sep="")
writeTable(overlaps.down.pe3, file.path(wd.de.data, "genes_DE+SRC_p1e-3_n464_down.txt"), colnames=F, rownames=F, sep="")

### 30/08/21
## Overlapping with early IZ genes
overlaps.up.pe3.iz.e <- intersect(overlaps.up.pe3, rownames(tpm.gene.log2.m.rfd.ctr.iz.e))





###
## Overlaps with early CD genes
genes.iz.e.cd <- readTable(file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_e_cd_n53.txt"), header=T, rownames=T, sep="")
genes.iz.e.ho <- readTable(file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_e_ho_n53.txt"), header=T, rownames=T, sep="")
genes.iz.e <- c(genes.iz.e.cd, genes.iz.e.ho)
overlaps.iz.e.cd <- intersect(rownames(genes.iz.e.cd), overlaps)
overlaps.iz.e <- intersect(c(rownames(genes.iz.e.cd), rownames(genes.iz.e.ho)), overlaps)

de.tpm.gene.iz.e.cd <- de.tpm.gene[overlaps.iz.e.cd,]
de.tpm.gene.iz.e.cd <- de.tpm.gene.iz.e.cd[order(de.tpm.gene.iz.e.cd$P),]

de.tpm.gene.iz.e <- de.tpm.gene[overlaps.iz.e,]
de.tpm.gene.iz.e <- de.tpm.gene.iz.e[order(de.tpm.gene.iz.e$P),]

##
overlaps.up.pe3.iz.e <- intersect(rownames(subset(subset(de.tpm.gene.iz.e, P < 1E-3), LOG2_FC > 0)), overlaps.pe3)
overlaps.down.pe3.iz.e <- intersect(rownames(subset(subset(de.tpm.gene.iz.e, P < 1E-3), LOG2_FC < 0)), overlaps.pe3)

writeTable(overlaps.up.pe3.iz.e, file.path(wd.de.data, "genes_DE+SRC+E+IZ_p1e-3_n22_up.txt"), colnames=F, rownames=F, sep="")
writeTable(overlaps.down.pe3.iz.e, file.path(wd.de.data, "genes_DE+SRC+E+IZ_p1e-3_n67_down.txt"), colnames=F, rownames=F, sep="")

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE+SRC+E+IZ_p1e-3_ICT1")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high-risk, proliferative, early IZ genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene[overlaps.iz.e,], 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3.iz.e, overlaps.down.pe3.iz.e)

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC+E_p1e-3_22q11.22")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL early IZ genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene[overlaps.iz.e,], 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3.iz.e, overlaps.down.pe3.iz.e)







##
overlaps.up.pe3.iz.e.cd <- rowname(subset(subset(de.tpm.gene.iz.e.cd, P < 1E-3), LOG2_FC > 0))
overlaps.down.pe3.iz.e.cd <- rowname(subset(subset(de.tpm.gene.iz.e.cd, P < 1E-3), LOG2_FC < 0))

writeTable(overlaps.up.pe3.iz.e.cd, file.path(wd.de.data, "genes_DE+SRC+CD_p1e-3_n11_up.txt"), colnames=F, rownames=F, sep="")
writeTable(overlaps.down.pe3.iz.e.cd, file.path(wd.de.data, "genes_DE+SRC+CD_p1e-3_n36_down.txt"), colnames=F, rownames=F, sep="")

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC+CD_p1e-3")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL early IZ, CD genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene[overlaps.iz.e.cd,], 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e.cd)

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC+CD_p1e-3")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL early IZ, CD genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene[overlaps.iz.e.cd,], 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps.up.pe3.iz.e.cd)

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC+CD_p1e-3_22q11.22")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL early IZ, CD genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene[overlaps.iz.e.cd,], 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3.iz.e.cd, overlaps.down.pe3.iz.e.cd)

#overlaps.r5p47 <- intersect(rownames(de.tpm.gene.iz.e.cd), rownames(tpm.gene.log2))
#xlab.text <- "High/Low risk fold change [log2]"
#plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_r5p47_DE+SRC+CD_p1e-3")
#genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
#file.de <- paste0(plot.de, ".pdf")
#file.main <- c("NBL early IZ, CD genes", "Low (L) vs. High (H)")
#plotVolcano(de.tpm.gene.iz.e.cd[overlaps.r5p47,], 0.001, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=28, 17 LR vs 10 TERT)
# Last Modified: 13/09/21
# -----------------------------------------------------------------------------
samples.tert       <- rbind(subset(samples.atrx, group == 0), subset(samples.atrx, group == 2))
tpm.gene.log2.tert <- tpm.gene.log2[, rownames(samples.tert)]

## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : TERT (n=10) vs LR (n=17) as factor
argv      <- data.frame(predictor="risk", predictor.wt="low", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0+1_TERT_wilcox_q_n27")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2.tert, samples.tert, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

## NBL expressed genes
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_TERT_wilcox_q_n27.RData")

xlab.text <- "TERT/LR fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE-TERT_p1e-2")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL TERT-associated genes", "LR (n=17) vs. TERT (n=10)")
plotVolcano2(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, NA, NA, "606 LR < TERT")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=26, 17 LR vs 9 MNA)
# Last Modified: 13/09/21
# -----------------------------------------------------------------------------
samples.mna       <- rbind(subset(samples.atrx, group == 0), subset(samples.atrx, group == 3))
tpm.gene.log2.mna <- tpm.gene.log2[, rownames(samples.mna)]

## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : MNA (n=9) vs LR (n=17) as factor
argv      <- data.frame(predictor="risk", predictor.wt="low", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0+1_MNA_wilcox_q_n26")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2.mna, samples.mna, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

## NBL expressed genes
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_MNA_wilcox_q_n26.RData")

xlab.text <- "MNA/LR fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE-MNA_p1e-2")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL MNA-associated genes", "LR (n=17) vs. MNA (n=9)")
plotVolcano2(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, NA, NA, "1,524 LR < MNA")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=24, 17 LR vs 7 ATRX)
# Last Modified: 13/09/21
# -----------------------------------------------------------------------------
samples.atrx2       <- rbind(subset(samples.atrx, group == 0), subset(samples.atrx, group == 1))
tpm.gene.log2.atrx <- tpm.gene.log2[, rownames(samples.atrx2)]

## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : High (n=36) vs Low (n=17) as factor
argv      <- data.frame(predictor="risk", predictor.wt="low", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0+1_ATRX_wilcox_q_n24")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2.atrx, samples.atrx2, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

## NBL expressed genes
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_ATRX_wilcox_q_n24.RData")

xlab.text <- "ATRX/LR fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE-ATRX_p1e-2")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL ATRX-associated genes", "LR (n=17) vs. ATRX (n=7)")
plotVolcano2(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, NA, NA, "347 LR < ATRX")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=27, 17 LR vs 10 HR)
# Last Modified: 10/09/21
# -----------------------------------------------------------------------------
samples.hr       <- rbind(subset(samples.atrx, group == 0), subset(samples.atrx, group == 4))
tpm.gene.log2.hr <- tpm.gene.log2[, rownames(samples.hr)]

## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : High (n=10) vs Low (n=17) as factor
argv      <- data.frame(predictor="risk", predictor.wt="low", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-median0+1_HR_wilcox_q_n27")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2.hr, samples.hr, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
nrow(subset(de.tpm.gene, P <= 0.001))

## NBL expressed genes
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_HR_wilcox_q_n27.RData")

xlab.text <- "HR/LR fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE-HR_p1e-2")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL HR-associated genes", "LR (n=17) vs. HR (n=10)")
plotVolcano2(de.tpm.gene, 0.01, genes, file.de, file.main, xlab.text, NA, NA, "485 LR < HR")

# -----------------------------------------------------------------------------
# 
# Last Modified: 11/08/21
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_TERT_wilcox_q_n27.RData")
tert.pe3 <- rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC >= 0))

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_MNA_wilcox_q_n26.RData")
mna.pe3  <- rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC >= 0))

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_ATRX_wilcox_q_n24.RData")
atrx.pe3 <- rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC >= 0))

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_HR_wilcox_q_n27.RData")
hr.pe3 <- rownames(subset(subset(de.tpm.gene, P < 1E-3), LOG2_FC >= 0))

overlaps4.pe3 <- intersect(intersect(intersect(tert.pe3, mna.pe3), atrx.pe3), hr.pe3)

###
##
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_TERT_wilcox_q_n27.RData")
tert.pe2 <- rownames(subset(subset(de.tpm.gene, P < 1E-2), LOG2_FC >= 0))

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_MNA_wilcox_q_n26.RData")
mna.pe2  <- rownames(subset(subset(de.tpm.gene, P < 1E-2), LOG2_FC >= 0))

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_ATRX_wilcox_q_n24.RData")
atrx.pe2 <- rownames(subset(subset(de.tpm.gene, P < 1E-2), LOG2_FC >= 0))

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-median0+1_HR_wilcox_q_n27.RData")
hr.pe2 <- rownames(subset(subset(de.tpm.gene, P < 1E-2), LOG2_FC >= 0))

overlaps4.pe2 <- intersect(intersect(intersect(tert.pe2, mna.pe2), atrx.pe2), hr.pe2)
overlaps4.pe2.up <- intersect(overlaps.up.pe3, overlaps4.pe2)

de <- de.tpm.gene[overlaps4.pe2.up,]
de <- de[order(de$P),]

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0+1_DE+SRC+MNA_p1e-3_FKBP3")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high-risk, proliferative, universal genes", "Low (L) vs. High (H)")
plotVolcano2(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, overlaps.up.pe3, overlaps4.pe2.up, "")


# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 11/08/21
# -----------------------------------------------------------------------------
plotVolcano <- function(de, pvalue, genes, file.de, file.main, xlab.text, ymax=0) {
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   xmin <- min(de$LOG2_FC)
   if (ymax ==0) ymax <- max(de$log10P)
   
   pdf(file.de, height=6, width=6)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(xmin, xmax), ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab="P-value significance [-log10]", col="lightgray", main=file.main[1], cex=1.4, cex.axis=1.2, cex.lab=1.25, cex.main=1.3)
   #abline(v=c(-log2(1.5), log2(1.5)), lty=5, col="darkgray")
   
   abline(h=c(-log10(pvalue)), lty=5)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col=red, cex=1.4)
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col=blue, cex=1.4)
 
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
   
   #axis(side=1, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1), cex.axis=1.1)
   axis(side=1, at=seq(-12, 12, by=4), labels=c(-12, -8, -4, 0, 4, 8, 12), cex.axis=1.2)
   mtext(file.main[2], cex=1.25, line=0.15)
   legend("topleft", legend=c("Positively", "Negatively"), col=c(red, blue), pch=19, cex=1.25)
   dev.off()
}

## NBL expressed genes
xlab.text <- "S/G1-like fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_SRC_p1e-3_TERT")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c(expression(bold("NBL proliferation-correlated genes")), expression("Expression vs."~italic('in silico')~"sorting"))
plotVolcano(src.tpm.gene, 0.001, genes, file.de, file.main, xlab.text)

## NBL expressed genes
xlab.text <- "High/Low risk fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_DE+SRC_p1e-3_TERT")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL high risk, proliferation-associated genes", "Low (L) vs. High (H)")
plotVolcano(de.tpm.gene[overlaps,], 0.001, genes, file.de, file.main, xlab.text)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# Last Modified: 11/08/21
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_median0")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC_p1e-3_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC_p1e-3_down")

setTable <- function(wd.de.reactome) {
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
}

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC_p1e-3_up")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[6, 2] <- "Respiratory electron transport, ATP synthesis by chemiosmotic coupling"

reactome.up <- subset(reactome, Entities.pValue <= 1e-8)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "158 genes")
file.de <- file.path(wd.de.reactome, "genes_DE+SRC_p1e-3_n158_up.pdf")

pdf(file.de, height=4, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 12), xaxt="n", names.arg=reactome.up$Pathway.name, col=flowjo.red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 12, by=4))
mtext(main.text[2], line=0.3)
abline(v=8, lty=5)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC_p1e-3_down")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "464 genes")
file.de <- file.path(wd.de.reactome, "genes_DE+SRC_p1e-3_n464_down.pdf")

pdf(file.de, height=2.5, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col=flowjo.blue, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-8, lty=5)

axis(side=1, at=seq(-12, 0, by=4))
mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# Last Modified: 13/08/21
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_median0")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_SRC_p1e-3_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_SRC_p1e-3_down")

setTable <- function(wd.de.reactome) {
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
}

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_SRC_p1e-3_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[6, 2] <- "Respiratory electron transport, ATP synthesis by chemiosmotic coupling"

reactome.up <- subset(reactome, Entities.pValue <= 1e-15)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "426 genes")
file.de <- file.path(wd.de.reactome, "genes_SRC_p1e-3_n426_up.pdf")

pdf(file.de, height=4, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 16), xaxt="n", names.arg=reactome.up$Pathway.name, col=red.light, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 16, by=4))
mtext(main.text[2], line=0.3)
abline(v=15, lty=5)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_SRC_p1e-3_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 4e-4)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "1,055 genes")
file.de <- file.path(wd.de.reactome, "genes_SRC_p1e-3_n1055_down.pdf")

pdf(file.de, height=2.5, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-16, 0), xaxt="n", names.arg="", col=blue.light, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-15, lty=5)

axis(side=1, at=seq(-16, 0, by=4))
mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# Last Modified: 11/08/21
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_median0")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC+E_p1e-3_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC+E_p1e-3_down")

setTable <- function(wd.de.reactome) {
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
}

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC+E_p1e-3_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
#reactome[1, 2] <- "TP53 Regulates Transcription Involved in Cytochrome C Release"

reactome.up <- subset(reactome, Entities.pValue <= 2e-3)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated", "11 genes")
file.de <- file.path(wd.de.reactome, "genes_DE+SRC+CD_p1e-3_n11_up.pdf")

pdf(file.de, height=4.5, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 4, by=2))
mtext(main.text[2], line=0.3)
abline(v=3, lty=5)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC_p1e-3_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "464 genes")
file.de <- file.path(wd.de.reactome, "genes_DE+SRC_p1e-3_n464_down.pdf")

pdf(file.de, height=2.5, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col=blue.light, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-8, lty=5)

axis(side=1, at=seq(-12, 0, by=4))
mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# 
# Last Modified: 27/08/21
# -----------------------------------------------------------------------------
genes <- c("TERT", "ANAPC11", "TRIAP1", "ICT1")
genes <- c("FKBP3", "MRPL11", "PSMA6")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, flowjo.red, "bottomright")
}

genes <- c("ALK", "ATRX", "MYCN")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotSRC(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, flowjo.blue, "bottomright")
}

###
##
genes <- c("TERT", "ALK", "ATRX", "MYCN", "ANAPC11", "TRIAP1", "ICT1")
genes <- c("FKBP3", "MRPL11", "PSMA6")
genes <- c("KLF15", "RAD51C", "PALB2", "NEDD8")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   
   file.name <- paste0("boxplot0_NBL_L+H_", genes[g])
   plotBox02(wd.de.plots, file.name, tpm.gene.log2[id, rownames(subset(samples, risk == "low"))], tpm.gene.log2[id, rownames(subset(samples, risk == "high"))], main=genes[g], names=c("Low", "High"), cols=c(blue.lightest, red.lightest))
   #file.name <- paste0("boxplot0_NBL_L+H_test_", genes[g])
   #plotBox02(wd.de.plots, file.name, tpm.gene.log2[id, rownames(subset(samples, risk == "low"))], tpm.gene.log2[id, rownames(subset(samples, risk == "high"))], main=genes[g], names=c("Low", "High"), cols=c(green.lighter, orange.lighter))
}

###
## 09/09/21
yims <- c(min(phenos$COR), max(phenos$COR))

file.name <- paste0("boxplot0_NBL_L+H")
plotBox020(wd.de.plots, file.name, subset(phenos, risk == "low")$COR, subset(phenos, risk == "high")$COR, main=paste0("High median = ", round0(median(subset(phenos, risk == "high")$COR), digits=4)), names=c("Low", "High"), cols=c(blue.lightest, red.lightest), yims)

file.name <- paste0("boxplot0_NBL_LR+MNA")
plotBox020(wd.de.plots, file.name, subset(phenos, risk == "low")$COR, subset(phenos, MYCN_amp == 1)$COR, main=paste0("MNA median = ", round0(median(subset(phenos, MYCN_amp == 1)$COR), digits=4)), names=c("LR", "MNA"), cols=c(blue.lightest, red.lightest), yims)

file.name <- paste0("boxplot0_NBL_LR+TERT")
plotBox020(wd.de.plots, file.name, subset(phenos, risk == "low")$COR, subset(phenos, TERT_Rearrangement == 1)$COR, main=paste0("TERT median = ", round0(median(subset(phenos, TERT_Rearrangement == 1)$COR), digits=4)), names=c("LR", "TERT"), cols=c(blue.lightest, red.lightest), yims)

file.name <- paste0("boxplot0_NBL_LR+HR")
hr <- subset(subset(subset(phenos, risk == "high"), MYCN_amp == 0), TERT_Rearrangement == 0)
plotBox020(wd.de.plots, file.name, subset(phenos, risk == "low")$COR, hr$COR, main=paste0("HR median = ", round0(median(hr$COR), digits=4)), names=c("LR", "HR"), cols=c(blue.lightest, red.lightest), yims)

# -----------------------------------------------------------------------------
# plotStripchart (LR, MNA, TERT, HR)
# Last Modified: 13/09/21
# -----------------------------------------------------------------------------
#phenos$group <- 0
#phenos[rownames(subset(phenos, risk == "low")),]$group <- 0
#phenos[rownames(subset(phenos, TERT_Rearrangement == 1)),]$group <- 2
#phenos[rownames(subset(phenos, MYCN_amp == 1)),]$group <- 1
#hr <- rownames(subset(subset(subset(phenos, risk == "high"), group != 1), group != 2))
#phenos[hr,]$group <- 3

#phenos.group <- phenos[order(phenos$group),]

#file.name <- paste0("boxplot+stripchart_NBL_LR+HR_Q4")
#main <- "Cell cycle status by NBL groups"
#plotStripchart(wd.de.plots, file.name, phenos.group, main, names=c("LR", "MNA", "TERT", "HR"), cols=c(blue, blue.light, red.light, red), height=6, width=10)

# -----------------------------------------------------------------------------
# plotStripchart (LR, MNA, TERT, ATRX, HR)
# Last Modified: 13/09/21
# -----------------------------------------------------------------------------
phenos.atrx <- phenos
phenos.atrx$group <- 0
phenos.atrx[rownames(subset(phenos.atrx, risk == "low")),]$group <- 0
phenos.atrx[rownames(subset(phenos.atrx, TERT_Rearrangement == 1)),]$group <- 2
phenos.atrx[rownames(subset(phenos.atrx, MYCN_amp == 1)),]$group <- 3
hr <- rownames(subset(subset(subset(phenos.atrx, risk == "high"), group != 2), group != 3))
phenos.atrx[rownames(subset(phenos.atrx[hr,], ATRX != "0")),]$group <- 1
phenos.atrx[rownames(subset(phenos.atrx[hr,], ATRX == "0")),]$group <- 4

phenos.atrx <- phenos.atrx[order(phenos.atrx$group),]

file.name <- paste0("boxplot+stripchart_NBL_LR+HR+ATRX_Q4")
main <- "Cell cycle status by NBL groups"
plotStripchartATRX(wd.de.plots, file.name, phenos.atrx, main, names=c("LR", "ATRX", "TERT", "MNA", "HR"), cols=c(blue, blue.light, red.light, red), height=6, width=10)





# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# Last Modified: 11/08/21
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway_median0")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC+CD_p1e-3_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC+CD_p1e-3_down")

setTable <- function(wd.de.reactome) {
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
}

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC+CD_p1e-3_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[1, 2] <- "TP53 Regulates Transcription Involved in Cytochrome C Release"

reactome.up <- subset(reactome, Entities.pValue <= 3e-3)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Up-regulated", "11 genes")
file.de <- file.path(wd.de.reactome, "genes_DE+SRC+CD_p1e-3_n11_up.pdf")

pdf(file.de, height=3.5, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=reactome.up$Pathway.name, col=red, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 4, by=2))
mtext(main.text[2], line=0.3)
abline(v=3, lty=5)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_DE+SRC_p1e-3_down")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-4)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "464 genes")
file.de <- file.path(wd.de.reactome, "genes_DE+SRC_p1e-3_n464_down.pdf")

pdf(file.de, height=2.5, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-12, 0), xaxt="n", names.arg="", col=blue.light, xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-8, lty=5)

axis(side=1, at=seq(-12, 0, by=4))
mtext(main.text[2], line=0.3)
dev.off()











# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.nbl   ## 17/05/20 WHY nrds.RT.NRFD.sclc.nl? ## MUY MUY IMPORTANTE!!

tpm.gene.log2.m <- tpm.gene.log2.m
nrow(tpm.gene.log2.m)
# [1] 34908

tpm.gene.log2.m <- tpm.gene.log2.m[rownames(de.tpm.gene),]
nrow(tpm.gene.log2.m)
# [1] 22899

tpm.gene.log2.m <- tpm.gene.log2.m[intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m)),]
nrow(tpm.gene.log2.m)
# [1] 18764

tpm.gene.log2.m.rfd <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD)
tpm.gene.log2.m.rfd$length <- abs(tpm.gene.log2.m.rfd$end_position - tpm.gene.log2.m.rfd$start_position)
nrow(tpm.gene.log2.m.rfd)
# [1] 31703
# [1] 21419
# [1] 17654
nrow(subset(tpm.gene.log2.m.rfd, TRC == 0))
# [1] 3
# [1] 2
# [1] 2

###
## TTR and CTR (IZ +TZ)
save(tpm.gene.log2.m.rfd, file=file.path(wd.de.data, "tpm_gene_log2+1_m_rfd0.RData"))

load(file=file.path(wd.de.data, "tpm_gene_log2+1_m_rfd0.RData"))
file.name <- file.path(wd.de.data, "tpm_gene_log2_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm_gene_median0+1_log2_m_rfd.RData")
#file.name <- file.path(wd.de.data, "tpm_gene_r5p47_log2_m_rfd.RData")

###
## ADD: 29/08/21
load(file=file.path(wd.de.data, "tpm_gene_log2+1_m_rfd0.RData"))
file.name <- file.path(wd.de.data, "tpm_gene_log2+1_m_rfd.RData")
setTRC(tpm.gene.log2.m.rfd, rfd=0.9, file.name)
load(file.name)

## CTR
length(which(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD < 0))
# [1] 125
# [1] 111
# [1] 104
length(which(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD > 0))
# [1] 147
# [1] 125
# [1] 113

# -----------------------------------------------------------------------------
# Shared genes (All)
# Last Modified: 10/09/20
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/tpm_gene_log2+1_m_rfd.RData")

gene.nbl   <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz),   rownames(tpm.gene.log2.m.rfd.ctr.tz),   rownames(tpm.gene.log2.m.rfd.ttr))
gene.nbl.e <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz.e), rownames(tpm.gene.log2.m.rfd.ctr.tz.e), rownames(tpm.gene.log2.m.rfd.ttr.e))
gene.nbl.l <- list(rownames(tpm.gene.log2.m.rfd.ctr.iz.l), rownames(tpm.gene.log2.m.rfd.ctr.tz.l), rownames(tpm.gene.log2.m.rfd.ttr.l))


# -----------------------------------------------------------------------------
# Density (All)
# Last Modified: 17/03/21
# -----------------------------------------------------------------------------
nrow(tpm.gene.log2.m.rfd.ctr.iz)/nrow(nrds.RT.NRFD.nbl.ctr.iz)*1000
# [1] 15.9889
nrow(tpm.gene.log2.m.rfd.ctr.tz)/nrow(nrds.RT.NRFD.nbl.ctr.tz)*1000
# [1] 12.17577
nrow(tpm.gene.log2.m.rfd.ttr)/nrow(nrds.RT.NRFD.nbl.ttr)*1000
# [1] 11.30821

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/nrow(nrds.RT.NRFD.nbl.ctr.iz.e)*1000
# [1] 20.23631
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/nrow(nrds.RT.NRFD.nbl.ctr.iz.l)*1000
# [1] 8.101177

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/nrow(nrds.RT.NRFD.nbl.ctr.tz.e)*1000
# [1] 16.13775
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/nrow(nrds.RT.NRFD.nbl.ctr.tz.l)*1000
# [1] 8.227748

nrow(tpm.gene.log2.m.rfd.ttr.e)/nrow(nrds.RT.NRFD.nbl.ttr.e)*1000
# [1] 14.8517
nrow(tpm.gene.log2.m.rfd.ttr.l)/nrow(nrds.RT.NRFD.nbl.ttr.l)*1000
# [1] 7.297637

# -----------------------------------------------------------------------------
# E vs. L (All)
# Last Modified: 10/09/20
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/tpm_gene_log2+1_m_rfd.RData")

nrow(tpm.gene.log2.m.rfd.ttr)/31703
# [1] 0.7274075
nrow(tpm.gene.log2.m.rfd.ctr.iz)/31703
# [1] 0.155632
nrow(tpm.gene.log2.m.rfd.ctr.tz)/31703
# [1] 0.1166767

nrow(tpm.gene.log2.m.rfd.ttr.e)/31703
# [1] 0.5072075
nrow(tpm.gene.log2.m.rfd.ttr.l)/31703
# [1] 0.2202

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/31703
# [1] 0.128032
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/31703
# [1] 0.02759991

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/31703
# [1] 0.07718512
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/31703
# [1] 0.03949153

###
## median0
nrow(tpm.gene.log2.m.rfd.ttr)/
 # [1] 
 nrow(tpm.gene.log2.m.rfd.ctr.iz)/
 # [1] 
 nrow(tpm.gene.log2.m.rfd.ctr.tz)/
 # [1] 
 
 nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/
 # [1] 
 nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/
 # [1] 
 
 nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/
 # [1] 
 nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/
 # [1] 
 
 ###
 ## r5p47
 nrow(tpm.gene.log2.m.rfd.ttr)/17654
# [1] 0.7132095
nrow(tpm.gene.log2.m.rfd.ctr.iz)/17654
# [1] 0.1759941
nrow(tpm.gene.log2.m.rfd.ctr.tz)/17654
# [1] 0.1105698

nrow(tpm.gene.log2.m.rfd.ctr.iz.e)/17654
# [1] 
nrow(tpm.gene.log2.m.rfd.ctr.iz.l)/17654
# [1] 

nrow(tpm.gene.log2.m.rfd.ctr.tz.e)/17654
# [1] 
nrow(tpm.gene.log2.m.rfd.ctr.tz.l)/17654
# [1] 

# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 27/11/20; 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
file.name <- file.path(wd.de.data, "tpm_gene_log2+1_m_rfd.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 12)
# > min(tpm.gene.log2.m.rfd$MEDIAN)
# [1] 0
# > max(tpm.gene.log2.m.rfd$MEDIAN)
# [1] 12.40235

file.name <- paste0("boxplot3_nbl_tpm.gene+1_RFD")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="NBL expression", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim)

file.name <- paste0("boxplot4_nbl_tpm.gene+1_RFD")
plotBox4(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL expression (CTR)", names=c("L", "E", "L", "E"), cols=c(blue, blue, red, red), ylim)

file.name <- paste0("boxplot6_nbl_tpm.gene+1_RFD")
plotBox6(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL expression", names=c("L", "E", "L", "E", "L", "E"), cols=c("black", "black", blue, blue, red, red), ylim)

###
## 30/08/21
file.name <- file.path(wd.de.data, "tpm_gene_median0_log2+1_m_rfd.RData")
load(file.name)

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 12)
# > min(tpm.gene.log2.m.rfd$MEDIAN)
# [1] 0
# > max(tpm.gene.log2.m.rfd$MEDIAN)
# [1] 12.74563

file.name <- paste0("boxplot3_nbl_tpm.gene.median0+1_RFD")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="NBL expressed genes", names=c("TTR", "TZ", "IZ"), cols=c("black", blue, red), ylim)

file.name <- paste0("boxplot4_nbl_tpm.gene.median0+1_RFD")
plotBox4(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL expressed genes (CTR)", names=c("L", "E", "L", "E"), cols=c(blue, blue, red, red), ylim)

file.name <- paste0("boxplot6_nbl_tpm.gene.median0+1_RFD")
plotBox6(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL expressed genes", names=c("L", "E", "L", "E", "L", "E"), cols=c("black", "black", blue, blue, red, red), ylim)




# -----------------------------------------------------------------------------
# RFD vs. TPM (All)
# Last Modified: 01/09/20; 29/05/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), max(tpm.gene.log2.m.rfd$MEDIAN))

file.name <- paste0("boxplot0_nbl_tpm.gene_rfd_TTR+IZ_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes", names=c("TTR", "IZ"), cols=c("black", red), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene_rfd_TTR+TZ_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="NBL genes", names=c("TTR", "TZ"), cols=c("black", blue), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene_rfd_TZ+IZ_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes (CTR)", names=c("TZ", "IZ"), cols=c(blue, red), ylim)

##
file.name <- paste0("boxplot0_nbl_tpm.gene_rfd_IZ_E+L_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL IZ genes", names=c("Late", "Early"), cols=c(red, red), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene_rfd_TZ_E+L_nasa.blue")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="NBL TZ genes", names=c("Late", "Early"), cols=c(blue, blue), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene_rfd_TTR_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr.l, tpm.gene.log2.m.rfd.ttr.e, main="NBL TTR genes", names=c("Late", "Early"), cols=c("black", "black"), ylim)

# -----------------------------------------------------------------------------
# RFD vs. TPM (r5p47)
# Last Modified: 29/08/21; 02/09/20; 25/08/20; 29/05/20
# -----------------------------------------------------------------------------
# load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/tpm_gene_median0_log2_m_rfd.RData")
# > nrow(tpm.gene.log2.m.rfd.ctr.iz.e)
# [1] 3231
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), max(tpm.gene.log2.m.rfd$MEDIAN))

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_TTR+IZ_Figure1")
plotBox00(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_IZ_E+L_Figure1")
plotBox00(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL IZ genes", names=c("Late", "Early"), cols=c("red", "red"), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_IZ_E_HO+CD_Figure1")
plotBox00(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="NBL early IZ genes", names=c("HO", "CD"), cols=c("red", "red"), ylim)



#max(tpm.gene.log2.m.rfd$MEDIAN)
#[1] 12.40235
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 12.5)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_TTR+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)
file.name <- paste0("boxplot_nbl_tpm.gene.median0+1_rfd_TTR+IZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_TTR+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="NBL genes", names=c("TTR", "TZ"), cols=c("black", "blue"), ylim)

#file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_IZ+TZ")
#plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="NBL genes (CTR)", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_TZ+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes (CTR)", names=c("TZ", "IZ"), cols=c("blue", "red"), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_TZ+IZ_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz, tpm.gene.log2.m.rfd.ctr.iz, main="NBL genes", names=c("TZ", "IZ"), cols=c("blue", "red"), ylim)

##
file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_IZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL IZ genes", names=c("Late", "Early"), cols=c("red", "red"), ylim)
file.name <- paste0("boxplot_nbl_tpm.gene.median0+1_rfd_IZ_E+L")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.l, tpm.gene.log2.m.rfd.ctr.iz.e, main="NBL IZ genes", names=c("Late", "Early"), cols=c("red", "red"), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_IZ_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="NBL early IZ genes", names=c("HO", "CD"), cols=c("red", "red"), ylim)
file.name <- paste0("boxplot_nbl_tpm.gene.median0+1_rfd_IZ_E_HO+CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.iz.e.cd, main="NBL early IZ genes", names=c("HO", "CD"), cols=c("red", "red"), ylim)

##
file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_TZ_E+L")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="NBL TZ genes", names=c("Late", "Early"), cols=c("blue", "blue"), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_TZ_E+L_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.l, tpm.gene.log2.m.rfd.ctr.tz.e, main="NBL TZ genes", names=c("Late", "Early"), cols=c("blue", "blue"), ylim)

file.name <- paste0("boxplot0_nbl_tpm.gene.median0_rfd_TZ_E_HO+CD")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="NBL early TZ genes", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)
#file.name <- paste0("boxplot_sclc_tpm.gene.median0_rfd_TZ_E_HO+CD_long")
#plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, main="NBL early TZ genes", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)

# -----------------------------------------------------------------------------
# GENE_NRFD plotBoxTSSNRFD (median0)
# Link: https://www.r-bloggers.com/whisker-of-boxplot/
# Last Modified: 31/08/20
# -----------------------------------------------------------------------------
IQR <- 0.011192776 - 0.001195048
upper <- 0.011192776 + 1.5*IQR
IQR <- 0.007344794 - 0.001042600
lower <- 0.001042600 - 1.5 * IQR
#ylim <- c(-0.008410691, 0.02618937)
ylim <- c(-0.01, 0.028)

plotBoxNRFD(base, BASE, ylim, tpm.gene.log2.m.rfd.ctr.iz.e.cd, tpm.gene.log2.m.rfd.ctr.iz.e.ho, tpm.gene.log2.m.rfd.ctr.tz.e.cd, tpm.gene.log2.m.rfd.ctr.tz.e.ho)

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/05/20
# -----------------------------------------------------------------------------
genes <- c("RPL38")
for (g in 1:length(genes)) {
 id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
 plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}

##
#plotNLSCLC <- function(id, gene) {
#   median(as.numeric(tpm.gene.lung.log2[id,]))
#   median(as.numeric(tpm.gene.log2[id,]))
#   testU(as.numeric(tpm.gene.lung.log2[id,]), as.numeric(tpm.gene.log2[id,]))
#   file.name <- paste0("boxplot_sclc_tpm.gene.r5p47_gene_", gene, "_")
#   plotBox2(wd.de.plots, file.name, tpm.gene.lung.log2[id,], tpm.gene.log2[id,], main=gene, names=c("Lung", "SCLC"))
#}

#genes <- c("GTF3C2", "SUPT7L", "RAD9A", "E2F3", "ALK")
#genes <- c("MARS", "KIF18B", "BRCA2")
#genes <- c("BLM", "POLE")
#for (g in 1:length(genes)) {
#   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
#   plotNLSCLC(id, genes[g])
#}

# -----------------------------------------------------------------------------
# SCLC-specific IZ genes
# Last Modified: 31/08/20
# -----------------------------------------------------------------------------
## Expressed IZ-E, CD genes
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.e.cd))
de.tpm.gene.rfd.iz.e.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.e.cd$Q <- qvalue(de.tpm.gene.rfd.iz.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.e.cd, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_e_cd_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.e.ho))
de.tpm.gene.rfd.iz.e.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.e.ho$Q <- qvalue(de.tpm.gene.rfd.iz.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.e.ho, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_e_ho_n53.txt"), colnames=T, rownames=F, sep="\t")

## Expressed TZ-E, CD genes
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.e.cd))
de.tpm.gene.rfd.tz.e.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.e.cd$Q <- qvalue(de.tpm.gene.rfd.tz.e.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.e.cd, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_tz_e_cd_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.e.ho))
de.tpm.gene.rfd.tz.e.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.e.ho$Q <- qvalue(de.tpm.gene.rfd.tz.e.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.e.ho, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_tz_e_ho_n53.txt"), colnames=T, rownames=F, sep="\t")




























## NBL ALL genes
xlab.text <- "NBL S/G1 fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-3_all_Helicases")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL expressed genes (n=22,899)", "Expression vs. In-silico sorting")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=11)

xlab.text <- "NBL S/G1 fold change [log2]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-3_all_TFBS")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.de <- paste0(plot.de, ".pdf")
file.main <- c("NBL expressed genes (n=22,899)", "Expression vs. In-silico sorting")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main, xlab.text, ymax=11)

#overlaps <- intersect(tfbs, de.tpm.gene$external_gene_name)
#setdiff(tfbs, overlaps)



## NBL IZ-E, CD genes
xlab.text <- "NBL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-5_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL early IZ, CD genes (n=1,579)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.00001, genes, file.de, file.main, xlab.text, ymax=8.5)

## SCLC TZ-E, HO genes
xlab.text <- "NBL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-5_tz_e_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL early TZ, HO genes (n=951)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.tz.e.ho, 0.00001, genes, file.de, file.main, xlab.text, ymax=8.5)





## SCLC IZ-E, CD genes
xlab.text <- "NBL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-5_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL early IZ, CD genes (n=1,579)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.00001, genes, file.de, file.main, xlab.text, ymax=8.5)

## SCLC IZ-E, CD genes
xlab.text <- "NBL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-2_iz_e_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL early IZ, CD genes (n=1,579)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.e.cd, 0.01, genes, file.de, file.main, xlab.text, ymax=8.5)

## SCLC TZ-E, HO genes
xlab.text <- "NBL S/G1 [log2 fold change]"
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_median0_rfd_p1e-2_tz_e_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL early TZ, HO genes (n=951)", "Correlation between expression and in-silico sorting")
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.tz.e.ho, 0.01, genes, file.de, file.main, xlab.text, ymax=8.5)











# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.pathway <- file.path(wd.de, "pathway")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_down")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_up")
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_down")

setTable <- function(wd.de.reactome) {
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
}

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-5_up")
setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[2, 2] <- "TP53 Regulates Transcription of Cytochrome C Release"
#reactome[6, 2] <- "TP53 regulates transcription of genes involved in cytochrome c release"
 
reactome.up <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "8 genes")
file.de <- file.path(wd.de.reactome, "genes_IZ-E-CD_p1e-5_n8_up.pdf")

pdf(file.de, height=5.5, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 8), xaxt="n", names.arg=reactome.up$Pathway.name, col="gold", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 8, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_IZ-E-CD_p1e-2_down")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "176 genes")
file.de <- file.path(wd.de.reactome, "genes_IS-E-CD_p1e-3_n176_down.pdf")

pdf(file.de, height=1.87, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-8, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-3, lty=5)

axis(side=1, at=seq(-8, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()

###
## Link: https://www.statmethods.net/graphs/bar.html
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_up")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")

reactome.up <- subset(reactome, Entities.pValue <= 1e-3)[, 2:7]
reactome.up$log10P <- -log10(reactome.up$Entities.pValue)
reactome.up <- reactome.up[order(reactome.up$log10P),]

main.text <- c("Positively-correlated", "29 genes")
file.de <- file.path(wd.de.reactome, "genes_TZ-E-HO_p1e-2_n29_up.pdf")

pdf(file.de, height=2.95, width=7.5)
par(mar=c(4,27.5,4,1))   # increase y-axis margin.
barplot(reactome.up$log10P, main=main.text[1], las=1, horiz=T, xlim=c(0, 8), xaxt="n", names.arg=reactome.up$Pathway.name, col="gold", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3

axis(side=1, at=seq(0, 8, by=1))
mtext(main.text[2], line=0.3)
dev.off()

##
wd.de.reactome <- file.path(wd.de.pathway, "reactome_TZ-E-HO_p1e-2_down")
#setTable(wd.de.reactome)
reactome <- readTable(file.path(wd.de.reactome, "result.tsv"), header=T, rownames=F, sep="\t")
reactome[1, 2] <- "ATM-mediated phosphorylation of repair proteins at DNA DSBs"
reactome[5, 2] <- "Regulation of gene expression in branching morphogenesis"

reactome.down <- subset(reactome, Entities.pValue <= 1e-2)[, 2:7]
reactome.down$log10P <- -log10(reactome.down$Entities.pValue)
reactome.down <- reactome.down[order(reactome.down$log10P),]
reactome.down$log10P <- reactome.down$log10P * -1

main.text <- c("Negatively-correlated", "99 genes")
file.de <- file.path(wd.de.reactome, "genes_TS-E-HO_p1e-2_n99_down.pdf")

pdf(file.de, height=3.3, width=7.5)
par(mar=c(4,1,4,27.5))   # increase y-axis margin.
posbar <- barplot(reactome.down$log10P, main=main.text[1], las=1, horiz=T, xlim=c(-8, 0), xaxt="n", names.arg="", col="steelblue1", xlab="-log10(p-value)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=reactome.down$Pathway.name, xpd=NA)
abline(v=-3, lty=5)

axis(side=1, at=seq(-8, 0, by=1))
mtext(main.text[2], line=0.3)
dev.off()














## All
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_p3e04")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC expressed genes", paste0("n=18,004"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene, 0.0003, genes, file.de, file.main, xlab.text, ymax=6.8)

##
xlab.text <- "SCLC M2/M1 [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.1")

## IZ-S
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_p3e04_iz_sp_cd")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC expressed, IZ-S, CD genes", paste0("n=642"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp.cd, 0.0003, genes, file.de, file.main, xlab.text, ymax=4.5)

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08_iz_sp_ho")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC IZ-S, HO genes", paste0("n=657"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.iz.sp.ho, 0.08, genes, file.de, file.main, xlab.text)

## TOR1
de.tpm.gene.rfd.ctr.sp <- rbind(de.tpm.gene.rfd.iz.sp, de.tpm.gene.rfd.tz.sp)
de.tpm.gene.rfd.ctr.sp$Q <- qvalue(de.tpm.gene.rfd.ctr.sp$P)$qvalue

plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_rfd_fdr0.08_ctr_sp")
genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC CTR-S genes", paste0("n=2,173"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de.tpm.gene.rfd.ctr.sp, 0.08, genes, file.de, file.main, xlab.text)






















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
# D.E.
# Last Modified: 08/01/20
# -----------------------------------------------------------------------------
#tpm.gene.log2.0      <- tpm.gene.log2[rownames(de.tpm.gene.iz),]
#tpm.gene.lung.log2.0 <- tpm.gene.lung.log2[rownames(de.tpm.gene.iz),]
overlaps <- intersect(intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.cd), rownames(tpm.gene.log2)), rownames(tpm.gene.lung.log2))
length(overlaps)
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

de2.iz.sp <- de2[intersect(rownames(de2), rownames(tpm.gene.log2.m.rfd.ctr.iz.cd)),]
de2.iz.sp$FDR <- p.adjust(de2.iz.sp$P, method="BH", n=length(de2.iz.sp$P))
#de2.iz.sp$Q <- qvalue(de2.iz.sp$P)$qvalue
de2.iz.sp <- de2.iz.sp[order(de2.iz.sp$P),]

## FDR
#library(qvalue)
#de2$Q   <- qvalue(de2$P)$qvalue
de2$FDR <- p.adjust(de2$P, method="BH", n=length(de2$P))
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de2 <- de2[order(de2$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de2.tpm.gene <- cbind(annot[rownames(de2.iz.sp),], de2.iz.sp)   ## BE EXTRA CAREFUL!!

save(de2.tpm.gene, samples, file=file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_cd.RData"))
writeTable(de2.tpm.gene, file.path(wd.de.data, "de2_sclc_tpm-gene-r5p47_lung_U_q_n70+41_iz_sp_cd.txt"), colnames=T, rownames=F, sep="\t")
# nrow(de2.tpm.gene)
# [1] 616

## Volcano
xlab.text <- "SCLC/Lung [log2FC]"
plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_r5p47_iz_sp_cd_bh1e-16_lung")

genes <- readTable(paste0(plot.de, ".tab"), header=T, rownames=F, sep="\t")
file.main <- c("SCLC, IZ-S and CD genes", paste0("n=616"))
file.de <- paste0(plot.de, ".pdf")
plotVolcano(de2.tpm.gene, 1E-16, genes, file.de, file.main, xlab.text)









overlaps <- intersect(rownames(de2.tpm.gene), rownames(tpm.gene))
length(overlaps)
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
