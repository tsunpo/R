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
BASE <- "SCLC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-wgs-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

wd.rt      <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data <- file.path(wd.rt, "data")

samples <- readTable(file.path(wd.wgs, "sclc_wgs_n101.cm2"), header=T, rownames=T, sep="")[,-1]

## Expression level (11/04/19)
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
expressed <- rownames(tpm.gene)
# > length(expressed)
# [1] 19131

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; CM2 vs CM1)
# Last Modified: 22/05/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : RB1_MUT (1) vs RB1_WT (0) as factor
file.name <- paste0("de_", base, "_wgs-gene_cm2_wilcox_q_n101")

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
   colnames(de.chr)[3:4] <- c("CM2_WT", "CM2")
   de <- rbind(de, de.chr)
}
de$FDR <- testFDR(de$P, argv$test.fdr)
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.wgs.gene <- cbind(annot[rownames(de),], de)

save(de.wgs.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.wgs.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

wgs.gene.log2 <- log2(wgs.gene + 0.01)
save(wgs.gene, file=file.path(wd.de.data, paste0(base, "_wgs.gene.RData")))

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure 1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, pvalue, genes, file.de, file.main) {
   #pvalue <- fdrToP(fdr, de)
   fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #xmax <- 1
   ymax <- max(de$log10P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="log2FC(SCLC CM2/CM1)", ylab="-log10(p-value)", col="darkgray", main=file.main[1])
 
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
#plot.main <- "Differential read depth between CM2 and CM1 in SCLC"
#plot.de <- file.path(wd.de.plots, "volcanoplot_sclc_rt_p1e-6_cm2_n70_")
plot.main <- "Differential expression between CM2 and CM1 in SCLC"
plot.de <- file.path(wd.de.plots, "volcanoplot-r5p47_sclc_p1e-6_cm2_n70_")

## Chr2
genes <- readTable(paste0(plot.de, "CM2.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.wgs.gene$external_gene_name),]

file.main <- c(plot.main, "")
#file.de <- paste0(plot.de, "_chr2_log2FC1.pdf")
file.de <- paste0(plot.de, "CM2.pdf")
#plotVolcano(de.wgs.gene, 1.00E-06, genes, file.de, file.main)
plotVolcano(de.tpm.gene, 1.00E-06, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_sclc_wgs-gene_cm2_wilcox_q_n101")
writeRNKformat(de.wgs.gene, wd.de.gsea, file.name)

## Tirosh et al 2016 
load(file.path(wd.src.ref, "cycle.RData"))                           ## See guide-to-the/cycle.R
core.G1S <- intersect(core.G1S, rownames(wgs.gene.log2))             ## 41/43/43 are expressed in the dataset
core.G2M <- intersect(core.G2M, rownames(wgs.gene.log2))             ## 54/54/55
core.Stemness <- intersect(core.Stemness, rownames(wgs.gene.log2))   ## 54/58/63

writeGRPformat(core.G1S, wd.de.gsea, "core.G1-S")
writeGRPformat(core.G2M, wd.de.gsea, "core.G2-M")
writeGRPformat(core.Stemness, wd.de.gsea, "core.Stemness")

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(wgs.gene.log2))   ## 262/304
periodic.G2M <- intersect(periodic.G2M, rownames(wgs.gene.log2))   ## 792/876

writeGRPformat(periodic.G1S, wd.de.gsea, "G1-S")
writeGRPformat(periodic.G2M, wd.de.gsea, "G2-M")

# -----------------------------------------------------------------------------
# SRC
# Last Modified: 23/05/19
# -----------------------------------------------------------------------------
overlaps  <- intersect(colnames(wgs.gene), colnames(tpm.gene))
expressed <- intersect(rownames(wgs.gene), rownames(tpm.gene))
wgs.gene.o <- wgs.gene[expressed, overlaps]
tpm.gene.o <- tpm.gene[expressed, overlaps]
samples.o  <- samples[overlaps,]

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; CM2 vs CM1)
# Last Modified: 23/05/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : RB1_MUT (1) vs RB1_WT (0) as factor
file.name <- paste0("de_", base, "_wgs-gene-o_cm2_wilcox_q_n70")
wgs.gene.o.log2 <- log2(wgs.gene.o + 0.01)

de <- c()
for (c in 1:22) {
   chr <- chrs[c]
   samples.o[, chr] <- as.factor(samples.o[, chr])
   
   ensGene.chr <- subset(ensGene, chromosome_name == chr)
   wgs.gene.o.log2.chr <- wgs.gene.o.log2[intersect(rownames(wgs.gene.o.log2), rownames(ensGene.chr)),]

   argv <- data.frame(predictor=chr, predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
   de.chr <- differentialAnalysis(wgs.gene.o.log2.chr, samples.o, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
   colnames(de.chr)[3:4] <- c("CM2_WT", "CM2")
   de <- rbind(de, de.chr)
}
de$FDR <- testFDR(de$P, argv$test.fdr)
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.wgs.gene <- cbind(annot[rownames(de),], de)

save(de.wgs.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.wgs.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; CM2 vs CM1)
# Last Modified: 23/05/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : RB1_MUT (1) vs RB1_WT (0) as factor
file.name <- paste0("de_", base, "_tpm-gene-o_cm2_wilcox_q_n70")
tpm.gene.o.log2 <- log2(tpm.gene.o + 0.01)

de <- c()
for (c in 1:22) {
   chr <- chrs[c]
   samples.o[, chr] <- as.factor(samples.o[, chr])
 
   ensGene.chr <- subset(ensGene, chromosome_name == chr)
   tpm.gene.o.log2.chr <- tpm.gene.o.log2[intersect(rownames(tpm.gene.o.log2), rownames(ensGene.chr)),]
 
   argv <- data.frame(predictor=chr, predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
   de.chr <- differentialAnalysis(tpm.gene.o.log2.chr, samples.o, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
   colnames(de.chr)[3:4] <- c("CM2_WT", "CM2")
   de <- rbind(de, de.chr)
}
de$FDR <- testFDR(de$P, argv$test.fdr)
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# SRC between WGS and TPM
# Last Modified: 23/05/19
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "PCC_R", "PCC_P", "PCC_Q", "LOG2FC_WGS", "LOG2FC_TPM")
de <- toTable(0, length(colnames), nrow(tpm.gene.o.log2), colnames)
rownames(de) <- rownames(tpm.gene.o.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.o.log2), function(x) cor.test(as.numeric(tpm.gene.o.log2[x,]), as.numeric(wgs.gene.o.log2[x,]), method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.o.log2), function(x) cor.test(as.numeric(tpm.gene.o.log2[x,]), as.numeric(wgs.gene.o.log2[x,]), method="spearman", exact=F)[[3]])

## Pearson
de$PCC_R <- mapply(x = 1:nrow(tpm.gene.o.log2), function(x) cor.test(as.numeric(tpm.gene.o.log2[x,]), as.numeric(wgs.gene.o.log2[x,]), method="pearson")$estimate)
de$PCC_P <- mapply(x = 1:nrow(tpm.gene.o.log2), function(x) cor.test(as.numeric(tpm.gene.o.log2[x,]), as.numeric(wgs.gene.o.log2[x,]), method="pearson")$p.value)

## Log2 fold change
de$LOG2FC_WGS <- de.wgs.gene[rownames(de),]$LOG2_FC
de$LOG2FC_TPM <- de.tpm.gene[rownames(de),]$LOG2_FC
de$KEEP <- 0
de$KEEP <- de$LOG2FC_WGS * de$LOG2FC_TPM

## FDR
library(qvalue)
de$Q     <- qvalue(de$P)$qvalue
de$PCC_Q <- qvalue(de$PCC_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.wgs.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.wgs.tpm.gene, samples, file=file.path(wd.de.data, "de_sclc_wgs-tpm-gene-o_rt_src_pcc_q_n70.RData"))
writeTable(de.wgs.tpm.gene, file.path(wd.de.data, "de_sclc_wgs-tpm-gene-o_rt_src_pcc_q_n70.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Final report
# Last Modified: 23/05/19
# -----------------------------------------------------------------------------
de.wgs.tpm.gene.out <- subset(de.wgs.tpm.gene, KEEP > 0)
#de.wgs.tpm.gene.out <- subset(de.wgs.tpm.gene.out, P <= 1E-03)

overlaps <- intersect(rownames(de.wgs.tpm.gene.out), rownames(subset(de.tpm.gene, P <= 1E-03)))
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
plotBoxTPM <- function(gene, wd.de.plots, tpm.gene.log2, de.tpm.gene, samples, ylim=NULL) {
   ids <- subset(ensGene, external_gene_name == gene)   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
   chr <- ids$chromosome_name[1]
   
   for (i in 1:nrow(ids)) {
      id <- ids$ensembl_gene_id[i]

      ##
      file.name <- paste0("boxplot_tpm.gene.o.log2_", gene)
      if (nrow(ids) != 1)
         file.name <- paste0(file.name, "_", id)
      
      gene.tpm <- cbind(t(tpm.gene.o.log2[id, rownames(samples.o)]), samples.o)
      colnames(gene.tpm)[1] <- "MEDIAN"
  
      pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.5)
      boxplot(gene.tpm$MEDIAN ~ gene.tpm[, chr], outline=T, names=c("CM1", "CM2"), ylab="log2(TPM+0.01)", main=paste0(gene, " (", id, ")"))
      
      mtext(paste0("p-value = ", scientific(de.tpm.gene[id,]$P)), cex=1.2, line=0.3)
      dev.off()
   }
}

plotBox <- function(gene, wd.de.plots, wgs.gene.o.log2, tpm.gene.o.log2, samples.o, ylim=NULL) {
 ids <- subset(ensGene, external_gene_name == gene)   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
 chr <- ids$chromosome_name[1]
 
 for (i in 1:nrow(ids)) {
  id <- ids$ensembl_gene_id[i]
  
  ## 
  file.name <- paste0("boxplot_tpm.gene.o.log2_", gene)
  if (nrow(ids) != 1)
   file.name <- paste0(file.name, "_", id)
  
  gene.tpm <- cbind(t(tpm.gene.o.log2[id, rownames(samples.o)]), samples.o)
  colnames(gene.tpm)[1] <- "MEDIAN"
  
  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.5)
  boxplot(gene.tpm$MEDIAN ~ gene.tpm[, chr], outline=T, names=c("CM1", "CM2"), ylab="log2(TPM+0.01)", main=paste0(gene, " (", id, ")"))
  dev.off()
  
  ##
  file.name <- paste0("boxplot_wgs.gene.o.log2_", gene)
  if (nrow(ids) != 1)
   file.name <- paste0(file.name, "_", id)
  
  gene.wgs <- cbind(t(wgs.gene.o.log2[id, rownames(samples.o)]), samples.o)
  colnames(gene.wgs)[1] <- "MEDIAN"
  
  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.5)
  boxplot(gene.wgs$MEDIAN ~ gene.tpm[, chr], outline=T, names=c("CM1", "CM2"), ylab="log2(RD+0.01)", main=paste0(gene, " (", id, ")"))
  mtext(paste0("p-value = ", scientifc(gene.wgs$)), cex=1.2, line=0.3)
  dev.off()
  
  ##
  file.name <- paste0("plot_wgs.tpm.gene.o.log2_", gene)
  if (nrow(ids) != 1)
   file.name <- paste0(file.name, "_", id)
  
  pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=6)
  plot(gene.tpm$MEDIAN ~ gene.wgs$MEDIAN, xlab="log2(RD+0.01)", ylab="log2(TPM+0.01)", main=paste0(gene, " (", id, ")"))
  dev.off()
 }
}


##
genes <- c("RP11-141C7.3", "BRD9", "MAP9", "RAD9A", "POLE", "MCM10", "DNA2", "CCAR1", "SMARCD1", "E2F3", "ERCC8")
for (g in 1:length(genes)) {
   plotBoxTPM(genes[g], wd.de.plots, tpm.gene.o.log2, de.tpm.gene.sclc, samples.o)
   #plotBoxWGS(genes[g], wd.de.plots, wgs.gene.o.log2, de.wgs.tpm.gene, samples.o)
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
