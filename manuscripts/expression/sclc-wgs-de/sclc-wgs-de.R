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
handbooks  <- c("Commons.R", "DifferentialExpression.R")
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
tpm.gene.log2 <- tpm.gene.log2[,rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!

colnames <- c("RHO", "P", "Q", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$RB1_RATE, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$RB1_RATE, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## Log2 fold change
de$LUAD  <- median00(tpm.gene.log2, rownames(subset(samples, Cancer_Type == 0)))
de$LCNEC <- median00(tpm.gene.log2, rownames(subset(samples, Cancer_Type == 1)))
de$SCLC  <- median00(tpm.gene.log2, rownames(subset(samples, Cancer_Type == 2)))
de$LCNEC_LUAD <- de$LCNEC - de$LUAD
de$SCLC_LUAD  <- de$SCLC  - de$LUAD

## FDR
library(qvalue)
de$Q   <- qvalue(de$P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_luad+lcnec+sclc_tpm-gene-r5p47_rb1_src_q_n198.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "de_luad+lcnec+sclc_tpm-gene-r5p47_rb1_src_q_n198.txt"), colnames=T, rownames=F, sep="\t")








# -----------------------------------------------------------------------------
# Cannoli/Kilonova plot
# -----------------------------------------------------------------------------
plotKilonovaFishers <- function(de, genes, file.de, file.main, p, isIg) {
   de.sig.positive <- subset(subset(subset(de, P_Fishers <= p), Effect1 > 0), Effect2 > 0)
   de.sig.negative <- subset(subset(subset(de, P_Fishers <= p), Effect1 < 0), Effect2 < 0)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect1, de$Effect2, pch=16, xlab="Median fold change (log2 SCLC/LUAD)", ylab="Median fold change in CCLE", col="darkgray", main=file.main)
 
   points(de.sig.positive$Effect1, de.sig.positive$Effect2, pch=16, col="red")
   points(de.sig.negative$Effect1, de.sig.negative$Effect2, pch=16, col="dodgerblue")
   if (isIg) {   ## ADD 02/11/17: Immunoglobulin (Ig) variable chain and T-cell receptor (TCR) genes
      de.ig <- de[grep("IG", de$gene_biotype),]
      points(de.ig$Effect1, de.ig$Effect2, pch=16, col="gold")
      de.tr <- de[grep("TR", de$gene_biotype),]
      points(de.tr$Effect1, de.tr$Effect2, pch=16, col="forestgreen")
   }
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect1, gene$Effect2, pch=1, col="black")  
   
         if (genes[g] == "PMAIP1" || genes[g] == "XRCC3")
            text(gene$Effect1, gene$Effect2, genes[g], col="black", adj=c(-0.07, 1.2), cex=0.75)  
         else if (genes[g] == "XRCC2")
            text(gene$Effect1, gene$Effect2, genes[g], col="black", adj=-0.1, cex=0.75)
         else
            if (gene$Effect1 > 0)
               text(gene$Effect1, gene$Effect2, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect1, gene$Effect2, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   if (isIg)
      legend("topleft", legend=c("Upregulated in both", "Downregulated in both", "Ig VJC genes (no D)", "TCR VC genes (no JD)"), col=c("red", "dodgerblue", "gold", "forestgreen"), pch=19)
   else
      legend("topleft", legend=c("Upregulated in both", "Downregulated in both"), col=c("red", "dodgerblue"), pch=19)
 
   abline(h=0, lty=5)
   abline(v=0, lty=5)
   dev.off()
}

## Kilonova plot
de.kilo <- de.fishers
de.kilo$Effect1 <- de.kilo$George_Effect
de.kilo$Effect2 <- de.kilo$CCLE_Effect

plot.main <- "Mutual fold change between primary tumors and cancer cell lines"
plot.de <- paste0(wd.ccle, "kilonovaplot_all-ccle_fishers")









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
