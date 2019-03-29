# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss additive effect on gene expression among LUAD, LCNEC and SCLC
# Figure(s)    : Figure 1 (D and E), plus Figure 1 (C)
# Name         : manuscripts/expression/luad+lcnec+sclc-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 09/01/19
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
BASE <- "LUAD+LCNEC+SCLC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "pheno.all_n198.txt"), header=T, rownames=T, sep="")
samples$RB1[which(is.na(samples$RB1))] <- "NA"
samples$Cancer_Type <- as.factor(samples$Cancer_Type)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
#  among LUAD, LCNEC and SCLC
# Last Modified: 08/01/19; 17/08/17
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
# Volcano plots of RB1-loss AE genes amongst lung cancers (LUAD, LCNEC and SCLC)
# Figure(s)    : Figure 1 (D)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
rhoToP <- function(rho, de) {
   de.rho <- rbind(subset(de, RHO >= rho), subset(de, RHO <= rho*-1))
 
   return(max(de.rho$P))
}

plotVolcano <- function(de, rho, genes, file.de, file.main) {
   pvalue <- rhoToP(rho, de)  
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   ymax <- max(de$log10P)

   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="SCLC/LUAD log2 fold change", ylab="-log10(p-value)", col="darkgray", main=file.main[1])
 
   abline(h=c(-log10(pvalue)), lty=5)
   text(xmax*-1 + 2*xmax/40, -log10(pvalue) + ymax/42, paste0("rho=±", rho), cex=0.85)

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
 
   mtext(paste0("(", file.main[2], ")"), cex=1.2, line=0.3)
   legend("topleft", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "RB1-loss additive expression among LUAD, LCNEC and SCLC"
plot.de <- file.path(wd.de.plots, "volcanoplot_luad+lcnec+sclc_rb1_rho0.7")
de.tpm.gene$LOG2_FC <- de.tpm.gene$SCLC_LUAD

## Cell cycle regulation
genes <- readTable(file.path(wd.de.plots, "volcanoplot_luad+lcnec+sclc_rb1_rho0.7_cycle.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Cell cycle regulation")
file.de <- paste0(plot.de, "_cycle.pdf")
plotVolcano(de.tpm.gene, 0.7, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Principal component analysis (PCA) of LUAD, LCNEC and SCLC samples on RB1-loss AE genes
# Figure(s)    : Figure 1 (E)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
genes.rb1.rho0.7 <- rownames(subset(de.tpm.gene, P <= rhoToP(0.7, de.tpm.gene)))
# > length(genes.rb1.rho0.7)
# [1] 203

## Figure 1 (E)
test <- tpm.gene[genes.rb1.rho0.7, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$RB1
trait[which(samples$Cancer_Type == 0)] <- "LUAD"
trait[which(samples$Cancer_Type == 2)] <- "SCLC"
trait[which(trait == "RB1")] <- "LCNEC (RB1)"
trait[which(trait == "WT")]  <- "LCNEC (WT)"
trait[which(trait == "NA")]  <- "LCNEC (NA)"

file.main <- c("LUAD, LCNEC and SCLC samples on top 203 AE genes", "RB1-loss additive effect; absolute rho > 0.7")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_luad+lcenc+sclc_rb1_rho0.7_ae203", size=6.5, file.main, "topleft", c("lightgray", "red", "dodgerblue", "yellowgreen", "orange"), NULL, flip.x=-1, flip.y=-1)

##
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

file.main <- c("LUAD, LCNEC and SCLC samples on top 145 DE genes", "RB1-loss differential effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_luad+lcenc+sclc_rb1_q0.05_de145", size=6.5, file.main, "topright", c("lightgray", "red", "dodgerblue", "yellowgreen", "orange"), NULL, flip.x=1, flip.y=-1)

# -----------------------------------------------------------------------------
# PCA of LUAD, LCNEC and SCLC samples on RB1-loss DE genes (FDR < 0.05)
# Figure(s)    : Figure 1 (C)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
## > length(genes.rb1.q0.05)
## [1] 145   ## Original 145 from lcnec-tpm-de-pca.R

## Keep only LCNEC and SCLC samples
samples <- subset(samples, Cancer_Type != 0)   ## MOVE HERE 08/01/19

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$RB1
trait[which(samples$Cancer_Type == 2)] <- "SCLC"
trait[which(trait == "RB1")] <- "LCNEC (RB1)"
trait[which(trait == "WT")]  <- "LCNEC (WT)"
trait[which(trait == "NA")]  <- "LCNEC (NA)"

## Figure 1 (E)
file.main <- c("LCNEC and SCLC samples on top 145 DE genes", "RB1-loss differential effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcenc+sclc_rb1_q0.05_de145", size=6.5, file.main, "bottomright", c("lightgray", "red", "dodgerblue", "orange"), NULL, flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on SCLC/LUAD ranked gene lists
# Figure(s)    : Figure S1 (C and D)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_luad+lcnec+sclc_tpm-gene-r5p47_rb1_src_q_n198")
#de.tpm.gene$LOG2_FC <- de.tpm.gene$SCLC_LUAD
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
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 279/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 834/876

writeGRPformat(periodic.G1S, wd.de.gsea, "G1-S")
writeGRPformat(periodic.G2M, wd.de.gsea, "G2-M")

# -----------------------------------------------------------------------------
# Expression levels amongst lung cancers (LUAD, LCNEC and SCLC)
# Last Modified: 10/12/18
# -----------------------------------------------------------------------------
plotBox <- function(gene, wd.de.plots, tpm.gene.log2, pheno.all, ylim=NULL) {
   ids <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
   
   for (i in 1:length(ids)) {
      id <- ids[i]

      gene.tpm <- cbind(t(tpm.gene.log2)[rownames(pheno.all), id], pheno.all)
      colnames(gene.tpm)[1] <- "MEDIAN"
      file.name <- paste0("boxplot_tpm.gene.log2_", gene)
      if (is.null(ylim))
         ylim.id <- c(min(gene.tpm$MEDIAN), max(gene.tpm$MEDIAN))
      else
         ylim.id <- ylim
      
      if (length(ids) != 1)
         file.name <- paste0(file.name, "_", id)
      
      pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=4)
      boxplot(MEDIAN ~ Cancer_Type, data=gene.tpm, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=ylim.id, ylab="log2(TPM+0.01)", main=paste0(gene, " (", id, ")"))
      dev.off()
   }
}

##
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

#genes <- c("B2M", "BAX", "NUCB1", "RBL1")
#genes <- c("PATL2", "RBL1")
#genes <- c("MCM2", "CDKN2A", "BARD1")
genes <- c("SPG11", "TRIM69", "DHDH", "NUCB1")
for (g in 1:length(genes))
   plotBox(genes[g], wd.de.plots, tpm.gene.log2, samples)

# -----------------------------------------------------------------------------
# Expression levels amongst lung cancers (LUAD, LCNEC and SCLC)
# Last Modified: 10/12/18
# -----------------------------------------------------------------------------
plotBoxWithNormal <- function(gene, wd.de.plots, tpm.gene.log2, pheno.all, ylim=NULL) {
   ids <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
 
   for (i in 1:length(ids)) {
      id <- ids[i]
  
      gene.tpm <- cbind(t(tpm.gene.log2)[rownames(pheno.all), id], pheno.all)
      colnames(gene.tpm)[1] <- "MEDIAN"
      file.name <- paste0("boxplot_tpm.gene.log2_", gene)
      if (is.null(ylim))
         ylim.id <- c(min(gene.tpm$MEDIAN), max(gene.tpm$MEDIAN))
      else
         ylim.id <- ylim
  
      if (length(ids) != 1)
         file.name <- paste0(file.name, "_", id)
  
      pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=4.5)
      boxplot(MEDIAN ~ Cancer_Type, data=gene.tpm, outline=T, names=c("Lung", "LUAD", "LCNEC", "SCLC"), ylim=ylim.id, ylab="log2(TPM+0.01)", main=paste0(gene, " (", id, ")"))
      dev.off()
   }
}

load("/Users/tpyang/Work/uni-koeln/tyang2/LUSC/analysis/expression/kallisto/normal-tpm-de/data/lusc_kallisto_0.43.1_tpm.gene_r5p47.RData")
# > dim(tpm.gene)
# [1] 19872    41
tpm.gene.normal <- tpm.gene
tpm.gene.normal.log2 <- log2(tpm.gene.normal + 0.01)

genes <- c("KIAA1456")
for (g in 1:length(genes))
   plotBoxWithNormal(genes[g], wd.de.plots, tpm.gene.log2, samples)

