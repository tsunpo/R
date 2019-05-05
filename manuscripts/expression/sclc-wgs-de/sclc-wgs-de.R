# =============================================================================
# Manuscript   : 
# Chapter I    : 
# Figure(s)    : 
# Name         : manuscripts/expression/cll-wgs-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/03/19
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

samples <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
samples$SAMPLE_ID <- paste0(samples$SAMPLE_ID, "_T")
rownames(samples) <- samples$SAMPLE_ID
samples$RT <- as.factor(samples$RT)

## Expression level (11/04/19)
wd.tpm      <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.tpm.data <- file.path(wd.tpm, "data")
load(file.path(wd.tpm.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
expressed <- rownames(tpm.gene)
# > length(expressed)
# [1] 19131

load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r30p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01
# > dim(tpm.gene.log2)
# [1] 33357   101

overlaps <- intersect(rownames(tpm.gene.log2), expressed)
# > length(overlaps)
# [1] 18917
tpm.gene.log2 <- tpm.gene.log2[overlaps,]

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=69-15NA, 20 RB1 vs 34 WT)
# Last Modified: 22/05/18
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : RB1_MUT (1) vs RB1_WT (0) as factor
argv      <- data.frame(predictor="RT", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101")
file.main <- paste0("RT (n=50) vs WT (n=51) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

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
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   #xmax <- 1
   ymax <- max(de$log10P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="SCLC M1/M0 log2 fold change", ylab="-log10(p-value)", col="darkgray", main=file.main[1])#, xaxt="n")
   #axis(side=1, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1))

   abline(h=c(-log10(pvalue)), lty=5)
   #text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.05", cex=0.85)
   text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, "FDR=0.075", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

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
   legend("topright", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "Differential read depth between SCLC T29 and T33"
plot.de <- file.path(wd.de.plots, "volcanoplot-r30p47-r5p47_sclc_rt_p1e-6")

## Chr2
genes <- readTable(paste0(plot.de, "_chr2.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "chr2:74.3-85.9Mb")
#file.de <- paste0(plot.de, "_chr2_log2FC1.pdf")
file.de <- paste0(plot.de, "_chr2.pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main)

## Natural killer cell receptor phenotypes
genes <- readTable(paste0(plot.de, "_nk.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "Natural killer cell receptors")
#file.de <- paste0(plot.de, "_nk_log2FC1.pdf")
file.de <- paste0(plot.de, "_nk.pdf")
plotVolcano(de.tpm.gene, 0.001, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_sclc_tpm-gene-r30p47-r5p47_rt_wilcox_q_n101")
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
