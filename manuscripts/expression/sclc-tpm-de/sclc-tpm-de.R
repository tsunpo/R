# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/lcnec-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 08/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
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
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "SCLC"
base <- tolower(BASE)
wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata/Tirosh 2016")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
wd.de.gsea  <- file.path(wd.de, "gsea")

samples <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
samples$RT <- as.factor(samples$RT)

##
samples.rna <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")
overlaps <- intersect(rownames(samples.rna), rownames(samples))
samples <- samples[overlaps,]
# > length(which(samples$RT == 1))
# [1] 36
# > length(which(samples$RT == 0))
# [1] 34

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# D.E. using non-parametric test (44 RT vs 10 WT; n=54)
# Last Modified: 04/04/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RT (44) vs WT (10) as factor
argv      <- data.frame(predictor="RT", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47_rt_wilcox_q_n70")
file.main <- paste0("RT (n=36) vs WT (n=34) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_sclc_tpm-gene-r5p47_rt_wilcox_q_n70")
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
# Meta-analysis on wgs.de first
# Last Modified: 24/04/19; 12/04/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-tpm-de/data/de_sclc_tpm-gene-r5p47_rt_wilcox_q_n70.RData")
de.tpm.gene.sclc <- de.tpm.gene

load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/de_nbl_tpm-gene-r5p47_rt_wilcox_q_n53.RData")
de.tpm.gene.nbl <- de.tpm.gene

#load("/Users/tpyang/Work/uni-koeln/tyang2/CLL/analysis/expression/kallisto/cll-wgs-de/data/de_cll_tpm-gene-r30p47-r5p47_rt_wilcox_q_n93.RData")
#de.tpm.gene.cll <- de.tpm.gene

overlaps <- intersect(rownames(de.tpm.gene.sclc), rownames(de.tpm.gene.nbl))
de.tpm.gene.sclc.o <- de.tpm.gene.sclc[overlaps,]
de.tpm.gene.nbl.o  <- de.tpm.gene.nbl[overlaps,]
# > length(overlaps)
# [1] 17958

idx1 <- which(de.tpm.gene.sclc.o$LOG2_FC * de.tpm.gene.nbl.o$LOG2_FC > 0)
de.tpm.gene.sclc.o <- de.tpm.gene.sclc.o[idx1,]
de.tpm.gene.nbl.o  <- de.tpm.gene.nbl.o[idx1,]
# > length(idx1)
# [1] 9947

library(qvalue)
de.tpm.gene.sclc.o$FISHERS_P2   <- fishers(de.tpm.gene.sclc.o$P, de.tpm.gene.nbl.o$P)
de.tpm.gene.sclc.o$FISHERS_FDR2 <- qvalue(de.tpm.gene.sclc.o$FISHERS_P2)$qvalue

de.tpm.gene.sclc.o <- de.tpm.gene.sclc.o[order(de.tpm.gene.sclc.o$FISHERS_P2),]
file.name <- paste0("de_", base, "_tpm-gene-r5p47_rt_wilcox_q_fishers_sclc+nbl")
save(de.tpm.gene.sclc.o, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene.sclc.o, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")







overlaps <- intersect(intersect(rownames(de.tpm.gene.sclc), rownames(de.tpm.gene.nbl)), rownames(de.tpm.gene.wgs.cll))
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

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 15/08/18
# -----------------------------------------------------------------------------
plotVolcano <- function(de, fdr, genes, file.de, file.main) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="RB1/WT fold change (log2)", ylab="P-value (-log10)", col="darkgray", main=file.main)

   abline(h=c(-log10(fdrToP(fdr, de))), lty=5)
   text(xmax*-1 + 2*xmax/36, -log10(fdrToP(fdr, de)) + ymax/42, "FDR=0.05", cex=0.85)
   abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

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
   
   legend("topleft", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "RB1-loss differential expression in LCNEC"
plot.de <- file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05")

## Figure 1 (Cell cycle)
genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05_Cycle.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "", "Cell cycle pathway", "")
file.de <- paste0(plot.de, "_Cycle.pdf")
plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)

## Figure 2 (Hintone)
genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05_Histone.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "", "Histone modifications", "")
file.de <- paste0(plot.de, "_Histone.pdf")
plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)

## Figure 3 (DDR)
genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05_Repair.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "", "DNA repair pathway", "")
file.de <- paste0(plot.de, "_Repair.pdf")
plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)

## Figure ? (UP/TCD/HO/Q4)
genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05_UP_TCD_HO_Q4.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "", "Upregulated/TCD/Head-on/Q4", "")
file.de <- paste0(plot.de, "_UP_TCD_HO_Q4.pdf")
plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Gene sets (please refer to guide-to-the/cycle.R) that are expressed in the dataset
# -----------------------------------------------------------------------------
load(file.path(wd.src.ref, "cycle.RData"))

## Tirosh et al 2016 
core.G1S <- intersect(core.G1S, rownames(tpm.gene.log2))   ## 41/43/43
core.G2M <- intersect(core.G2M, rownames(tpm.gene.log2))   ## 54/54/55
core.Stemness  <- intersect(core.Stemness, rownames(tpm.gene.log2))    ## 54/58/63

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 262/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 792/876










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

# -----------------------------------------------------------------------------
# GSEA
# -----------------------------------------------------------------------------
file.name <- paste0("de_", base, "_tpm-gene-r5p47_rb1_wilcox_q_n54")
writeRNKformat(de.tpm.gene, wd.de.data, file.name)

## GRP gene sets
writeGRPformat(genes.rb1.q0.05, wd.de.data, "genes.rb1.q0.05")
writeGRPformat(genes.rb1.q0.1, wd.de.data, "genes.rb1.q0.1")

## Tirosh 2016
writeGRPformat(core.G1S, wd.de.data, "core.G1-S")
writeGRPformat(core.G2M, wd.de.data, "core.G2-M")
writeGRPformat(core.SC, wd.de.data, "core.Stemness")

## Dominguez 2016
genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.log2))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.log2))
# > length(genes.G1S)
# [1] 277   ## Original 279/304 from Dominguez et al.
# > length(genes.G2M)
# [1] 830   ## Original 838/876 from Dominguez et al.
writeGRPformat(genes.G1S, wd.de.data, "genes.G1S")
writeGRPformat(genes.G2M, wd.de.data, "genes.G2M")

# -----------------------------------------------------------------------------
# Gene length
# Last Modified: 22/08/18
# -----------------------------------------------------------------------------
plotBox <- function(gene, wd.de, expr.pheno.log2, pheno.all) {
   ensembl_gene_id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id[1]   ## To avoid ENSG00000269846 (one of RBL1)
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/boxplot/boxplot_tpm.gene.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   dev.off()
}

ens.genes.rb1.q0.1  <- initLength(genes.rb1.q0.1, 0)
ens.genes.rb1.q0.05 <- initLength(genes.rb1.q0.05, 1)
ens.core.G1S <- initLength(core.G1S, 2)
ens.core.G2M <- initLength(core.G2M, 3)
ens.core.SC  <- initLength(core.SC, 4)
ens.genes <- rbind(ens.genes.rb1.q0.1, ens.genes.rb1.q0.05, ens.core.G1S, ens.core.G2M, ens.core.SC)
ens.genes$Group <- as.factor(ens.genes$Group)
ens.genes$Length <- log10(ens.genes$Length)

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_RB1-G1S-G2M-SC_length.pdf"))
pdf(file.name, height=6, width=4.5)
ymin <- min(ens.genes$Length)
ymax <- max(ens.genes$Length)
boxplot(Length ~ Group, data=ens.genes, outline=T, names=c("RB1", "RB1", "G1-S", "G2-M", "SC"), col=c("white", "white", "lightgray", "lightgray", "lightgray"), ylim=c(ymin, ymax), ylab="Gene length (log10)", main=c("LCNEC RB1 D.E. and Tirosh 2016", "gene lists"))
dev.off()

## G1-S vs G2-M
testW(ens.core.G1S$Length, ens.core.SC$Length)
# [1] 0.7494677
testW(ens.core.G2M$Length, ens.core.SC$Length)
# [1] 0.8034764
testW(ens.core.G1S$Length, ens.core.G2M$Length)
# [1] 0.9490751

## RB1 (q<0.1) vs G1-S and G2-M
testW(ens.genes.rb1.q0.1$Length, ens.core.SC$Length)
# [1] 0.8274065
testW(ens.genes.rb1.q0.1$Length, ens.core.G1S$Length)
# [1] 0.8365886
testW(ens.genes.rb1.q0.1$Length, ens.core.G2M$Length)
# [1] 0.4585982

## RB1 (q<0.05) vs G1-S and G2-M
testW(ens.genes.rb1.q0.05$Length, ens.core.SC$Length)
# [1] 0.8064687
testW(ens.genes.rb1.q0.05$Length, ens.core.G1S$Length)
# [1] 0.7674659
testW(ens.genes.rb1.q0.05$Length, ens.core.G2M$Length)
# [1] 0.9107349
