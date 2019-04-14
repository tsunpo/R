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
BASE <- "CLL"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-wgs-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "cll_rna_n93.txt"), header=T, rownames=T, sep="")
samples <- subset(samples, RT != 0)
samples$RT <- as.factor(samples$RT)

## Expression level (07/04/19)
wd.tpm    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.tpm.data  <- file.path(wd.tpm, "data")
load(file.path(wd.tpm.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
expressed <- rownames(tpm.gene)
# > length(expressed)
# [1] 18502

## Read depth (07/04/19)
load(file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r30p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01
# > dim(tpm.gene.log2)
# [1] 32974    62
# [1] 28364    62   ## min_reads=30, min_prop=1.00

overlaps <- intersect(rownames(tpm.gene.log2), expressed)
# > length(overlaps)
# [1] 18163
# [1] 16784   ## min_reads=30, min_prop=1.00
tpm.gene.log2 <- tpm.gene.log2[overlaps,]

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; n=69-15NA, 20 RB1 vs 34 WT)
# Last Modified: 22/05/18
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : RB1_MUT (1) vs RB1_WT (0) as factor
argv      <- data.frame(predictor="RT", predictor.wt=-1, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
#file.name <- paste0("de_", base, "_tpm-gene-r30p100_rt_wilcox_q_n62")
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n62")
file.main <- paste0("RT (n=29) vs WT (n=33) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure 1 (A)
# Last Modified: 08/04/19; 07/01/19
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, pvalue, genes, file.de, file.main) {
   #pvalue <- fdrToP(fdr, de)
   de.sig <- subset(de, FISHERS_P <= pvalue)
   de.sig$log10P <- -log10(de.sig$FISHERS_P)
 
   de$log10P <- -log10(de$FISHERS_P)
   #xmax <- max(de$LOG2_FC)
   xmax <- 1
   ymax <- max(de$log10P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="CLL T29/T33 log2 fold change", ylab="-log10(Fishers' combined p-value)", col="darkgray", main=file.main[1], xaxt="n")
   axis(side=1, at=seq(-1, 1, by=0.5), labels=c(-1, -0.5, 0, 0.5, 1))
   
   abline(h=c(-log10(pvalue)), lty=5)
   #text(xmax*-1 + 2*xmax/28, -log10(pvalue) + ymax/42, "FDR=0.002", cex=0.85)
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

   mtext(file.main[2], cex=1.2, line=0.3)
   legend("topright", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "CLL and SCLC combined differential read depth analyses"
plot.de <- file.path(wd.de.plots, "volcanoplot-r30p47-r5p47_cll_rt_fishers_p1e4")

## Chr2
genes <- readTable(paste0(plot.de, "_chr2.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "chr2:74.3-85.9Mb")
file.de <- paste0(plot.de, "_chr2_log2FC1.pdf")
#file.de <- paste0(plot.de, "_chr2.pdf")
plotVolcano(de.tpm.gene, 0.0001, genes, file.de, file.main)

## Natural killer cell receptor phenotypes
genes <- readTable(paste0(plot.de, "_nk.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "Natural killer cell receptors")
file.de <- paste0(plot.de, "_nk_log2FC1.pdf")
#file.de <- paste0(plot.de, "_nk.pdf")
plotVolcano(de.tpm.gene, 0.0001, genes, file.de, file.main)

##
plot.main <- "CLL and SCLC combined differential read depth analyses"
plot.de <- file.path(wd.de.plots, "volcanoplot-r30p47-r5p47_cll_rt_fishers_p1e6")

## EEF1A1
genes <- readTable(paste0(plot.de, "_EEF1A1.tab"), header=T, rownames=F, sep="\t")
rownames(genes) <- genes$GENE
genes <- genes[intersect(genes$GENE, de.tpm.gene$external_gene_name),]

file.main <- c(plot.main, "")
file.de <- paste0(plot.de, "_EEF1A1_log2FC1.pdf")
#file.de <- paste0(plot.de, "_EEF1A1.pdf")
plotVolcano(de.tpm.gene, 0.000001, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_cll_tpm-gene-r30p47-r5p47_rt_wilcox_q_n62")
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
# Test
# Last Modified: 12/04/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/CLL/analysis/expression/kallisto/cll-wgs-de/data/de_cll_tpm-gene-r30p47-r5p47_rt_wilcox_q_n62.RData")
de.tpm.gene.cll <- de.tpm.gene

load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-wgs-de/data/de_sclc_tpm-gene-r30p47-r5p47_rt_wilcox_q_n77.RData")
de.tpm.gene.sclc <- de.tpm.gene

de.tpm.gene.cll.down  <- subset(de.tpm.gene.cll, LOG2_FC < 0)
de.tpm.gene.sclc.down <- subset(de.tpm.gene.sclc, LOG2_FC < 0)
de.tpm.gene.cll.down.p  <- subset(de.tpm.gene.cll.down, P <= 0.001)
de.tpm.gene.sclc.down.p <- subset(de.tpm.gene.sclc.down, P <= 0.001)
intersect(rownames(de.tpm.gene.cll.down.p), rownames(de.tpm.gene.sclc.down.p))
# [1] "ENSG00000121351" "ENSG00000139116" "ENSG00000139218" "ENSG00000156508"
# > de.tpm.gene.cll.down.p[c("ENSG00000121351", "ENSG00000139116", "ENSG00000139218", "ENSG00000156508"),]
# ensembl_gene_id external_gene_name chromosome_name strand start_position end_position   gene_biotype            P        FDR    RT_WT       RT     LOG2_FC
# ENSG00000121351 ENSG00000121351               IAPP           chr12      1       21507893     21532912 protein_coding 5.463067e-05 0.01285717 5.015109 4.786841 -0.22826819
# ENSG00000139116 ENSG00000139116             KIF21A           chr12     -1       39687030     39837192 protein_coding 1.653357e-04 0.01320173 6.432993 6.103671 -0.32932176
# ENSG00000139218 ENSG00000139218             SCAF11           chr12     -1       46312914     46385903 protein_coding 2.312490e-04 0.01320173 7.805038 7.742476 -0.06256168
# ENSG00000156508 ENSG00000156508             EEF1A1            chr6     -1       74225473     74233520 protein_coding 6.738102e-04 0.01320173 5.926873 5.836459 -0.09041426

# -----------------------------------------------------------------------------
# Test
# Last Modified: 12/04/19
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene.cll), rownames(de.tpm.gene.sclc))
de.tpm.gene.cll.o  <- de.tpm.gene.cll[overlaps,]
de.tpm.gene.sclc.o <- de.tpm.gene.sclc[overlaps,]

de.tpm.gene.cll.o$FISHERS_P <- fishers(de.tpm.gene.cll.o$P, de.tpm.gene.sclc.o$P)
library(qvalue)
de.tpm.gene.cll.o$FISHERS_FDR <- qvalue(de.tpm.gene.cll.o$FISHERS_P)$qvalue

de.tpm.gene.cll.o <- de.tpm.gene.cll.o[order(de.tpm.gene.cll.o$FISHERS_P),]
de.tpm.gene <- de.tpm.gene.cll.o
file.name <- paste0("de_", base, "_tpm-gene-r30p47-r5p47_rt_wilcox_q_n62_fishers")
save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")





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
