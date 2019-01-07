# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/lcnec-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/01/19
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
BASE <- "LCNEC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.path  <- file.path(wd.de, "pathway")

samples <- readTable(file.path(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "AVG_FRAGMENT_LENGTH", "MAX_INSERT_SIZE", "RB1_MUT")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## Use pseudocount=0.01

# -----------------------------------------------------------------------------
# D.E. using non-parametric test (20 RB1 vs 34 WT; n=69-15NA)
# Last Modified: 22/05/18
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RB1_MUT (1) vs RB1_WT(0) as factor
argv      <- data.frame(predictor="RB1_MUT", predictor.wt=0, test="Wilcox", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47_rb1_wilcox_q_n54")
file.main <- paste0("RB1 MUT (n=20) vs WT (n=34) in ", BASE)

de.tpm.gene <- pipeDE(tpm.gene.log2, samples, argv, ensGene)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
#load(file.path("/Users/tpyang/Work/uni-koeln/tyang2", "LCNEC", "analysis/expression/kallisto", paste0("lcnec", "-tpm-de/data/", "de_", "lcnec", "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
## > length(genes.rb1.q0.05)
## [1] 145

## RB1 status on D.E genes
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- as.numeric(samples[,"RB1_MUT"])
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

##
file.main <- c("LCNEC samples on top 145 D.E. genes", "RB1-loss effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_rb1_q0.05_145de", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)

# -----------------------------------------------------------------------------
# Volcano plots
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
plotVolcano <- function(de, fdr, genes, file.de, file.main) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="RB1/WT fold change [log2]", ylab="P-value [-log10]", col="darkgray", main=file.main[1])

   abline(h=c(-log10(fdrToP(fdr, de))), lty=5)
   text(xmax*-1 + 2*xmax/36, -log10(fdrToP(fdr, de)) + ymax/42, "FDR=0.05", cex=0.85)
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
   legend("topleft", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "RB1-loss differential expression in LCNEC"
plot.de <- file.path(wd.de.plots, "volcanoplot_lcnec_rb1_q0.05")

## Figure 1 (Cell cycle)
genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_rb1_q0.05_cycle.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Cell cycle regulation")
file.de <- paste0(plot.de, "_cycle.pdf")
plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)










## Figure 1b (Hintone)
#genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_rb1_q0.05_histone.tab"), header=T, rownames=F, sep="\t")
#file.main <- c(plot.main, "Histone modifications")
#file.de <- paste0(plot.de, "_histone.pdf")
#plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)

## Figure 1c (Repair)
#genes <- readTable(file.path(wd.de.plots, "volcanoplot_lcnec_rb1_q0.05_repair.tab"), header=T, rownames=F, sep="\t")
#file.main <- c(plot.main, "DNA repair pathway")
#file.de <- paste0(plot.de, "_repair.pdf")
#plotVolcano(de.tpm.gene, 0.05, genes, file.de, file.main)

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
