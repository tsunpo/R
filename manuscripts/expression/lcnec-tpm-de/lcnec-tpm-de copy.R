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
handbooks  <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "LCNEC"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "AVG_FRAGMENT_LENGTH", "MAX_INSERT_SIZE", "RB1_MUT")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

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
# Volcano plots Lite
# Last Modified: 25/10/17
# -----------------------------------------------------------------------------
plotVolcanoLite <- function(de, fdr, genes, file.de, file.main) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$Effect)
   ymax <- max(de$log10P)
   p <- max(de.sig$P)
 
   pdf(file.de, height=7, width=7)
   plot(de$Effect, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="Median fold change (log2 RB1/WT)", ylab="Significance (-log10 P-value)", col="darkgray", main=file.main)
   abline(h=c(-log10(fdrToP(0.05, de))), lty=5, col="darkgray")
   text(xmax*-1 + 2*xmax/29, -log10(fdrToP(0.05, de)) + ymax/42, "FDR=5%", col="darkgray")
 
   de.up <- subset(de.sig, Effect > 0)
   points(de.up$Effect, de.up$log10P, pch=16, col="red")
 
   de.down <- subset(de.sig, Effect < 0)
   points(de.down$Effect, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])
      if (nrow(gene) > 0) {
         points(gene$Effect, gene$log10P, pch=1, col="black")
    
         if (genes[g] == "XRCC2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.10, cex=0.75)
         else if (genes[g] == "RAD52")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.11, cex=0.75)
         else if (genes[g] == "CDK6" || genes[g] == "E2F4" || genes[g] == "MDM2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=1.13, cex=0.75)
         else if (genes[g] == "RNASEH2A")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.05, cex=0.75)
         else if (genes[g] == "MAD2L2")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         #else if (genes[g] == "RBL2")
         #   text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.09, cex=0.75)
         else if (genes[g] == "PARP1" || genes[g] == "XRCC4" || genes[g] == "RNF138" || genes[g] == "BRCA2" || genes[g] == "ASF1A" || genes[g] == "HMGB1" || genes[g] == "CHAF1A")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.1, cex=0.75)
         else if (genes[g] == "PSIP1" || genes[g] == "RBBP4")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.12, cex=0.75)
         else if (genes[g] == "E2F8")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.14, cex=0.75)
         else if (genes[g] == "LIG1")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=-0.16, cex=0.75)
         else if (genes[g] == "HMGB2" )
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.07, 1.3), cex=0.75)
         else if (genes[g] == "CREBBP")
            text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(-0.05, 1.3), cex=0.75)
         else
            if (gene$Effect > 0)
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$Effect, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
 
   abline(h=c(-log10(p)), lty=5)
   text(xmax*-1 + 2*xmax/24, -log10(p) + ymax/42, "FDR=10%")
   legend("topleft", legend=c("Up-regulation", "Down-regulation"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

plot.main <- "RB1-loss differential expression in LCNEC"
plot.de    <- file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05")
plot.genes <- file.path(wd.de.plots, "volcanoplot_lcnec_RB1_Q0.05_Cycle.tab")

###
## Figure 1 (Cell cycle)
genes <- c("TP53", "MDM2", "RB1", "E2F1", "E2F7", "E2F8", "E2F4", "E2F5", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN1B", "CDKN2C", "BARD1", "CDC7", "CDK2", "CCNE2", "CCND1", "POLD4", "TOPBP1", "CHEK2", "ATR", "ATM")
genes <- c("STMN1", "FBXO5", "TUBA1A", "TMPO", "MAD2L2", genes)
file.main <- c(plot.main, "", "Cell cycle control", "")
file.de <- paste0(plot.de, "_Figure_1_Cycle2.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)

## Figure 1
genes <- c("XRCC1", "XRCC2", "XRCC3", "XRCC4", "NHEJ1", "FEN1", "LIG1", "LIG4", "TP53", "MDM2", "TP53BP1", "PARP1", "RB1", "RBL1", "RBL2", "E2F1", "E2F7", "E2F8", "E2F4", "E2F5", "CDK6", "CDK4", "MCM2", "RNASEH2A", "CDKN2A", "CDKN1A", "CDKN1B", "BARD1", "CDC7", "CDK2", "CCNE2", "CCND1", "POLD4", "TOPBP1", "CHEK2", "BRCA1", "BRCA2", "PSIP1", "RBBP8", "RNF168", "BARD1")
file.main <- c(plot.main, "", "DSB repair pathway choice & Cell cycle control", "")
file.de <- paste0(plot.de, "_Figure_1.pdf")
plotVolcanoLite(de.lcnec, 0.1, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results (Up- and Down-regulation)
# -----------------------------------------------------------------------------
wd.de.data.reactome <- file.path(wd.de.data, "pathway_q0.05_up")
#wd.de.data.reactome <- file.path(wd.de.data, "pathway_q0.05_down")

list <- ensGene[,c("ensembl_gene_id",	"external_gene_name")]

reactome <- read.csv(file.path(wd.de.data.reactome, "result.csv"))
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
writeTable(reactome, file.path(wd.de.data.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")
