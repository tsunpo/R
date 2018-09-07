# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/hela-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 09/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Asymmetry.R", "DifferentialExpression.R", "ReplicationOrigin.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Human Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "HeLa"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata/Dominguez 2016")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
wd.de.path  <- file.path(wd.de, "pathway")

samples <- readTable(file.path(wd.rna, "hela_rna_n14.txt"), header=T, rownames=T, sep="\t")
rownames(samples) <- gsub("D", "T", rownames(samples))
samples$SAMPLE_ID <- gsub("D", "T", samples$SAMPLE_ID)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
colnames(tpm.gene) <-  gsub("D", "T", colnames(tpm.gene))
tpm.gene <- tpm.gene[,rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# PCA (on G1-S/G2-M genes)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "genes_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "genes_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.log2))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.log2))
# > length(genes.G1S)
# [1] 279   ## Original 304 from Dominguez et al.
# > length(genes.G2M)
# [1] 838   ## Original 876 from Dominguez et al.

test <- tpm.gene[genes.G1S,]
#test <- tpm.gene[genes.G2M,]
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

##
file.main <- "PCA of HeLa (n=14) on 279/304 G1-S genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_279G1-S", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=1)

file.main <- "PCA of HeLa (n=14) on 838/876 G2-M genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_838G2-M", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# PCA (on LCNEC RB1 D.E. genes; Q < 0.05/0.1)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.q0.1  <- intersect(genes.rb1.q0.1, rownames(tpm.gene))   ## Do NOT use filtered tpm.gene
genes.rb1.q0.05 <- intersect(genes.rb1.q0.05, rownames(tpm.gene))
## > length(genes.rb1.q0.1)
## [1] 578   ## Original 639 from lcnec-tpm-de-pca.R
## > length(genes.rb1.q0.05)
## [1] 135   ## Original 145 from lcnec-tpm-de-pca.R

## HeLa on RB1 D.E genes
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene[genes.rb1.q0.1, rownames(samples)]
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

##
file.main <- "HeLa (n=14) on 135/145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.05_135DE", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=1)

file.main <- "HeLa (n=14) on 578/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.1_578DE", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=-1)

# -----------------------------------------------------------------------------
# Gene length
# Last Modified: 10/08/18
# -----------------------------------------------------------------------------
initLength <- function(genes, group) {
   ens.genes <- ensGene[genes,]
   ens.genes$Length <- mapply(x = 1:nrow(ens.genes), function(x) getLengthTx(genes[x]))
   ens.genes$Group  <- group
   
   return(ens.genes)
}

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
ens.genes.ALL <- initLength(rownames(tpm.gene), 2)   ## ADD 17/08/18
ens.genes.G1S <- initLength(genes.G1S, 3)
ens.genes.G2M <- initLength(genes.G2M, 4)
ens.genes <- rbind(ens.genes.rb1.q0.1, ens.genes.rb1.q0.05, ens.genes.ALL, ens.genes.G1S, ens.genes.G2M)
ens.genes$Group <- as.factor(ens.genes$Group)
ens.genes$Length <- log10(ens.genes$Length)

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_RB1-G1S-G2M-ALL_length.pdf"))
pdf(file.name, height=6, width=4.5)
ymin <- min(ens.genes$Length)
ymax <- max(ens.genes$Length)
boxplot(Length ~ Group, data=ens.genes, outline=T, names=c("RB1", "RB1", "All", "G1-S", "G2-M"), col=c("white", "white", "gray", "gray", "gray"), ylim=c(ymin, ymax), ylab="Gene length (log10)", main=c("LCNEC RB1 D.E. and HeLa cell-cycle", "gene lists"))
dev.off()

## G1-S vs G2-M
testW(ens.genes.G1S$Length, ens.genes.ALL$Length)
# [1] 0.831461
testW(ens.genes.G2M$Length, ens.genes.ALL$Length)
# [1] 9.266777e-32
testW(ens.genes.G1S$Length, ens.genes.G2M$Length)
# [1] 1.012973e-11

## RB1 (q<0.1) vs G1-S and G2-M
testW(ens.genes.rb1.q0.1$Length, ens.genes.ALL$Length)
# [1] 9.892606e-05
testW(ens.genes.rb1.q0.1$Length, ens.genes.G1S$Length)
# [1] 3.853642e-03
testW(ens.genes.rb1.q0.1$Length, ens.genes.G2M$Length)
# [1] 5.253057e-08

## RB1 (q<0.05) vs G1-S and G2-M
testW(ens.genes.rb1.q0.05$Length, ens.genes.ALL$Length)
# [1] 0.53314
testW(ens.genes.rb1.q0.05$Length, ens.genes.G1S$Length)
# [1] 0.3546051
testW(ens.genes.rb1.q0.05$Length, ens.genes.G2M$Length)
# [1] 4.708969e-06

# -----------------------------------------------------------------------------
# D.E. using non-parametric test: 5 First cycle vs 4 Second cycle (on D.E. LCNCE RB1 genes)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RB1_MUT (1) vs RB1_WT(0) as factor
samples$CYCLE_DE[13] <- NA   ## Preclude T13 which is similar to G1-S

argv      <- data.frame(predictor="CYCLE_DE", predictor.wt=0, test="Wilcox", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-135rb1_cycle2_wilcox_q_n9")
#file.name <- paste0("de_", base, "_tpm-gene-578rb1_cycle2_wilcox_q_n9")

de.tpm.gene <- pipeDE(tpm.gene.log2[genes.rb1.q0.05, rownames(samples)], samples, argv, ensGene)
#de.tpm.gene <- pipeDE(tpm.gene.log2[genes.rb1.q0.1, rownames(samples)], samples, argv, ensGene)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
