# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/hela-tpm-de-pca.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 09/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R")
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

samples <- readTable(file.path(wd.rna, "hela_rna_n14.txt"), header=T, rownames=T, sep="\t")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene <- tpm.gene[,rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# PCA (on G1-S/G2-M genes)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "genes_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "genes_G2-M.list"), header=F, rownames=F, sep="")
# > length(genes.G1S)
# [1] 304
# > length(genes.G2M)
# [1] 876

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.log2))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.log2))
# > length(genes.G1S)
# [1] 279
# > length(genes.G2M)
# [1] 838

test <- tpm.gene.log2[genes.G1S,]
#test <- tpm.gene.log2[genes.G2M,]
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

file.main <- "PCA of HeLa (n=14) on 279/304 G1-S genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_279G1-S_PC1-PC2", file.main, "bottomright", c("purple3", "forestgreen", "gold"), rownames(samples))

file.main <- "PCA of HeLa (n=14) on 838/876 G2-M genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_838G2-M_PC1-PC2", file.main, "bottomright", c("purple3", "forestgreen", "gold"), rownames(samples))


# -----------------------------------------------------------------------------
# PCA (on LCNEC RB1 D.E. genes; Q < 0.1)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm_gene_r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1 <- rownames(subset(de.tpm.gene, FDR <= 0.1))
#genes.rb1 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.overlaps <- intersect(genes.rb1, rownames(tpm.gene.log2))
## > length(genes.rb1.overlaps)
## [1] 578
## > length(genes.rb1.overlaps)
## [1] 135

test <- tpm.gene.log2[genes.rb1.overlaps, rownames(samples)]
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

file.main <- "HeLa (n=14) on 578/639 D.E. (LCNEC RB1 vs WT; Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_RB1_Q0.1_578DE_PC1-PC2", file.main, "bottomright", c("purple3", "forestgreen", "gold"), rownames(samples))

file.main <- "HeLa (n=14) on 135/145 D.E. (RB1 vs WT; Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_RB1_Q0.05_135DE_PC1-PC2", file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples))

# -----------------------------------------------------------------------------
# D.E. using non-parametric test: 5 First cycle vs 4 Second cycle (on only 578 D.E. in RB1-loss LCNEC)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm_gene_r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1 <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.overlaps <- intersect(genes.rb1, rownames(tpm.gene.log2))
## > length(genes.rb1.overlaps)
## [1] 578

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene <- tpm.gene[genes.rb1.overlaps, rownames(samples)]
tpm.gene.log2 <- log2(tpm.gene + 0.01)

## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RB1_MUT (1) vs RB1_WT(0) as factor
argv      <- data.frame(predictor="CYCLE_DE", predictor.wt=0, test="Wilcox", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-578de_cycle2_wilcox_q_n10")
file.main <- paste0("First cycle (n=5) vs Second cycle (n=5) in ", BASE)

de.tpm.gene <- pipeDE(tpm.gene.log2, samples, argv, ensGene)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")












# -----------------------------------------------------------------------------
# PCA (on LCNEC RB1 D.E. genes; Q < 0.05)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
genes.rb1 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.overlaps <- intersect(genes.rb1, rownames(tpm.gene.log2))
## > length(genes.rb1.overlaps)
## [1] 135

test <- tpm.gene.log2[genes.rb1.overlaps, rownames(samples)]
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

