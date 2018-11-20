# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/luad-lcnec-sclc-tpm-de-pca.R
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
BASE <- "LUAD-LCNEC-SCLC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.de    <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "pheno.all_n198.txt"), header=T, rownames=T, sep="")
samples$RB1[which(is.na(samples$RB1))] <- "NA"
samples$Cancer_Type <- as.factor(samples$Cancer_Type)
samples <- subset(samples, Cancer_Type != 0)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)

# -----------------------------------------------------------------------------
# PCA (on LCNEC RB1 D.E. genes; Q < 0.1)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.q0.1  <- intersect(genes.rb1.q0.1, rownames(tpm.gene))   ## Do NOT use filtered tpm.gene
genes.rb1.q0.05 <- intersect(genes.rb1.q0.05, rownames(tpm.gene))
## > length(genes.rb1.q0.1)
## [1] 632   ## Original 639 from lcnec-tpm-de-pca.R
## > length(genes.rb1.q0.05)
## [1] 145   ## Original 145 from lcnec-tpm-de-pca.R

## ALL on RB1 D.E genes
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene[genes.rb1.q0.1, rownames(samples)]
pca.de <- getPCA(t(test))

trait <- samples$RB1
trait[which(samples$Cancer_Type == 2)] <- "SCLC"
trait[which(trait == "RB1")] <- "LCNEC (RB1)"
trait[which(trait == "WT")]  <- "LCNEC (WT)"

##
file.main <- "NETs on 145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcenc-sclc_RB1_Q0.05_145DE", size=6.5, file.main, "bottomright", c("red", "dodgerblue", "gray", "orange"), NULL, flip.x=1, flip.y=1)

file.main <- "NETs on 632/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcnec-sclc_RB1_Q0.1_632DE", size=6.5, file.main, "bottomright", c("red", "dodgerblue", "gray", "orange"), NULL, flip.x=1, flip.y=1)
