# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/luad-lcnec-sclc-tpm-de-pca.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 08/08/18
# =============================================================================
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
tpm.gene.log2 <- log2(tpm.gene + 0.01)




# -----------------------------------------------------------------------------
# PCA (on LCNEC RB1 D.E. genes; FDR < 0.05)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.q0.05 <- intersect(genes.rb1.q0.05, rownames(tpm.gene))
## > length(genes.rb1.q0.05)
## [1] 145   ## Original 145 from lcnec-tpm-de-pca.R

## ALL on RB1 D.E genes
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene.log2 <- log2(tpm.gene + 0.01)   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$RB1
trait[which(samples$Cancer_Type == 2)] <- "SCLC"
trait[which(trait == "RB1")] <- "LCNEC (RB1)"
trait[which(trait == "WT")]  <- "LCNEC (WT)"

##
file.main <- c("LCNEC + SCLC samples on top 145 D.E. genes", "RB1-loss effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcenc+sclc_rb1_q0.05_145de", size=6.5, file.main, "bottomright", c("red", "dodgerblue", "gray", "orange"), NULL, flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# Expression levels between LUAD, LCNEC and SCLC
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
genes <- c("PATL2", "RBL1")
for (g in 1:length(genes))
   plotBox(genes[g], wd.de.plots, tpm.gene.log2, pheno.all)
