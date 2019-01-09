# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss additive effect on gene expression in lung cancers
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
wd.de.gsea  <- file.path(wd.de, "gsea")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "pheno.all_n198.txt"), header=T, rownames=T, sep="")
samples$RB1[which(is.na(samples$RB1))] <- "NA"
samples$Cancer_Type <- as.factor(samples$Cancer_Type)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# Find genes with "additive" effects between SCLC, LUAD and LCNEC
# Last Modified: 08/01/19; 17/08/17
# -----------------------------------------------------------------------------
tpm.gene.log2 <- tpm.gene.log2[,rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!

colnames <- c("SRC_RHO", "SRC_P", "SRC_Q", "LUAD", "LCNEC", "SCLC", "LCNEC_LUAD", "SCLC_LUAD")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$SRC_RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$RB1_RATE, method="spearman", exact=F)[[4]])
de$SRC_P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$RB1_RATE, method="spearman", exact=F)[[3]])
## ANOVA
#samples$Cancer_Type <- as.factor(pheno.all$Cancer_Type)
#de$ANOVA_P <- mapply(x = 1:nrow(expr.pheno.log2), function(x) testANOVA(x, expr.pheno.log2, pheno.all))

## FC
de$LUAD  <- median00(tpm.gene.log2, rownames(subset(samples, Cancer_Type == 0)))
de$LCNEC <- median00(tpm.gene.log2, rownames(subset(samples, Cancer_Type == 1)))
de$SCLC  <- median00(tpm.gene.log2, rownames(subset(samples, Cancer_Type == 2)))
de$LCNEC_LUAD <- de$LCNEC - de$LUAD
de$SCLC_LUAD  <- de$SCLC  - de$LUAD

library(qvalue)
de$SRC_Q   <- qvalue(de$SRC_P)$qvalue
#de$ANOVA_Q <- qvalue(de$ANOVA_P)$qvalue
de <- de[order(de$SRC_P),]

##
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_luad+lcnec+sclc_tpm-gene-r5p47_rb1_src_q_n198.RData"))
writeTable(de.tpm.gene, file.path(wd.de.data, "de_luad+lcnec+sclc_tpm-gene-r5p47_rb1_src_q_n198.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Principal component analysis (PCA)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
rhoToP <- function(rho, de) {
   de.rho <- rbind(subset(de, SRC_RHO >= rho), subset(de, SRC_RHO <= rho*-1))
 
   return(max(de.rho$SRC_P))
}

genes.rb1.rho0.7 <- rownames(subset(de.tpm.gene, SRC_P <= rhoToP(0.7, de.tpm.gene)))
# > length(genes.rb1.rho0.7)
# [1] 203

##
test <- tpm.gene[genes.rb1.rho0.7, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$RB1
trait[which(samples$Cancer_Type == 0)] <- "LUAD"
trait[which(samples$Cancer_Type == 2)] <- "SCLC"
trait[which(trait == "RB1")] <- "LCNEC (RB1)"
trait[which(trait == "WT")]  <- "LCNEC (WT)"
trait[which(trait == "NA")]  <- "LCNEC (NA)"

file.main <- c("LUAD, LCNEC and SCLC on top 203 A.E. genes", "RB1-loss additive effect; absolute rho > 0.7")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_luad+lcenc+sclc_rb1_rho0.7_ae203", size=6.5, file.main, "topleft", c("lightgray", "red", "dodgerblue", "yellowgreen", "orange"), NULL, flip.x=-1, flip.y=-1)

##
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

file.main <- c("LUAD, LCNEC and SCLC on top 145 D.E. genes", "RB1-loss differential effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_luad+lcenc+sclc_rb1_q0.05_de145", size=6.5, file.main, "topright", c("lightgray", "red", "dodgerblue", "yellowgreen", "orange"), NULL, flip.x=1, flip.y=-1)

# -----------------------------------------------------------------------------
# GSEA
# -----------------------------------------------------------------------------
file.name <- paste0("de_luad+lcnec+sclc_tpm-gene-r5p47_rb1_src_q_n198")
de.tpm.gene$LOG2_FC <- de.tpm.gene$SCLC_LUAD
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
# PCA (on RB1-loss differential effect genes; FDR < 0.05)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
## > length(genes.rb1.q0.05)
## [1] 145   ## Original 145 from lcnec-tpm-de-pca.R

## LCNEC+SCLC on RB1 D.E genes
samples <- subset(samples, Cancer_Type != 0)   ## MOVE HERE 08/01/19

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- samples$RB1
trait[which(samples$Cancer_Type == 2)] <- "SCLC"
trait[which(trait == "RB1")] <- "LCNEC (RB1)"
trait[which(trait == "WT")]  <- "LCNEC (WT)"
trait[which(trait == "NA")]  <- "LCNEC (NA)"

##
file.main <- c("LCNEC and SCLC on top 145 D.E. genes", "RB1-loss differential effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcenc+sclc_rb1_q0.05_de145", size=6.5, file.main, "bottomright", c("lightgray", "red", "dodgerblue", "orange"), NULL, flip.x=1, flip.y=1)

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
