# =============================================================================
# Manuscript   : 
# Chapter I    : 
# Name         : manuscripts/expression/sclc-tpm-cycle.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/09/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Asymmetry.R", "CellCycle.R", "ReplicationOrigin.R")
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
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "AVG_FRAGMENT_LENGTH", "MAX_INSERT_SIZE")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=F, proteinCodingOnly=T, proteinCodingNonRedundantOnly=F)
tpm.gene.log2 <- getLog2andMedian(tpm.gene)

#plotDensityCount(tpm.gene.log2$MEDIAN, nrow(tpm.gene.log2), file.path(wd.de.data, paste0(base, "_kallisto_0.43.1_tpm.gene_r5p47_pcg.pdf")), NULL)

# -----------------------------------------------------------------------------
# Gene sets that are expressed in the dataset
# -----------------------------------------------------------------------------
load(file.path(wd.src.ref, "cycle.RData"))   ## See guide-to-the/cycle.R

## Tirosh et al 2016 
core.G1S <- intersect(core.G1S, rownames(tpm.gene.log2))   ## 41/43/43
core.G2M <- intersect(core.G2M, rownames(tpm.gene.log2))   ## 54/54/55

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 279/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 834/876

# -----------------------------------------------------------------------------
# Distribution of cell cycle genes in quantiles
# -----------------------------------------------------------------------------
tx.q4 <- getTxQ4(tpm.gene.log2, NA)
for (q in 1:4)
   print(length(tx.q4[[q]]))
# [1] 4783
# [1] 4783
# [1] 4782
# [1] 4783

tx.q4.cycle <- getTxQ4Cycle(tx.q4, core.G1S, core.G2M, periodic.G1S, periodic.G2M)
save(tx.q4.cycle, file=file.path(wd.de.plots, paste0(base, "_genes_tx_q4_cycle.RData")))

barplotTxQ4Cycle(wd.de.plots, base, BASE, tx.q4.cycle, beside=T)

# -----------------------------------------------------------------------------
# Gene length in quantiles
# -----------------------------------------------------------------------------
tx.q4.length <- getTxQ4Length(tx.q4)
save(tx.q4.length, file=file.path(wd.de.plots, paste0(base, "_genes_tx_q4_length.RData")))

boxplotTxQ4Length(wd.de.plots, base, BASE, tx.q4.length)

# -----------------------------------------------------------------------------
# Gene-to-gene minimum distance in quantiles
# -----------------------------------------------------------------------------
tx.q4.g2g <- getTxQ4G2G(tx.q4)
save(tx.q4.g2g, file=file.path(wd.de.plots, paste0(base, "_genes_tx_q4_g2g.RData")))

boxplotTxQ4G2G(wd.de.plots, base, BASE, tx.q4.g2g, samples)
plotTxQ4G2GDensity(wd.de.plots, base, BASE, tx.q4.g2g, samples, count=T, ymax=NULL)

testW(c(tx.q4.g2g[[1]], tx.q4.g2g[[2]], tx.q4.g2g[[3]]), tx.q4.g2g[[4]])
# [1] 3.72577e-60

# -----------------------------------------------------------------------------
# Gene length in ALL, G1/S and G2/M genes
# Last Modified: 09/09/18
# -----------------------------------------------------------------------------
genes.G1S <- unique(c(core.G1S, periodic.G1S))
genes.G2M <- unique(c(core.G2M, periodic.G2M))

ens.genes.ALL <- initLength(rownames(tpm.gene), 0)
ens.genes.G1S <- initLength(genes.G1S, 1)
ens.genes.G2M <- initLength(genes.G2M, 2)
ens.genes <- rbind(ens.genes.ALL, ens.genes.G1S, ens.genes.G2M)
ens.genes$Group <- as.factor(ens.genes$Group)
ens.genes$Length <- log10(ens.genes$Length)
save(ens.genes, file=file.path(wd.de.plots, paste0(base, "_genes_ALL-G1S-G2M_length.RData")))

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_ALL-G1S-G2M_length.pdf"))
pdf(file.name, height=6, width=3)
ymin <- min(ens.genes$Length)
ymax <- max(ens.genes$Length)
boxplot(Length ~ Group, data=ens.genes, outline=T, names=c("All", "G1/S", "G2/M"), col=c("white", "cornsilk", "darkgray"), ylim=c(ymin, ymax), ylab="Gene length (log10)", main=BASE)
dev.off()

## ALL vs G1-S vs G2-M
testW(ens.genes.G1S$Length, ens.genes.ALL$Length)
# [1] 0.5767282
testW(ens.genes.G2M$Length, ens.genes.ALL$Length)
# [1] 1.063058e-33
testW(ens.genes.G1S$Length, ens.genes.G2M$Length)
# [1] 2.836401e-11

# -----------------------------------------------------------------------------
# Gene-to-gene min distance in ALL, G1/S and G2/M genes
# Last Modified: 09/09/18
# -----------------------------------------------------------------------------
initG2G <- function(genes, group) {
   ens.genes <- ensGene[genes,]
   ens.genes$G2G <- as.numeric(getG2G(genes))
   ens.genes$Group  <- group

   return(ens.genes)
}

ens.genes.ALL <- initG2G(rownames(tpm.gene), 0)
ens.genes.G1S <- initG2G(genes.G1S, 1)
ens.genes.G2M <- initG2G(genes.G2M, 2)
ens.genes <- rbind(ens.genes.ALL, ens.genes.G1S, ens.genes.G2M)
ens.genes$Group <- as.factor(ens.genes$Group)
save(ens.genes, file=file.path(wd.de.plots, paste0(base, "_genes_ALL-G1S-G2M_g2g.RData")))

ens.genes <- subset(ens.genes, G2G != 0)  ## ADD 09/09/18
ens.genes$G2G <- log10(ens.genes$G2G)

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_ALL-G1S-G2M_g2g.pdf"))
pdf(file.name, height=6, width=3)
ymin <- min(ens.genes$G2G)
ymax <- max(ens.genes$G2G)
boxplot(G2G ~ Group, data=ens.genes, outline=T, names=c("All", "G1/S", "G2/M"), col=c("white", "cornsilk", "darkgray"), ylim=c(ymin, ymax), ylab="Gene-to-gene min dist. (log10)", main=BASE)
dev.off()

## ALL vs G1-S vs G2-M
testW(ens.genes.G1S$G2G, ens.genes.ALL$G2G)
# [1] 1.263306e-159
testW(ens.genes.G2M$G2G, ens.genes.ALL$G2G)
# [1] 0
testW(ens.genes.G1S$G2G, ens.genes.G2M$G2G)
# [1] 2.43581e-12
