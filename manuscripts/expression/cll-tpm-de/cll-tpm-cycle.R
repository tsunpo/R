# =============================================================================
# Manuscript   : 
# Chapter I    : 
# Name         : manuscripts/expression/cll-tpm-cycle.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/09/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Asymmetry.R", "CellCycle.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "cycle.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "CLL"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "cll_rna_n74.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "MAX_INSERT_SIZE")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene)

# -----------------------------------------------------------------------------
# Gene sets that are expressed in the dataset
# -----------------------------------------------------------------------------
load(file.path(wd.src.ref, "cycle.RData"))   ## See guide-to-the/cycle.R

## Tirosh et al 2016 
core.G1S <- intersect(core.G1S, rownames(tpm.gene.log2))   ## 40/43/43
core.G2M <- intersect(core.G2M, rownames(tpm.gene.log2))   ## 50/54/55

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 262/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 792/876

# -----------------------------------------------------------------------------
# Distribution of cell cycle genes in quantiles
# -----------------------------------------------------------------------------
tx.q4 <- getTxQ4(tpm.gene.log2, NA)
for (q in 1:4)
   print(length(tx.q4[[q]]))
# [1] 4626
# [1] 4625
# [1] 4625
# [1] 4626

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

boxplotTxQ4G2G(wd.de.plots, base, BASE, tx.q4.g2g)
