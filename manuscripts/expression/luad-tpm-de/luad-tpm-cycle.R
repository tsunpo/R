# =============================================================================
# Manuscript   : 
# Chapter I    : 
# Name         : manuscripts/expression/luad-tpm-cycle.R
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

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "LUAD"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "luad_rna_n49-1.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "AVG_FRAGMENT_LENGTH", "MAX_INSERT_SIZE")

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene)

# -----------------------------------------------------------------------------
# Gene sets (please refer to guide-to-the/cycle.R) that are expressed in the dataset
# -----------------------------------------------------------------------------
load(file.path(wd.src.ref, "cycle.RData"))

## Tirosh et al 2016 
core.G1S <- intersect(core.G1S, rownames(tpm.gene.log2))   ## 41/43/43
core.G2M <- intersect(core.G2M, rownames(tpm.gene.log2))   ## 54/54/55

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 280/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 826/876

# -----------------------------------------------------------------------------
# Distribution of cell cycle genes in four quantiles
# -----------------------------------------------------------------------------
tx.q4 <- getTxQ4(NA, tpm.gene.log2)
for (q in 1:4)
   print(length(tx.q4[[q]]))
# [1] 4669
# [1] 4668
# [1] 4668
# [1] 4669

tx.q4.cycle <- getTxQ4Cycle(tx.q4, core.G1S, core.G2M, periodic.G1S, periodic.G2M)
save(tx.q4.cycle, file=file.path(wd.de.plots, paste0(base, "_genes_tx_q4_cycle.Rdata")))

barplotTxQ4Cycle(wd.de.plots, base, BASE, tx.q4.cycle, beside=T)

# -----------------------------------------------------------------------------
# Gene length in four quantiles
# -----------------------------------------------------------------------------
tx.q4.length <- getTxQ4Length(tx.q4)
save(tx.q4.length, file=file.path(wd.de.plots, paste0(base, "_genes_tx_q4_length.Rdata")))

boxplotTxQ4Length(wd.de.plots, base, BASE, tx.q4.length)
