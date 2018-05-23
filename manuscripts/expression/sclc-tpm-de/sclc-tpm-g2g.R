# =============================================================================
# Manuscript   : The dangeous case of DNA replication
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/sclc-tpm-g2g.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Asymmetry.R", "ReplicationOrigin.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
base <- tolower(BASE)

wd.rna <- file.path(wd, BASE, "ngs/RNA")

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx"))
wd.asym.plots <- file.path(wd.asym, "plots")
setwd(wd.asym)

samples <- readTable(file.path(wd.rna, "sclc_rna_n81.list"), header=F, rownames=T, sep="")

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2651
# [1] 2651
# [1] 2651
# [1] 2651

g2g.q4 <- getG2GQ4(tx.q4)
p3 <- testT(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testT(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testT(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 0.0005114360 0.0006123511 0.0143906231
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.001534308 0.001837053 0.043171869

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", nrow(samples), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

# =============================================================================
# Step 2: Are most expressed genes overlapped?
# Last Modified: 14/05/18
# =============================================================================
genes.q4.sclc  <- tx.q4[[4]]
genes.q4.lcnec <- tx.q4[[4]]
genes.q4.lcnec.rb1    <- tx.q4[[4]]
genes.q4.lcnec.rb1.wt <- tx.q4[[4]]
genes.q4.luad <- tx.q4[[4]]
genes.q4.nbl  <- tx.q4[[4]]

# > length(genes.q4.sclc)
# [1] 2651
# > length(genes.q4.lcnec)
# [1] 2633
# > length(genes.q4.luad)
# [1] 2588
# > length(genes.q4.nbl)
# [1] 2582

# > length(intersect(genes.q4.lcnec.rb1, genes.q4.lcnec.rb1.wt))
# [1] 2388
# > length(genes.q4.lcnec.rb1)
# [1] 2633
# > length(genes.q4.lcnec.rb1.wt)
# [1] 2633
# > 2388/2633
# [1] 0.9069502

# > length(intersect(genes.q4.sclc, genes.q4.lcnec))
# [1] 2336
# > length(intersect(genes.q4.sclc, genes.q4.luad))
# [1] 2046
# > length(intersect(genes.q4.sclc, genes.q4.nbl))
# [1] 1954
# > length(intersect(genes.q4.lcnec, genes.q4.nbl))
# [1] 1874
# > length(intersect(genes.q4.luad, genes.q4.nbl))
# [1] 1740
# > length(intersect(intersect(intersect(genes.q4.sclc, genes.q4.lcnec), genes.q4.luad), genes.q4.nbl))
# [1] 1631
