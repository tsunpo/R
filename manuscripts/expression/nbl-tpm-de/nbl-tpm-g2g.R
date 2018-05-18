# =============================================================================
# Manuscript   : The dangeous case of DNA replication 
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/nbl-tpm-g2g.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "Asymmetry.R", "ReplicationOrigin.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
BASE <- "NBL"
#wd     <- paste0("/ngs/cangen/tyang2/", BASE, "/analysis/")                   ## tyang2@gauss
#wd.ngs <- paste0("/ngs/cangen/tyang2/", BASE, "/ngs/WGS/")
wd     <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/")   ## tpyang@localhost
wd.ngs <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/ngs/WGS/")

wd.asym       <- paste0(wd, "asymmetries/", tolower(BASE), "-asym-tx/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "nbl_wgs_n57.list"), header=F, rownames=F, sep="")

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load(paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/expression/kallisto/", tolower(BASE), "-tpm-de/data/", tolower(BASE), "_kallisto_0.43.1_tpm.gene_r5_p47.RData"))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2582
# [1] 2582
# [1] 2582
# [1] 2582

g2g.q4 <- getG2GQ4(tx.q4)
p3 <- testT(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testT(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testT(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 0.002381438 0.074382250 0.312440900
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.007144314 0.223146750 0.937322700

plotG2GQ4(wd.asym.plots, BASE, length(samples), g2g.q4)
