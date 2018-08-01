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
BASE <- "HeLa"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.plots <- file.path(wd.asym, "plots")

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%       25%       50%       75%      100% 
# -5.9033046  0.6216422  3.5646662  5.2402453 12.8536280 

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2349
# [1] 2349
# [1] 2348
# [1] 2349

g2g.q4 <- getG2GQ4(tx.q4)
testOnewayANOVA(g2g.q4)
# [1] 8.924556e-03
testW(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]]), g2g.q4[[4]])
# [1] 3.525678e-26
median(as.numeric(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]])))
# [1] 229176.5
median(as.numeric(g2g.q4[[4]]))
# [1] 142399

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 1.538768e-17 7.415199e-23 4.107799e-14
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 4.616303e-17 2.224560e-22 1.232340e-13

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.30103, 7.103921))
# > log10(max(distances))
# [1] 6.997238
# > log10(min(distances))
# [1] 0.30103

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=1599.855)
# > max(getDensity(g2g.q4[[3]], count)$y)
# [1] 1561.846










# =============================================================================
# Step 3: Look at the little bumps (1KB)
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)

g2g.q4 <- get1KBG2GQ4(tx.q4, 1000)
p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 1.042204e-02 8.234912e-06 1.196299e-01
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 3.126612e-02 2.470474e-05 3.588896e-01

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=NULL)

# =============================================================================
# Step 4: Look at the little bumps (G1-S, S-G2, and G2-M)
# Last Modified: 23/05/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
cycle <- "G2-M"

samples.G1toS <- rownames(subset(samples, CELL_CYCLE == cycle))
tpm.gene.input <- tpm.gene.input[,samples.G1toS]
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[,samples.G1toS])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)

g2g.q4 <- getG2GQ4(tx.q4)
p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 1.042204e-02 8.234912e-06 1.196299e-01
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 3.126612e-02 2.470474e-05 3.588896e-01

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_", cycle, "_g2g.pdf"))
file.main <- paste0(BASE, " (", cycle, "; n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_", cycle, "_g2g_d.pdf"))
file.main <- paste0(BASE, " (", cycle, "; n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, NULL)






###
##
g2g.q4 <- getQ1G2GQ4(tx.q4)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 1.042204e-02 8.234912e-06 1.196299e-01
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 3.126612e-02 2.470474e-05 3.588896e-01

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_", cycle, "_q1_g2g.pdf"))
file.main <- paste0(BASE, " (", cycle, "; n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_", cycle, "_q1_g2g_d.pdf"))
file.main <- paste0(BASE, " (", cycle, "; n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=NULL)













# =============================================================================
# Step 3: Look at the little bumps (Q1)
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getQ1G2GQ4(tx.q4)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 3.577538e-48 3.060852e-42 1.242962e-26
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.073261e-47 9.182557e-42 3.728887e-26

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_q1_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_q1_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)

