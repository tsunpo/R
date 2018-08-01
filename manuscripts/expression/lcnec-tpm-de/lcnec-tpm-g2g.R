# =============================================================================
# Manuscript   : The dangeous case of DNA replication
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/lcnec-tpm-g2g.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/05/18
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
BASE <- "LCNEC"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.plots <- file.path(wd.asym, "plots")

wd.rna <- file.path(wd, BASE, "ngs/RNA")
samples <- readTable(file.path(wd.rna, "lcnec_rna_n69.list"), header=F, rownames=T, sep="")

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%       25%       50%       75%      100% 
# -5.423718  1.144500  3.351758  4.806042 12.469278

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2638
# [1] 2638
# [1] 2637
# [1] 2638

g2g.q4 <- getG2GQ4(tx.q4)
p <- testOnewayANOVA(g2g.q4)
# [1] 1.762388e-05

p <- testW(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]]), g2g.q4[[4]])
# [1] 8.175331e-34

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 2.594506e-24 1.067951e-27 7.484947e-18
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 7.783518e-24 3.203852e-27 2.245484e-17

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.30103, 7.103921))
# > log10(max(distances))
# [1] 6.947882
# > log10(min(distances))
# [1] 0.30103

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=1881.203)
# > max(getDensity(g2g.q4[[2]], count)$y)
# [1] 1881.203

# =============================================================================
# Step 1: Gene-to-gene minmum distance (RB1)
# Last Modified: 18/05/18
# =============================================================================
samples <- subset(samples, V5 == 1)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene[,samples$V1], ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%       25%       50%       75%      100% 
# -5.654260  1.160195  3.417934  4.849166 12.365296

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
testOnewayANOVA(g2g.q4)
# [1] 1.310906e-05
testW(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]]), g2g.q4[[4]])
# [1] 1.002245e-26
median(as.numeric(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]])))
# [1] 210403
median(as.numeric(g2g.q4[[4]]))
# [1] 139858

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_RB1.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.30103, 7.103921))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d_RB1.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=1881.203)

# =============================================================================
# Step 1: Gene-to-gene minmum distance (RB1-WT)
# Last Modified: 18/05/18
# =============================================================================
samples <- subset(samples, V5 == 0)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene[,samples$V1], ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%       25%       50%       75%      100% 
# -5.616245  1.236817  3.387824  4.813272 12.651321

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
testOnewayANOVA(g2g.q4)
# [1] 3.73474e-05
testW(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]]), g2g.q4[[4]])
# [1] 1.017538e-34
median(as.numeric(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]])))
# [1] 210403
median(as.numeric(g2g.q4[[4]]))
# [1] 134036

testW(g2g.q4.lenec.rb1, g2g.q4.lenec.wt)
# [1] 0.4365798

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_RB1-WT.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.30103, 7.103921))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d_RB1-WT.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=1881.203)













###
## LCENC (RB1-WT)
c(p1, p2, p3)
# [1] 9.650667e-28 1.316295e-23 1.017159e-17
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 2.895200e-27 3.948886e-23 3.051478e-17

ymax <- log10(max(as.numeric(g2g.q4[[2]])))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_RB1-WT.pdf"))
file.main <- paste0(BASE, " (RB1-WT; n=", nrow(subset(samples, V5 == 0)), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.3, ymax))

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", nrow(samples), ")")
plotDensity(g2g.q4, file.name, file.main, T)

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d_RB1.pdf"))
file.main <- paste0(BASE, " (RB1; n=", nrow(subset(samples, V5 == 1)), ")")
plotDensity(g2g.q4, file.name, file.main, T)

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d_RB1-WT.pdf"))
file.main <- paste0(BASE, " (RB1-WT; n=", nrow(subset(samples, V5 == 0)), ")")
plotDensity(g2g.q4, file.name, file.main, T)

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
# [1] 0.01463031 0.35342132 0.13099764
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.04389094 1.00000000 0.39299291

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)

###
## RB1
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input <- tpm.gene.input[,rownames(subset(samples, V5 == 1))]
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- get1KBG2GQ4(tx.q4, 1000)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 0.008760673 0.144840721 0.030408631
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.02628202 0.43452216 0.09122589

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_RB1.pdf"))
file.main <- paste0(BASE, " (RB1; n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d_RB1.pdf"))
file.main <- paste0(BASE, " (RB1; n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)

###
## RB1-WT
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input <- tpm.gene.input[,rownames(subset(samples, V5 == 0))]
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- get1KBG2GQ4(tx.q4, 1000)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 0.04344371 0.11290129 0.18574353
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.1303311 0.3387039 0.5572306

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_RB1-WT.pdf"))
file.main <- paste0(BASE, " (RB1-WT; n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d_RB1-WT.pdf"))
file.main <- paste0(BASE, " (RB1-WT; n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)
