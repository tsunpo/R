# =============================================================================
# Manuscript   : The dangeous case of DNA replication
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/luad-tpm-g2g.R
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
BASE <- "LUAD"
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
# -5.623874  1.257033  3.372138  4.770867 13.149860

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2588
# [1] 2588
# [1] 2588
# [1] 2588

g2g.q4 <- getG2GQ4(tx.q4)
testOnewayANOVA(g2g.q4)
# [1] 2.068892e-05
testW(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]]), g2g.q4[[4]])
# [1] 1.262442e-31
median(as.numeric(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]])))
# [1] 204480.5
median(as.numeric(g2g.q4[[4]]))
# [1] 133284

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 2.299704e-29 5.478248e-19 1.787152e-17
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 6.899112e-29 1.643474e-18 5.361456e-17

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.30103, 7.103921))
# > log10(max(distances))
# [1] 6.979327
# > log10(min(distances))
# [1] 0.30103

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/histdens.html
# Last Modified: 23/05/18
# =============================================================================
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=1881.203)
# > max(getDensity(g2g.q4[[3]], count)$y)
# [1] 1653.476











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
# [1] 3.449370e-45 1.008533e-65 2.228323e-30
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.034811e-44 3.025599e-65 6.684969e-30

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_q1_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_q1_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=NULL)









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
# [1] 0.0003957046 0.0142010600 0.1295683109
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.001187114 0.042603180 0.388704933

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)










# =============================================================================
# Step: Smoking status
# Last Modified: 08/06/18
# =============================================================================
wd.ngs <- file.path(wd, "LUAD/ngs/WGS")
wd.rna <- file.path(wd, "LUAD/ngs/RNA")
samples.broad <- readTable(file.path(wd.ngs, "luad_wgs_Broad_n15.list"), header=F, rownames=F, sep="")
samples.dkfz  <- readTable(file.path(wd.ngs, "luad_wgs_DKFZ_n16.list"), header=F, rownames=F, sep="")
samples.pan   <- readTable(file.path(wd.ngs, "luad_wgs_PanNeg_n7.list"), header=F, rownames=F, sep="")
samples.ngx   <- readTable(file.path(wd.ngs, "luad_wgs_NGXBIO_n1.list"), header=F, rownames=F, sep="")

samples.wgs <- c(samples.broad, samples.dkfz, samples.pan, samples.ngx)
samples.rna <- readTable(file.path(wd.rna, "luad_rna_n49-1.list"), header=F, rownames=T, sep="")[,1]

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main)
