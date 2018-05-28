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

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx"))
wd.asym.plots <- file.path(wd.asym, "plots")
setwd(wd.asym)

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%       25%       50%       75%      100% 
# -5.737860  1.390236  3.487460  4.898448 12.365721

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
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.plots <- file.path(wd.asym, "plots")
setwd(wd.asym)

tpm.gene.input <- pipeTPM(BASE)
g2g.q4 <- pipeRO(tpm.gene.input)

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)









getDensity <- function(distances) {
   return(density(log10(as.numeric(distances))))
}

plotDensity <- function(g2g.q4, file.name, file.main) {
   ymax <- max(getDensity(g2g.q4[[1]])$y)
   for (q in 2:4)
      if (max(getDensity(g2g.q4[[q]])$y) > ymax)
         ymax <- max(getDensity(g2g.q4[[q]])$y)
   
   cols <- c("blue", "deepskyblue", "salmon", "red")
   pdf(file.name, height=6, width=6)
   plot(getDensity(g2g.q4[[1]]), col="blue", ylab="Frequency", xlab="Gene-to-gene min dist. (log10)", ylim=c(0, ymax), main=file.main)
   lines(getDensity(g2g.q4[[2]]), col="deepskyblue")
   lines(getDensity(g2g.q4[[3]]), col="salmon")
   lines(getDensity(g2g.q4[[4]]), col="red")
   
   legend("topleft", legend=c("25%", "50%", "75%", "100%"), levels(cols), fill=cols) 
   dev.off()
}

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d_test.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main)





distances <- c()
medians <- c()
for (q in 1:4) {
   distances <- c(distances, as.numeric(g2g.q4[[q]]))
   medians <- c(medians, tpm.gene.input.log2[tx.q4[[q]],]$MEDIAN)
}

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
pdf(file.name, height=6, width=6)
plot(density(log10(distances)), ylab="Frequency", xlab="Gene-to-gene min dist. (log10)")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_test.pdf"))
pdf(file.name, height=6, width=6)
plot(log10(distances), medians, ylab="Expression medians", xlab="Gene-to-gene min dist. (log10)")
dev.off()



# =============================================================================
# Step 2: Are most expressed genes overlapped?
# Last Modified: 14/05/18
# =============================================================================
genes.q4.sclc  <- tx.q4[[4]]
genes.q4.lcnec <- tx.q4[[4]]
genes.q4.lcnec.rb1    <- tx.q4[[4]]
genes.q4.lcnec.rb1.wt <- tx.q4[[4]]
genes.q4.luad <- tx.q4[[4]]
genes.q4.nbl <- tx.q4[[4]]
genes.q4.cll <- tx.q4[[4]]

# > length(genes.q4.sclc)
# [1] 2651
# > length(genes.q4.lcnec)
# [1] 2633
# > length(genes.q4.luad)
# [1] 2588
# > length(genes.q4.nbl)
# [1] 2582
# > length(genes.q4.cll)
# [1] 2326

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
overlaps.lung <- intersect(intersect(genes.q4.sclc, genes.q4.lcnec), genes.q4.luad)
overlaps <- intersect(intersect(intersect(genes.q4.sclc, genes.q4.lcnec), genes.q4.luad), genes.q4.nbl)
# > length(intersect(intersect(intersect(genes.q4.sclc, genes.q4.lcnec), genes.q4.luad), genes.q4.nbl))
# [1] 1631

# length(intersect(overlaps, genes.q4.cll))
# > 1148/2304
# [1] 0.4982639

# =============================================================================
# Step 2: Are most expressed genes overlapped?
# Last Modified: 14/05/18
# =============================================================================
# > length(intersect(genes.q4.sclc, genes.q4.lcnec.rb1))
# [1] 2361
# > length(intersect(genes.q4.sclc, genes.q3.lcnec.rb1))
# [1] 261
# > length(intersect(genes.q4.sclc, genes.q2.lcnec.rb1))
# [1] 22

# > length(intersect(genes.q4.sclc, genes.q4.lcnec.rb1.wt))
# [1] 2285
# > length(intersect(genes.q4.sclc, genes.q3.lcnec.rb1.wt))
# [1] 326
# > length(intersect(genes.q4.sclc, genes.q2.lcnec.rb1.wt))
# [1] 26

# > length(intersect(genes.q4.sclc, genes.q4.luad))
# [1] 2046
# > length(intersect(genes.q4.sclc, genes.q3.luad))
# [1] 482
# > length(intersect(genes.q4.sclc, genes.q2.luad))
# [1] 86

# > length(intersect(genes.q4.sclc, genes.q4.nbl))
# [1] 1954
# > length(intersect(genes.q4.sclc, genes.q3.nbl))
# [1] 569
# > length(intersect(genes.q4.sclc, genes.q2.nbl))
# [1] 87

# > length(intersect(genes.q4.sclc, genes.q4.cll))
# [1] 1477
# > length(intersect(genes.q4.sclc, genes.q3.cll))
# [1] 700
# > length(intersect(genes.q4.sclc, genes.q2.cll))
# [1] 252
