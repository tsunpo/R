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
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.plots <- file.path(wd.asym, "plots")

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_tpm0.RData")))
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%       25%       50%       75%      100% 
# -4.547471  2.075755  3.755650  5.033693 12.365721

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2406
# [1] 2406
# [1] 2406
# [1] 2406

g2g.q4 <- getG2GQ4(tx.q4)
p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 3.593696e-23 5.269174e-22 7.505582e-14
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.078109e-22 1.580752e-21 2.251675e-13

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)

# =============================================================================
# Step 3: Leading strands (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
genes.cd <- c()
for (q in 1:4)
   genes.cd <- c(genes.cd, intersect(ens.tx.rt.input, tx.q4.rt.input[[q]]))
# > length(genes.cd)
# [1] 4232
save(genes.cd, genes.ho, file=file.path(wd.asym.data, paste0(base,"_asym_tx_q4_rt_genes-cd+ho.RData")))

tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.cd,])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
# > log10(min(as.numeric(g2g.q4[[4]])))
# [1] 3.428135
# > ymax
# [1] 763.311

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 5.046751e-06 3.286694e-08 4.209444e-06
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.514025e-05 9.860083e-08 1.262833e-05

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_cd_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.428135, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_cd_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)

# =============================================================================
# Step 4: Lagging strands (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
genes.ho <- c()
for (q in 1:4)
 genes.ho <- c(genes.ho, intersect(ens.tx.rt.input, tx.q4.rt.input[[q]]))
# > length(genes.ho)
# [1] 4418

tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.ho,])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
# > log10(max(as.numeric(g2g.q4[[2]])))
# [1] 7.468407

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 6.362036e-08 1.076820e-05 3.289446e-09
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.908611e-07 3.230460e-05 9.868337e-09

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.428135, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)











# =============================================================================
# Step 5: Right-lagging strands (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
# > nrow(subset(ensGene.tx.rt.input[genes.ho,], RT == 1))
# [1] 2149
# > nrow(subset(ensGene.tx.rt.input[genes.ho,], RT == -1))
# [1] 2269
genes.ho.right <- rownames(subset(ensGene.tx.rt.input[genes.ho,], RT == 1))
genes.ho.left  <- rownames(subset(ensGene.tx.rt.input[genes.ho,], RT == -1))

tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.ho.right,])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
# > log10(max(as.numeric(g2g.q4[[2]])))
# [1] 7.468407

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 0.0084963829 0.0006914775 0.0026034907
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.025489149 0.002074433 0.007810472

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho-right_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.428135, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho-right_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)

# =============================================================================
# Step 6: Left-lagging strands (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.ho.left,])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
# > log10(max(as.numeric(g2g.q4[[2]])))
# [1] 7.468407

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 0.02278788 0.46311374 0.00922025
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.06836365 1.00000000 0.02766075

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho-left_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.428135, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho-left_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)
















# =============================================================================
# Step 5: All genes (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[c(genes.cd, genes.ho),])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)
# > log10(min(as.numeric(g2g.q4[[2]])))
# [1] 0.60206
# > ymax
# [1] 1408.151

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 5.046751e-06 3.286694e-08 4.209444e-06
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.514025e-05 9.860083e-08 1.262833e-05

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_cd_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.60206, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_cd_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)





# =============================================================================
# Step 5: Right replication (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================


genes.ho <- c()
for (q in 1:4)
 genes.ho <- c(genes.ho, intersect(ens.tx.rt.input, tx.q4.rt.input[[q]]))
# > length(genes.ho)
# [1] 4418

tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.ho,])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 6.362036e-08 1.076820e-05 3.289446e-09
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.908611e-07 3.230460e-05 9.868337e-09

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.428135, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)

# =============================================================================
# Step 4: Leading strands (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
genes.cd <- c()
for (q in 1:4)
 genes.cd <- c(genes.cd, intersect(ens.tx.rt.input, tx.q4.rt.input[[q]]))
# > length(genes.cd)
# [1] 4232

tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.cd,])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 5.046751e-06 3.286694e-08 4.209444e-06
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.514025e-05 9.860083e-08 1.262833e-05

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_cd_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.428135, 7.468407))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_cd_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=763.311)
















# =============================================================================
# Step 5: Leading+lagging strands (Following Step 6.2 from sclc-asym-tx-rt.R)
# Last Modified: 14/06/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[c(genes.cd, genes.ho),])
tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
g2g.q4 <- getG2GQ4(tx.q4)

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 3.331408e-20 7.748811e-20 1.176142e-12
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 9.994224e-20 2.324643e-19 3.528426e-12

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho+cd_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ho+cd_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
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
# [1] 3.449370e-45 1.008533e-65 2.228323e-30
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.034811e-44 3.025599e-65 6.684969e-30

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ro_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ro_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)









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
# [1] 0.0007606195 0.0647591740 0.0391856377
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.002281858 0.194277522 0.117556913

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d.pdf"))
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
