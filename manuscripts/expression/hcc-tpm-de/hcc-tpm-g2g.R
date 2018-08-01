# =============================================================================
# Manuscript   : The dangeous case of DNA replication 
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/cll-tpm-g2g.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 23/05/18
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
BASE <- "HCC"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.plots <- file.path(wd.asym, "plots")

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%        25%        50%        75%       100% 
# -5.6001991  0.5970103  2.7232855  4.2414861 14.5334664 

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2489
# [1] 2488
# [1] 2488
# [1] 2489

g2g.q4 <- getG2GQ4(tx.q4)
save(g2g.q4, file=file.path(wd.asym.data, paste0(base, "_asym_tx_g2g.RData")))
testOnewayANOVA(g2g.q4)
# [1] 0.0003166389
testW(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]]), g2g.q4[[4]])
# [1] 7.983115e-35





median(as.numeric(c(g2g.q4[[1]], g2g.q4[[2]], g2g.q4[[3]])))
# [1] 
median(as.numeric(g2g.q4[[4]]))
# [1] 

p3 <- testW(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testW(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testW(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 1.484115e-37 4.173536e-18 4.373463e-18
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 4.452346e-37 1.252061e-17 1.312039e-17

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(0.30103, 7.103921))
# > log10(max(distances))
# [1] 6.936843
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
# [1] 1613.451










# =============================================================================
# Test: All Tx(+)
# Last Modified: 15/06/18
# =============================================================================
tpm.gene.input <- pipeTPM(wd, BASE)
genes.plus  <- rownames(subset(ensGene[rownames(tpm.gene.input),], strand == 1))
genes.minus <- rownames(subset(ensGene[rownames(tpm.gene.input),], strand == -1))
# > length(genes.plus)
# [1] 4488
# > length(genes.minus)
# [1] 4309

g2g.q4 <- pipeG2G(wd, BASE, genes.plus)
# [1] 6.193839e-05 7.985932e-02 1.434225e-01
# [1] 0.0001858152 0.2395779548 0.4302673601
# > ymax
# [1] 754.7507

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_plus_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.446848, 7.350468))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_plus_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=702.1426)

##
g2g.q4.plus <- g2g.q4[[4]]

# =============================================================================
# Test: All Tx(-)
# Last Modified: 15/06/18
# =============================================================================
g2g.q4 <- pipeG2G(wd, BASE, genes.minus)
# [1] 1.402395e-12 9.014271e-06 1.476496e-05
# [1] 4.207184e-12 2.704281e-05 4.429489e-05
# > ymax
# [1] 702.1426

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_minus_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(3.446848, 7.350468))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_minus_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T, ymax=702.1426)

##
g2g.q4.minus <- g2g.q4[[4]]
testW(g2g.q4.plus, g2g.q4.minus)
# [1] 0.2006883












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
# [1] 0.2379321 0.1344062 0.9925154
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 0.7137963 0.4032187 1.0000000

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=c(log10(2), 3))

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_1kb_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)










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
# [1] 6.648437e-31 7.239728e-01 5.664283e-01
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.994531e-30 1.000000e+00 1.000000e+00

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ro_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

##
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_ro_g2g_d.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotDensity(g2g.q4, file.name, file.main, count=T)















# -----------------------------------------------------------------------------
# Max insert size
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

samples <- readTable(file.path(wd.rna, "cll_rna_n74.list"), header=F, rownames=T, sep="")
colnames(samples) <- c("SAMPLE_ID", "FILE_NAME", "MAX_INSERT_SIZE")

samples$MEDIAN <- mapply(x = 1:nrow(samples), function(x) median(as.numeric(tpm.gene.log2[,x])))

pdf(file.path(wd.de.plots, "max_insert_size.pdf"))
plot(samples$MEDIAN, samples$MAX_INSERT_SIZE)
lm.fit <- lm(samples$MEDIAN ~ samples$MAX_INSERT_SIZE)
dev.off()

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
pca.de <- getPCA(t(tpm.gene.log2))   ## BUG FIX 13/02/17: Perform PCA using normalised data

###
## MAX_INSERT_SIZE
trait <- as.numeric(samples[,"MAX_INSERT_SIZE"])
greater <- which(trait > 300)
lesser  <- which(trait <= 300)
trait[greater] <- ">300"
trait[lesser]  <- "<300"

file.name <- "pca_MAX_INSERT_SIZE"
file.main <- "Max insert size in CLL"
plotPCA(1, 2, pca.de, trait, wd.de.plots, file.name, file.main, NA, NA, c("dodgerblue", "red"))









# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
plotDensity <- function(g2g.q4, file.name, file.main) {
   distances <- c()
   for (q in 1:4)
      distances <- c(distances, as.numeric(g2g.q4[[q]]))
   
   pdf(file.name, height=6, width=6)
   plot(density(distances), ylab="Frequency", xlab="Gene-to-gene min dist. (log10)")
   
   lines(density(as.numeric(g2g.q4[[4]])), col="red")
   lines(density(as.numeric(g2g.q4[[3]])), col="salmon")
   lines(density(as.numeric(g2g.q4[[2]])), col="deepskyblue")
   lines(density(as.numeric(g2g.q4[[1]])), col="blue")
   
   dev.off()
}

plotDensity <- function(g2g.q4, file.name, file.main) {
   pdf(file.name, height=6, width=6)
   plot(density(log10(as.numeric(g2g.q4[[1]])), bw=0.15), col="blue", ylab="Frequency", xlab="Gene-to-gene min dist. (log10)")
   lines(density(log10(as.numeric(g2g.q4[[2]])), bw=0.15), col="deepskyblue")
   lines(density(log10(as.numeric(g2g.q4[[3]])), bw=0.15), col="salmon")
   lines(density(log10(as.numeric(g2g.q4[[4]])), bw=0.15), col="red")
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

d <- density(mtcars$mpg)
plot(d) 


