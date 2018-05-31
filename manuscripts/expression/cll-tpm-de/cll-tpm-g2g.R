# =============================================================================
# Manuscript   : The dangeous case of DNA replication 
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/cll-tpm-g2g.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 23/05/18
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
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "CLL"
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
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)
# > quantile(tpm.gene.input.log2$MEDIAN)
# 0%        25%        50%        75%       100% 
# -5.7550939  0.4394228  3.2257788  4.9586013 13.9241600 

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2326
# [1] 2325
# [1] 2325
# [1] 2326

g2g.q4 <- getG2GQ4(tx.q4)
p3 <- testT(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testT(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testT(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 2.251716e-06 1.655007e-01 2.680425e-01
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 6.755149e-06 4.965022e-01 8.041276e-01

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g.pdf"))
file.main <- paste0(BASE, " (n=", ncol(tpm.gene.input), ")")
plotG2GQ4(g2g.q4, file.name, file.main, ylim=NULL)

# =============================================================================
# Step 2: Density plots
# https://www.statmethods.net/graphs/density.html
# Last Modified: 23/05/18
# =============================================================================
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "CLL"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.plots <- file.path(wd.asym, "plots")

tpm.gene.input <- pipeTPM(BASE)
g2g.q4 <- pipeRO("CLL", tpm.gene.input)

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_g2g_d.pdf"))
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


