# =============================================================================
# Manuscript   : 
# Chapter I    : 
# Name         : manuscripts/expression/sclc-tpm-cycle.R
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

BASE <- "LUSC"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples <- readTable(file.path(wd.rna, "lusc_rna_n423.txt"), header=F, rownames=T, sep="")
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
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 842/876

# -----------------------------------------------------------------------------
# Distribution of cell cycle genes in four quantiles
# -----------------------------------------------------------------------------
tx.q4 <- getTxQ4(NA, tpm.gene.log2)
for (q in 1:4)
   print(length(tx.q4[[q]]))
# [1] 4939
# [1] 4938
# [1] 4938
# [1] 4938

tx.q4.cycle <- getTxQ4Cycle(tx.q4, core.G1S, core.G2M, periodic.G1S, periodic.G2M)
save(tx.q4.cycle, file=file.path(wd.de.plots, paste0(base, "_genes_tx_q4_cycle.Rdata")))

barplotTxQ4Cycle(wd.de.plots, base, BASE, tx.q4.cycle, beside=T)






# -----------------------------------------------------------------------------
# Gene length in four quantiles
# -----------------------------------------------------------------------------
tx.q4.length <- initLength(tx.q4[[1]], 1)[0,]
for (q in 1:4) {
   #genes <- intersect(ens.tx.snv.input, tx.q4[[q]])
   genes <- tx.q4[[q]]
 
   tx.q4.length <- rbind(tx.q4.length, initLength(genes, q))
 
   genes.g1s <- intersect(genes, unique(c(core.G1S, genes.G1S)))
   genes.g2m <- intersect(genes, unique(c(core.G2M, genes.G2M)))
   tx.q4.cycle[[q]][[1]] <- length(genes.g1s)
   tx.q4.cycle[[q]][[2]] <- length(genes.g2m)
   tx.q4.cycle[[q]][[3]] <- length(setdiff(genes, c(genes.g1s, genes.g2m)))
}

##
tx.q4.length$Group <- as.factor(tx.q4.length$Group)
tx.q4.length$Length <- log10(tx.q4.length$Length)
colnames <- c("Q1", "Q2", "Q3", "Q4")
rownames <- c("G1-S", "G2-M", "Others")

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_length.pdf"))
pdf(file.name, height=6, width=4)
ymin <- min(tx.q4.length$Length)
ymax <- max(tx.q4.length$Length)
boxplot(Length ~ Group, data=tx.q4.length, outline=T, names=colnames, ylim=c(ymin, ymax), ylab="Gene length (log10)", xlab="Expression", main="SCLC")
dev.off()

##
data <- toTable(0, length(colnames), 3, colnames)
rownames(data) <- rownames
for (q in 1:4)
   for (r in 1:3)
      data[r, q] <- tx.q4.cycle[[q]][[r]]
writeTable(data, file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_length.txt")), colnames=T, rownames=T, sep="\t")
data <- as.matrix(data)

##
data <- toTable(0, length(colnames), 3, colnames)
for (q in 1:4) {
   sum <- tx.q4.cycle[[q]][[1]]
   for (r in 2:3)
      sum <- sum + tx.q4.cycle[[q]][[r]]
   
   for (r in 1:3)
      data[r, q] <- tx.q4.cycle[[q]][[r]] / sum * 100
}
writeTable(data, file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_percent.txt")), colnames=T, rownames=T, sep="\t")
data <- as.matrix(data)

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_percent.pdf"))
cols <- c("white", "darkgray")
pdf(file.name, height=6, width=4)
barplot(data[-3,], ylab="Cell cycle genes (%)", xlab="Expression", col=cols, main="SCLC", beside=T)
legend("topleft", legend=rownames[-3], fill=cols)
dev.off()
