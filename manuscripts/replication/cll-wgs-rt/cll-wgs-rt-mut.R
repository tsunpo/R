# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/cll-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/02/19
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Common.R", "Mutation.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE1 <- "CLL"
PAIR1 <- "T"
BASE0 <- "CLL"
PAIR0 <- "T"
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

wd.anlys <- file.path(wd, BASE1, "analysis")
wd.rt       <- file.path(wd.anlys, "replication", paste0(base1, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs <- file.path(wd, BASE1, "ngs/WGS")
#samples1 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
samples  <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
samples1 <- readTable(file.path(wd.ngs, "cll_rna_n93-rt29.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "cll_rna_n93-wt33.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# 
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
phenos <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/CLL/ngs/WGS/peiflyne/cll_cn_summary.txt", header=F, rownames=T, sep="")
colnames(phenos) <- c("V1", "purity", "ploidy", "V4", "sex", "V6", "V7")

samples.pos <- samples1
samples.neg <- samples0
samples.var <- setdiff(samples, c(samples1, samples0))
# > 29+33+34
# [1] 96

# -----------------------------------------------------------------------------
# Mutation burdens between T29 and T33
# Last Modified: 15/03/19
# -----------------------------------------------------------------------------
muts <- c()
for (s in 1:length(samples)) {
   sample <- samples[s]
   txt <- read.peiflyne.muts.txt(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_muts.txt")), type="SNM")
 
   muts <- unique(c(muts, txt$Gene_Hugo))
}

samples.rts <- list(samples1, samples0, setdiff(samples, c(samples1, samples0)))
txts <- toTable(0, length(samples), length(muts), c(samples1, samples0, setdiff(samples, c(samples1, samples0))))
rownames(txts) <- muts
for (r in 1:length(samples.rts)) {
   samples.rt <- samples.rts[[r]]

   for (s in 1:length(samples.rt)) {
      sample <- samples.rt[s]
      txt <- read.peiflyne.muts.txt(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_muts.txt")), type="SNM")
      mut <- data.frame(table(txt$Gene_Hugo))
      rownames(mut) <- mut$Var1

      txts[rownames(mut), sample] <- mut$Freq
   }
}

txts.pos <- txts[,samples.rts[[1]]]
txts.pos$T29 <- 0
txts.pos$T29 <- mapply(x = 1:nrow(txts.pos), function(x) sum(txts.pos[x,]))
writeTable(txts.pos, file.path(wd.rt.data, "muts_pos29.txt"), colnames=T, rownames=T, sep="\t")

txts.neg <- txts[,samples.rts[[2]]]
txts.neg$T33 <- 0
txts.neg$T33 <- mapply(x = 1:nrow(txts.neg), function(x) sum(txts.neg[x,]))
writeTable(txts.neg, file.path(wd.rt.data, "muts_neg33.txt"), colnames=T, rownames=T, sep="\t")

txts.var <- txts[,samples.rts[[3]]]
txts.var$T34 <- 0
txts.var$T34 <- mapply(x = 1:nrow(txts.var), function(x) sum(txts.var[x,]))
writeTable(txts.var, file.path(wd.rt.data, "muts_var34.txt"), colnames=T, rownames=T, sep="\t")

