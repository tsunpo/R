# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/10/19; 05/03/19; 25/02/19; 30/01/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_WB.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

samples.sclc <- setSamplesQ4(wd.rt.data, samples1)

nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

test.sclc <- nrds.T.chr.d.all

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "NBL"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

samples.nbl <- setSamplesQ4(wd.rt.data, samples1)

nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

test.nbl <- nrds.T.chr.d.all

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "CLL"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

samples.cll <- setSamplesQ4(wd.rt.data, samples1)

nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]

   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
   
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

test.cll <- nrds.T.chr.d.all

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-normal"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

#wd.ngs.data <- file.path(wd.ngs, "data")
wd.ngs.data <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"), "data")
samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
## samples1 are modified to the normals in this script

samples.sclc.nl <- setSamplesQ4(wd.rt.data, samples1)

nrds.T.chr.d.all <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

test.sclc.nl <- nrds.T.chr.d.all

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "NBL"
PAIR1 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-normal"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

#wd.ngs.data <- file.path(wd.ngs, "data")
wd.ngs.data <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"), "data")
samples1 <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
## samples1 are modified to the normals in this script

samples.nbl.wb <- setSamplesQ4(wd.rt.data, samples1)

nrds.T.chr.d.all <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

test.nbl.wb <- nrds.T.chr.d.all

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "CLL"
PAIR1 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt-normal"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

#wd.ngs.data <- file.path(wd.ngs, "data")
wd.ngs.data <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"), "data")
samples1 <- readTable(file.path(wd.ngs, "cll_wgs_n96.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
## samples1 are modified to the normals in this script

samples.cll.wb <- setSamplesQ4(wd.rt.data, samples1)

nrds.T.chr.d.all <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]

   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.ngs.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
   
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

test.cll.wb <- nrds.T.chr.d.all
rm(nrds.T.chr.d.all)

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")
samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n92_NL.list"), header=F, rownames=F, sep="")
#samples0 <- readTable(file.path(wd.ngs, "sclc_wgs_n9_WB.list"), header=F, rownames=F, sep="")
n1 <- length(samples1)
n0 <- length(samples0)

overlaps <- intersect(intersect(intersect(intersect(intersect(rownames(test.sclc), rownames(test.nbl)), rownames(test.cll)), rownames(test.sclc.nl)), rownames(test.nbl.wb)), rownames(test.cll.wb))
test.sclc <- test.sclc[overlaps, -1]
test.nbl  <- test.nbl[overlaps, -1]
test.cll  <- test.cll[overlaps, -1]
test.sclc.nl <- test.sclc.nl[overlaps, -1]
test.nbl.wb <- test.nbl.wb[overlaps, -1]
test.cll.wb <- test.cll.wb[overlaps, -1]

test <- cbind(test.sclc, test.nbl, test.cll, test.sclc.nl, test.nbl.wb, test.cll.wb)
rm(test.sclc)
rm(test.nbl)
rm(test.cll)
rm(test.sclc.nl)
rm(test.nbl.wb)
rm(test.cll.wb)

pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc+nbl+cll+sclc.nl+nbl.wb+cll.wb_pca.de.RData")))
