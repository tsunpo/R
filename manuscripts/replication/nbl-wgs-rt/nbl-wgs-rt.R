# =============================================================================
# Manuscript   : Journey to the centre of the cancer cells
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/sclc-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "ReplicationTiming.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Calculate normalised read counts (RPKM)
# Last Modified: 01/05/17
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "NBL"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.rt    <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")
wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.data  <- file.path(wd.asym,  "data")
wd.asym.plots <- file.path(wd.asym,  "plots")

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.list"), header=F, rownames=F, sep="")

#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_tpm0.RData")))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)

# -----------------------------------------------------------------------------
# Insert size
# Last Modified: 10/09/18
# -----------------------------------------------------------------------------
wd.meta <- file.path(wd, BASE, "metadata/Peifer 2015")
table <- readTable(file.path(wd.meta, "nature14980-s1_Table1.txt"), header=T, rownames=T, sep="\t")

table1 <- table[,c("Insert.Size", "Mean.Coverage")]
table1$Group <- 1
table0 <- table[,c("Insert.Size.1", "Mean.Coverage.1")]
colnames(table0) <- c("Insert.Size", "Mean.Coverage")
table0$Group <- 0
table <- rbind(table1, table0)
table$Group <- as.factor(table$Group)

file.name <- file.path(wd.rt.plots, paste0("boxplot_", base, "_wgs_insert-size.pdf"))
pdf(file.name, height=6, width=3)
boxplot(Insert.Size ~ Group, data=table, outline=F, names=c("Normal", "Tumour"), col=c("dodgerblue", "red"), ylab="Insert size", main=BASE)
dev.off()

median(table1$Insert.Size)
# [1] 323.245
median(table0$Insert.Size)
# [1] 318.325
sd(table1$Insert.Size)
# [1] 58.56802
sd(table0$Insert.Size)
# [1] 58.36447

## Tumour vs. Normal
testW(table1$Insert.Size, table0$Insert.Size)
# [1] 0.5744207






# -----------------------------------------------------------------------------
# Step 6.1: Define replicaiton timing direction for expressed genes (Following Step 4 in "asym-sclc-tx.R" and Step 5 from rt-sclc-wgs.R)
# Link(s):  http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#           https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
BASE  <- "NBL"
PAIR1 <- "T"
PAIR0 <- "N"
PAIR  <- paste0(PAIR1, "-", PAIR0)
#CHR   <- 2
CUTOFF <- 0.2

###
##
bed.gc <- bed[which(bed$GC > 0),]   ## Only keep partitions (in the BED file) with a GC content
ensGene.tx <- ensGene[rownames(tpm.gene.input),]

ensGene.tx.rt <- ensGene.tx[1,]
ensGene.tx.rt$SLOPE_START <- 0
ensGene.tx.rt$SLOPE_END <- 0
ensGene.tx.rt <- ensGene.tx.rt[-1,]
for (c in 1:22) {
   #chr <- chrs[CHR]
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Replication timing
   rpkms.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR, "_n", length(samples), ".txt.gz")), header=T, rownames=T, sep="\t") 
 
   ##
   rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN > -CUTOFF),]
   rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN < CUTOFF),]
   bed.gc.chr <- bed.gc.chr[rownames(rpkms.chr.rt),]
 
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 5000000, 25000000, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 60000000, 80000000, rpkms.chr.rt$MEDIAN, bed.gc.chr,PAIR1, PAIR0, "png")
 
   ## Determin replication direction for each expressed gene
   slopes <- diff(smooth.spline(rpkms.chr.rt$MEDIAN)$y)/diff((bed.gc.chr$START)/1E7)   ## WHY?
 
   ensGene.tx.chr <- subset(ensGene.tx, chromosome_name == chr)
   ensGene.tx.chr$SLOPE_START <- NA
   ensGene.tx.chr$SLOPE_START <- NA
   for (g in 1:nrow(ensGene.tx.chr)) {
      gene <- ensGene.tx.chr[g,]
      bed.s <- getEnsGeneBED(gene$start_position, bed.gc.chr)
      bed.e <- getEnsGeneBED(gene$end_position, bed.gc.chr)
  
      if (length(bed.s) != 0) ensGene.tx.chr$SLOPE_START[g] <- slopes[which(rownames(bed.gc.chr) == bed.s[1])]
      if (length(bed.e) != 0) ensGene.tx.chr$SLOPE_END[g] <- slopes[which(rownames(bed.gc.chr) == bed.e[1])]
   }
   ensGene.tx.rt <- rbind(ensGene.tx.rt, ensGene.tx.chr)
}
save(ensGene.tx.rt, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt.RData")))
# > nrow(ensGene.tx.rt)
# [1] ?   ## tpm.gene_r5_p47
# [1] 9718   ## tpm.gene_tpm0










# -----------------------------------------------------------------------------
# Finalise and plot replication timing
# Name: 
# Last Modified: 05/07/17
# -----------------------------------------------------------------------------
BASE  <- "NBL"
PAIR1 <- "T"
PAIR0 <- "N"
PAIR  <- paste0(PAIR1, "-", PAIR0)
CUTOFF <- 0

###
##
bed.gc <- bed[which(bed$GC > 0),]   ## Only keep partitions (in the BED file) with a GC content
#ensGene.tx <- ensGene[rownames(tpm.gene.input),]

#ensGene.tx.rt <- ensGene.tx[1,]
#ensGene.tx.rt$SLOPE_START <- 0
#ensGene.tx.rt$SLOPE_END <- 0
#ensGene.tx.rt <- ensGene.tx.rt[-1,]
for (c in 1:22) {
   #chr <- chrs[CHR]
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   ## Replication timing
   rpkms.chr <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR, "_n", length(samples), ".txt.gz")), header=T, rownames=T, sep="\t") 
 
   ##
   #rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN > -CUTOFF),]
   #rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN < CUTOFF),]
   bed.gc.chr <- bed.gc.chr[rownames(rpkms.chr),]
 
   plotRT0(wd.rt.plots, BASE, chr, length(samples), NA, NA, rpkms.chr, bed.gc.chr, PAIR1, PAIR0, "png")
}




# -----------------------------------------------------------------------------
# Finalise and plot replication timing
# Name: 
# Last Modified: 05/07/17
# -----------------------------------------------------------------------------
BASE <- "SCLC"
pair1 <- "T"
pair0 <- "N"

for (c in 1:22) {
   chr <- chrs[c]

   rpkms.rt.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv_", chr, "_", pair1, "-", pair0, "_n", length(samples), ".txt.gz"))[,samples]

}





# -----------------------------------------------------------------------------
# Replication timing in NBL
# Name: 
# Last Modified: 30/08/18; 14/06/17
# -----------------------------------------------------------------------------
for (c in 2:2) {
   chr <- chrs[c]

   ## RT (no overlapping)
   rpkms.T.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "/nbl_rpkm.corr.gc.d_", chr, "_T_n", length(samples), ".txt.gz"))[,samples]
   rpkms.N.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "/nbl_rpkm.corr.gc.d_", chr, "_N_n", length(samples), ".txt.gz"))[,samples]

   rpkms.T.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.T.chr), function(x) median(as.numeric(rpkms.T.chr[x,])))
   rpkms.N.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.N.chr), function(x) median(as.numeric(rpkms.N.chr[x,])))
   
   rpkms.T.chr$MEDIAN_FC_LOG2 <- log2(rpkms.T.chr$MEDIAN_FC + 0.01)
   rpkms.N.chr$MEDIAN_FC_LOG2 <- log2(rpkms.N.chr$MEDIAN_FC + 0.01)
   rpkms.T.chr.cut <- rpkms.T.chr[which(rpkms.T.chr$MEDIAN_FC_LOG2 < 9.5),]
   rpkms.T.chr     <- rpkms.T.chr.cut[which(rpkms.T.chr.cut$MEDIAN_FC_LOG2 > 7.5),]
   rpkms.N.chr.cut <- rpkms.N.chr[which(rpkms.N.chr$MEDIAN_FC_LOG2 < 9.5),]
   rpkms.N.chr     <- rpkms.N.chr.cut[which(rpkms.N.chr.cut$MEDIAN_FC_LOG2 > 7.5),]
   
   bed.gc$BED <- rownames(bed.gc)             ## ADD 07/07/17
   bed.gc.chr <- subset(bed.gc, CHR == chr)   ## ADD 07/07/17

   bed.gc.chr.rd  <- bed.gc.chr[rownames(rpkms.T.chr),]
   plotRD(wd.rt.plots, "NBL tumour", "_rpkm.corr.gc.d.rd_", chr, "T", rpkms.T.chr$MEDIAN_FC_LOG2, bed.gc.chr.rd, "png")
   
   bed.gc.chr.rd  <- bed.gc.chr[rownames(rpkms.N.chr),]
   plotRD(wd.rt.plots, "NBL normal", "_rpkm.corr.gc.d.rd_", chr, "N", rpkms.N.chr$MEDIAN_FC_LOG2, bed.gc.chr.rd, "png")
   
   
   
   ## RT (overlaps; 29/08/18)
   rpkms.T.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "/nbl_rpkm.corr.gc.d_", chr, "_T_n", length(samples), ".txt.gz"))[,samples]
   rpkms.N.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "/nbl_rpkm.corr.gc.d_", chr, "_N_n", length(samples), ".txt.gz"))[,samples]
   overlaps <- intersect(rownames(rpkms.T.chr), rownames(rpkms.N.chr))
   rpkms.T.chr <- rpkms.T.chr[overlaps,]
   rpkms.N.chr <- rpkms.N.chr[overlaps,]
   
   ## RT
   rpkms.chr <- toTable(NA, length(samples), length(overlaps), samples)
   rownames(rpkms.chr) <- overlaps
   for (s in 1:length(samples))
      rpkms.chr[,s] <- mapply(x = 1:nrow(rpkms.chr), function(x) getLog2RDRatio(rpkms.T.chr[x, s], rpkms.N.chr[x, s]))
   ## 
   rpkms.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.chr), function(x) median(as.numeric(rpkms.chr[x,])))
   rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN_FC >= -0.25),]
   rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN_FC <= 0.25),]
   bed.gc.chr <- subset(bed.gc, CHR == chr)               ## ADD 07/07/17
   bed.gc.chr.rt <- bed.gc.chr[rownames(rpkms.chr.rt),]   ## ADD 07/07/17
   
   plotRT(paste0(wd.rt.plots, "NBL_RT"), "NBL", chr, length(samples), NA, NA, rpkms.chr.rt$MEDIAN_FC, bed.gc.chr.rt, "T", "N", "png")
   
   ###
   ## RE-VISIT 29/08/18
   rpkms.T.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.T.chr), function(x) median(as.numeric(rpkms.T.chr[x,])))
   rpkms.N.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.N.chr), function(x) median(as.numeric(rpkms.N.chr[x,])))
   bed.gc$BED <- rownames(bed.gc)             ## ADD 07/07/17
   bed.gc.chr <- subset(bed.gc, CHR == chr)   ## ADD 07/07/17
   bed.gc.chr.rd <- bed.gc.chr[overlaps,]     ## ADD 07/07/17
   
   plotRD(wd.rt.plots, "NBL tumour", "_rpkm.corr.gc.d_", chr, "T", log2(rpkms.T.chr$MEDIAN_FC + 0.01), bed.gc.chr.rd, "png")
   plotRD(wd.rt.plots, "NBL normal", "_rpkm.corr.gc.d_", chr, "N", log2(rpkms.N.chr$MEDIAN_FC + 0.01), bed.gc.chr.rd, "png")
   
   ##
   rpkms.T.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.T.chr), function(x) median(as.numeric(rpkms.T.chr[x,])))
   rpkms.N.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.N.chr), function(x) median(as.numeric(rpkms.N.chr[x,])))
   bed.gc$BED <- rownames(bed.gc)             ## ADD 07/07/17
   bed.gc.chr <- subset(bed.gc, CHR == chr)   ## ADD 07/07/17
   
   rpkms.T.chr$MEDIAN_FC_LOG2 <- log2(rpkms.T.chr$MEDIAN_FC + 0.01)
   rpkms.T.chr.cut <- rpkms.T.chr[which(rpkms.T.chr$MEDIAN_FC_LOG2 < 9.5),]
   rpkms.T.chr.cut <- rpkms.T.chr.cut[which(rpkms.T.chr.cut$MEDIAN_FC_LOG2 > 7.5),]
   overlaps.T <- intersect(rownames(rpkms.T.chr.cut), bed.gc.chr$BED) 
   plotRD(wd.rt.plots, "NBL tumour", "_rpkm.corr.gc.d_", chr, "T", rpkms.T.chr.cut$MEDIAN_FC_LOG2, bed.gc.chr[overlaps.T,], "png")      
   
   rpkms.N.chr$MEDIAN_FC_LOG2 <- log2(rpkms.N.chr$MEDIAN_FC + 0.01)
   rpkms.N.chr.cut <- rpkms.N.chr[which(rpkms.N.chr$MEDIAN_FC_LOG2 < 9.5),]
   rpkms.N.chr.cut <- rpkms.N.chr.cut[which(rpkms.N.chr.cut$MEDIAN_FC_LOG2 > 7.5),]
   overlaps.N <- intersect(rownames(rpkms.N.chr.cut), bed.gc.chr$BED) 
   plotRD(wd.rt.plots, "NBL normal", "_rpkm.corr.gc.d_", chr, "N", rpkms.N.chr.cut$MEDIAN_FC_LOG2, bed.gc.chr[overlaps.N,], "png")      
   
   
   
   
   
   
   
   
   
   
   
   for (s in 1:length(samples)) {
      sample <- samples[s]
      plotRD(wd.rt.plots, sample, "_rpkm.corr.gc.d_", chr, "T-N", rpkms.chr[,sample], bed.gc.chr, "png")
   }
   
   ## 
   rpkms.T.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.T.chr), function(x) median(as.numeric(rpkms.T.chr[x,])))
   rpkms.N.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.N.chr), function(x) median(as.numeric(rpkms.N.chr[x,])))
   bed.gc.chr <- bed.gc[rownames(rpkms.T.chr),]
   bed.gc.chr$BED <- rownames(bed.gc.chr)

   #bed.gc.chr <- bed.gc[rownames(rpkms.T.chr),]
   plotRD(wd.rt.plots, "SCLC tumour", "_rpkm.corr.gc.d_", chr, "T", log2(rpkms.T.chr$MEDIAN_FC + 0.01), bed.gc.chr, "png")   
   plotRD(wd.rt.plots, "SCLC normal", "_rpkm.corr.gc.d_", chr, "N", log2(rpkms.N.chr$MEDIAN_FC + 0.01), bed.gc.chr, "png")      
   


   
   ##
   rpkms.chr <- toTable(NA, length(samples), nrow(rpkms.T.chr), samples)
   rownames(rpkms.chr) <- overlaps
   for (s in 1:length(samples))
      rpkms.chr[,s] <- mapply(x = 1:nrow(rpkms.chr), function(x) getLog2RDRatio(rpkms.T.chr[x, s], rpkms.N.chr[x, s]))
   
   ## 
   rpkms.chr$MEDIAN_FC <- mapply(x = 1:nrow(rpkms.chr), function(x) median(as.numeric(rpkms.chr[x,])))
   #rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN_FC > -0.2),]
   #rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN_FC < 0.2),]
   rpkms.chr.rt <- rpkms.chr
   
   bed.gc.chr <- bed.gc[rownames(rpkms.chr.rt),]
   plotRT(paste0(wd.rt.plots, "SCLC_RT"), "SCLC", chr, length(samples), NA, NA, rpkms.chr.rt$MEDIAN_FC, bed.gc.chr, "png")
   
   ## 
   rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN_FC > -0.2),]
   rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN_FC < 0.2),]
   
   bed.gc.chr <- bed.gc[rownames(rpkms.chr.rt),]
   plotRT(paste0(wd.rt.plots, "SCLC_RT"), "SCLC", chr, length(samples), NA, NA, rpkms.chr.rt$MEDIAN_FC, bed.gc.chr, "png")
   
   ## TEST
   sample <- "S00038"
   bed.gc.chr <- bed.gc[overlaps,]
   plotRD(wd.rt.plots, sample, "_rpkm.corr.gc.q_", chr, "T-N", rpkms.chr[,sample], bed.gc.chr, "png")
  
   ## TEST
   sample <- "S00944"
   q <- quantile(rpkm.chr$RPKM_CORR_GC, c(0.001, 0.999))
   #rpkm.chr.q <- rpkm.chr[which(rpkm.chr$RPKM_CORR_GC > q[1]),]
   rpkm.chr.q <- rpkm.chr[which(rpkm.chr$RPKM_CORR_GC < q[2]),]
   bed.gc.chr.q <- bed.gc.chr[rpkm.chr.q$BED,]
   plotRD(wd.rt.plots, sample, "_rpkm.corr.gc.d.q_", chr, "T", rpkms.T.chr.q[,sample]$RPKM_CORR_GC, bed.gc.chr.q, "png")
   
   
}

pipePlotRD <- function(rpkms.chr, bed.gc.chr) {
   q <- quantile(rpkms.chr$RPKM_CORR_GC, c(0.001, 0.999))
   #rpkm.chr.q <- rpkms.chr[which(rpkms.chr$RPKM_CORR_GC > q[1]),]
   rpkms.chr.q <- rpkms.chr[which(rpkms.chr$RPKM_CORR_GC < q[2]),]
   bed.gc.chr.q <- bed.gc.chr[rpkms.chr.q$BED,]
   plotRD(wd.rt.plots, sample, "_rpkm.corr.gc.d.q_", chr, "T", rpkms.T.chr.q[,sample]$RPKM_CORR_GC, bed.gc.chr.q, "png")
 
}



load(file=paste0(wd.rt.data, "sclc_rpkm.corr.gc.d_T_n", length(samples), ".RData"))
rpkms.T <- rpkms.d
load(file=paste0(wd.rt.data, "sclc_rpkm.corr.gc.d_N_n", length(samples), ".RData"))
rpkms.N <- rpkms.d
rm(rpkms.d)

overlaps <- intersect(rownames(rpkms.T), rownames(rpkms.N))
rpkms.T <- rpkms.T[overlaps, samples]
rpkms.N <- rpkms.N[overlaps, samples]

for (c in 1:length(chrs)) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   bed.gc.chr <- intersectbed.gc.chr
   
   rpkms.T.chr <- 
   rpkms.N.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.d_", chr, "_N.txt.gz"))
   overlaps <- intersect(rownames(rpkms.T.chr), rownames(rpkms.N.chr))
   rpkms.T.chr <- rpkms.T.chr[overlaps, samples]
   rpkms.N.chr <- rpkms.N.chr[overlaps, samples]
 
   ##
   rpkms.chr <- rpkms.N.chr
   
   
   bed.gc.chr <- bed.gc[overlaps,]
   rpkms.T.chr$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.T.chr), function(x) median(as.numeric(rpkms.T.chr[x,])))
   rpkms.N.chr$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.N.chr), function(x) median(as.numeric(rpkms.N.chr[x,])))
 
   rt.chr <- log2(rpkms.T.chr$MEDIAN_RPKM + 0.01) - log2(rpkms.N.chr$MEDIAN_RPKM + 0.01)
 
   plotRT(paste0(wd.rt.plots, "SCLC_RT"), "SCLC", chr, NA, NA, rt.chr, bed.gc.chr, "png")
}











# -----------------------------------------------------------------------------
# Most variable replication timing profiles
# Last Modified: 19/06/17
# -----------------------------------------------------------------------------
BASE <- "SCLC"
CUTOFF <- 0.01
PAIR1 <- "T"
PAIR0 <- "N"
PAIR <- paste0(PAIR1, "-", PAIR0)

for (c in 1:22) {
   chr <- chrs[c]
 
   rt.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt_", chr, "_", PAIR, "_n", length(samples), ".txt.gz"))
   if (c == 1)
      rt <- rt.chr
   else
      rt <- rbind(rt, rt.chr)
}
save(rt, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt_", PAIR, "_n", length(samples), ".RData"))

rt$MEAN <- mapply(x = 1:nrow(rt), function(x) mean(as.numeric(rt[x, samples])))
rt$CV2 <- rt$VAR / rt$MEAN^2
rt <- rt[,c(samples, "MEDIAN", "MEAN", "VAR", "CV2")]
save(rt, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv_", PAIR, "_n", length(samples), ".RData"))

# -----------------------------------------------------------------------------
# Most variable replication timing profiles
# Link: http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 19/06/17
# -----------------------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite(c("DESeq", "statmod", "pcaMethods", "fastICA"))
library(DESeq); library(statmod); library(pcaMethods); library(fastICA)

##
ed <- rt
means <- rt$MEAN^2
cv2 <- rt$CV2

pdf(paste0(wd.rt.plots, "scatter_mean2-cv2.pdf"))
par(mar=c(3.5, 3.5, 1,1), mgp=c(2, 0.65, 0), cex=0.9)
smoothScatter(log(means^2), log(cv2))
dev.off()

## Fit a regression line based on the controls:
require(statmod)
minMeanForFit <- unname( quantile( means[ which(cv2 > .3) ], .95) )
useForFit <- means >= minMeanForFit # & spikeins
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

## Add the fit and the 95% confidence interval to our plot:
pdf(paste0(wd.rt.plots, "scatter_mean2-cv2_int.pdf"))
# repeat previous plot
par(mar=c(3.5,3.5,1,1), mgp=c(2,0.65,0), cex=0.9); smoothScatter(log(means), log(cv2));
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- nrow(ed) - 1
# add confidence interval
lines(log(xg), log(vfit * qchisq(0.975,df)/df), lty=2, col="black")
lines(log(xg), log(vfit * qchisq(0.025,df)/df), lty=2, col="black")
dev.off()

## Rank genes by the significance of deviation from the fit
afit <- a1/means+a0
varFitRatio <- rt$VAR/(afit*rt$MEAN^2)
varorder <- order(varFitRatio, decreasing=T)
oed <- ed[varorder,]

# repeat previous plot
pdf(paste0(wd.rt.plots, "scatter_mean2-cv2_int_sig_top132809.pdf"))
par(mar=c(3.5,3.5,1,1), mgp=c(2,0.65,0), cex=0.9); smoothScatter(log(means), log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg), log(vfit * qchisq(0.975,df)/df), lty=2, col="black"); lines(log(xg), log(vfit * qchisq(0.025,df)/df), lty=2, col="black");
# add top 5% (132809 genes)
points(log(means[varorder[1:132809]]), log(cv2[varorder[1:132809]]), col=2)
dev.off()

### 
## Difference with the var() top 5% (Which one is better?)
pdf(paste0(wd.rt.plots, "scatter_mean2-cv2_int_sig_var5_top132809.pdf"))
par(mar=c(3.5,3.5,1,1), mgp=c(2,0.65,0), cex=0.9); smoothScatter(log(means), log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg), log(vfit * qchisq(0.975,df)/df), lty=2, col="black"); lines(log(xg), log(vfit * qchisq(0.025,df)/df), lty=2, col="black");
# add top 5% (132809 genes)
points(log(rt.var5$MEAN^2), log(rt.var5$CV2), col="blue")
dev.off()




# -----------------------------------------------------------------------------
# Thin down top 5% most variable replication timing profiles
# Last Modified: 19/06/17
# -----------------------------------------------------------------------------
## Reports
colnames <- c("N", "N_CV2")
report.rt.var <- toTable(NA, length(colnames), 22, colnames)
rownames(report.rt.var) <- chrs[1:22]

rt.cv2 <- rt[varorder[1:132809],]
save(rt.cv2, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv2.top5_", PAIR, "_n", length(samples), ".RData"))
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   
   rt.chr     <- rt[intersect(rownames(bed.gc.chr), rownames(rt)),]
   rt.var.chr <- rt.var[intersect(rownames(bed.gc.chr), rownames(rt.var)),]
   
   report.rt.var$N[c]     <- nrow(rt.chr)
   report.rt.var$N_CV2[c] <- nrow(rt.var.chr)
}
report.rt.var$RATIO <- report.rt.var$N_CV2 / report.rt.var$N
   
   
## We can also evaluate statistical significance of the deviation
pval <- pchisq(varFitRatio*df, df=df, lower.tail=F)
adj.pval <- p.adjust(pval, "fdr")
sigVariedGenes <- adj.pval<1e-3;
table(sigVariedGenes)
# sigVariedGenes
# FALSE    TRUE 
# 938949 1717221 

#m <- oed[1:50,]
#heatmap(m/apply(m,1,max),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",ColSideColors=ifelse(grepl("ES",colnames(m)),"red","blue"))

## Apply Winsorization procedure
winsorize <- function (x, fraction=0.05) {
   if(length(fraction) != 1 || fraction < 0 ||
      fraction > 0.5) {
      stop("bad value for 'fraction'")
   }
   lim <- quantile(x, probs=c(fraction, 1-fraction))
   x[ x < lim[1] ] <- lim[1]
   x[ x > lim[2] ] <- lim[2]
   x
}

# winsorize to remove 2 most extreme cells (from each side)
wed <- t(apply(ed, 1, winsorize, fraction=2/ncol(ed)))

# now let's recalculate the most variable genes with the winsorized matrix (wed)
means = rowMeans(wed); vars = apply(wed, 1, var); cv2 <- rt$VAR/rt$MEAN^2
useForFit <- means >= unname( quantile( means[ which( cv2 > .3 ) ], .95 ) ) 
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
afit <- fit$coef["a1tilde"]/means+fit$coef["a0"]
vfit <- fit$coef["a1tilde"]/xg+fit$coef["a0"]
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio, decreasing=T)
oed <- wed[varorder,]
# save for the next exercise
save(oed,file="oed_win.RData")

xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=2)





top10 <- quantile(rt$VAR, 0.9)
rt.var <- subset(rt, VAR >= top10)
save(rt.var, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv.var10_", PAIR, "_n", length(samples), ".RData"))

top5 <- quantile(rt$VAR, 0.95)
rt.var5 <- subset(rt, VAR >= top5)
save(rt.var, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv.var5_", PAIR, "_n", length(samples), ".RData"))

top1 <- quantile(rt$VAR, 0.99)
rt.var <- subset(rt, VAR >= top1)
save(rt.var, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv.var1_", PAIR, "_n", length(samples), ".RData"))

# > length(intersect(rownames(rt.var), rownames(rt.var5)))
# [1] 112261
# > nrow(rt.var)
# [1] 132809
# > nrow(rt.var5)
# [1] 132809
# > 112261/132809
# [1] 0.8452816
##
#top5 <- quantile(rt$VAR, 0.05)
#rt.var <- subset(rt, VAR <= top5)
save(rt.var, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv.var5_", PAIR, "_n", length(samples), ".RData"))

# -----------------------------------------------------------------------------
# Sex-specific replication timing variation
# Last Modified: 19/06/17
# -----------------------------------------------------------------------------
pipeDE <- function(expr, pheno.expr, dao, file.de, file.main) {
   ## Differential expression
   de <- differentialAnalysis(expr, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
 
   writeTable(de, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
   writeTable(onlyVariable(de, dao$effect), paste0(file.de, "_Effect-", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
 
   ## Volcano plot
   plotVolcano(de, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)
 
   return(de)
}

filenameDE <- function(wd.de, title, dao) {
   return(paste0(wd.de, title, "_", dao$predictor, "_", dao$test, "_", dao$test.fdr))
}

## Data access object for test parameters
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR:  Q/BH
## DE:   MUT vs WT(female) as factor
dao <- data.frame(test="Wilcox", test.fdr="Q", fdr=0.05, effect=0.4, predictor="sex", predictor.wt="female", stringsAsFactors=F)
load(file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv.var5_", PAIR, "_n", length(samples), ".RData"))

###
##
file.de <- filenameDE(wd.de, "SCLC_DE_rpkm-corr-gc-d-rt-varorder5", dao)
file.main <- c("Replication Timing Variation in SCLC", "Male (n=63) vs Female (n=38)")

##
rt.var <- rt[varorder[1:132809],]
#rt.var <- rt[varorder[1:26562],]
expr <- rt.var[,samples]
pheno.expr <- pheno[samples,]

de.rt.var <- pipeDE(expr, pheno.expr, dao, file.de, file.main)
#plotVolcano(de.rt.var, dao$fdr, dao$effect, file.de, file.main, legend.x=NA)

save(de.rt.var, file=paste0(wd.de, "sclc_de_rpkm.corr.gc.d.rt.var1_sex.RData"))

# -----------------------------------------------------------------------------
# rt-eQTL
# Last Modified: 20/06/17
# -----------------------------------------------------------------------------
library(qvalue)
BASE <- "SCLC"
PAIR <- "T-N"

rt.eqtl <- NA
for (c in 1:22) {
   chr <- chrs[c]

   qtl <- read.rt.eqtl.txt.gz(paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt_eqtl_", chr, "_", PAIR, ".txt.gz"))
   if (is.na(rt.eqtl))
      rt.eqtl <- qtl
   else
      rt.eqtl <- rbind(rt.eqtl, qtl)
}
rt.eqtl$FDR <- qvalue(rt.eqtl$P)$qvalue
rt.eqtl <- rt.eqtl[order(rt.eqtl$P),]

writeTable(rt.eqtl, gzfile(paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt_eqtl_", PAIR, "_n70.txt.gz")), colnames=T, rownames=F, sep="\t")
save(rt.eqtl, file=paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt.cv2.top5_eqtl_", PAIR, "_n70.RData"))




library(qvalue)
rt$VAR_FDR <- qvalue(rt$VAR)$qvalue
rt <- subset(rt$VAR_FDR >= CUTOFF)

## Reports
colnames <- c("N", "N_VAR_FDR")
rt.chr.var <- toTable(NA, length(colnames), 22, colnames)
rownames(rt.chr.var) <- chrs[1:22]






# -----------------------------------------------------------------------------
# Replication timing? or just read depth for each chomosomes in each samples
# Name: ngs/WGS/data/cmd-sclc-rt-?.R (as in command line version)
# Last Modified: 20/05/17
# -----------------------------------------------------------------------------
pipeRT <- function(expr, pheno.expr, dao, file.de, file.main, annot) {
   ## Differential expression
   de <- differentialAnalysis(expr, pheno.expr, dao$predictor, dao$predictor.wt, dao$test, dao$test.fdr)
   de <- cbind(annot[rownames(de),], de)
 
   writeTable(de, paste0(file.de, ".txt"), colnames=T, rownames=T, sep="\t")
   writeTable(onlyVariable(de, dao$effect), paste0(file.de, "_Effect>", dao$effect, ".txt"), colnames=T, rownames=T, sep="\t")   ## Only varible genes
 
   ## Volcano plot

   return(de)
}

getFinalWindows <- function(rpkms.chr.q.sd) {   ## Remove windows if FPKM = NA in any of the samples 
   return(mapply(x = 1:nrow(rpkms.chr.q.sd), function(x) !any(is.na(rpkms.chr.q.sd[x, -1]))))
}

###
## Plot replication timing per chrom per sampe
dao <- data.frame(lower=0.001, upper=0.999, ymin=2^5, ymax=2^10, xmin=0, xmax=0)   ## Data access object for cutoff parameters
# > quantile(rpkm.chr$RPKM_CORR_GC, c(.001,.01,.05,.25,.5,.75,.95,.99,.999))
# 0.1%          1%          5%         25%         50%         75%         95%         99%       99.9% 
# 0.000000    7.737862  216.762235  305.117795  363.932227  423.044198  501.127115  564.204090 1239.837149 
#
# > mean + sd*3
# [1] 641.2868
# > mean + sd*2
# [1] 548.8635
#
# > mean - sd*2
# [1] 179.1704
# > mean - sd*3
# [1] 86.74711

plotRD(wd.rt.plots, sample, "_rpkm.corr.gc_", chr, PAIR, rpkm.chr$RPKM_CORR_GC, bed.gc.chr, "png")

##
rpkm.chr.d <- rpkm.chr[getDetected(rpkm.chr),]   ## ADD 13/06/17


q <- quantile(rpkm.chr$RPKM_CORR_GC, c(dao$lower, dao$upper))
rpkm.chr.q <- rpkm.chr[which(rpkm.chr$RPKM_CORR_GC > q[1]),]
rpkm.chr.q <- rpkm.chr.q[which(rpkm.chr.q$RPKM_CORR_GC < q[2]),]
bed.gc.chr.q <- bed.gc.chr[rpkm.chr.q$BED,]
plotRD(wd.rt.plots, sample, "_rpkm.corr.gc.q_", chr, PAIR, rpkm.chr.q$RPKM_CORR_GC, bed.gc.chr.q, "png")


for (c in 1:length(chrs)) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   dao$xmax <- subset(chromInfo, chrom == chr)$size
 
   rpkms.chr <- readTable(paste0(wd.ngs.data, "sclc_rpkm.corr.gc_", chr, ".txt.gz"), header=T, rownames=T, sep="")
   rpkms.chr.q.sd <- rpkms.chr
   for (s in 1:length(samples)) {
      sample <- samples[s]
  
      rpkm.chr <- rpkms.chr[,c("BED", sample)]
      colnames(rpkm.chr) <- c("BED", "RPKM_CORR_GC")
      q <- quantile(rpkm.chr$RPKM_CORR_GC, c(dao$lower, dao$upper))
      dao$ymin <- q[1]
      dao$ymax <- q[2]
  
      rpkm.chr.q <- rpkm.chr[which(rpkm.chr$RPKM_CORR_GC > q[1]),]
      rpkm.chr.q <- rpkm.chr.q[which(rpkm.chr.q$RPKM_CORR_GC < q[2]),]
      bed.gc.chr.q <- bed.gc.chr[rpkm.chr.q$BED,]
  
      ##
      mean <- mean(rpkm.chr.q$RPKM_CORR_GC)
      sd <- sd(rpkm.chr.q$RPKM_CORR_GC)
      rpkm.chr.q.sd <- rpkm.chr.q[which(rpkm.chr.q$RPKM_CORR_GC < mean + sd*3),]
      rpkm.chr.q.sd <- rpkm.chr.q.sd[which(rpkm.chr.q.sd$RPKM_CORR_GC > mean - sd*3),]
      #bed.gc.chr.q.sd <- bed.gc.chr.q[rpkm.chr.q.sd$BED,]

      #dao$ymax <- max(rpkm.chr.q.sd$RPKM_CORR_GC)
      #dao$ymin <- min(rpkm.chr.q.sd$RPKM_CORR_GC)
      #plotRTPerChromPerSample(wd.rt.plots, chr, sample, rpkm.chr.q.sd$RPKM_CORR_GC, bed.gc.chr.q.sd, dao, "png")
      
      rpkms.chr.q.sd[setdiff(rpkms.chr.q.sd$BED, rpkm.chr.q.sd$BED), sample] <- NA
   }
   
   ##
   #rpkms.chr.q.sd$KEEP <- getFinalWindows(rpkms.chr.q.sd)
   #rpkms.chr.q.sd <- subset(rpkms.chr.q.sd, KEEP == T)[,-c(ncol(rpkms.chr.q.sd))]
   rpkms.chr.q.sd <- removeNAinAnyRows(rpkms.chr.q.sd)
   writeTable(rpkms.chr.q.sd, paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chr, ".txt.gz"), colnames=T, rownames=F, sep="\t")
      
   ##
   bed.gc.chr.q.sd <- bed.gc.chr.q.sd[rpkms.chr.q.sd$BED,]
   rpkms.chr.q.sd$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.chr.q.sd), function(x) median(as.numeric(rpkms.chr.q.sd[x, -1])))
   
   dao$ymin <- log2(min(rpkms.chr.q.sd$MEDIAN_RPKM) + 0.01)
   dao$ymax <- log2(max(rpkms.chr.q.sd$MEDIAN_RPKM) + 0.01)
   plotRT(wd.rt.plots, chr, sample, log2(rpkms.chr.q.sd$MEDIAN_RPKM + 0.01), bed.gc.chr.q.sd, cytoBand.chr, dao, "png")
}

# -----------------------------------------------------------------------------
# Gather filtered mean read depth (in both normal and tumour)
# Last Modified: 02/06/17
# -----------------------------------------------------------------------------
rpkms <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", "chr1", ".txt.gz"))
for (c in 2:22)
   rpkms <- rbind(rpkms, read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chrs[c], ".txt.gz")))
save(rpkms, file=paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_T-FF_n108.RData"))

##
rpkms <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", "chr1", ".txt.gz"))
na <- c()
for (c in 2:22) {
   rpkms.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chrs[c], ".txt.gz"))
   n1 <- nrow(rpkms.chr)
   rpkms2.chr <- removeNAinAnyRows(rpkms.chr)
   n2 <- nrow(rpkms2.chr)
   if(n1 != n2)
      na <- c(na, chrs[c])
   
   rpkms <- rbind(rpkms, rpkms2.chr)
}
save(rpkms, file=paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_T-FF_n108.RData"))

## 
rpkms <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", "chr1", "_N.txt.gz"))
na <- c()
for (c in 2:22) {
   rpkms.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chrs[c], "_N.txt.gz"))
   n1 <- nrow(rpkms.chr)
   rpkms2.chr <- removeNAinAnyRows(rpkms.chr)
   n2 <- nrow(rpkms2.chr)
   if(n1 != n2)
      na <- c(na, chrs[c])
 
   rpkms <- rbind(rpkms, rpkms2.chr)
}
save(rpkms, file=paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_N-FF_n99.RData"))

# -----------------------------------------------------------------------------
# Replication timing? in SCLC (TEST)
# Name: 
# Last Modified: 07/06/17
# -----------------------------------------------------------------------------
for (c in 1:length(chrs)) {
   chr <- chrs[c]

   rpkms.T.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chr, ".txt.gz"))
   rpkms.N.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chr, "_N.txt.gz"))
   overlaps <- intersect(rownames(rpkms.T.chr), rownames(rpkms.N.chr))
   rpkms.T.chr <- rpkms.T.chr[overlaps,]
   rpkms.N.chr <- rpkms.N.chr[overlaps,]
   
   ##
   bed.gc.chr <- bed.gc[overlaps,]
   rpkms.T.chr$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.T.chr), function(x) median(as.numeric(rpkms.T.chr[x,])))
   rpkms.N.chr$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.N.chr), function(x) median(as.numeric(rpkms.N.chr[x,])))
   
   rt.chr <- log2(rpkms.T.chr$MEDIAN_RPKM + 0.01) - log2(rpkms.N.chr$MEDIAN_RPKM + 0.01)
   
   plotRT(paste0(wd.rt.plots, "SCLC_RT"), "SCLC", chr, NA, NA, rt.chr, bed.gc.chr, "png")
}


overlaps <- intersect(rpkms.T, rpkms.T)
rpkms.T <- rpkms.T[]

# -----------------------------------------------------------------------------
# Correlations
# Last Modified: 02/06/17
# -----------------------------------------------------------------------------
rpkms$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms), function(x) median(as.numeric(rpkms[x,])))

cor(rpkms[,1], rpkms$MEDIAN_RPKM, method="kendall") 

cor(rpkms.N[,1], rpkms.N$MEDIAN_RPKM, method="kendall") 

# -----------------------------------------------------------------------------
# PCA (Filtered mean read depth in both normal and tumour)
# Last Modified: 03/06/17
# -----------------------------------------------------------------------------
summaryProportionofVariance <- function(summaries, pc) {
   return(round(summaries[2, pc]*100, 1))
}

plotPCA <- function(x, y, scores, summaries, trait, wd.pca, variable, file.main, legend.x, legend.y, cols, isID, isFlip) {
   trait[is.na(trait)] <- "NA"
   trait.v <- sort(unique(trait))
 
   cols <- cols[1:length(trait.v)]
   trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes
 
   cols[which(trait.v == "NA")] <- "lightgrey"
   trait.col[which(trait == "NA")] <- "lightgrey"
 
   xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
   ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")
 
   if (isFlip)
      scores[,y] <- -scores[,y]
   if (isID)
      pdf(paste0(wd.pca, "pca_", variable, "_", names(scores)[x], "-", names(scores)[y], "_ID.pdf"))
   else
      pdf(paste0(wd.pca, "pca_", variable, "_", names(scores)[x], "-", names(scores)[y], ".pdf"))
   plot(scores[,x], scores[,y], col=trait.col, pch=1, cex=1.5, main=file.main, xlab=xlab, ylab=ylab)
   if (isID)
      text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)

   if (is.na(legend.x))
      legend.x <- min(scores[,x])
   if (is.na(legend.y))
      legend.y <- max(scores[,y])
   legend(legend.x, legend.y, trait.v, col=cols, pch=1, cex=1)   ##bty="n")
   dev.off()
}

###
##
load(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_T-FF_n108.RData"))
pca.rt <- getPCA(t(rpkms))
scores <- pcaScores(pca.rt)
summaries <- pcaSummary(pca.rt)
save(scores, summaries, file=paste0(wd.rt.data, "sclc_pca.rt_T-FF_n108.RData"))

rpkms.C <- rpkms[,setdiff(colnames(rpkms), c("S02397", "S02398", "S02400", "S02401", "S02402", "S02403", "S02404"))]
pca.rt <- getPCA(t(rpkms.C))
scores <- pcaScores(pca.rt)
summaries <- pcaSummary(pca.rt)
save(scores, summaries, file=paste0(wd.rt.data, "sclc_pca.rt_T-FF_n101.RData"))

##
load(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_N-FF_n99.RData"))
pca.rt <- getPCA(t(rpkms))
scores <- pcaScores(pca.rt)
summaries <- pcaSummary(pca.rt)
save(scores, summaries, file=paste0(wd.rt.data, "sclc_pca.rt_N-FF_n99.RData"))

rpkms.C <- rpkms[,setdiff(colnames(rpkms), c("S02397", "S02398", "S02400", "S02401", "S02402", "S02403", "S02404"))]
pca.rt <- getPCA(t(rpkms.C))
scores <- pcaScores(pca.rt)
summaries <- pcaSummary(pca.rt)
save(scores, summaries, file=paste0(wd.rt.data, "sclc_pca.rt_N-FF_n92.RData"))

###
## TUMOURS
pheno <- readTable(paste0(wd.meta, "George 2015/nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
pheno <- pheno[intersect(samples, rownames(pheno)),]
pheno1 <- readTable(paste0(wd.rt.data, "rpkm~size_sclc.txt"), header=T, rownames=T, sep="")[,-1]
pheno <- cbind(pheno1[rownames(pheno),], pheno[,-1])

load(file=paste0(wd.rt.data, "sclc_pca.rt_T-FF_n108.RData"))
file.main <- "Read depth in SCLC (n=108)"
trait <- pheno[,"ethnicity"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T108_Ethnicity", file.main, NA, -440, c("blue", "black"), F, F)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T108_Ethnicity", file.main, NA, -440, c("blue", "black"), T, F)

##
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")

pdf(paste0(wd.rt.plots, "pca_T108_", names(scores)[x], "-", names(scores)[y], ".pdf"))
plot(scores[,x], scores[,y], pch=1, cex=1.5, main="Read depth in SCLC (n=108)", xlab=xlab, ylab=ylab)
#text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

###
##
pheno <- readTable(paste0(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
pheno <- pheno[intersect(samples, rownames(pheno)),]
pheno1 <- readTable(paste0(wd.rt.data, "rpkm~size_sclc.txt"), header=T, rownames=T, sep="")[,-1]
pheno <- cbind(pheno1[rownames(pheno),], pheno[,-1])

load(file=paste0(wd.rt.data, "sclc_pca.rt_T-FF_n101.RData"))
file.main <- "Read depth in SCLC (n=101)"
pheno <- pheno[rownames(scores),]
trait <- pheno[,"ethnicity"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T101_Ethnicity", file.main, 2600, 1550, c("blue", "black"), F, T)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T101_Ethnicity", file.main, 2600, 1550, c("blue", "black"), T, T)

##
trait <- pheno[,"sex"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T101_Sex", file.main, 2600, 1550, c("red", "blue"), F, T)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T101_Sex", file.main, 2600, 1550, c("red", "blue"), T, T)

##
trait <- pheno[,"smoking_status"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T101_Smoking", file.main, 2000, 1550, c("red", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "blue"), F, T)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "T101_Smoking", file.main, 2000, 1550, c("red", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "blue"), T, T)

##
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")

pdf(paste0(wd.rt.plots, "pca_T101_", names(scores)[x], "-", names(scores)[y], ".pdf"))
plot(scores[,x], -scores[,y], pch=1, cex=1.5, main="Read depth in SCLC (n=101)", xlab=xlab, ylab=ylab)
#text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

###
## NORMALS
pheno <- readTable(paste0(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
pheno <- pheno[intersect(samples, rownames(pheno)),]
pheno1 <- readTable(paste0(wd.rt.data, "rpkm~size_sclc.txt"), header=T, rownames=T, sep="")[,-1]
pheno <- cbind(pheno1[rownames(pheno),], pheno[,-1])

load(file=paste0(wd.rt.data, "sclc_pca.rt_N-FF_n99.RData"))
file.main <- "Read depth in SCLC normals (n=99)"
pheno <- pheno[rownames(scores),]
trait <- pheno[,"ethnicity"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "N99_Ethnicity", file.main, NA, -450, c("blue", "black"), F, F)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "N99_Ethnicity", file.main, NA, -450, c("blue", "black"), T, F)

##
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")

pdf(paste0(wd.rt.plots, "pca_N99_", names(scores)[x], "-", names(scores)[y], ".pdf"))
plot(scores[,x], scores[,y], pch=1, cex=1.5, main="Read depth in SCLC normals (n=99)", xlab=xlab, ylab=ylab)
#text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

###
##
pheno <- readTable(paste0(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
pheno <- pheno[intersect(samples, rownames(pheno)),]
pheno1 <- readTable(paste0(wd.rt.data, "rpkm~size_sclc.txt"), header=T, rownames=T, sep="")[,-1]
pheno <- cbind(pheno1[rownames(pheno),], pheno[,-1])

load(file=paste0(wd.rt.data, "sclc_pca.rt_N-FF_n92.RData"))
pheno <- pheno[rownames(scores),]
trait <- pheno[,"ethnicity"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "N92_Ethnicity", file.main, NA, NA, c("blue", "black"), F, T)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "N92_Ethnicity", file.main, NA, NA, c("blue", "black"), T, T)

##
trait <- pheno[,"sex"]
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "N92_Sex", file.main, NA, NA, c("red", "blue"), F, T)
plotPCA(1, 2, scores, summaries, trait, wd.rt.plots, "N92_Sex", file.main, NA, NA, c("red", "blue"), T, T)

##
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")

pdf(paste0(wd.rt.plots, "pca_N92_", names(scores)[x], "-", names(scores)[y], ".pdf"))
plot(scores[,x], -scores[,y], pch=1, cex=1.5, main="Read depth in SCLC normals (n=92)", xlab=xlab, ylab=ylab)
#text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

# -----------------------------------------------------------------------------
# Find if there is any association between RPKMs and insert sizes
# Last Modified: 03/06/17
# -----------------------------------------------------------------------------
pheno <- readTable(paste0(wd.meta, "George 2015/nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
pheno <- pheno[intersect(samples, rownames(pheno)),]
pheno1 <- readTable(paste0(wd.rt.data, "rpkm~size_sclc.txt"), header=T, rownames=T, sep="")[,-1]
pheno <- cbind(pheno1[rownames(pheno),], pheno[,-1])

trait <- pheno[,"ethnicity"]
trait[is.na(trait)] <- "NA"
trait.v <- sort(unique(trait))

cols <- c("blue", "black")
cols <- cols[1:length(trait.v)]
trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes

cols[which(trait.v == "NA")] <- "lightgrey"
trait.col[which(trait == "NA")] <- "lightgrey"

pdf(paste0(wd.rt.data, "rpkm~size_sclc_Ethnicity.pdf"))
plot(MEDIAN_RPKM ~ INSERT_SIZE, data=pheno, col=trait.col,  pch=1, cex=1.5, xlab="Insert size", ylab="Median RPKM", main="SCLC WGS (n=108)")
#text(pheno$INSERT_SIZE, pheno$MEDIAN_RPKM, rownames(pheno), cex=0.6, col="black", pos=3)
lm <- lm(pheno$MEDIAN_RPKM ~ pheno$INSERT_SIZE)
abline(lm[[1]][[1]], lm[[1]][[2]])

legend(395, 411.5, trait.v, col=cols, pch=1, cex=1)   ##bty="n")
#text(412, max(pheno$MEDIAN_RPKM) - 3.7, "p-value=0.7633", cex=0.6, pos=3) 
#text(min(pheno$INSERT_SIZE) + 35, max(pheno$MEDIAN_RPKM) - 2.5, "p-value=9.704e-07", cex=0.6, pos=3) 
dev.off()

##
pdf(paste0(wd.rt.data, "rpkm~size_sclc.pdf"))
plot(MEDIAN_RPKM ~ INSERT_SIZE, data=pheno, pch=1, cex=1.5, xlab="Insert size", ylab="Median RPKM", main="SCLC WGS (n=108)")
text(pheno$INSERT_SIZE, pheno$MEDIAN_RPKM, rownames(pheno), cex=0.6, col="black", pos=3)
lm <- lm(pheno$MEDIAN_RPKM ~ pheno$INSERT_SIZE)
abline(lm[[1]][[1]], lm[[1]][[2]])

#text(412, max(pheno$MEDIAN_RPKM) - 3.7, "p-value=0.7633", cex=0.6, pos=3) 
#text(min(pheno$INSERT_SIZE) + 35, max(pheno$MEDIAN_RPKM) - 2.5, "p-value=9.704e-07", cex=0.6, pos=3) 
dev.off()

###
##
trait <- pheno[,"sex"]
trait[is.na(trait)] <- "NA"
trait.v <- sort(unique(trait))

cols <- c("red", "blue")
cols <- cols[1:length(trait.v)]
trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes

cols[which(trait.v == "NA")] <- "lightgrey"
trait.col[which(trait == "NA")] <- "lightgrey"

pdf(paste0(wd.rt.data, "rpkm~size_sclc_Sex.pdf"))
plot(MEDIAN_RPKM ~ INSERT_SIZE, data=pheno, col=trait.col,  pch=1, cex=1.5, xlab="Insert size", ylab="Median RPKM", main="SCLC WGS (n=108)")
#text(pheno$INSERT_SIZE, pheno$MEDIAN_RPKM, rownames(pheno), cex=0.6, col="black", pos=3)
lm <- lm(pheno$MEDIAN_RPKM ~ pheno$INSERT_SIZE)
abline(lm[[1]][[1]], lm[[1]][[2]])

legend(395, 411.5, trait.v, col=cols, pch=1, cex=1)   ##bty="n")
#text(412, max(pheno$MEDIAN_RPKM) - 3.7, "p-value=0.7633", cex=0.6, pos=3) 
#text(min(pheno$INSERT_SIZE) + 35, max(pheno$MEDIAN_RPKM) - 2.5, "p-value=9.704e-07", cex=0.6, pos=3) 
dev.off()

# -----------------------------------------------------------------------------
# Koren et al., 2012
# Last Modified: 06/06/17
# -----------------------------------------------------------------------------
plotRTKoren <- function(wd.rt.plots.filename, title, chr, xmin, xmax, koren.chr, ext) {
   main.text <- " read depth"
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Read depth"
 
   ymin <- min(koren.chr$Fraction)
   ymax <- max(koren.chr$Fraction)
 
   wd.rt.plots.filename <- paste0(wd.rt.plots.filename, chr)
   png(paste0(wd.rt.plots.filename, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
 
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=paste0(title, main.text))
   points(koren.chr$Pos/1E6, koren.chr$Fraction, col="red", cex=0.3)
   lines(koren.chr$Pos/1E6, smooth.spline(koren.chr$Fraction)$y)
 
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.5, col="lightgrey")
 
   dev.off()
}

###
##
koren <- readTable(paste0(wd.meta, "Koren 2012/Koren et al Table S2.txt"), header=F, rownames=F, sep="")
colnames(koren) <- c("Chr", "Pos", "Fraction")
koren$Pos <- as.numeric(koren$Pos)

koren.chr <- subset(koren, Chr == "2")
koren.chr <- koren.chr[!is.na(koren.chr$Fraction),]
plotRTKoren(paste0(wd.rt.plots, "png/koren_table-s2_"), "Koren et al., 2012", "chr2", 0, 249250621, koren.chr, "png")









##
load(file=paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_N-FF_n92.RData"))


x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", summaryProportionofVariance(summaries, x), "%)")
ylab <- paste0("Principal component ", y, " (", summaryProportionofVariance(summaries, y), "%)")

pdf(paste0(filenameRT(wd.rt.plots), "pca_", names(scores)[x], "-", names(scores)[y], "_T-FF_n108.pdf"))
plot(scores[,x], scores[,y], pch=1, cex=1.5, main="Read depth in SCLC (n=108)", xlab=xlab, ylab=ylab)
#text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

## PCA
#pca.rt <- getPCA(t(rpkms))
scores <- pcaScores(pca.rt)
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", pcaProportionofVariance(pca, x), "%)")
ylab <- paste0("Principal component ", y, " (", pcaProportionofVariance(pca, y), "%)")


pdf(paste0(filenameRT(wd.rt.plots), "pca_", names(scores)[x], "-", names(scores)[y], ".pdf"))
plot(scores[,x], scores[,y], pch=1, cex=1.5, main="PCA of read depth in SCLC", xlab=xlab, ylab=ylab)
text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

pdf(paste0(filenameRT(wd.rt.plots), "pca_", names(scores)[x], "-", names(scores)[y], "_N.pdf"))
plot(scores[,x], scores[,y], pch=1, cex=1.5, main="PCA on read depth in SCLC normals (n=99)", xlab=xlab, ylab=ylab)
text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

###
## PCA (Removed 7 Koreans)
rpkms.C <- rpkms[,setdiff(colnames(rpkms), c("S02397", "S02398", "S02400", "S02401", "S02402", "S02403", "S02404"))]

pca.rt <- getPCA(t(rpkms.C))
scores <- pcaScores(pca.rt)
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", pcaProportionofVariance(pca, x), "%)")
ylab <- paste0("Principal component ", y, " (", pcaProportionofVariance(pca, y), "%)")

pdf(paste0(filenameRT(wd.rt.plots), "pca_", names(scores)[x], "-", names(scores)[y], "_T101.pdf"))
plot(scores[,x], -scores[,y], pch=1, cex=1.5, main="PCA on read depth in SCLC (n=101)", xlab=xlab, ylab=ylab)
text(scores[,x], -scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()

###
## PCA (Removed 7 Koreans)
rpkms.C <- rpkms[,setdiff(colnames(rpkms), c("S02397", "S02398", "S02400", "S02401", "S02402", "S02403", "S02404"))]

pca <- getPCA(t(rpkms.C))
scores <- pcaScores(pca)
x <- 1
y <- 2
xlab <- paste0("Principal component ", x, " (", pcaProportionofVariance(pca, x), "%)")
ylab <- paste0("Principal component ", y, " (", pcaProportionofVariance(pca, y), "%)")

pdf(paste0(filenameRT(wd.rt.plots), "pca_", names(scores)[x], "-", names(scores)[y], "_N92.pdf"))
plot(scores[,x], -scores[,y], pch=1, cex=1.5, main="PCA on read depth in SCLC normals (n=92)", xlab=xlab, ylab=ylab)
text(scores[,x], -scores[,y], rownames(scores), cex=0.6, col="black", pos=3)
dev.off()




# -----------------------------------------------------------------------------
# Replication timing for each chomosomes (Filtered mean read depth in both normal and tumour)
# Last Modified: 02/06/17
# -----------------------------------------------------------------------------
filenameRT <- function(wd.rt.plots) {
   return(paste0(wd.rt.plots, "png/sclc_rpkm.corr.gc.q_"))
}

for (c in 1:length(chrs)) {
   chr <- chrs[c]

   ## Read filtered (*.q.sd) coverage files
   rpkms.chr <- read.rpkm.corr.gc.txt.gz(paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chr, ".txt.gz"))

   ## PCA
   pca <- getPCA(t(rpkms.chr))
   scores <- pcaScores(pca)
   xlab <- paste0("Principal component ", x, " (", pcaProportionofVariance(pca, x), "%)")
   ylab <- paste0("Principal component ", y, " (", pcaProportionofVariance(pca, y), "%)")
   
   pdf(paste0(filenameRT(wd.rt.plots), "pca_", chr, "_", names(scores)[x], "-", names(scores)[y], ".pdf"))
   plot(scores[,x], scores[,y], pch=16, cex=1.5, main="PCA of read depth (Chromosome 1)", xlab=xlab, ylab=ylab)
   text(scores[,x], scores[,y], rownames(scores), cex=0.6, col="black", pos=3) 
   dev.off()
   
   ##
   #bed.gc.chr.q.sd <- bed.gc.chr[rpkms.chr.q.sd$BED,]
   rpkms.chr$MEDIAN_RPKM   <- mapply(x = 1:nrow(rpkms.chr),   function(x) median(as.numeric(rpkms.chr[x,])))

   #dao$xmax <- subset(chromInfo, chrom == chr)$size
   plotRT(filenameRT(wd.rt.plots), chr, NA, NA, sample, rpkms.chr, "png")
}



plotRT(filenameRT(wd.rt.plots), "SCLC median", chr, 34600000, 44100000, rpkms.chr, bed.gc, cytoBand, "png")

#plotRT(paste0(wd.rt.plots, "png/sclc_rpkm.corr.gc.q_"), "chr2", 50000000, 100000000, "SCLC median", rpkms.chr, bed.gc, cytoBand, "png")

 
# -----------------------------------------------------------------------------
# Molecular classification / PCA
# Last Modified: 02/06/17
# -----------------------------------------------------------------------------
pca.de      <- getPCA(t(tpm.gene.pcg.exp2.log2[rownames(significantAndVariable(de.tpm.gene.pcg.exp2,      dao$effect, dao$fdr)),]))   ## BUG FIX 13/02/17: Perform PCA using normalised data

###
## RB1 status on D.E genes
pheno.expr <- pheno[samples.expr,]

trait <- pheno.expr[,"RB1NEW"]
trait[which(trait == 0)] <- "WT"
trait[which(trait == 1)] <- "RB1"

file.main <- "LCNEC RB1 Status on 54 D.E. Genes"
plotPCA(1, 2, pca.de, trait, wd.de, "RB1_54DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))

file.main <- c("LCNEC RB1 Status on 20 D.E. Genes", "Two-Subtype Residual")
plotPCA(1, 2, pca.de.res2, trait, wd.de, "RB1_20DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))

file.main <- c("LCNEC RB1 Status on 7 D.E. Genes", "Four-Subtype Residual")
plotPCA(1, 2, pca.de.res, trait, wd.de, "RB1_7DE", file.main, NA, NA, c("purple", "red", "deepskyblue", "blue"))

## Find outliers
scores <- pcaScores(pca.de)

outliers <- subset(scores[,c("PC1", "PC2")], PC1 < -2.5)
subset(cbind(outliers, pheno.expr[rownames(outliers),]), RB1NEW == 0)

outliers <- subset(scores[,c("PC1", "PC2")], PC1 > -2.5)
subset(cbind(outliers, pheno.expr[rownames(outliers),]), RB1NEW == 1)







 
for (c in 1:length(chrs)) {
   chr <- chrs[c]
   dao$xmax <- subset(chromInfo, chrom == chr)$size
   
   rpkms.chr <- readTable(paste0(wd.ngs.data, "sclc_rpkm.corr.gc_", chr, ".txt.gz"), header=T, rownames=T, sep="")
   rpkms.chr.q <- rpkms.chr
   #rpkms.chr.q.sd <- rpkms.chr
   #rpkms.chr.q.sd2 <- rpkms.chr
   for (s in 1:length(samples)) {
      sample <- samples[s]
      
      rpkm.chr <- rpkms.chr[,c("BED", sample)]
      colnames(rpkm.chr) <- c("BED", "RPKM_CORR_GC")
      q <- quantile(rpkm.chr$RPKM_CORR_GC, c(dao$lower, dao$upper))
      dao$ymin <- q[1]
      dao$ymax <- q[2]
            
      rpkm.chr.q <- rpkm.chr[which(rpkm.chr$RPKM_CORR_GC > q[1]),]
      rpkm.chr.q <- rpkm.chr.q[which(rpkm.chr.q$RPKM_CORR_GC < q[2]),]
      bed.gc.chr.q <- bed.gc.chr[rpkm.chr.q$BED,]
      
      #plotRTPerChromPerSample(wd.rt.plots, chr, sample, rpkm.chr.q$RPKM_CORR_GC, bed.gc.chr.q, dao, "pdf")
      plotRTPerChromPerSample(wd.rt.plots, chr, sample, rpkm.chr.q$RPKM_CORR_GC, bed.gc.chr.q, dao, "png")
      
      rpkms.chr.q[setdiff(rpkms.chr.q$BED, rpkm.chr.q$BED), sample] <- NA
      
      ##
      #mean <- mean(rpkm.chr.q$RPKM_CORR_GC)
      #sd <- sd(rpkm.chr.q$RPKM_CORR_GC)
      #rpkm.chr.q.sd <- rpkm.chr.q[which(rpkm.chr.q$RPKM_CORR_GC < mean + sd*3),]
      #bed.gc.chr.q.sd <- bed.gc.chr.q[rpkm.chr.q.sd$BED,]
      
      #dao$ymax <- max(rpkm.chr.q.sd$RPKM_CORR_GC)
      #plotRTPerChromPerSample(wd.rt.plots, chr, sample, rpkm.chr.q.sd$RPKM_CORR_GC, bed.gc.chr.q.sd, dao, "png")
      
      ##
      #rpkm.chr.q.sd2 <- rpkm.chr.q.sd[which(rpkm.chr.q.sd$RPKM_CORR_GC > q[2]),]
      #bed.gc.chr.q.sd2 <- bed.gc.chr.q.sd[rpkm.chr.q.sd2$BED,]

      #dao$ymin <- min(rpkm.chr.q.sd2$RPKM_CORR_GC)
      #plotRTPerChromPerSample2(wd.rt.plots, chr, sample, rpkm.chr.q.sd2$RPKM_CORR_GC, bed.gc.chr.q.sd2, dao, "png")
      
      ##
      #rpkms.chr.q.sd[setdiff(rpkms.chr.q.sd$BED, rpkm.chr.q.sd$BED), sample] <- NA
      #rpkms.chr.q.sd2[setdiff(rpkms.chr.q.sd2$BED, rpkm.chr.q.sd2$BED), sample] <- NA
   }
   
   ##
   rpkms.chr.q$KEEP <- getFinalWindows(rpkms.chr.q)
   rpkms.chr.q <- subset(rpkms.chr.q, KEEP == T)[,-c(ncol(rpkms.chr.q))]
   writeTable(rpkms.chr.q, paste0(wd.rt.data, "sclc_rpkm.corr.gc.q_", chr, ".txt.gz"), colnames=T, rownames=F, sep="\t")
   
   ##
   bed.gc.chr.q <- bed.gc.chr.q[rpkms.chr.q$BED,]
   dao$ymin <- min(rpkms.chr.q$MEDIAN_RPKM)
   dao$ymax <- max(rpkms.chr.q$MEDIAN_RPKM)
   
   rpkms.chr.q$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.chr.q), function(x) median(as.numeric(rpkms.chr.q[x, -1])))
   plotRT(wd.rt.plots, chr, sample, rpkms.chr.q$MEDIAN_RPKM, bed.gc.chr.q, dao, "png")
}   
   ##
   rpkms.chr.q.sd$KEEP <- getFinalWindows(rpkms.chr.q.sd)
   rpkms.chr.q.sd <- subset(rpkms.chr.q.sd, KEEP == T)[,-c(ncol(rpkms.chr.q.sd))]
   writeTable(rpkms.chr.q.sd, paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd_", chr, ".txt.gz"), colnames=T, rownames=F, sep="\t")

   ##
   bed.gc.chr.q.sd <- bed.gc.chr.q.sd[rpkms.chr.q.sd$BED,]
   dao$ymin <- min(rpkms.chr.q.sd$MEDIAN_RPKM)
   dao$ymax <- max(rpkms.chr.q.sd$MEDIAN_RPKM)
   
   rpkms.chr.q.sd$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.chr.q.sd), function(x) median(as.numeric(rpkms.chr.q.sd[x, -1])))
   plotRT(wd.rt.plots, chr, sample, rpkms.chr.q.sd$MEDIAN_RPKM, bed.gc.chr.q.sd, dao, "png")
   
   ####
   ##
   rpkms.chr.q.sd2$KEEP <- getFinalWindows(rpkms.chr.q.sd2)
   rpkms.chr.q.sd2 <- subset(rpkms.chr.q.sd2, KEEP == T)[,-c(ncol(rpkms.chr.q.sd2))]
   writeTable(rpkms.chr.q.sd2, paste0(wd.rt.data, "sclc_rpkm.corr.gc.q.sd2_", chr, ".txt.gz"), colnames=T, rownames=F, sep="\t")
   
   ##
   bed.gc.chr.q.sd2 <- bed.gc.chr.q.sd2[rpkms.chr.q.sd2$BED,]
   dao$ymin <- min(rpkms.chr.q.sd2$MEDIAN_RPKM)
   dao$ymax <- max(rpkms.chr.q.sd2$MEDIAN_RPKM)
   
   rpkms.chr.q.sd2$MEDIAN_RPKM <- mapply(x = 1:nrow(rpkms.chr.q.sd2), function(x) median(as.numeric(rpkms.chr.q.sd2[x, -1])))
   plotRT2(wd.rt.plots, chr, sample, rpkms.chr.q.sd2$MEDIAN_RPKM, bed.gc.chr.q.sd2, dao, "png")
}

# -----------------------------------------------------------------------------
# Replication timing for each chomosomes (mean read depth)
# Last Modified: 18/05/17
# -----------------------------------------------------------------------------
for (c in 1:length(chrs)) {
   chr <- chrs[c]
   
   ## Read coverage from PeifLyne
   rt.chr <- readTable(paste0(wd.ngs.data, "sclc_rpkm.corr.gc_", chr, ".txt.gz"), header=T, rownames=T, sep="")[,-1]
   #plotReplicationTime(wd.rt.plots, chr, rt.chr, "SCLC", pdf")
   plotReplicationTime(wd.rt.plots, chr, rt.chr, "SCLC", "png")
   
   ## Read coverage from bedtools
   rt2.chr <- readTable(paste0(wd.ngs.data2, "sclc_rpkm.corr.gc_", chr, ".txt.gz"), header=T, rownames=T, sep="")[,-1]
   #plotReplicationTime(wd.rt.plot2, chr, rt2.chr, "SCLC", "pdf")
   plotReplicationTime(wd.rt.plots2, chr, rt2.chr, "SCLC", "png")
   
   ###
   ## TEST
   for (s in 1:length(samples)) {
      sample <- samples[s]
      
      #png(paste0(wd.rt.data, "plot_rpkm.corr.gc_PL-Bedtools_", sample, ".png"))
      #plot(log2(rt.chr[,s] + 0.01)~log2(rt2.chr[,s] + 0.01), xlab="peiflyne", ylab="bedtools multicov q1", main=sample)
      #dev.off()
      
      plotRTPerChromPerSample(wd.rt.plots, chr, sample, rpkm.T.chr.q$RPKM_CORR_GC, bed.gc.chr.q, chromInfo, cutoff, "png")
   }

}






# -----------------------------------------------------------------------------
# Find if there is any association between RPKMs and repeat masks
# Name: ngs/WGS/data/cmd-sclc-rt.R (as in command line version)
# Last Modified: 22/05/17
# -----------------------------------------------------------------------------
getRPKMForRMSK <- function(bed.gc.chr, rmsk.low.chr) {
 bed.gc.chr.end <- bed.gc.chr[which(bed.gc.chr$END > rmsk.low.chr$genoStart),]   ## Important! Not >= seg$START)
 bed.gc.chr.end.start <- bed.gc.chr.end[which(bed.gc.chr.end$START <= rmsk.low.chr$genoEnd),]
 
 return(max(bed.gc.chr.end.start[,6]))
}

#cbind(bed.gc[rpkm.chr[16270:16530,]$BED,], rpkm.chr[16270:16530,])
rmsk <- subset(rmsk, repClass=="Low_complexity")

for (c in 1:length(chrs)) {
 chr <- chrs[c]
 rmsk.chr <- subset(rmsk, genoName == chr)
 
 rpkms.chr <- readTable(paste0(wd.ngs.data, "sclc_rpkm.corr.gc_", chr, ".txt.gz"), header=T, rownames=T, sep="")
 for (s in 1:length(samples)) {
  sample <- samples[s]
  rpkms.chr.sample <- rpkms.chr[,c("BED", sample)]
  bed.gc.chr <- cbind(bed.gc[rpkms.chr.sample$BED,], rpkms.chr.sample)
  
  rmsk.chr$RPKM_CORR_GC <- mapply(r = 1:nrow(rmsk.chr), function(r) getRPKMForRMSK(bed.gc.chr, rmsk.chr[r,]))
  rmsk.chr <- rmsk.low.chr[order(rmsk.low.chr$RPKM_CORR_GC, decreasing=T),]
  rmsk.chr$repName <- gsub("-rich", "", rmsk.low.chr$repName)
  rmsk.chr$repName <- gsub("_rich", "", rmsk.low.chr$repName)
  
  png(paste0(wd.rt.plots, "repClass~rpkm.corr.gc_", chr, "_", sample, "_repClass.png"), height=4, width=21, units="in", res=300)
  boxplot(log2(RPKM_CORR_GC + 0.01) ~ repClass, data=rmsk.chr, ylab="repClass", xlab="Log2(RPKM_CORR_GC)")
  dev.off()
  
  png(paste0(wd.rt.plots, "repFamily~rpkm.corr.gc_", chr, "_", sample, "_repFamily.png"), height=4, width=50, units="in", res=300)
  boxplot(log2(RPKM_CORR_GC + 0.01) ~ repFamily, data=rmsk.chr, ylab="repFamily", xlab="Log2(RPKM_CORR_GC)")
  dev.off()
 } 
}

# -----------------------------------------------------------------------------
# Multithreads
# Last Modified: 10/05/17
# -----------------------------------------------------------------------------
for (s in 1:length(samples)) {
 print(paste0("sbatch --cpus-per-task=2 --mem=2000mb --time=02:00:00 --account=UniKoeln R CMD BATCH --no-save --no-restore '--args ", s, "' sclc-cmd-rpkm-corr.R logs/", samples[s], ".out"))
 print("wait")
}

# -----------------------------------------------------------------------------
# GC contents
# Last Modified: 12/05/17
# -----------------------------------------------------------------------------
load(paste0(wd.rt.data, "human-genome.1kb-grid.bed.gc.RData"))
gc.mean <- mean(bed.gc$GC)

samples <- readTable(paste0(wd.ngs, "sclc_wgs_n101+7.list"), header=F, rownames=F, sep="")
chrs <- paste0("chr", c(1:22, "X", "Y"))
#writeTable(chrs, paste0(wd.ngs, "chrs.list"), colnames=F, rownames=F, sep="")
for (c in 1:length(chrs)) {
   chr <- chrs[c]
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   pdf(paste0(wd.rt.plots, "gc.smooth_", chr, ".pdf"), height=4, width=10)
   plot(NULL, ylim=c(0, 1), xlim=c(0, subset(chromInfo, chrom == chr)$size/1E6), xlab=xlab.text)
   #points(bed.gc.chr$START/1E6, rpkm.T.chr, col="red", cex=0.3)
   lines(bed.gc.chr$START/1E6, smooth.spline(bed.gc.chr$GC)$y)
   dev.off()
 
   pdf(paste0(wd.rt.plots, "gc_", chr, ".pdf"), height=4, width=10)
   plot(NULL, ylim=c(0, 1), xlim=c(0, subset(chromInfo, chrom == chr)$size/1E6), xlab=xlab.text)
   #points(bed.gc.chr$START/1E6, rpkm.T.chr, col="red", cex=0.3)
   lines(bed.gc.chr$START/1E6, bed.gc.chr$GC)
   dev.off()
}






# -----------------------------------------------------------------------------
# Read counts ratio?
# Last Modified: 11/05/17
# -----------------------------------------------------------------------------
samples <- readTable(paste0(wd.ngs, "sclc_wgs_n101+7.list"), header=F, rownames=F, sep="")
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   rpkm.T <- readTable(paste0(wd.ngs.data, sample, "/", sample, "_T.bam.rpkm.corr.gc.txt.gz"), header=T, rownames=T, sep="")
   rpkm.N <- readTable(paste0(wd.ngs.data, sample, "/", sample, "_N.bam.rpkm.txt.gz"), header=T, rownames=T, sep="")   
   rpkm.N <- rpkm.N[rownames(rpkm.T),]
   
   rpkm.T <- subset(rpkm.T, RPKM != 0)
   rpkm.N <- subset(rpkm.N, RPKM != 0)
   overlaps <- intersect(rownames(rpkm.T), rownames(rpkm.N))
   rpkm.T <- rpkm.T[overlaps,]
   rpkm.N <- rpkm.N[overlaps,]
   
   rpkm.T$RPKM_LOG2         <- log2(rpkm.T$RPKM + 0.01)
   rpkm.T$RPKM_CORR_LOG2    <- log2(rpkm.T$RPKM_CORR + 0.01)
   rpkm.T$RPKM_CORR_GC_LOG2 <- log2(rpkm.T$RPKM_CORR_GC + 0.01)
   rpkm.N$RPKM_LOG2         <- log2(rpkm.N$RPKM + 0.01)
   
   ###
   ##
   chr1 <- rpkm.T[rownames(subset(bed.gc, CHR == "chr1")),]$RPKM_CORR_GC
   chr7 <- rpkm.T[rownames(subset(bed.gc, CHR == "chr7")),]$RPKM_CORR_GC
   
   pdf(paste0(wd.rt.data, "chr7_rpkm.cln.log2_", sample, "_T.corr.gc.pdf"))
   plot(NULL,ylim=c(2^5,2^10),xlim=c(0,length(chr7)))
   points(chr7,col="red",cex=0.3)
   lines(smooth.spline(chr7))
   dev.off()
   
   ###
   ##
   rpkm.log.fc <- rpkm.T$RPKM_CORR_LOG2- rpkm.N$RPKM_LOG2
   pdf(paste0(wd.rt.data, "density_rpkm.cln.log2_", sample, "_T.corr-N.pdf"))
   plot(density(rpkm.log.fc), xlab="Effect size (log2 FC)", ylab="Frequency", main="Tumour$RPKM_CORR - Normal$RPKM")
   dev.off()

   rpkm.log.fc <- rpkm.T$RPKM_LOG2- rpkm.N$RPKM_LOG2   
   pdf(paste0(wd.rt.data, "density_rpkm.cln.log2_", sample, "_T-N.pdf"))
   plot(density(rpkm.log.fc), xlab="Effect size (log2 FC)", ylab="Frequency", main="Tumour$RPKM - Normal$RPKM")
   dev.off()
   
   ###
   ##
   idx <- which(rpkm.log.fc >= 0.5)
   
   pdf(paste0(wd.rt.data, "hist_rpkm.corr.gc.log2_", sample, "_T.pdf"))
   hist(rpkm.T$RPKM_CORR_GC_LOG2, xlab="log2(RPKM)", ylab="Frequency", main="Tumour$RPKM_CORR_GC")
   dev.off()
   
   pdf(paste0(wd.rt.data, "hist_rpkm.corr.gc.log2_", sample, "_T.RT.pdf"))
   hist(rpkm.T$RPKM_CORR_GC_LOG2[idx], xlab="log2(RPKM)", ylab="Frequency", main="RT Tumour$RPKM_CORR_GC")
   dev.off()
   
   ##
   pdf(paste0(wd.rt.data, "density_rpkm.corr.gc.log2_", sample, "_T.pdf"))
   plot(density(rpkm.T$RPKM_CORR_GC_LOG2), xlab="log2(RPKM)", ylab="Frequency", main="Tumour$RPKM_CORR_GC")
   dev.off()
   
   pdf(paste0(wd.rt.data, "density_rpkm.corr.gc.log2_", sample, "_T.RT.pdf"))
   plot(density(rpkm.T$RPKM_CORR_GC_LOG2[idx]), xlab="log2(RPKM)", ylab="Frequency", main="RT Tumour$RPKM_CORR_GC")
   dev.off()
   
   pdf(paste0(wd.rt.data, "density_rpkm.cln.log2.idx_", sample, "_T.corr-N.pdf"))
   plot(density(rpkm.log.fc[idx]), xlab="Effect size (log2 FC)", ylab="Frequency", main="Tumour$RPKM_CORR - Normal$RPKM")
   dev.off()
}




   rpkm.corr.gc.ratio$RPKM_CORR_LOG2 <- log2(rpkm.corr.gc.ratio$RPKM_CORR + 0.01)
   rpkm.corr.gc.ratio$RPKM_CORR_GC_LOG2 <- log2(rpkm.corr.gc.ratio$RPKM_CORR_GC + 0.01)

   ##
   rpkm2.corr.gc.ratio <- rpkm.corr.gc.ratio
   
   rpkm.corr.gc.ratio <- readTable(paste0(wd.ngs.data, sample, "/", sample, "_T.bam.rpkm.corr.gc.txt.gz"), header=T, rownames=F, sep="")
   rpkm.corr.gc.ratio$RPKM_CORR_LOG2 <- log2(rpkm.corr.gc.ratio$RPKM_CORR + 0.01)
   rpkm.corr.gc.ratio$RPKM_CORR_GC_LOG2 <- log2(rpkm.corr.gc.ratio$RPKM_CORR_GC + 0.01)
   
   ###
   ##
   pdf(paste0(wd.rt.data, "hist_rpkm.corr.gc_CORR_", sample, ".pdf"))
   hist(rpkm.corr.gc.ratio$RPKM_CORR_LOG2, xlab="log2(RPKM_CORR)", ylab="Frequency", main="Copy Number-corrected RPKM")
   #plot(density(rpkm.corr.gc.ratio$RPKM_CORR_LOG2), xlab="log2(RPKM_CORR)", ylab="Frequency", main="Copy Number-corrected RPKM")
   dev.off()
   
   pdf(paste0(wd.rt.data, "hist_rpkm.corr.gc_CORR_GC_", sample, ".pdf"))
   hist(rpkm.corr.gc.ratio$RPKM_CORR_GC_LOG2, xlab="log2(RPKM_CORR_GC)", ylab="Frequency", main="Copy Number-, GC-corrected RPKM")
   dev.off()
   
   ##
   rpkm.corr.gc.ratio.log2 <- subset(rpkm.corr.gc.ratio, RPKM_CORR_GC_LOG2 >= 5)
   rpkm.corr.gc.ratio.log2 <- subset(rpkm.corr.gc.ratio.log2, RPKM_CORR_GC_LOG2 <= 10)
   
   pdf(paste0(wd.rt.data, "hist_rpkm.corr.gc.log2_CORR_", sample, ".pdf"))
   hist(rpkm.corr.gc.ratio.log2$RPKM_CORR_LOG2, xlab="log2(RPKM_CORR)", ylab="Frequency", main="Copy Number-corrected RPKM")
   dev.off()
   
   pdf(paste0(wd.rt.data, "hist_rpkm.corr.gc.log2_CORR_GC_", sample, ".pdf"))
   hist(rpkm.corr.gc.ratio.log2$RPKM_CORR_GC_LOG2, xlab="log2(RPKM_CORR_GC)", ylab="Frequency", main="Copy Number-, GC-corrected RPKM")
   dev.off()
   
   pdf(paste0(wd.rt.data, "density_rpkm.corr.gc.log2_CORR_", sample, ".pdf"))
   plot(density(rpkm.corr.gc.ratio.log2$RPKM_CORR_LOG2), xlab="log2(RPKM_CORR)", ylab="Frequency", main="Copy Number-corrected RPKM")
   dev.off()
   
   pdf(paste0(wd.rt.data, "density_rpkm.corr.gc.log2_CORR_GC_", sample, ".pdf"))
   plot(density(rpkm.corr.gc.ratio.log2$RPKM_CORR_GC_LOG2), xlab="log2(RPKM_CORR_GC)", ylab="Frequency", main="Copy Number-, GC-corrected RPKM")
   dev.off()
#}

# -----------------------------------------------------------------------------
# Multithreads
# Last Modified: 10/05/17
# -----------------------------------------------------------------------------
samples <- readTable(paste0(wd.ngs, "sclc_wgs_n101+7.list"), header=F, rownames=F, sep="")
for (s in 1:length(samples))
   print(paste0("sbatch --cpus-per-task=2 --mem=5000mb --time=72:00:00 --account=UniKoeln R CMD BATCH --no-save --no-restore '--args ", s, "' sclc-cmd-rpkm-corr.R logs/", samples[s], ".out"))






# -----------------------------------------------------------------------------
# Test
# Last Modified: 10/05/17
# -----------------------------------------------------------------------------
samples <- readTable(paste0(wd.ngs, "sclc_wgs_n101+7.list"), header=F, rownames=F, sep="")

colnames <- c("SAMPLE", "OLD", "NEW")
segments <- toTable(0, length(colnames), length(samples), colnames)
segments$SAMPLE <- samples
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   ## Old and new SEGs
   segs.old <- read.peiflyne.cn.seg(paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_cn.seg"))
   segs.new <- read.peiflyne.cn.seg(paste0(wd.ngs, "03_coverage/", sample, "/", sample, "_ANALYSIS/", sample, "_cn.seg"))
   
   segments$OLD[s] <- nrow(segs.old)
   segments$NEW[s] <- nrow(segs.new)   
}

pdf(paste0(wd.rt.data, "plot_cn.seg_OLD-NEW.pdf"))
plot(segments$OLD~segments$NEW)
dev.off()




## Multithreads
#   runs <- round(nrow(segs)/100) - 1
#   start <- 0
#   end <- 0
#   for (r in 1:runs) {
#      start <- 100 * (r-1) + 1
#      if (r != runs)
#         end <- 100 * r
#      else
#         end <- nrow(segs)
#      print(paste0("sbatch --cpus-per-task=2 --mem=2000mb --time=02:00:00 --account=UniKoeln R CMD BATCH --no-save --no-restore '--args ", s, " ", start, " ", end, "' sclc-rpkm-corr.R logs/", sample, "_", start, "-", end, ".out"))
#   }
#}






rt.rr.cln <- subset(rt.rr, BAM_TUMOR != 0)
rt.rr.cln <- subset(rt.rr.cln, BAM_NORMAL != 0)

##
pdf(paste0(wd, "hist_rt_ratio.pdf"))
hist(rr$RPM_RATIO, xlab="Ratio", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

pdf(paste0(wd, "hist_rt_ratio_log2.pdf"))
hist(log2(rr$RPM_RATIO), xlab="log2(Ratio)", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

pdf(paste0(wd, "density_rt_ratio_log2.pdf"))
plot(density(log2(rt.rr$RPM_RATIO)), xlab="log2(Ratio)", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

##
pdf(paste0(wd, "hist_rt.rr_RPM.log2.pdf"))
hist(rt.rr$RPM_LOG2, xlab="RPM (log2 FC)", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

pdf(paste0(wd, "hist_rt.rr.cln_RPM.log2.pdf"))
hist(rt.rr.cln$RPM_LOG2, xlab="RPM (log2 FC)", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

##
pdf(paste0(wd, "density_rt.rr_RPM.log2.pdf"))
plot(density(rt.rr$RPM_LOG2), xlab="RPM (log2 FC)", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

pdf(paste0(wd, "density_rt.rr.cln_RPM.log2.pdf2"))
plot(density(rt.rr$RPM_LOG2), xlab="RPM (log2 FC)", ylab="Frequency", main="RPM between Tumor and Normal")
dev.off()

# -----------------------------------------------------------------------------
# Density of breakpoint length
# Last Modified: 25/04/17
# ----------------------------------------------------------------------------- 
p1 <- sv.bp.tpm$Length_Pair_1
p2 <- sv.bp.tpm$Length_Pair_2

xr <- range(density(p1)$x)
yr <- range(density(p1)$y)

pdf(paste0(wd, "density_sv.bp.tpm_length-pair.pdf"))
plot(density(p1), xlim=xr, ylim=yr, main="Length (Max - Min_Pos_Pair) in Intra + Inter")
lines(density(p2), xlim=xr, ylim=yr, col="darkgray")
legend(x=max(xr)-136, y=max(yr), legend=c("Pair 1", "Pair 2"), col=c("black", "darkgray"), lty=1)
dev.off()



de.all.chr1[,-1]








##
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/rearrangements/"
wd <- "/re/home/tyang2/SCLC/analysis/rearrangements/"
wd.ngs <- "/ngs/cangen/tyang2/ngs/SCLC/WGS/"
setwd(wd)
chrs.length <- readTable("/re/home/tyang2/local/R/guide-to-the/chromosome_length.txt", header=F, rownames=T, sep="\t")
colnames(chrs.length) <- c("Chromosome", "Length")

samples <- readTable("sclc_wgs_sv_n106.list", header=F, rownames=F, sep="\t")
chrs <- paste0("chr", c(1, 2, 8, 3, 10, 11, 19, 13, 12, 4, 6, 5, 17, 9, 18, 20, "X", "Y", 7, 14, 15, 16, 21)) 

# -----------------------------------------------------------------------------
# Samples
# -----------------------------------------------------------------------------
svReport <- function(sv) {
   r1 <- nrow(subset(sv, Size < 10000))
   r2 <- nrow(subset(subset(sv, Size >= 10000), Size < 20000))
   r3 <- nrow(subset(sv, Size >= 20000))
 
   return(c(r1, r2, r3))
}

###
##
reports <- 0
for (s in 1:length(samples)) {
   sample <- samples[s]   
   sv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt")
   sv <- readTable(sv.file, header=T, rownames=F, sep="\t")
 
   #sv.inv <- subset(sv, Inversion != "n")
   #sv.inv.inter <- subset(sv.inv, Type1 == "Inter")
 
   #sv.inv.intra <- subset(sv.inv, Type1 == "Intra")
   #sv.inv.intra$Size <- abs(sv.inv.intra$Exact_Pos_Pair_1 - sv.inv.intra$Exact_Pos_Pair_2)
   #sv.inv.intra.qc <- subset(sv.inv.intra, Total_Coverage >= 4)
   #sv.inv.intra.qc <- subset(sv.inv.intra.qc, Not_In_Normals == 1)   
 
   ##
   #writeTable(sv.inv, paste(path, "data/sv/inv/", sample, "_unsnarl_inv.txt", sep=""), colnames=T, rownames=F, sep="\t")
 
   #sv.inv.10k <- rbind(sv.inv.inter, subset(sv.inv.intra, Size >= 10000)[,-ncol(sv.inv.intra)])
   #writeTable(sv.inv.10k, paste(path, "data/sv/inv/", sample, "_unsnarl_inv_10k.txt", sep=""), colnames=T, rownames=T, sep="\t")
 
   ## Report
   #report <- svReport(inv.intra.qc)
   #report <- c(sample, report)
   report <- c(sample, nrow(sv))
 
   reports <- rbind(reports, report)
}
writeTable(reports, "sclc_wgs_sv_n106.txt", colnames=F, rownames=F, sep="\t")




# -----------------------------------------------------------------------------
# Method: Find clustered SV/INV (PNG/PDF)   750, 405
# Usage: plotSV()
# Last Modified: 09/11/17
# -----------------------------------------------------------------------------
plotINV <- function(s, sample, chr, sv, cna, device, min, max) {
   sv.intra <- subset(sv, Type1 == "Intra")
   sv.inter <- subset(sv, Type1 == "Inter")
   sv.inter.out <- subset(sv.inter, Chr_Pair1 == chr)
   sv.inter.in  <- subset(sv.inter, Chr_Pair2 == chr)
   if (min == 0 && max == 0) {
      min <- min(cna$Start)
      max <- max(cna$End)
   }

   file <- paste0("plots/cna_inv_by-chrs/", chr, "/", s, "_cna_", chr, "_", sample, "_", min, "-", max)
   if (device == "png") {
      png(paste(file, device, sep="."), height=405, width=750)
   } else if (device == "pdf") {
      pdf(paste(file, device, sep="."), height=4, width=10)
   }
   #par(mfrow=c(5,1), mar=c(3.5,4,0,0), cex=0.2)
   attach(mtcars)
   par(mfrow=c(5,1), mar=c(3.5,4,0,0), cex=0.2)
   
   ###
   ## Intra:NO
   plot(min, min, xlim=c(min, max), ylim=c(-1, 0), xaxt='n', yaxt="n", col="white")
   plotSV(subset(sv.intra, Inversion == "n"), min, max, "darkgray")
   
   ## Intra:H2H
   plot(min, min, xlim=c(min, max), ylim=c(-1, 0), xaxt='n', yaxt="n", col="white")
   plotSV(subset(sv.intra, Inversion == "h"), min, max, "orange")

   ## Intra:T2T
   plot(min, min, xlim=c(min, max), ylim=c(-1, 0), xaxt='n', yaxt="n", col="white")
   plotSV(subset(sv.intra, Inversion == "t"), min, max, "blue")

   ## CNA
   col1 <- ifelse (cna$iCN >= 2, "red", "blue")
   col2 <- ifelse (cna$iCN == 2, "black", col1)
   plot(cna$Start, cna$iCN, xlim=c(min, max), xaxt='n', ann=F, cex.axis=3, col="white")
   for (c in 1:nrow(cna)) {
      segments(cna$Start, cna$iCN, cna$End, cna$iCN, col=col2)
   }
   #abline(h=10, lty=2, lwd=1, col="dimgrey")
   
   #for (c in 1:nrow(cytoband)) {
   #   abline(v=cytoband$End, lty=5, lwd=0.5, col="lightgrey")
   #}
   
   ## Inter
   plot(min, 0, xlim=c(min, max), ylim=c(-1, 0), xaxt='n', yaxt="n", col="white")
   for (c in 1:nrow(sv.inter.out)) {
      abline(v=sv.inter.out$Exact_Pos_Pair_1[c], lty=1, lwd=0.5, col="green")
   }
   for (c in 1:nrow(sv.inter.in)) {
      abline(v=sv.inter.in$Exact_Pos_Pair_2[c], lty=1, lwd=0.5, col="darkgreen")
   }
   axis(side=1, at=seq(min, max, by=2500000), las=1, cex.axis=3, tck=-0.04, mgp=c(0,2.5,0))
  
   ## Log R
   #plot(snv$Position, snv$Log.R, xlim=c(min, max), xaxt='n', xaxt='n', ann=F, cex.axis=3)
   
   ## BAF
   #plot(snv$Position, snv$BAF, xlim=c(min, max), xaxt='n', ann=F, cex.axis=3)
   #axis(side=1, at=seq(min, max, by=2500000), las=1, cex.axis=3, tck=-0.04, mgp=c(0,2.5,0))

   dev.off()
}

plotSV <- function(sv, min, max, color) {
   if (nrow(sv) != 0) {
      for (c in 1:nrow(sv)) {
         breakpointMean <- (sv$Exact_Pos_Pair_1[c] + sv$Exact_Pos_Pair_2[c]) / 2
         
         arc <- function(x) -(x - breakpointMean)^2 / (sv$Exact_Pos_Pair_2[c] - breakpointMean)^2
         curve(arc, sv$Exact_Pos_Pair_1[c], sv$Exact_Pos_Pair_2[c], xlim=c(min, max), col=color, add=T)
      }
   }
}

binWindows <- function(sample, breaks.chr, window, chr.length) {  ##, report.invs.sample, report.invs.max) {
   report.bin.chr <- data.frame(matrix(NA, 1, 0))
   rownames(report.bin.chr) <- sample
   	  
   start <- 0
   while (start < chr.length$Length) {
      end <- start + window
      report.bin.bit <- data.frame(matrix(NA, 1, 1))
   	  rownames(report.bin.bit) <- sample
      colnames(report.bin.bit) <- paste(gsub("chr", "", chr), ":", paste(start/1000000, end/1000000, sep="-"), "Mb", sep="")
         
      breaks.bin <- breaks.chr[apply(breaks.chr, 1, function(x) as.numeric(x[2]) >= start && as.numeric(x[2]) <= end) == T,]
      report.bin.bit[1,1] <- nrow(breaks.bin) ##/ report.invs.sample / window * 1000000 * report.invs.max   ## BAD IDEA 01/04/17: TPM ?
   	  
   	  report.bin.chr <- cbind(report.bin.chr, report.bin.bit)
   	  start <- start + window/2
   }
   return(report.bin.chr)
}

###
## TO-DO: Now its only "nrow(sv)" entries, not breakpoints!!
report <- c(1:length(samples))
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   ## INPUT: INV
   #sv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt")   ## All rearrangements
   sv.file <- paste0(wd.ngs, "02_rearrangements/", sample, "_unsnarl.txt")   ## All rearrangements
   sv <- readTable(sv.file, header=T, rownames=F, sep="\t")
   ##
   #sv <- subset(sv, Type1 != "Inter")   ## BUG FIX: 01/04/17
   #sv$Size <- abs(sv$Exact_Pos_Pair_1 - sv$Exact_Pos_Pair_2)
   #sv <- subset(sv, Size >= 10000)
   
   report[s] <- nrow(sv)
}

#window <- 10000000
window <- 5000000

report.bin <- NA
#report.intra <- c(1:length(samples))   ## TO-DO
for (s in 1:length(samples)) {
   sample <- samples[s]
   
   ## INPUT: INV
   #sv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt")   ## All rearrangements
   sv.file <- paste0(wd.ngs, "02_rearrangements/", sample, "_unsnarl.txt")   ## All rearrangements
   sv <- readTable(sv.file, header=T, rownames=F, sep="\t")
   ##
   sv.inter <- subset(sv, Type1 == "Inter")   ## ADD: 20/04/17
   sv.intra <- subset(sv, Type1 == "Intra")   ## ADD: 20/04/17
   sv.intra$Size <- abs(sv.intra$Exact_Pos_Pair_1 - sv.intra$Exact_Pos_Pair_2)
   sv.intra <- subset(sv.intra, Size >= 10000)
   sv.intra <- sv.intra[,-27]   ## TO-DO
   sv <- rbind(sv.intra, sv.inter)
   
   break1 <- sv[,c("Chr_Pair1", "Exact_Pos_Pair_1")]
   colnames(break1) <- c("Break_Chr", "Break_Exact_Pos")
   break2 <- sv[,c("Chr_Pair2", "Exact_Pos_Pair_2")]
   colnames(break2) <- c("Break_Chr", "Break_Exact_Pos")
   breaks <- rbind(break1, break2)
   
   ## INPUT: CNA
   cna.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_iCN.seg")
   cna <- readTable(cna.file, header=F, rownames=F, sep="\t")
   colnames(cna) <- c("Sample", "Chr", "Start", "End", "V5", "iCN")
   
   ###
   ## Bin the genome into windows
   report.bin.sample <- data.frame(matrix(NA, 1, 0))
   rownames(report.bin.sample) <- sample
   
   for (c in 1:length(chrs)) {
   #for (c in 1:1) {
      chr <- chrs[c]
         
   	  sv.chr <- subset(sv, Chr_Pair1 == chr)
   	  sv.tmp <- subset(sv, Chr_Pair1 != chr)
      sv.chr <- rbind(sv.chr, subset(sv.tmp, Chr_Pair2 == chr))
      
      breaks.chr <- subset(breaks, Break_Chr == chr) 
      cna.chr <- subset(cna, Chr == chr)
      chr.length <- subset(chrs.length, Chromosome == chr)

      ## OUTPUT: BFB regional plot
      #if (nrow(sv.chr) != 0)
      #   plotINV(s, sample, chr, sv.chr, cna.chr, "png", 0, 0)

   	  ## Bin the genome into windows
   	  report.bin.sample <- cbind(report.bin.sample, binWindows(sample, breaks.chr, window, chr.length))   ##, report.invs[s], max(report.invs)))
   }
      
   if (is.na(report.bin)) {
      report.bin <- report.bin.sample
   } else {
      report.bin <- rbind(report.bin, report.bin.sample)
   }
}

###
## OUTPUT: BFB report
#report.5mb <- report.bin
#report.10mb <- report.bin

#save(report.10mb, file=paste(path, "analysis/report_sv_10mb.RData", sep=""))
writeTable(report.bin, "report_sv_5mb_10k.txt", colnames=T, rownames=T, sep="\t")

##
#report.bin <- report.5mb
#report.bin <- report.bin[samples,]

##
#window <- 5000000
cutoff <- 5
report.bin.p75 <- NA
for (s in 1:nrow(report.bin)) {
   report.bin.sample <- report.bin[s,]
   #report.bin.sample <- report.bin.sample / window * 1000000

   percentile75 <- as.numeric(quantile(report.bin.sample[report.bin.sample > 0])[4])
   idx1 <- which(report.bin.sample >= percentile75)
   idx2 <- which(report.bin.sample >= cutoff)
   #idx <- idx1   
   idx <- intersect(idx1, idx2)
   #report.bin.sample[,idx]
   if (length(idx) == 0)
      report.bin.sample[1,] <- 0
   else {
      report.bin.sample[1, -idx] <- 0
   }
   
   if (is.na(report.bin.p75)) {
      report.bin.p75 <- report.bin.sample
   } else {
      report.bin.p75 <- rbind(report.bin.p75, report.bin.sample)
   }
}

##
rm <- c()
for (b in 1:ncol(report.bin.p75)) {
   if (sum(report.bin.p75[,b]) == 0) {
   	  rm <- c(rm, b)
   }
}
report.bin.p75 <- report.bin.p75[,-rm]

writeTable(report.bin.p75, "report_sv_10k_5mb_p75_c5.txt", colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# RNA-seq samples
# -----------------------------------------------------------------------------
rna <- readTable("/re/home/tyang2/ngs/SCLC/RNA/sclc_rna_n81-8.list", header=F, rownames=F, sep="\t")
rna.list <- rna[seq(1, length(rna), 2)]
rna.list <- gsub("_1.fq.gz", "", rna.list)
rna.list <- cbind(rna.list, rna.list)
rna.list[,1] <- gsub("T", "", rna.list[,1])
rna.list[,1] <- gsub("AB", "", rna.list[,1])
writeTable(rna.list, "/re/home/tyang2/ngs/SCLC/RNA/sclc_rna_n81-8-1.list", colnames=F, rownames=F, sep="\t")







report.bin.chr <- data.frame(matrix(NA, 24, 2))
report.bin.chr[,1] <- chrs
bins <- paste("chr", colnames(report.bin.p75), sep="")
for (c in 1:length(chrs)) {
   chr <- report.bin.p75[, grep(paste(chrs[c], ":", sep=""), bins)]

   sum <- 0
   if (is.null(ncol(chr))) {
      sum <- sum(chr)
   } else {
   	  if (ncol(chr) != 0) {
         for (b in 1:ncol(chr)) {
   	        sum <- sum + sum(chr[,b])
         }
      }
   }
   report.bin.chr[c, 2] <- sum
}
writeTable(report.bin.chr, paste0(wd, "report_sv_inv_5mb_p75_chrs.txt"), colnames=F, rownames=F, sep="\t")

##
rm <- c()
for (b in 1:ncol(report.bin.p75)) {
   if (sum(report.bin.p75[,b]) == 0) {
   	  rm <- c(rm, b)
   }

   report.bin.p75[which(report.bin.p75[,b] == 0), b] <- ""
}
report.bin.p75 <- report.bin.p75[,-rm]
report.5mb.p75.n10 <- report.bin.p75
#report.5mb.p75.c1.5 <- report.bin.p75
#report.5mb.p75.c2 <- report.bin.p75
#report.5mb.p75.c3 <- report.bin.p75
#report.5mb.p75.c3.5 <- report.bin.p75
#report.5mb.p75.c4 <- report.bin.p75
#report.5mb.p75.c5 <- report.bin.p75

#report.10mb.p75.c1.5 <- report.bin.p75
#report.10mb.p75.c2 <- report.bin.p75

##
## 5 MB (Upper quantile, cutoff 2)
> dim(report.5mb)            ## report.5mb
[1]  106 1253
> dim(report.5mb.p75.c1.5)   ## report.5mb.p75.c1.5
[1] 106 778
> dim(report.5mb.p75.c2)     ## report.5mb.p75.c2
[1] 106 609
> dim(report.5mb.p75.c3)     ## report.5mb.p75.c3
[1] 106 304
> dim(report.5mb.p75.c3.5)     ## report.5mb.p75.c3.5
[1] 106 204
> dim(report.5mb.p75.c4)     ## report.5mb.p75.c4
[1] 106 161
> dim(report.5mb.p75.c5)     ## report.5mb.p75.c5
[1] 106 103

> dim(report.5mb.p75.n10)
[1] 106 140

save(report.5mb, report.5mb.p75.c1.5, report.5mb.p75.c2, report.5mb.p75.c3, report.5mb.p75.c3.5, report.5mb.p75.c4, report.5mb.p75.c5, file=paste(path, "analysis/report_sv_5mb.RData", sep=""))
writeTable(report.5mb.p75.n10, paste(path, "analysis/report_sv_5mb_p75_n10_20161122.txt", sep=""), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Method: Size of INVs
# Usage: 
# Last Modified: 24/11/16
# -----------------------------------------------------------------------------
path <- "/Users/tpyang/Work/uni-koeln/tyang2/Link2Analysis/SCLC/"
setwd(path)
source("../the_daily_package.R")
samples <- readTable("analysis/sv_inv_bfb_sclc_c3_H19-L87.list", header=T, rownames=T, sep="\t")

report.size <- NA
for (s in 1:nrow(samples)) {
   sample <- samples[s, 1]
   
   ## INPUT: INV
   sv.file <- paste(path, "ngs/WGS/", sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt", sep="")   ## All rearrangements
   sv <- readTable(sv.file, header=T, rownames=F, sep="\t")
   sv.intra <- subset(sv, Type1 == "Intra")
   sv.intra$Size <- abs(sv.intra$Exact_Pos_Pair_1 - sv.intra$Exact_Pos_Pair_2)
   sv.intra$iCN_Max <- NA
   
   ## INPUT: CNA
   cna.file <- paste(path, "ngs/WGS/", sample, "/", sample, "_ANALYSIS/", sample, "_iCN.seg", sep="")
   cna <- readTable(cna.file, header=F, rownames=F, sep="\t")
   colnames(cna) <- c("Sample", "Chr", "Start", "End", "V5", "iCN")
   rownames(cna) <- paste(cna$Chr, cna$Start, cna$End, sep="-")
   
   for (r in 1:nrow(sv.intra)) {
   	  inv <- sv.intra[r,]
   	  
   	  ## Breakpoint1
   	  cna.inv <- subset(cna, Chr == inv$Chr_Pair1)
   	  cna.inv.b1.s <- subset(cna.inv, Start >= inv$Exact_Pos_Pair_1)
   	  cna.inv.b1.e <- subset(cna.inv, End >= inv$Exact_Pos_Pair_1)
   	  overlaps1 <- intersect(rownames(cna.inv.b1.s), rownames(cna.inv.b1.e))
   	  
   	  cna.inv.b1 <- rbind(cna.inv.b1.s[overlaps1,], cna.inv.b1.s[setdiff(rownames(cna.inv.b1.s), overlaps1),])
   	  cna.inv.b1 <- rbind(cna.inv.b1, cna.inv.b1.e[setdiff(rownames(cna.inv.b1.e), overlaps1),])

   	  ## Breakpoint2
   	  cna.inv.b2.s <- subset(cna.inv, Start >= inv$Exact_Pos_Pair_2)
   	  cna.inv.b2.e <- subset(cna.inv, End >= inv$Exact_Pos_Pair_2)
   	  overlaps2 <- intersect(rownames(cna.inv.b2.s), rownames(cna.inv.b2.e))
   	  
   	  cna.inv.b2 <- rbind(cna.inv.b2.s[overlaps2,], cna.inv.b2.s[setdiff(rownames(cna.inv.b2.s), overlaps2),])
   	  cna.inv.b2 <- rbind(cna.inv.b2, cna.inv.b2.e[setdiff(rownames(cna.inv.b2.e), overlaps2),])
   	  
   	  ## Joint data set
   	  overlaps12 <- intersect(rownames(cna.inv.b1), rownames(cna.inv.b2))
   	  
   	  cna.inv.b12 <- rbind(cna.inv.b1[overlaps12,], cna.inv.b1[setdiff(rownames(cna.inv.b1), overlaps12),])
   	  cna.inv.b12 <- rbind(cna.inv.b12, cna.inv.b2[setdiff(rownames(cna.inv.b2), overlaps12),])
   	  
   	  ##
   	  sv.intra[r,]$iCN_Max <- max(na.omit(cna.inv.b12$iCN))
   }
   
   if (is.na(report.size)) {
      report.size <- sv.intra[, c("Sample", "Chr_Pair1", "Exact_Pos_Pair_1", "Exact_Pos_Pair_2", "Inversion", "Size", "iCN_Max")]
   } else {
      report.size <- rbind(report.size, sv.intra[, c("Sample", "Chr_Pair1", "Exact_Pos_Pair_1", "Exact_Pos_Pair_2", "Inversion", "Size", "iCN_Max")])
   }
}
#writeTable(report.size, paste(path, "analysis/report_inv_Size_iCN_Max_20161124.txt", sep=""), colnames=T, rownames=T, sep="\t")

##
#pdf("report_inv_size-iCN.pdf")
#plot(report.size$Size, report.size$iCN_Max)
#dev.off()

## 100kb
report.size.table <- toTable(NA, 9, 10, c("Size", "Freq", "iCN<2", "iCN=2", "2<iCN<=5", "5<iCN<=10", "10<iCN<=25", "25<iCN", "iCN_Max"))
start <- 0
bp <- 100000
for (p in 1:10) {
   size <- p * bp
   t <- subset(subset(report.size, Size <= size), Size >= start)
   
   report.size.table[p, 1] <- size
   report.size.table[p, 2] <- nrow(t)
   if (nrow(t) != 0) {
      report.size.table[p, 3] <- nrow(subset(t, iCN_Max < 2))
      report.size.table[p, 4] <- nrow(subset(t, iCN_Max == 2))
      report.size.table[p, 5] <- nrow(subset(subset(t, iCN_Max > 2), iCN_Max <= 5))
      report.size.table[p, 6] <- nrow(subset(subset(t, iCN_Max > 5), iCN_Max <= 10))
      report.size.table[p, 7] <- nrow(subset(subset(t, iCN_Max > 10), iCN_Max <= 25))
      report.size.table[p, 8] <- nrow(subset(t, iCN_Max > 25))
      report.size.table[p, 9] <- max(na.omit(t$iCN_Max))
   } else {
   	  report.size.table[p, 3:8] <- ""
   }
   rownames(report.size.table)[p] <- paste(start/1000, "kb-", size/1000, "kb", sep="")
      
   start <- size
}
writeTable(report.size.table, paste(path, "analysis/report_inv_Size_100kb_20161124.txt", sep=""), colnames=T, rownames=T, sep="\t")

## 1Mb
report.size.table <- toTable(NA, 9, 10, c("Size", "Freq", "iCN<2", "iCN=2", "2<iCN<=5", "5<iCN<=10", "10<iCN<=25", "25<iCN", "iCN_Max"))
start <- 0
bp <- 1000000
for (p in 1:226) {
   size <- p * bp
   t <- subset(subset(report.size, Size <= size), Size >= start)
   
   report.size.table[p, 1] <- size
   report.size.table[p, 2] <- nrow(t)
   if (nrow(t) != 0) {
      report.size.table[p, 3] <- nrow(subset(t, iCN_Max < 2))
      report.size.table[p, 4] <- nrow(subset(t, iCN_Max == 2))
      report.size.table[p, 5] <- nrow(subset(subset(t, iCN_Max > 2), iCN_Max <= 5))
      report.size.table[p, 6] <- nrow(subset(subset(t, iCN_Max > 5), iCN_Max <= 10))
      report.size.table[p, 7] <- nrow(subset(subset(t, iCN_Max > 10), iCN_Max <= 25))
      report.size.table[p, 8] <- nrow(subset(t, iCN_Max > 25))
      report.size.table[p, 9] <- max(na.omit(t$iCN_Max)) 
   } else {
   	  report.size.table[p, 3:8] <- ""
   }
   rownames(report.size.table)[p] <- paste(start/1000000, "Mb-", size/1000000, "Mb", sep="")
      
   start <- size
}
#writeTable(report.size.table, paste(path, "analysis/report_inv_Size_1Mb_20161124.txt", sep=""), colnames=T, rownames=T, sep="\t")

#pdf("report_inv_size_1-50.pdf")
#plot(report.size.table[1:50,])
#dev.off()

report.size.table <- data.frame(table(report.size$Size))
report.size.table <- report.size.table[order(report.size.table$Freq, decreasing=T),]
report.size.table$Var1 <- as.vector(report.size.table$Var1)
report.size.table$Var1 <- as.numeric(report.size.table$Var1)

#pdf("report_inv_size.pdf")
#plot(report.size.table$Var1, report.size.table$Freq)
#dev.off()

# -----------------------------------------------------------------------------
# Method: Two different y axes on the same plot
# Usage: 
# Last Modified: 29/11/16
# -----------------------------------------------------------------------------
report.size.table.tmp <- report.size.table

report.size.table <- report.size.table[1:150,]
report.size.table <- report.size.table.tmp

## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)

## Plot first set of data and draw its axis
plot(report.size.table$Size/1000000, report.size.table$Freq, pch=16, axes=FALSE, xlab="", ylab="", type="b",col="black", main="Distribution of Rearrangement")
axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
mtext("Frequency",side=2,line=2.5)
box()

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(report.size.table$Size/1000000, report.size.table$iCN_Max, pch=15, xlab="", ylab="", axes=FALSE, type="b", col="red")
## a little farther out (line=4) to make room for labels
mtext("Max iCN",side=4,col="red",line=4) 
axis(4, ylim=c(0,max(report.size.table$iCN_Max)), col="red",col.axis="red",las=1)

## Draw the time axis
axis(1,pretty(range(report.size.table$Size/1000000),10))
mtext("Size (Mb)", side=1, col="black", line=2.5)  

## Add Legend
legend("topright",legend=c("Frequency","Max iCN"), text.col=c("black","red"), pch=c(16,15), col=c("black","red"))








##
## 10 MB (Upper quantile, cutoff 2)
> dim(report.10mb)            ## report.5mb
[1]  106 633
> dim(report.10mb.p75.c1.5)   ## report.5mb.p75.c1.5
[1] 106 521
> dim(report.10mb.p75.c2)     ## report.5mb.p75.c2
[1] 106 452

save(report.10mb, report.10mb.p75.c1.5, report.10mb.p75.c2, file=paste(path, "analysis/report_sv_10mb.RData", sep=""))

##
## 5 MB (ESAC)
> dim(report.bin)            ## report.5mb
[1]  119 1253
> dim(report.bin.p75.c1.5)   ## report.5mb.p75.c1.5
[1] 119 235
> dim(report.bin.p75.c2)        ## report.5mb.p75.c2
[1] 119 177

writeTable(report.bin.p75, "/Users/tpy20/Work/mrc-cu/tpy20/ICGC/sv/production/SV/INV/BFB/plots/5mb_q75/report_inv_4r0c_10kb_5mb_q75_c1.5.txt", colnames=T, rownames=T, sep="\t")
writeTable(report.5mb.p75.c3.5, "/Users/tpy20/Work/mrc-cu/tpy20/ICGC/sv/production/SV/INV/BFB/plots/5mb_q75/report_inv_4r0c_10kb_5mb_q75_c3.5.txt", colnames=T, rownames=T, sep="\t")
writeTable(report.bin.p75, "/Users/tpy20/Work/mrc-cu/tpy20/ICGC/sv/production/SV/INV/BFB/plots/5mb_q75/report_inv_bin_4r0c_10kb_5mb_q75_c2.txt", colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Method: Fisher's exact
# Usage: 
# Last Modified: 10/11/16
# -----------------------------------------------------------------------------
###
## Main
#samples.old <- readTable("analysis/sv_inv_bfb_sclc_c3_H54-L52.list", header=T, rownames=T, sep="\t")
#samples <- readTable("analysis/sv_inv_bfb_sclc_c3_H55-L51.list", header=T, rownames=T, sep="\t")
#samples <- readTable("analysis/sv_inv_bfb_sclc_c3_H47-L59.list", header=T, rownames=T, sep="\t")
samples <- readTable("analysis/sv_inv_bfb_sclc_c3_H19-L87.list", header=T, rownames=T, sep="\t")

pheno1 <- readTable("metadata/nature14664-s1_ST1.txt", header=T, rownames=T, sep="\t")
pheno2 <- readTable("metadata/nature14664-s1_ST2.txt", header=T, rownames=T, sep="\t")
phenos <- cbind(samples, pheno1[samples$Sample,], pheno2[samples$Sample,])

## SCLC (n=106; 54 High vs 52 Low)
## Significant mutated genes
muts <- readTable("metadata/nature14664-s1_ST3.txt", header=T, rownames=F, sep="\t")
muts <- subset(muts, PAT_ID != "S00841")
muts <- subset(muts, PAT_ID != "S00935")
muts <- subset(muts, PAT_ID != "S02297")
muts <- subset(muts, PAT_ID != "S02353")
muts.non <- subset(muts, Type_1 != "silent")

high <- phenos[phenos$Group == "H",]$Sample
low  <- phenos[phenos$Group == "L",]$Sample
muts.non.high <- NA
for (s in 1:length(high)) {
   sample <- high[s]
   muts.non.high <- rbind(muts.non.high, subset(muts.non, PAT_ID == sample))
}
muts.non.high <- muts.non.high[-1,]

muts.non.low <- NA
for (s in 1:length(low)) {
   sample <- low[s]
   muts.non.low <- rbind(muts.non.low, subset(muts.non, PAT_ID == sample))
}
muts.non.low <- muts.non.low[-1,]

genes <- c("TP53", "RB1", "KIAA1211", "COL22A1", "RGS7", "FPR1", "EP300", "CREBBP", "ASPM", "ALMS1", "PDE4DIP", "XRN1", "PTGFRN", "TP73", "RBL1", "RBL2", "FMN2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "PTEN", "KIT", "PIK3CA", "BRAF")
num <- c()
mut <- toTable(0, 6, length(genes), c("High", "Low", "High_cDNA", "Low_cDNA", "High_Strand", "Low_Strand"))
rownames(mut) <- genes
for (g in 1:length(genes)) {
   gene <- genes[g]
   num <- c(num, length(unique(subset(muts.non, Gene_Hugo == gene)$PAT_ID)))
   
   mut[g, 1] <- length(unique(subset(muts.non.high, Gene_Hugo == gene)$PAT_ID))
   mut[g, 2] <- length(unique(subset(muts.non.low, Gene_Hugo == gene)$PAT_ID))
   
   mut[g, 3] <- subset(muts.non.high, Gene_Hugo == gene)$Change_cDNA[1]
}

###
##
genes <- c("NRG1", "XRN2", "ATRX")
mut <- toTable(0, 1, length(genes), c("Samples"))
rownames(mut) <- genes
for (g in 1:length(genes)) {
   gene <- genes[g]
   num <- unique(subset(muts.non, Gene_Hugo == gene)$PAT_ID)
   
   mut[g, 1] <- length(unique(subset(muts.non.high, Gene_Hugo == gene)$PAT_ID))
   mut[g, 2] <- length(unique(subset(muts.non.low, Gene_Hugo == gene)$PAT_ID))
   
   mut[g, 3] <- subset(muts.non.high, Gene_Hugo == gene)$Change_cDNA[1]
}

###
## Fisher's test for each mutations
results.gene <- toTable(0, 2, length(genes), c("P", "BH"))
rownames(results.gene) <- genes
for (g in 1:length(genes)) {
   test <- toTable(0, 2, 2, c("WT", "MUT"))
   
   test[1, 2] <- mut[g, 1]
   test[1, 1] <- length(high) - mut[g, 1] 
   test[2, 2] <- mut[g, 2]
   test[2, 1] <- length(low) - mut[g, 2] 

   results.gene[g, 1] <- fisher.test(test)[[1]]
}
results.gene <- results.gene[order(results.gene$P),]
results.gene$BH <- p.adjust(results.gene$P, "BH")

results.gene <- cbind(mut[rownames(results.gene), 1:2], results.gene)
writeTable(results.gene, "report_sv_5mb_q75_c3_muts_gene_19H-87L.txt", colnames=T, rownames=T, sep="\t")

> results.gene
         High Low           P        BH
XRN1        8   1 0.009891179 0.2472795
KIAA1211   12   7 0.079399569 0.7938822
ASPM        8   4 0.127021152 0.7938822
ALMS1       8   4 0.127021152 0.7938822
NOTCH1     10   6 0.171160151 0.8558008
RBL2        1   5 0.223692203 0.9320508
NOTCH2      1   4 0.379459735 1.0000000
BRAF        1   0 0.443396226 1.0000000
NOTCH4      2   1 0.583158822 1.0000000
RB1        37  49 0.622884108 1.0000000
RBL1        1   3 0.627675841 1.0000000
PTGFRN      3   2 0.653227499 1.0000000
PTEN        2   4 0.690994801 1.0000000
CREBBP      5   5 0.747698848 1.0000000
NOTCH3      5   5 0.747698848 1.0000000
TP53       45  56 1.000000000 1.0000000
COL22A1     8  11 1.000000000 1.0000000
RGS7        5   6 1.000000000 1.0000000
FPR1        3   4 1.000000000 1.0000000
EP300       6   7 1.000000000 1.0000000
PDE4DIP     5   6 1.000000000 1.0000000
TP73        2   3 1.000000000 1.0000000
FMN2        9  11 1.000000000 1.0000000
KIT         3   3 1.000000000 1.0000000
PIK3CA      1   2 1.000000000 1.0000000















## Remove two Chromothripsis samples
samples <- setdiff(samples, c("S02297", "S02353"))

# -----------------------------------------------------------------------------
# Remove short reads in INV (?)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# SV/INV rates per Mb (24/10/16)
# -----------------------------------------------------------------------------
svReportChromosome <- function(chrs, sv) {
 breaks <- c(sv$Chr_Pair1, sv$Chr_Pair2)
 table <- data.frame(table(breaks))
 rownames(table) <- table$breaks
 
 diff <- setdiff(chrs, names(table))
 empty <- data.frame(matrix(0, length(diff), 2))
 empty[,1] <- diff
 rownames(empty) <- diff
 colnames(empty) <- c("breaks", "Freq")
 
 report <- rbind(table, empty)
 return(report[chrs, 2])
}

###
## Main
setwd(path)

reports.inv.chr <- toTable(NA, 24, 0, chrs)
for (x in 1:length(samples)) {
 sample <- samples[x]
 
 sv.file <- paste(path, "ngs/WGS/", sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt", sep="")   ## All rearrangements
 #sv.file <- paste(path, "data/sv/inv/", sample, "_unsnarl_inv.txt", sep="")                             ## INVs
 #sv.file <- paste(path, "data/sv/inv/", sample, "_unsnarl_inv_10k.txt", sep="")                         ## INVs greater than 10k
 sv <- readTable(sv.file, header=T, rownames=F, sep="\t")
 
 report.inv.chr <- svReportChromosome(chrs, sv)
 reports.inv.chr <- rbind(reports.inv.chr, report.inv.chr)
}
colnames(reports.inv.chr) <- chrs
rownames(reports.inv.chr) <- samples

## Sort by samples
total <- toTable(0, 1, length(samples), c("Total"))

reports.inv.chr <- cbind(reports.inv.chr, total)
for (x in 1:length(samples)) {
 reports.inv.chr$Total[x] <- sum(reports.inv.chr[x, 1:24])
}
reports.inv.chr <- reports.inv.chr[order(reports.inv.chr$Total, decreasing=T),]
reports.inv.chr <- reports.inv.chr[,-25]

## Sort by chrs
total <- toTable(0, 2, 24, c("Total", "Rate"))

reports.inv.chr.t <- t(reports.inv.chr)
reports.inv.chr.t <- cbind(reports.inv.chr.t, total)
for (x in 1:nrow(reports.inv.chr.t)) {
 reports.inv.chr.t$Total[x] <- sum(reports.inv.chr.t[x, 1:length(samples)])
 reports.inv.chr.t$Rate[x] <- reports.inv.chr.t$Total[x] / chrs.length[chrs[x],]$Length * 1000000
}
reports.inv.chr.t <- reports.inv.chr.t[order(reports.inv.chr.t$Rate, decreasing=T),]

samples <- colnames(reports.inv.chr.t)[1:(ncol(reports.inv.chr.t)-2)]
chrs <- rownames(reports.inv.chr.t)

writeTable(reports.inv.chr.t, paste(path, "analysis/sv_sorted.txt", sep=""), colnames=T, rownames=T, sep="\t")
#writeTable(reports.inv.chr.t, paste(path, "analysis/sv_inv_sorted.txt", sep=""), colnames=T, rownames=T, sep="\t")
#writeTable(reports.inv.chr.t, paste(path, "analysis/sv_inv_sorted_10k.txt", sep=""), colnames=T, rownames=T, sep="\t")

##
t <- reports.inv.chr
for (x in 1:nrow(t)) {
 t$Total[x] <- sum(reports.inv.chr[x, ])
}

t[new,]$Total












# -----------------------------------------------------------------------------
# Method: Differential expression (Read counts)
# Usage: 
# Last Modified: 25/11/16
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Samples
# -----------------------------------------------------------------------------
path <- "/Users/tpyang/Work/uni-koeln/tyang2/Link2Analysis/SCLC/"
setwd(path)
source("../the_daily_package.R")

getEnsemblGeneID <- function(ensembl_gene_id, ensembl) {
   return(getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters='ensembl_gene_id', values=ensembl_gene_id, mart=ensembl))
}

#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

> dim(exp.counts)
[1] 60504    81
> length(unique(rownames(exp.counts)))
[1] 60504
> dim(expr)
[1] 31640    69

##
maps.ensembl <- getEnsemblGeneID(rownames(expr)[1], ensembl)
for (e in 2:nrow(expr)) {
   maps.ensembl <- rbind(maps.ensembl, getEnsemblGeneID(rownames(expr)[e], ensembl))
}

# -----------------------------------------------------------------------------
# Method: Differential expression (Read counts)
# Usage: 19H vs 87L
# Last Modified: 28/11/16
# -----------------------------------------------------------------------------
## Read counts
load("/Users/tpyang/Work/uni-koeln/tyang2/Link2Analysis/SCLC/metadata/exp.counts.rda")

colnames <- gsub("AB", "", colnames(exp.counts))
colnames <- gsub("T", "", colnames)
colnames <- gsub("o", "", colnames)
colnames(exp.counts) <- colnames

overlaps <- intersect(samples$Sample, colnames(exp.counts))
expr <- exp.counts[,overlaps]

samples.high <- intersect(overlaps, high)
samples.low  <- intersect(overlaps, low)

> length(samples.high)
[1] 9
> length(samples.low)
[1] 58

## log2(1 + FPKM)
expr.log2 <- log2(1 + expr)

###
## Differential expression (t-Test)
ttest <- function(p, expr.log2, samples.high, samples.low) {
   expr.log2.high <- expr.log2[p, samples.high]
   expr.log2.low  <- expr.log2[p, samples.low]

   return(t.test(expr.log2.high, expr.log2.low)$p.value)
}

wtest <- function(p, expr.log2, phenos.expr) {
   return(wilcox.test(as.numeric(expr.log2[p,]) ~ phenos.expr$Group, exact=FALSE)[[3]])
}

cmedian <- function(p, expr.log2, samples.group) {
   return(median(as.numeric(expr.log2[p, samples.group])))
}

###
## Main (ALL)
expr.log2.tmp <- expr.log2

results.expr <- toTable(0, 2, nrow(expr.log2.tmp), c("P", "BH"))
rownames(results.expr) <- rownames(expr.log2.tmp)

## Effect size
results.expr$High <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.high))
results.expr$Low  <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.low))
results.expr$Effect <- results.expr$High - results.expr$Low

###
## Main (Only varible)
#hist(abs(results.expr$Effect))

results.expr.effect <- rbind(subset(results.expr, Effect >= 1.5), subset(results.expr, Effect <= -1.5))
expr.log2.effect <- expr.log2.tmp[rownames(results.expr.effect),]

##
expr.log2.tmp <- expr.log2.effect

##
maps.ensembl <- getEnsemblGeneID(rownames(expr.log2.tmp)[1], ensembl)
for (e in 2:nrow(expr.log2.tmp)) {
   maps.ensembl <- rbind(maps.ensembl, getEnsemblGeneID(rownames(expr.log2.tmp)[e], ensembl))
}
rownames(maps.ensembl) <- maps.ensembl$ensembl_gene_id

##
results.expr <- toTable(0, 2, nrow(expr.log2.tmp), c("P", "BH"))
rownames(results.expr) <- rownames(expr.log2.tmp)

## Effect size
results.expr$High <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.high))
results.expr$Low  <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.low))
results.expr$Effect <- results.expr$High - results.expr$Low

## P
results.expr$P <- mapply(x = 1:nrow(expr.log2.tmp), function(x) ttest(x, expr.log2.tmp, samples.high, samples.low))
results.expr <- results.expr[order(results.expr$P),]

## BH
results.expr <- results.expr[order(results.expr$P),]
results.expr$BH <- p.adjust(results.expr$P, "BH")

results.expr.maps <- cbind(results.expr, maps.ensembl[rownames(results.expr),])

## W
phenos.expr <- phenos[overlaps,]
phenos.expr$Group <- as.factor(phenos.expr$Group)

results.expr <- results.expr[rownames(expr.log2.tmp),]
results.expr$W <- mapply(x = 1:nrow(expr.log2.tmp), function(x) wtest(x, expr.log2.tmp, phenos.expr))

## WBH
results.expr <- results.expr[order(results.expr$W),]
results.expr$WBH <- p.adjust(results.expr$W, "BH")

results.expr.maps <- cbind(results.expr, maps.ensembl[rownames(results.expr),])

writeTable(results.expr.maps, "analysis/report_expr_counts_log2_19H_ttest_wilcox_effect1.5.txt", colnames=T, rownames=T, sep="\t")

sclc <- c("ENSG00000116990", "ENSG00000134323", "ENSG00000136997", "ENSG00000270141", "ENSG00000105173", "ENSG00000114127", "ENSG00000088930", "ENSG00000085224", "ENSG00000164362")

writeTable(results.expr[sclc,], "analysis/report_expr_counts_log2_19H_ttest_wilcox_sclc.txt", colnames=T, rownames=T, sep="\t")
save(maps.ensembl, file="analysis/maps.ensembl.RData")

# -----------------------------------------------------------------------------
# Method: Amplified genes from Chambers et al
# Last Modified: 07/12/16
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/Link2Analysis/SCLC/metadata/maps.ensembl.RData")
amps.all <- readTable("analysis/G4/metadata/Chambers et al/nbt.3295-ST8-amp.txt", header=T, rownames=F, sep="")

amps <- amps.all[-grep("p", amps.all$Oncogene),]
amps <- amps[-grep("q", amps$Oncogene),]
rownames(amps) <- amps$Oncogene

##
maps.ensembl.amps <- subset(maps.ensembl, hgnc_symbol == amps$Oncogene[1])
for (x in 2:nrow(amps)) {
   maps.ensembl.amp <- subset(maps.ensembl, hgnc_symbol == amps$Oncogene[x])
   if (nrow(maps.ensembl.amp) == 0) {
   	  maps.ensembl.amp <- maps.ensembl[grep(amps$Oncogene[x], maps.ensembl$hgnc_symbol),]	
   	  maps.ensembl.amps 
   }
   maps.ensembl.amps <- rbind(maps.ensembl.amps, maps.ensembl.amp)
}
rownames(maps.ensembl.amps) <- amps$Oncogene

maps.ensembl.amps <- cbind(maps.ensembl.amps, amps)

##
more <- c("MET", "CCND2", "ERBB3", "MYCN", "FGFR1", "AHR", "AHRR", "CYP1A1", "CYP1B1")
maps.ensembl.more <- subset(maps.ensembl, hgnc_symbol == more[1])
for (x in 2:length(more)) {
   maps.ensembl.more <- rbind(maps.ensembl.more, subset(maps.ensembl, hgnc_symbol == more[x]))
}
rownames(maps.ensembl.more) <- more

amps.more <- toTable("NA", ncol(amps), length(more), colnames(amps))
rownames(amps.more) <- more

maps.ensembl.more <- cbind(maps.ensembl.more, amps.more)

##
maps.ensembl.amps <- rbind(maps.ensembl.amps, maps.ensembl.more)

# -----------------------------------------------------------------------------
# Method: Expr ~ iCN
# Last Modified: 07/12/16
# -----------------------------------------------------------------------------
runLM_Mut <- function(ensembl_id, samples, expr.log2, results.src) {
   results.src.gene <- results.src[ensembl_id,]
   
   expr.log2.gene <- toTable(NA, 4, length(samples), c("Expr", "iCN_median", "iCN_max", "iCN_sum"))
   expr.log2.gene$Expr <- expr.log2[ensembl_id,]
   rownames(expr.log2.gene) <- samples
   
   for (s in 1:length(samples)) {
   	  sample <- samples[s]
   	  
   	  ## INPUT: CNA
      cna.file <- paste(path, "ngs/WGS/", sample, "/", sample, "_ANALYSIS/", sample, "_iCN.seg", sep="")
      cna <- readTable(cna.file, header=F, rownames=F, sep="\t")
      colnames(cna) <- c("Sample", "Chr", "Start", "End", "V5", "iCN")
      rownames(cna) <- paste(cna$Chr, cna$Start, cna$End, sep="-")
      
      ## Breakpoint1
      chr   <- paste("chr", results.src.gene$chromosome_name, sep="")
      start <- results.src.gene$start_position
      end   <- results.src.gene$end_position
            
   	  cna.gene <- subset(cna, Chr == chr)
   	  cna.gene <- subset(cna.gene, End >= start)
   	  cna.gene <- subset(cna.gene, Start <= end)   	  

      expr.log2.gene$iCN_median[s] <- median(cna.gene$iCN)    
      expr.log2.gene$iCN_max[s] <- max(cna.gene$iCN)         
      expr.log2.gene$iCN_sum[s] <- sum(cna.gene$iCN) 
   }
   
   pdf(paste("analysis/G4/plots/expr.log2.", results.src.gene$hgnc_symbol, ".pdf", sep=""))
   plot(expr.log2.gene$iCN_max, expr.log2.gene$Expr, xlab="Max iCN", ylab="log2(1+count)", pch=16, col="blue", main=results.src.gene$hgnc_symbol)
   abline(lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN_max)[[1]][[1]], lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN_max)[[1]][[2]])
   dev.off()
   
   writeTable(expr.log2.gene, paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), colnames=T, rownames=T, sep="\t")

   #return(summary(lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN))[[4]][[8]])
}

runLM <- function(ensembl_id, samples, expr.log2, results.src) {
   results.src.gene <- results.src[ensembl_id,]
   
   expr.log2.gene <- toTable(NA, 4, length(samples), c("Expr", "iCN_median", "iCN_max", "iCN_sum"))
   expr.log2.gene$Expr <- expr.log2[ensembl_id,]
   rownames(expr.log2.gene) <- samples
   
   for (s in 1:length(samples)) {
   	  sample <- samples[s]
   	  
   	  ## INPUT: CNA
      cna.file <- paste(path, "ngs/WGS/", sample, "/", sample, "_ANALYSIS/", sample, "_iCN.seg", sep="")
      cna <- readTable(cna.file, header=F, rownames=F, sep="\t")
      colnames(cna) <- c("Sample", "Chr", "Start", "End", "V5", "iCN")
      rownames(cna) <- paste(cna$Chr, cna$Start, cna$End, sep="-")
      
      ## Breakpoint1
      chr   <- paste("chr", results.src.gene$chromosome_name, sep="")
      start <- results.src.gene$start_position
      end   <- results.src.gene$end_position
            
   	  cna.gene <- subset(cna, Chr == chr)
   	  cna.gene <- subset(cna.gene, End >= start)
   	  cna.gene <- subset(cna.gene, Start <= end)   	  

      expr.log2.gene$iCN_median[s] <- median(cna.gene$iCN)    
      expr.log2.gene$iCN_max[s] <- max(cna.gene$iCN)         
      expr.log2.gene$iCN_sum[s] <- sum(cna.gene$iCN) 
   }
   
   pdf(paste("analysis/G4/plots/expr.log2.", results.src.gene$hgnc_symbol, ".pdf", sep=""))
   plot(expr.log2.gene$iCN_max, expr.log2.gene$Expr, xlab="Max iCN", ylab="log2(1+count)", pch=16, col="blue", main=results.src.gene$hgnc_symbol)
   abline(lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN_max)[[1]][[1]], lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN_max)[[1]][[2]])
   dev.off()
   
   writeTable(expr.log2.gene, paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), colnames=T, rownames=T, sep="\t")

   #return(summary(lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN))[[4]][[8]])
}

runLM_Slope <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(lm(expr.log2.gene$Expr ~ expr.log2.gene$iCN)[[1]][[2]])
}

iCN_median <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(median(expr.log2.gene$iCN_median))
}

iCN_max <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(max(expr.log2.gene$iCN_max))
}

iCN_sum <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(sum(expr.log2.gene$iCN_sum))
}

runSRC_iCNmedian <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN_median, method="spearman", alternative="g", exact=F)[3]))
}

runSRC_iCNmedian_rho <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN_median, method="spearman", alternative="g", exact=F)[4]))
}

runSRC_iCNmax <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN_max, method="spearman", alternative="g", exact=F)[3]))
}

runSRC_iCNmax_rho <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN_max, method="spearman", alternative="g", exact=F)[4]))
}

runSRC_iCNsum <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN_sum, method="spearman", alternative="g", exact=F)[3]))
}

runSRC_iCNsum_rho <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN_sum, method="spearman", alternative="g", exact=F)[4]))
}

runKendall <- function(ensembl_id) {
   results.src.gene <- results.src[ensembl_id,]
   expr.log2.gene <- readTable(paste("analysis/G4/data/expr.log2.", results.src.gene$hgnc_symbol, ".txt", sep=""), header=T, rownames=T, sep="\t")

   return(as.numeric(cor.test(expr.log2.gene$Expr, expr.log2.gene$iCN, method="kendall", alternative="g")[3]))
}

###
##
#sclc <- c("ENSG00000116990", "ENSG00000134323", "ENSG00000136997", "ENSG00000270141", "ENSG00000105173", "ENSG00000114127", "ENSG00000088930", "ENSG00000085224", "ENSG00000164362", "ENSG00000135446", "ENSG00000133703", "ENSG00000146648", "ENSG00000141736", "ENSG00000141510", "ENSG00000139687", "ENSG00000078900", "ENSG00000147889", "ENSG00000189283")
sclc <- maps.ensembl.amps$ensembl_gene_id

###
##
results.src <- toTable(NA, 10, length(sclc), c("Length", "SRC_median", "SRC_median_rho", "SRC_max", "SRC_max_rho", "SRC_sum", "SRC_sum_rho", "iCN_median", "iCN_max", "iCN_sum"))
rownames(results.src) <- sclc

##
maps.ensembl.sclc <- subset(maps.ensembl, ensembl_gene_id == sclc[1])
for (x in 2:length(sclc)) {
   maps.ensembl.sclc <- rbind(maps.ensembl.sclc, subset(maps.ensembl, ensembl_gene_id == sclc[x]))
}
results.src <- cbind(maps.ensembl.sclc, results.src)
rownames(results.src) <- sclc

samples <- names(expr.log2[1,])
results.src$Length <- results.src$end_position - results.src$start_position
mapply(x = 1:length(sclc), function(x) runLM(sclc[x], samples, expr.log2, results.src))
#results.src$LM_Slope <- mapply(x = 1:length(sclc), function(x) runLM_Slope(sclc[x]))
results.src$SRC_median <- mapply(x = 1:length(sclc), function(x) runSRC_iCNmedian(sclc[x]))
results.src$SRC_median_rho <- mapply(x = 1:length(sclc), function(x) runSRC_iCNmedian_rho(sclc[x]))
results.src$SRC_max <- mapply(x = 1:length(sclc), function(x) runSRC_iCNmax(sclc[x]))
results.src$SRC_max_rho <- mapply(x = 1:length(sclc), function(x) runSRC_iCNmax_rho(sclc[x]))
results.src$SRC_sum <- mapply(x = 1:length(sclc), function(x) runSRC_iCNsum(sclc[x]))
results.src$SRC_sum_rho <- mapply(x = 1:length(sclc), function(x) runSRC_iCNsum_rho(sclc[x]))

results.src$iCN_median <- mapply(x = 1:length(sclc), function(x) iCN_median(sclc[x]))
results.src$iCN_max <- mapply(x = 1:length(sclc), function(x) iCN_max(sclc[x]))
results.src$iCN_sum <- mapply(x = 1:length(sclc), function(x) iCN_sum(sclc[x]))

#results.src <- results.src[order(results.src$SRC),]
results.src$BH_median <- p.adjust(results.src$SRC_median, "BH")
results.src$BH_max <- p.adjust(results.src$SRC_max, "BH")
results.src$BH_sum <- p.adjust(results.src$SRC_sum, "BH")

#results.src <- data.frame(lapply(results.src, as.character), stringsAsFactors=FALSE)
writeTable(results.src, "analysis/G4/expr.log2.lm.src_sum.txt", colnames=T, rownames=F, sep="\t")

##
results.all <- cbind(results.src, maps.ensembl.amps)
results.all[27,]$hgnc_symbol <- "MYCL1"
results.all[27,]$Oncogene <- "MYCL1"
#results.all <- results.all[order(results.all$SRC),]
writeTable(results.all, "analysis/G4/expr.log2.lm.src.G4_sum.txt", colnames=T, rownames=F, sep="\t")

##
# Association between expression and G4
results.all.nona <- subset(results.all, Ratio.OQs.random != "NA")
#results.all.nona <- subset(results.all.nona, SRC_rho > 0)

minusLog10P <- log10(results.all.nona$SRC_max) * -1
G4 <- as.numeric(results.all.nona$Ratio.OQs.random)
#G4 <- as.numeric(results.all.nona$OQs)

pdf("analysis/G4/expr.log2_SRC-G4_max.pdf")
plot(G4, minusLog10P, xlab="Ratio.OQ.random from Chambers et al", ylab="Max iCN -log10(SRC P) from George et al", pch=16, col="blue", main="Transcription activity (Max iCN) vs G4")
abline(lm(G4 ~ minusLog10P)[[1]][[1]], lm(G4 ~ minusLog10P)[[1]][[2]])
text(G4, minusLog10P, results.all.nona$hgnc_symbol, cex=0.6, pos=3) 
dev.off()

cor.test(minusLog10P, G4, method="spearman", alternative="g", exact=F)
cor.test(minusLog10P, G4, method="spearman", alternative="g", exact=F)[3]

summary(lm(G4 ~ minusLog10P))
lm(iCN ~ minusLog10P)[[1]][[2]]

##
minusLog10P <- log10(results.all.nona$SRC_sum) * -1
G4 <- as.numeric(results.all.nona$Ratio.OQs.random)
#G4 <- as.numeric(results.all.nona$OQs)

pdf("analysis/G4/expr.log2_SRC-G4_sum.pdf")
plot(G4, minusLog10P, xlab="Ratio.OQ.random from Chambers et al", ylab="Total iCN -log10(SRC P) from George et al", pch=16, col="blue", main="Transcription activity (Total iCN) vs G4")
abline(lm(G4 ~ minusLog10P)[[1]][[1]], lm(G4 ~ minusLog10P)[[1]][[2]])
text(G4, minusLog10P, results.all.nona$hgnc_symbol, cex=0.6, pos=3) 
dev.off()

cor.test(minusLog10P, G4, method="spearman", alternative="g", exact=F)
cor.test(minusLog10P, G4, method="spearman", alternative="g", exact=F)[3]

##
minusLog10P <- log10(results.all.nona$SRC_max) * -1
iCN <- as.numeric(results.all.nona$iCN_max)

pdf("analysis/G4/expr.log2_SRC-iCN_max.pdf")
plot(iCN, minusLog10P, xlab="Max iCN from George et al", ylab="Max iCN -log10(SRC P) from George et al", pch=16, col="blue", main="Transcription activity vs Max iCN")
abline(lm(iCN ~ minusLog10P)[[1]][[1]], lm(iCN ~ minusLog10P)[[1]][[2]])
text(iCN, minusLog10P, results.all.nona$hgnc_symbol, cex=0.6, pos=3) 
dev.off()

cor.test(minusLog10P, iCN, method="spearman", alternative="g", exact=F)[2]
cor.test(minusLog10P, iCN, method="spearman", alternative="g", exact=F)[3]

##
minusLog10P <- log10(results.all.nona$SRC_sum) * -1
iCN <- as.numeric(results.all.nona$iCN_sum)

pdf("analysis/G4/expr.log2_SRC-iCN_sum.pdf")
plot(iCN, minusLog10P, xlab="Total iCN from George et al", ylab="Total iCN -log10(SRC P) from George et al", pch=16, col="blue", main="Transcription activity vs Total iCN")
abline(lm(iCN ~ minusLog10P)[[1]][[1]], lm(iCN ~ minusLog10P)[[1]][[2]])
text(iCN, minusLog10P, results.all.nona$hgnc_symbol, cex=0.6, pos=3) 
dev.off()

cor.test(minusLog10P, iCN, method="spearman", alternative="g", exact=F)
cor.test(minusLog10P, iCN, method="spearman", alternative="g", exact=F)[3]

##
##
minusLog10P_median <- log10(results.all$SRC_median) * -1
minusLog10P_max <- log10(results.all$SRC_max) * -1
minusLog10P_sum <- log10(results.all$SRC_sum) * -1

pdf("analysis/G4/expr.log2_SRC-iCN_max-median_AHRR.pdf")
plot(minusLog10P_median, minusLog10P_max, xlab="Median iCN -log10(SRC P)", ylab="Max iCN -log10(SRC P)", pch=16, col="blue", main="Transcription activity in Max iCN vs Median iCN")
abline(lm(minusLog10P_median ~ minusLog10P_max)[[1]][[1]], lm(minusLog10P_median ~ minusLog10P_max)[[1]][[2]])
text(minusLog10P_median, minusLog10P_max, results.all$hgnc_symbol, cex=0.6, pos=3) 
dev.off()

cor.test(minusLog10P_max, minusLog10P_sum, method="spearman", alternative="g", exact=F)[2]
cor.test(minusLog10P_max, minusLog10P_sum, method="spearman", alternative="g", exact=F)[3]

pdf("analysis/G4/expr.log2_SRC-iCN_max-sum_AHRR.pdf")
plot(minusLog10P_sum, minusLog10P_max, xlab="Total iCN -log10(SRC P)", ylab="Max iCN -log10(SRC P)", pch=16, col="blue", main="Transcription activity in Max iCN vs Total iCN")
abline(lm(minusLog10P_sum ~ minusLog10P_max)[[1]][[1]], lm(minusLog10P_sum ~ minusLog10P_max)[[1]][[2]])
text(minusLog10P_sum, minusLog10P_max, results.all$hgnc_symbol, cex=0.6, pos=3) 
dev.off()

cor.test(minusLog10P_max, minusLog10P_sum, method="spearman", alternative="g", exact=F)[2]
cor.test(minusLog10P_max, minusLog10P_sum, method="spearman", alternative="g", exact=F)[3]





# -----------------------------------------------------------------------------
# Method: Differential expression (Read counts)
# Usage: XRN1 (n=) vs XRN1 WT
# Last Modified: 08/12/16
# -----------------------------------------------------------------------------
## Read counts
load("/Users/tpyang/Work/uni-koeln/tyang2/Link2Analysis/SCLC/metadata/exp.counts.rda")

colnames <- gsub("AB", "", colnames(exp.counts))
colnames <- gsub("T", "", colnames)
colnames <- gsub("o", "", colnames)
colnames(exp.counts) <- colnames

overlaps <- intersect(samples$Sample, colnames(exp.counts))
expr <- exp.counts[,overlaps]

samples.high <- unique(c(subset(muts.non, Gene_Hugo == "XRN1")$PAT_ID, subset(muts.non, Gene_Hugo == "XRN2")$PAT_ID, subset(muts.non, Gene_Hugo == "ATRX")$PAT_ID))
samples.high <- intersect(overlaps, samples.high)
samples.low  <- setdiff(overlaps, samples.high)

> length(samples.high)
[1] 7
> length(samples.low)
[1] 60

> length(samples.high)
[1] 12
> length(samples.low)
[1] 55

## log2(1 + FPKM)
expr.log2 <- log2(1 + expr)

###
## Differential expression (t-Test)
ttest <- function(p, expr.log2, samples.high, samples.low) {
   expr.log2.high <- expr.log2[p, samples.high]
   expr.log2.low  <- expr.log2[p, samples.low]

   return(t.test(expr.log2.high, expr.log2.low)$p.value)
}

wtest <- function(p, expr.log2, phenos.expr) {
   return(wilcox.test(as.numeric(expr.log2[p,]) ~ phenos.expr$Group, exact=FALSE)[[3]])
}

cmedian <- function(p, expr.log2, samples.group) {
   return(median(as.numeric(expr.log2[p, samples.group])))
}

###
## Main (ALL)
expr.log2.tmp <- expr.log2

results.expr <- toTable(0, 2, nrow(expr.log2.tmp), c("P", "BH"))
rownames(results.expr) <- rownames(expr.log2.tmp)

## Effect size
results.expr$High <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.high))
results.expr$Low  <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.low))
results.expr$Effect <- results.expr$High - results.expr$Low

###
## Main (Only varible)
#hist(abs(results.expr$Effect))

results.expr.effect <- rbind(subset(results.expr, Effect >= 1.5), subset(results.expr, Effect <= -1.5))
expr.log2.effect <- expr.log2.tmp[rownames(results.expr.effect),]

##
expr.log2.tmp <- expr.log2.effect

##
maps.ensembl <- getEnsemblGeneID(rownames(expr.log2.tmp)[1], ensembl)
for (e in 2:nrow(expr.log2.tmp)) {
   maps.ensembl <- rbind(maps.ensembl, getEnsemblGeneID(rownames(expr.log2.tmp)[e], ensembl))
}
rownames(maps.ensembl) <- maps.ensembl$ensembl_gene_id

##
results.expr <- toTable(0, 2, nrow(expr.log2.tmp), c("P", "BH"))
rownames(results.expr) <- rownames(expr.log2.tmp)

## Effect size
results.expr$High <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.high))
results.expr$Low  <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.low))
results.expr$Effect <- results.expr$High - results.expr$Low

## P
results.expr$P <- mapply(x = 1:nrow(expr.log2.tmp), function(x) ttest(x, expr.log2.tmp, samples.high, samples.low))
results.expr <- results.expr[order(results.expr$P),]

## BH
results.expr <- results.expr[order(results.expr$P),]
results.expr$BH <- p.adjust(results.expr$P, "BH")

results.expr.maps <- cbind(results.expr, maps.ensembl[rownames(results.expr),])

## W
phenos.expr <- phenos[overlaps,]
phenos.expr$Group <- as.factor(phenos.expr$Group)

results.expr <- results.expr[rownames(expr.log2.tmp),]
results.expr$W <- mapply(x = 1:nrow(expr.log2.tmp), function(x) wtest(x, expr.log2.tmp, phenos.expr))

## WBH
results.expr <- results.expr[order(results.expr$W),]
results.expr$WBH <- p.adjust(results.expr$W, "BH")

results.expr.maps <- cbind(results.expr, maps.ensembl[rownames(results.expr),])

writeTable(results.expr.maps, "analysis/report_expr_counts_log2_XRN1_ttest_wilcox_effect1.5.txt", colnames=T, rownames=T, sep="\t")

sclc <- c("ENSG00000116990", "ENSG00000134323", "ENSG00000136997", "ENSG00000270141", "ENSG00000105173", "ENSG00000114127", "ENSG00000088930", "ENSG00000085224", "ENSG00000164362")

writeTable(results.expr[sclc,], "analysis/report_expr_counts_log2_19H_ttest_wilcox_sclc.txt", colnames=T, rownames=T, sep="\t")
save(maps.ensembl, file="analysis/maps.ensembl.RData")










# -----------------------------------------------------------------------------
# Method: Differential expression (Read counts)
# Usage: 
# Last Modified: 03/11/16
# -----------------------------------------------------------------------------
###
## Main
#samples <- readTable("analysis/sv_inv_bfb_sclc_c3.list", header=T, rownames=F, sep="\t")
#rownames(samples) <- samples$Sample

## FPKM
expr <- readTable("metadata/nature14664-s1_ST10.txt", header=T, rownames=2, sep="\t")
overlaps <- intersect(samples$Sample, colnames(expr))
maps <- expr[,1:2]
expr <- expr[,-c(1,2)]
expr <- expr[,overlaps]

## Read counts
colnames <- gsub("AB", "", colnames(exp.counts))
colnames <- gsub("T", "", colnames)
colnames <- gsub("o", "", colnames)
colnames(exp.counts) <- colnames

expr <- exp.counts
overlaps <- intersect(samples$Sample, colnames(exp.counts))
expr <- expr[,overlaps]

samples.high <- c("S02360", "S00022", "S00831", "S02273", "S02285", "S02291")
samples.high <- intersect(overlaps, samples.high)
samples.low  <- setdiff(overlaps, samples.high)

samples.high <- intersect(overlaps, subset(samples, Group == "H")$Sample)
samples.low  <- intersect(overlaps, subset(samples, Group == "L")$Sample)
> length(samples.high)
[1] 32
> length(samples.low)
[1] 37

## log2(1 + FPKM)
#expr[expr < 3] <- 0
expr.log2 <- log2(1 + expr)

## Remove IDs with missing data
remove0 <- function(p, expr.log2) {
   expr.log2.tmp <- as.numeric(expr.log2[p,])

   if (length(expr.log2.tmp[expr.log2.tmp == 0]) == 0)
      return(T)
   return (F)
}

remove00 <- function(p, expr.log2) {
   expr.log2.tmp <- as.numeric(expr.log2[p,])

   if (length(expr.log2.tmp[expr.log2.tmp == 0]) == 56)
      return(T)
   return (F)
}

remove50 <- function(p, expr) {
   expr.log2.tmp <- as.numeric(expr.log2[p,])

   if (length(expr.log2.tmp[expr.log2.tmp == 0]) == 56)
      return(T)
   return (F)
}

## Removed ID with all 0s
expr.log2$Missing <- NA
expr.log2$Missing <- mapply(x = 1:nrow(expr.log2), function(x) remove00(x, expr.log2))
expr.log2.no0 <- subset(expr.log2, Missing == F)
expr.log2.no0 <- expr.log2.no0[,-70]

> dim(expr.log2)
[1] 31640    69
> dim(expr.log2.no0)
[1] 31437    69

expr.log2$Missing  <- mapply(x = 1:nrow(expr.log2), function(x) remove0(x, expr.log2))
expr.log2.nona <- subset(expr.log2, Missing == T)
expr.log2.nona <- expr.log2.nona[,-70]

###
## Differential expression (t-Test)
ttest <- function(p, expr.log2, samples.high, samples.low) {
   expr.log2.high <- expr.log2[p, samples.high]
   expr.log2.low  <- expr.log2[p, samples.low]

   return(t.test(expr.log2.high, expr.log2.low)$p.value)
}

wtest <- function(p, expr.log2, phenos.expr) {
   return(wilcox.test(as.numeric(expr.log2[p,]) ~ phenos.expr$Group, exact=FALSE)[[3]])
}

cmedian <- function(p, expr.log2, samples.group) {
   return(median(as.numeric(expr.log2[p, samples.group])))
}

###
## Main (ALL)
expr.log2.tmp <- expr.log2
#expr.log2.tmp <- expr.log2.no0
#expr.log2.tmp <- expr.log2.nona

results.expr <- toTable(0, 2, nrow(expr.log2.tmp), c("P", "BH"))
rownames(results.expr) <- rownames(expr.log2.tmp)
results.expr$Gene <- maps[rownames(results.expr),]$gene

## Effect size
results.expr$High <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.high))
results.expr$Low  <- mapply(x = 1:nrow(expr.log2.tmp), function(x) cmedian(x, expr.log2.tmp, samples.low))
results.expr$Effect <- results.expr$High - results.expr$Low

## P
results.expr$P <- mapply(x = 1:nrow(expr.log2.tmp), function(x) ttest(x, expr.log2.tmp, samples.high, samples.low))
results.expr <- results.expr[order(results.expr$P),]

## BH
results.expr <- results.expr[order(results.expr$P),]
results.expr$BH <- p.adjust(results.expr$P, "BH")

## W
phenos.expr <- phenos[overlaps,]
phenos.expr$Group <- as.factor(phenos.expr$Group)

results.expr$W <- mapply(x = 1:nrow(expr.log2.tmp), function(x) wtest(x, expr.log2.tmp, phenos.expr))

## WBH
results.expr$WBH <- p.adjust(results.expr$W, "BH")

## Output
results.expr.no0 <- results.expr
results.expr.no0.effect <- results.expr

writeTable(results.expr, "analysis/report_sv_5mb_q75_c3_expr_no0_ttest_wilcox.txt", colnames=T, rownames=T, sep="\t")

## 
results.expr.no0 <- results.expr
results.expr.nona.effect <- results.expr

writeTable(results.expr, "analysis/report_sv_5mb_q75_c3_expr_no0_effect1_ttest_wilcox.txt", colnames=T, rownames=T, sep="\t")
results.expr.nona <- readTable("analysis/expression/report_sv_5mb_q75_c3_expr_nona_ttest_wilcox.txt", header=T, rownames=T, sep="\t")

###
## Main (Only varible)
results.expr.effect <- rbind(subset(results.expr, Effect >= 2), subset(results.expr, Effect <= -2))
expr.log2.effect <- expr.log2.tmp[rownames(results.expr.effect),]

expr.log2.tmp <- expr.log2.effect

# -----------------------------------------------------------------------------
# Method: Volvano plots
# Usage: 
# Last Modified: 05/11/16
# -----------------------------------------------------------------------------
results.tmp <- results.expr[which((!is.na(results.expr$W))),]
#results.tmp <- readTable("analysis/report_sv_5mb_q75_c3_expr_fpkm<3_ttest.txt", header=T, rownames=T, sep="\t")
xmax <- max(results.tmp$Effect)
ymax <- max(-log10(results.tmp$W))
p <- 2.622498e-04

pdf("report_sv_5mb_q75_n10_expr_19H_wilcox_volcano_effect1.5.pdf")
plot(results.tmp$Effect, -log10(results.tmp$W), xlim=c(-xmax, xmax), ylim=c(0,ymax), xlab="Effect size", ylab="-log10(P)", col="darkgray", main="Differential Expression of High (n=9) vs Low (n=58) in SCLC")

results.lightred <- subset(subset(results.tmp, W <= p), Effect < -1.5)
points(results.lightred$Effect, -log10(results.lightred$W))

results.red <- subset(results.lightred, Effect <= -1.5)
points(results.red$Effect, -log10(results.red$W), col="blue")

results.lightblue <- subset(subset(results.tmp, W <= p), Effect > 1.5)
points(results.lightblue$Effect, -log10(results.lightblue$W))

results.blue <- subset(results.lightblue, Effect >= 1.5)
points(results.blue$Effect, -log10(results.blue$W), col="red")

abline(h=c(-log10(p)), lty=5)
abline(v=c(1.5), col="red", lty=5)
abline(v=c(-1.5), col="blue", lty=5)
#legend(x=-xmax, y=ymax, legend=c("Up-regulation","Down"), fill=c("white","white"), border=c("red","blue"), pch=21)
#legend(x=-xmax, y=ymax, legend=c("Up-regulation","Down"), col=c("red","blue"), pch=21)
legend(x=-7.3, y=ymax, legend=c("Upregulation","Downregulation"), col=c("red","blue"), pch=21)
dev.off()








##
# Wilcox
probe <- "ENSG00000179546"
gene  <- "HTR1D"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 15), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=0, y=15, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000066248"
gene  <- "NGEF"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(1, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 12), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=1, y=12, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000126878"
gene  <- "AIF1L"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(6, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 22), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=6, y=22, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000265646"
gene  <- "TUFMP1"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 52), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=6.3, y=52, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000272626"
gene  <- "RP11-468N14.13 lincRNA"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 23), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=4.3, y=23, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

##
# SCLC
probe <- "ENSG00000116990"
gene  <- "MYCL"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(7, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 12), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=13.6, y=12, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000114127"
gene  <- "XRN1"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(9.5, 13), ylim=c(0, 24), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=9.5, y=24, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000088930"
gene  <- "XRN2"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(10, 14), ylim=c(0, 16), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=10, y=16, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000085224"
gene  <- "ATRX"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(11, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 17), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=11, y=17, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()

probe <- "ENSG00000164362"
gene  <- "TERT"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.counts.log2.19H.wilcox_", gene, ".pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(2, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 23), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+count) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=10, y=23, legend=c("High (n=9)", "Low (n=58)"), fill=c(red, blue), border=c(red, blue))
dev.off()










# -----------------------------------------------------------------------------
# Method: Volvano plots
# Usage: 
# Last Modified: 05/11/16
# -----------------------------------------------------------------------------
results.tmp <- results.expr
#results.tmp <- readTable("analysis/report_sv_5mb_q75_c3_expr_fpkm<3_ttest.txt", header=T, rownames=T, sep="\t")
xmax <- max(results.tmp$Effect)
ymax <- max(-log10(results.tmp$P))
p <- 0.001192877

pdf("report_sv_5mb_q75_c3_expr_nona_effect1_ttest1.2e-3_volcano.pdf")
plot(results.tmp$Effect, -log10(results.tmp$P), xlim=c(-xmax, xmax), ylim=c(0,ymax), xlab="Effect size", ylab="-log10(P)", col="darkgray", main="Differential Expression of High (n=32) vs Low (n=37) in SCLC")

results.lightred <- subset(subset(results.tmp, P <= p), Effect < -1)
points(results.lightred$Effect, -log10(results.lightred$P))

results.red <- subset(results.lightred, Effect <= -1)
points(results.red$Effect, -log10(results.red$P), col="blue")

results.lightblue <- subset(subset(results.tmp, P <= p), Effect > 1)
points(results.lightblue$Effect, -log10(results.lightblue$P))

results.blue <- subset(results.lightblue, Effect >= 1)
points(results.blue$Effect, -log10(results.blue$P), col="red")

abline(h=c(-log10(p)), lty=5)
abline(v=c(1), col="red", lty=5)
abline(v=c(-1), col="blue", lty=5)
#legend(x=-xmax, y=ymax, legend=c("Up-regulation","Down"), fill=c("white","white"), border=c("red","blue"), pch=21)
#legend(x=-xmax, y=ymax, legend=c("Up-regulation","Down"), col=c("red","blue"), pch=21)
legend(x=2, y=ymax, legend=c("Upregulation","Downregulation"), col=c("red","blue"), pch=21)
dev.off()

##
# t-test
probe <- "NM_003839"
gene  <- "TNFRSF11A"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.log2.ttest_", gene, "_overlapping.pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+FPKM) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=3.9, y=12, legend=c("High (n=32)", "Low (n=37)"), fill=c(red, blue), border=c(red, blue))
dev.off()

##
probe <- "NM_003012"
gene  <- "SFRP1"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.log2.ttest_", gene, "_overlapping.pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 12), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+FPKM) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=0, y=12, legend=c("High (n=32)", "Low (n=37)"), fill=c(red, blue), border=c(red, blue))
dev.off()

##
probe <- "NM_002825"
gene  <- "PTN"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.log2.ttest_", gene, "_overlapping.pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 12), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+FPKM) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=0, y=12, legend=c("High (n=32)", "Low (n=37)"), fill=c(red, blue), border=c(red, blue))
dev.off()

## 
# Wilcox
probe <- "NM_020630"
gene  <- "RET"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.log2.wilcox_", gene, "_overlapping.pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), ylim=c(0, 11), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+FPKM) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=3.2, y=11, legend=c("High (n=32)", "Low (n=37)"), fill=c(red, blue), border=c(red, blue))
dev.off()

##
probe <- "NM_001172410"
gene  <- "SFTPC"
red  <- rgb(1, 0, 0, 0.5)
blue <- rgb(0, 0, 1, 0.5)
pdf(paste("expr.log2.wilcox_", gene, "_overlapping.pdf", sep=""))
hist(as.numeric(expr.log2[probe, samples.high]), col=red, xlim=c(0, max(as.numeric(expr.log2[probe,]))), main=paste(gene, " (", probe, ")", sep=""), xlab="log2(1+FPKM) Expression")
hist(as.numeric(expr.log2[probe, samples.low]), col=blue, add=T)
box()
legend(x=8.6, y=16, legend=c("High (n=32)", "Low (n=37)"), fill=c(red, blue), border=c(red, blue))
dev.off()

# -----------------------------------------------------------------------------
# Method: XRN1 and ATRX
# Usage: 
# Last Modified: 22/11/16
# -----------------------------------------------------------------------------
## ATRX expression
> median(as.numeric(expr.log2.no0["NM_000489", atrx.mut]))
[1] 2.788777
> median(as.numeric(expr.log2.no0["NM_000489", atrx.wt]))
[1] 2.931437

> mean(as.numeric(expr.log2.no0["NM_000489", atrx.mut]))
[1] 2.750602
> mean(as.numeric(expr.log2.no0["NM_000489", atrx.wt]))
[1] 2.712417

## XRN2 expression
writeTable(t(expr.log2.no0[c("NM_019001", "NM_001042604", "NM_012255"),]), "analysis/report_sv_5mb_q75_c3_expr_no0_XRN1_2.txt", colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Method: RB1
# Usage: 
# Last Modified: 23/11/16
# -----------------------------------------------------------------------------
###
## Main
rb1 <- readTable("metadata/nature14664-s1_ST4.txt", header=T, rownames=F, sep="\t")

> rb1.rb1.intra <- rb1.rb1[which(rb1.rb1$Chr_Pair1 == rb1.rb1$Chr_Pair2),]
> rb1.rb1.intra$Size <- abs(as.numeric(rb1.rb1.intra$Exact_Pos_Pair_1) - as.numeric(rb1.rb1.intra$Exact_Pos_Pair_2))
> rb1.rb1.intra

writeTable(rb1.samples, "metadata/nature14664-s1_ST4_RB1.samples.txt", colnames=F, rownames=F, sep="")
writeTable(p73.samples, "metadata/nature14664-s1_ST4_TP73.samples.txt", colnames=F, rownames=F, sep="")

rb1.p73.intra <- rb1.p73[which(rb1.p73$Chr_Pair1 == rb1.p73$Chr_Pair2),]
rb1.p73.intra$Size <- abs(as.numeric(rb1.p73.intra$Exact_Pos_Pair_1) - as.numeric(rb1.p73.intra$Exact_Pos_Pair_2))
rb1.p73.intra






##
pdf("expr.log2_SFRP1_high.pdf")
hist(as.numeric(expr.log2["NM_003012", samples.high]), main="SFRP1 in High (n=32)")
dev.off()
pdf("expr.log2_SFRP1_low.pdf")
hist(as.numeric(expr.log2["NM_003012", samples.low]), main="SFRP1 in Low (n=37)")
dev.off()

##
pdf("expr.log2_TNFRSF11A_high.pdf")
hist(as.numeric(expr.log2["NM_003839", samples.high]), main="TNFRSF11A in High (n=32)")
dev.off()
pdf("expr.log2_TNFRSF11A_low.pdf")
hist(as.numeric(expr.log2["NM_003839", samples.low]), main="TNFRSF11A in Low (n=37)")
dev.off()







##
pdf("expr.log2_CA2_high.pdf")
hist(as.numeric(expr.log2["NM_000067", samples.high]), main="CA2 in High (n=38)")
dev.off()
pdf("expr.log2_CA2_low.pdf")
hist(as.numeric(expr.log2["NM_000067", samples.low]), main="CA2 in Low (n=31)")
dev.off()

pdf("expr.log2_GLIPR2_high.pdf")
hist(as.numeric(expr.log2["NM_022343", samples.high]), main="GLIPR2 in High (n=38)")
dev.off()
pdf("expr.log2_GLIPR2_low.pdf")
hist(as.numeric(expr.log2["NM_022343", samples.low]), main="GLIPR2 in Low (n=31)")
dev.off()

pdf("expr.log2_CTNNB1_high.pdf")
hist(as.numeric(expr.log2["NM_001904", samples.high]), main="CTNNB1 in High (n=38)")
dev.off()
pdf("expr.log2_CTNNB1_low.pdf")
hist(as.numeric(expr.log2["NM_001904", samples.low]), main="CTNNB1 in Low (n=31)")
dev.off()

## 
samples.myc <- c("S01512", "S02273", "S02375", "S02397", "S02065", "S00938", "S01023", "S02292", "S02403", "S01022", "S02299")
samples.myc.overlaps <- intersect(colnames(expr.log2.nona), samples.myc)

expr.log2.nona["NM_022343", samples.myc.overlaps]








# -----------------------------------------------------------------------------
# Method: Anti-correction genes
# Usage: 
# Last Modified: 07/11/16
# -----------------------------------------------------------------------------










###
## Removed not-variale probes
#expr.log2.nona$SD <- mapply(x = 1:nrow(expr.log2.nona), function(x) sd(as.numeric(expr.log2.nona[x,])))
#dim(expr.log2.nona)
#dim(subset(expr.log2.nona, SD >= 0.5))
#expr.log2.nona.var <- subset(expr.log2.nona, SD >= 0.8)

###
## PCA
expr.pca <- prcomp(expr.log2.nona, retx=TRUE, center=TRUE, scale=TRUE)
scores <- expr.pca$x
sd <- expr.pca$sdev
loadings <- expr.pca$rotation
summary <- summary(expr.pca) 

pca <- as.data.frame(scores)

## PCA outliers
trait.v <- samples[overlaps,]$Group

trait <- sort(unique(trait.v))
color <- rainbow(length(trait))
pca.point.color = function(model.id) {
   return(color[which(model.id == trait)])
}

pdf("pca_PC1-2_nona.pdf")
plot(pca$PC1, pca$PC2, col=sapply(trait.v, pca.point.color))
dev.off()
legend(600, 580, trait, bty="n", fill=color)

## Q
library(qvalue)
q <- qvalue(results.expr.nona$P)
results.expr.nona$Q <- q$qvalue

q <- qvalue(results.expr.nona$W)
results.expr.nona$WQ <- q$qvalue

###
## Fisher's test for anti-corr mutations
data.test <- mut

results <- toTable(0, 2, 0, c("Pairs", "P"))
for (x in 1:(nrow(data.test)-1)) {
   for (y in (x+1):nrow(data.test)) {
      #println(paste(x, y, sep="\t"))
      result <- toTable(0, 2, 1, c("Pairs", "P"))
      
      data.fisher <- rbind(data.test[x,], data.test[y,])
      result[1, 1] <- paste(rownames(data.test)[x], rownames(data.test)[y], sep=" vs ")
      result[1, 2] <- fisher.test(data.fisher)[[1]]
      results <- rbind(results, result)
   }
}
results <- results[order(results$P),]
results$BH <- p.adjust(results$P, "BH")

## 
results.all <- toTable(0, length(genes), length(samples), genes)
for (x in 1:nrow(samples)) {
   sample <- samples[x, 1]
   muts.non.sample <- subset(muts.non, PAT_ID == sample)
   
   for (y in 1:length(genes)) {
   	  gene <- genes[y]
      results.all[x, y] <- nrow(subset(muts.non.sample, Gene_Hugo == gene))
   }
}
writeTable(results.all, "report_5mb_q75_c3_muts.txt", colnames=T, rownames=T, sep="\t")
writeTable(mut, "report_5mb_q75_c3_mut.txt", colnames=T, rownames=T, sep="\t")
writeTable(results, "report_5mb_q75_c3_mut_pairs.txt", colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Method: Permutation
# Usage: 
# Last Modified: 07/11/16
# -----------------------------------------------------------------------------
permutationP <- function(x, expr.log2, obvP) {
   perms <- c()
   for (y in 1:10000) {  
      s1 <- sample(1:69, 38, replace=F)
      s2 <- setdiff(1:69, s1)
 
      expr.log2.high <- expr.log2[x, s1]
      expr.log2.low  <- expr.log2[x, s2]

      expP <- t.test(expr.log2.high, expr.log2.low)$p.value
      perms <- c(perms, expP)
   }
   
   return(length(perms[(perms - obvP) <= 0])/10000)
}

permutationW <- function(x, expr.log2, obvW) {
   perms <- c()
   for (y in 1:10000) {  
      s1 <- sample(1:69, 38, replace=F)
      s2 <- setdiff(1:69, s1)
 
      expr.log2.high <- as.numeric(expr.log2[x, s1])
      expr.log2.low  <- as.numeric(expr.log2[x, s2])

      expW <- wilcox.test(expr.log2.high, expr.log2.low, , exact=FALSE)[[3]]
      perms <- c(perms, expW)
   }
   
   return(length(perms[(perms - obvW) <= 0])/10000)
}

results.expr.nona$FDRP <- mapply(x = 1:12, function(x) permutation(x, expr.log2.nona, results.expr.nona$P[x]))






## Remove probes
check0 <- function(p, expr.log2, samples) {
   probe <- rownames(expr.log2)[p]
   expr.log2.tmp <- as.numeric(expr.log2[probe, samples])

   return(length(expr.log2.tmp[expr.log2.tmp == 0]))
}

## Check missing data
results.expr$High0 <- mapply(x = 1:nrow(expr.log2), function(x) check0(x, expr.log2, samples.high))
results.expr$Low0  <- mapply(x = 1:nrow(expr.log2), function(x) check0(x, expr.log2, samples.low))

pdf("expr.log2_0_high.pdf")
hist(results.expr$High0, main="Missing in High (n=38)")
dev.off()
pdf("expr.log2_0_low.pdf")
hist(results.expr$Low0, main="Missing in Low (n=31)")
dev.off()

## Effect size
results.expr$High <- mapply(x = 1:nrow(expr.log2), function(x) cmedian(x, expr.log2, samples.high))
results.expr$Low  <- mapply(x = 1:nrow(expr.log2), function(x) cmedian(x, expr.log2, samples.low))
results.expr$Effect <- results.expr$High - results.expr$Low

results.expr$High_M <- mapply(x = 1:nrow(expr.log2), function(x) cmean(x, expr.log2, samples.high))
results.expr$Low_M  <- mapply(x = 1:nrow(expr.log2), function(x) cmean(x, expr.log2, samples.low))
results.expr$Effect_M <- results.expr$High_M - results.expr$Low_M

results.expr$P <- mapply(x = 1:nrow(expr.log2), function(x) ttest(x, expr.log2, samples.high, samples.low))
results.expr.nona <- results.expr[!is.na(results.expr$P),]

##
results.expr$W <- mapply(x = 1:nrow(expr.log2), function(x) wtest(x, expr.log2, phenos.expr))
results.expr.nona <- results.expr[!is.na(results.expr$W),]

##
results.expr.nona <- results.expr.nona[order(results.expr.nona$W),]
#results.expr.nona$BH <- p.adjust(results.expr.nona$P, "BH")

library(qvalue)
q <- qvalue(results.expr.nona$W)
results.expr.nona$WQ <- q$qvalue

##
#results.expr.nona$Gene <- maps[rownames(results.expr.nona),]$gene
writeTable(results.expr.nona, "analysis_sv_report_5mb_q75_c3_expr_w.txt", colnames=T, rownames=T, sep="\t")

###
## Differential expression (Linear model)
linearm <- function(p, expr.log2, phenos.expr) {
   fm1 <- lm(as.numeric(expr.log2[p,]) ~ 1 + phenos.expr$Group)
   fm2 <- lm(as.numeric(expr.log2[p,]) ~ 1)
   a1 <- anova(fm1, fm2)

   return(a1[2, 6])
}

phenos.expr <- phenos[overlaps,]
phenos.expr$Group <- as.factor(phenos.expr$Group)

results.expr <- toTable(0, 2, nrow(expr.log2), c("P", "Q"))
rownames(results.expr) <- rownames(expr.log2)

results.expr$P <- mapply(x = 1:nrow(expr.log2), function(x) linearm(x, expr.log2, phenos.expr))
results.expr.nona <- results.expr[!is.na(results.expr$P),]
results.expr.nona <- results.expr.nona[order(results.expr.nona$P),]

library(qvalue)
q <- qvalue(results.expr.nona$P)
results.expr.nona$Q <- q$qvalue

##
results.expr.nona$Gene <- maps[rownames(results.expr.nona),]$gene
writeTable(results.expr.nona, "analysis/report_sv_5mb_q75_c3_expr_ttest.txt", colnames=T, rownames=T, sep="\t")










# -----------------------------------------------------------------------------
# Method: Mutational signatures
# Usage: 
# Last Modified: 03/11/16
# -----------------------------------------------------------------------------
###
## Main
sigs <- readTable("analysis/signature/SignaturesSCLC.txt", header=T, rownames=F, sep="\t")






## SCLC (n=110)
> length(unique(subset(muts, Gene_Hugo == "RB1")$PAT_ID))
[1] 88
> length(unique(subset(muts.non, Gene_Hugo == "RB1")$PAT_ID))
[1] 87
> length(unique(subset(muts, Gene_Hugo == "CREBBP")$PAT_ID))
[1] 12
> length(unique(subset(muts.non, Gene_Hugo == "CREBBP")$PAT_ID))
[1] 11
> length(unique(subset(muts, Gene_Hugo == "EP300")$PAT_ID))
[1] 13
> length(unique(subset(muts.non, Gene_Hugo == "EP300")$PAT_ID))
[1] 13









##
test <- toTable(0, 2, 2, c("H", "L"))

> mean(as.numeric(phenos[samples[samples$Group == "H",]$Sample,]$age))
[1] 65.28302
> mean(as.numeric(phenos[samples[samples$Group == "L",]$Sample,]$age))
[1] 64.13462

t <- phenos[samples[samples$Group == "H",]$Sample,]$age
sd(t[!is.na(t)])
sd(as.numeric(phenos[samples[samples$Group == "L",]$Sample,]$age))

t.test(t[!is.na(t)], as.numeric(phenos[samples[samples$Group == "L",]$Sample,]$age))

> table(phenos[samples[samples$Group == "H",]$Sample,]$sex)
female   male 
    18     36 
> table(phenos[samples[samples$Group == "L",]$Sample,]$sex)
female   male 
    21     31 





##
## Kaplan-Meier survival analysis
phenos.sub <- phenos[!is.na(phenos$overall_survival..months.),]

## PFS_censor
phenos.sub$Groups <- phenos.sub$Group
phenos.sub$Groups <- gsub("L", 0, phenos.sub$Groups)
phenos.sub$Groups <- gsub("H", 1, phenos.sub$Groups)

## OS_censor
phenos.sub$OS_censor <- phenos.sub$Status..at.time.of.last.follow.up.
phenos.sub$OS_censor <- gsub("dead", 0, phenos.sub$OS_censor)
phenos.sub$OS_censor <- gsub("alive", 1, phenos.sub$OS_censor)
phenos.sub$OS_censor <- as.numeric(phenos.sub$OS_censor)

phenos.sub$OS_month <- phenos.sub$overall_survival..months.

##
## Kaplan-Meier survival analysis
filename <- paste(path, "analysis/survfit_", sep="")

## OS~kmeans
pdf(paste(filename, "OS~Groups", ".pdf", sep=""))
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub) 
plot(fit, ylim=c(0, 1), xlab="Months", ylab="Survival", col=c(1, 2), main="OS (Breakspoints)")
legend("topright", legend=c("Low", "High"), lwd=1, col=1:2)
dev.off()
