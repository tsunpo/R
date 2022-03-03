# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/icgc/clle-es-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/02/22; 14/10/19; 26/02/19
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "CLLE-ES"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, "ICGC", BASE, "ngs/WGS")
wd.anlys <- file.path(wd, "ICGC", BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "clle-es_wgs_n95.list"), header=F, rownames=F, sep="")$V3
samples0 <- readTable(file.path(wd.ngs, "clle-es_wgs_n95.list"), header=F, rownames=F, sep="")$V3
n1 <- length(samples1)
n0 <- length(samples0)

# -----------------------------------------------------------------------------
# CLL T vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8604102
# > max(cors.samples[,-c(1:4)])
# [1] 0.6677247

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.clle.es <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.clle.es, file.path(wd.ngs, "clle-es_wgs_m2_n95.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7863882 -0.7634611 -0.7497204 -0.7409367  0.4914746 

writeTable(subset(samples.clle.es, Q4 %in% c(4,1)), file.path(wd.ngs, "clle-es_wgs_q4_n48.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
icgc.list <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/icgc-wgs.list", header=F, rownames=T, sep="")
icgc.list$N   <- NA
icgc.list$COR <- NA
for (i in 1:nrow(icgc.list)) {
   icgc <- icgc.list[i,]
   table.i <- subset(table, project_code == icgc$V1)
 
   samples <- readTable(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/", icgc$V1, "/ngs/WGS/", icgc$V2, "_wgs_m2_n", nrow(table.i), ".txt"), header=T, rownames=T, sep="")   
   icgc.list$N[i] <- nrow(samples)
   icgc.list$COR[i] <- median(samples$COR)
}

###
##
samples.1 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/LUAD-UK/ngs/WGS/luad-us_wgs_m2_n37.txt", header=T, rownames=T, sep="")
samples.2 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ESAD-UK/ngs/WGS/esad-uk_wgs_m2_n88.txt", header=T, rownames=T, sep="")
samples.3 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/BRCA-US/ngs/WGS/brca-us_wgs_m2_n89.txt", header=T, rownames=T, sep="")
samples.4 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/BRCA-UK/ngs/WGS/brca-uk_wgs_m2_n44.txt", header=T, rownames=T, sep="")
samples.5 <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/CLLE-ES/ngs/WGS/clle-es_wgs_m2_n95.txt", header=T, rownames=T, sep="")
n.1 <- nrow(samples.1)
n.2 <- nrow(samples.2)
n.3 <- nrow(samples.3)
n.4 <- nrow(samples.4)
n.5 <- nrow(samples.5)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.1] <- 0
samples$CANCER[(1+n.1):(n.1+n.2)] <- 1
samples$CANCER[(1+n.1+n.2):(n.1+n.2+n.3)] <- 2
samples$CANCER[(1+n.1+n.2+n.3):(n.1+n.2+n.3+n.4)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_brca+esad+clle_Spearman's.pdf"), height=6, width=6)
ymax <- 0.8   #max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("", "", ""), ylim=c(ymin, ymax), ylab="", main="Tumour read depth vs. RT", yaxt="n", cex.axis=1.8, cex.lab=1.9, cex.main=2.1)
abline(h=0, lty=5, lwd=2)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.9)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col=blue, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col=blue.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col=red.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col=red, pch=19, cex=1, add=T)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Spearman's rho", side=2, line=2.7, cex=1.9)
#mtext("", cex=1.2, line=0.3)
mtext(text=c("BRCA", "ESAD", "CLLE"), side=1, cex=1.9, line=1.3, at=c(1,2,3))
mtext(text=c("n=89", "n=88", "n=95"), side=1, cex=1.8, line=3, at=c(1,2,3))
dev.off()
