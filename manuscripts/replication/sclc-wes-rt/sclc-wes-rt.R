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
load(file.path(wd.src.ref, "hg19.bed.gc.peiflyne.wes.new.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WES_NEW")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wes-rt"))
wd.rt.data  <- file.path(wd.rt, "data_wes_new")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples1 <- readTable(file.path(wd.ngs, "sclc_wes_n21-4.list"), header=F, rownames=F, sep="\t")
n1 <- length(samples1)

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Refer to cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.m_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

##
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
nrow(test)
# [1] 188867
test.nona <- test[removeNA(test),]
nrow(test.nona)
# [1] 161273
test.no0 <- test[getExpressed(test),]
test.no0 <- test.no0[,-(ncol(test.no0))]
nrow(test.no0)
# [1] 161273
pca.de <- getPCA(t(test.no0))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc-wes_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_sclc_chrs.RData")))
file.main <- c("SCLC WES", "")
trait <- as.numeric(samples.wes$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_SCLC-WES", size=6, file.main, "bottomright", c(red, red.lighter, blue.lighter, blue), NULL, flip.x=1, flip.y=1, legend.title=NA)

scores.wes <- pcaScores(pca.de)
y <- samples.wes$COR
x <- scores$PC1
file.name <- file.path(wd.rt.plots, paste0("correlation_SORTING_vs_PC1"))
plotCorrelation(file.name, "SCLC WES", paste0("PC", 1, " (", pcaProportionofVariance(pca.de, 1), "%)"), expression(italic("In silico") ~ "sorting"), x, y, size=6)

overlaps <- intersect(rownames(scores.wgs), rownames(scores.wes))
y <- scores.wgs[overlaps,]$PC1
x <- scores.wes[overlaps,]$PC1
file.name <- file.path(wd.rt.plots, paste0("correlation_PC1_WGS-vs-WES"))
plotCorrelation(file.name, "SCLC", paste0("WES PC", 1, " (", pcaProportionofVariance(pca.wes, 1), "%)"), paste0("WGS PC", 1, " (", pcaProportionofVariance(pca.wgs, 1), "%)"), x, y, size=6)

###
##
wgds <- readTable(file.path(wd.ngs, "sclc_wgs_n110_WGD.txt"), header=T, rownames=T, sep="\t")
wgds[which(wgds$wgd_predicted == 0), ]$wgd_predicted <- "non-WGD"
wgds[which(wgds$wgd_predicted == 1), ]$wgd_predicted <- "WGD"

overlaps <- intersect(rownames(scores.wgs), rownames(scores.wes))
x <- scores.wes[overlaps,]$PC1
y <- wgds[overlaps,]$decision_value
file.name <- file.path(wd.rt.plots, paste0("correlation_WES-PC1-vs-WGD"))
plotCorrelation(file.name, "SCLC", paste0("WES PC", 1, " (", pcaProportionofVariance(pca.wes, 1), "%)"), "WGD", x, y, size=6)

overlaps <- intersect(rownames(scores.wgs), rownames(scores.wes))
x <- scores.wgs[overlaps,]$PC1
y <- wgds[overlaps,]$decision_value
file.name <- file.path(wd.rt.plots, paste0("correlation_WGS-PC1-vs-WGD"))
plotCorrelation(file.name, "SCLC", paste0("WGS PC", 1, " (", pcaProportionofVariance(pca.wgs, 1), "%)"), "WGD", x, y, size=6)



# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 19/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.wes <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.wes, file.path(wd.ngs, "sclc_wes_m2_n17.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.6019724 -0.5883589 -0.5664459 -0.5596107 -0.5185071 

# -----------------------------------------------------------------------------
# 
# Last Modified: 03/06/22
# -----------------------------------------------------------------------------
y <- samples.sclc[rownames(samples.wes.new),]$COR
x <- samples.wes.new$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wgs_vs_wes-new_n17_WES_NEW_SPLINE"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES (new partitions; no seg)", "WGS", x, y)

y <- samples.sclc[rownames(samples.wes.new.T),]$COR
x <- samples.wes.new.T$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wgs_vs_wes-new_n17_WES_NEW_T"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES (new partitions; no seg)", "WGS", x, y)







# -----------------------------------------------------------------------------
# 
# Last Modified: 03/06/22
# -----------------------------------------------------------------------------
x <- samples.wes.nocn[rownames(samples.wes),]$COR
y <- samples.wes$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes_vs_wes-nocn_n17"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES (without corr. for CN)", "WES", x, y)

x <- samples.wes.wgs.cn[rownames(samples.wes),]$COR
y <- samples.wes$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes_vs_wes-wgs-cn_n17"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES (Corr. for WGS CN)", "WES", x, y)

x <- samples.wgs[rownames(samples.wes),]$COR
y <- samples.wes$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes_vs_wgs_n17"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WGS", "WES", x, y)

y <- samples.wgs[rownames(samples.wes),]$COR
x <- samples.wes$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wgs_vs_wes_n17"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES", "WGS", x, y)

y <- samples.wgs[rownames(samples.wes.1k.rd200),]$COR
x <- samples.wes.1k.rd200$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wgs_vs_wes-1k-rd200_n17"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES (1 kb; >200 reads)", "WGS", x, y)

y <- samples.sclc[rownames(samples.wes.new),]$COR
x <- samples.wes.new$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wgs_vs_wes-new_n17"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=17)"), "WES (new partitions)", "WGS", x, y)





x <- samples.wes.wgs.cn$COR
y <- samples.wes.nocn$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes-wgs_vs_wes-nocn_n21"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=21)"), "WES (Corr. for WGS CN)", "WES (without corr. for CN)", x, y)

x <- samples.wes.nocn[rownames(samples.wes),]$COR
y <- samples.wes$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes_vs_wes_n5"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=5)"), "WES (without corr. for CN)", "WES", x, y)

#overlaps1 <- intersect(rownames(samples.wes.nocn), rownames(samples.wgs))
y <- samples.wgs[overlaps,]$COR
x <- samples.wes.nocn[overlaps,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes_vs_wgs_n21"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=21)"), "WES (without corr. for CN)", "WGS", x, y)

y <- samples.wgs[overlaps,]$COR
x <- samples.wes.wgs.cn$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes-wgs_vs_wgs_n21"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=21)"), "WES (Corr. for WGS CN)", "WGS", x, y)






samples.wgs.nooutlier <- subset(samples.wgs, COR < 0)

overlaps <- intersect(rownames(samples.wes.nocn), rownames(samples.wgs.nooutlier))
x <- samples.wgs.nooutlier[overlaps,]$COR
y <- samples.wes.nocn[overlaps,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_wes_vs_wgs_n18"))
plotCorrelation(file.name, expression(italic("In silico") ~ "sorting of SCLC samples (n=18)"), "WGS", "WES", x, y)

