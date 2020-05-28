# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/expression/nbl-tpm-rt-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 28/05/20
# =============================================================================
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Commons.R", "Transcription.R", "TranscriptionReplicationConflict.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.ensGene.bed.1kb.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "NBL"
base <- tolower(BASE)

wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")

samples.rna <- readTable(file.path(wd.rna, "nbl_rna_n54.list"), header=F, rownames=2, sep="\t")
samples.wgs <- readTable(file.path(wd.wgs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="\t")
overlaps <- intersect(samples.rna$V2, samples.wgs$SAMPLE_ID)
samples  <- cbind(samples.rna[overlaps,], samples.wgs[overlaps,])
samples$M2 <- as.factor(samples$M2)
rownames(samples) <- samples$V1

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
#load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene <- tpm.gene[, rownames(samples)]   ## VERY VERY VERY IMPORTANT!!!
tpm.gene.log2.m <- getLog2andMedian(tpm.gene, 0.01)
nrow(tpm.gene.log2.m)
# [1] 34908
# [1] 22793
# [1] 18764

tpm.gene <- tpm.gene[rownames(tpm.gene.log2.m),]
save(tpm.gene, file=file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
tpm.gene.log2 <- log2(tpm.gene + 0.01)

# -----------------------------------------------------------------------------
# Relationship between expression and in-silico sorting
# Last Modified: 10/01/20
# -----------------------------------------------------------------------------
colnames <- c("RHO", "P", "Q", "M1", "M2", "LOG2_FC")   ##"ANOVA_P", "ANOVA_Q", 
de <- toTable(0, length(colnames), nrow(tpm.gene.log2), colnames)
rownames(de) <- rownames(tpm.gene.log2)

## SRC
de$RHO <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[4]])
de$P   <- mapply(x = 1:nrow(tpm.gene.log2), function(x) cor.test(as.numeric(tpm.gene.log2[x,]), samples$COR, method="spearman", exact=F)[[3]])
de <- de[!is.na(de$P),]

## Log2 fold change
de$M1 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 0)))
de$M2 <- median00(tpm.gene.log2, rownames(subset(samples, M2 == 1)))
de$LOG2_FC <- de$M2 - de$M1

## FDR
library(qvalue)
de$Q <- qvalue(de$P)$qvalue
de <- de[order(de$P),]

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)   ## BE EXTRA CAREFUL!!

writeTable(de.tpm.gene, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_n53.txt"), colnames=T, rownames=F, sep="\t")
save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_nbl_tpm-gene_src_q_n53.RData"))
nrow(de.tpm.gene)
# [1] 22793

##
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
de.tpm.gene <- de.tpm.gene[intersect(rownames(de.tpm.gene), rownames(tpm.gene)),]
de.tpm.gene$Q <- qvalue(de.tpm.gene$P)$qvalue
de.tpm.gene <- de.tpm.gene[order(de.tpm.gene$P),]

writeTable(de.tpm.gene, file.path(wd.de.data, "de_nbl_tpm-gene-r5p47_src_q_n53.txt"), colnames=T, rownames=F, sep="\t")
save(de.tpm.gene, samples, file=file.path(wd.de.data, "de_nbl_tpm-gene-r5p47_src_q_n53.RData"))
nrow(de.tpm.gene)
# [1] ?

# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.nbl   ## MUY MUY IMPORTANTE!!
tpm.gene.log2.m <- tpm.gene.log2.m[rownames(de.tpm.gene),]
nrow(tpm.gene.log2.m)
# [1] 22793
# [1] ?

tpm.gene.log2.m.rfd <- getTRC(tpm.gene.log2.m, nrds.RT.NRFD)
nrow(tpm.gene.log2.m.rfd)
# [1] 21316
# [1] ?
nrow(subset(tpm.gene.log2.m.rfd, TRC == 0))
# [1] 2
# [1] ?

###
## TTR and CTR (IZ +TZ)
file.name <- file.path(wd.de.data, "tpm-gene-median0_log2_m_rfd.RData")
setTRC(tpm.gene.log2.m.rfd, rfd=0.9, file.name)
load(file.name)

length(which(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD < 0))
# [1] 111
# [1] 

length(which(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD > 0))
# [1] 125
# [1] 

## TTR
testU(tpm.gene.log2.m.rfd.ttr$MEDIAN, tpm.gene.log2.m.rfd.ctr$MEDIAN)
# [1] 9.962099e-07
testU(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 3.994546e-05
testU(tpm.gene.log2.m.rfd.ttr$MEDIAN, tpm.gene.log2.m.rfd.ctr.iz$MEDIAN)
# [1] 1.961802e-10
testU(tpm.gene.log2.m.rfd.ttr$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 0.6925411

## TTR (r5p47)
testU(tpm.gene.log2.m.rfd.ttr$MEDIAN, tpm.gene.log2.m.rfd.ctr$MEDIAN)
# [1] 
testU(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 
testU(tpm.gene.log2.m.rfd.ttr$MEDIAN, tpm.gene.log2.m.rfd.ctr.iz$MEDIAN)
# [1] 
testU(tpm.gene.log2.m.rfd.ttr$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 

# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd))
de.tpm.gene.rfd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd$Q <- qvalue(de.tpm.gene.rfd$P)$qvalue
writeTable(de.tpm.gene.rfd, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_n53.txt"), colnames=T, rownames=F, sep="\t")

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz))
de.tpm.gene.rfd.iz <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz$Q <- qvalue(de.tpm.gene.rfd.iz$P)$qvalue
writeTable(de.tpm.gene.rfd.iz, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.cd))
de.tpm.gene.rfd.iz.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.cd$Q <- qvalue(de.tpm.gene.rfd.iz.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.cd, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_cd_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.ho))
de.tpm.gene.rfd.iz.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.ho$Q <- qvalue(de.tpm.gene.rfd.iz.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.ho, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_iz_ho_n53.txt"), colnames=T, rownames=F, sep="\t")

## TZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz))
de.tpm.gene.rfd.tz <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz$Q <- qvalue(de.tpm.gene.rfd.tz$P)$qvalue
writeTable(de.tpm.gene.rfd.tz, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_tz_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.cd))
de.tpm.gene.rfd.tz.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.cd$Q <- qvalue(de.tpm.gene.rfd.tz.cd$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.cd, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_tz_cd_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.ho))
de.tpm.gene.rfd.tz.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.ho$Q <- qvalue(de.tpm.gene.rfd.tz.ho$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.ho, file.path(wd.de.data, "de_nbl_tpm-gene-median0_src_q_rfd_tz_ho_n53.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# RFD vs. TPM (Figures)
# Last Modified: 06/01/20
# -----------------------------------------------------------------------------
ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 14.5)

file.name <- paste0("boxplot_nbl_tpm.gene.median0_rfd_TTR+IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, main="NBL expression", names=c("TTR", "IZ"), cols=c("black", "red"), ylim)

file.name <- paste0("boxplot_nbl_tpm.gene.median0_rfd_TTR+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.tz, main="NBL expression", names=c("TTR", "TZ"), cols=c("black", "blue"), ylim)

file.name <- paste0("boxplot_nbl_tpm.gene.median0_rfd_IZ+TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="NBL expression", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

##
file.name <- paste0("boxplot_nbl_tpm.gene.median0_rfd_HO+CD_IZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.iz.cd, main="NBL IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

file.name <- paste0("boxplot_nbl_tpm.gene.median0_rfd_HO+CD_TZ")
plotBox0(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.tz.cd, main="NBL TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)











# -----------------------------------------------------------------------------
# Figures
# Last Modified: 06/01/20
# -----------------------------------------------------------------------------
plotBox3 <- function(wd.de.plots, file.name, tpm.1, tpm.2, tpm.3, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- c(trait, rep(2, nrow(tpm.3)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN, tpm.3$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=5, width=3.2)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main, boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
   text(2, 15, "****", col="black", cex=2.5)
 
   ##
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   axis(side=1, at=seq(1, 3, by=1), labels=names, font=2, cex.axis=1.2)
   axis(side=1, at=2, labels="Total n=30,078", line=1.3, col=NA, cex.axis=1.2)
 
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

ylim <- c(min(tpm.gene.log2.m.rfd$MEDIAN), 15)
file.name <- paste0("boxplot_nbl_tpm.gene.r5p47.rfd_TTR+IZ+TZ_3.2*5total")
plotBox3(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ttr, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="NBL expression", names=c("TTR", "IZ", "TZ"), cols=c("black", "red", "blue"), ylim)

plotBox <- function(wd.de.plots, file.name, tpm.1, tpm.2, main, names, cols, ylim) {
   trait <- rep(0, nrow(tpm.1))
   trait <- c(trait, rep(1, nrow(tpm.2)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(tpm.1$MEDIAN, tpm.2$MEDIAN))
 
   pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3)
   boxplot(expr ~ trait, outline=T, xaxt="n", ylab=paste0("log2(TPM + 0.01)"), main=main,boxcol=cols, whiskcol=cols, outcol=cols, medcol=cols, staplecol=cols, ylim=ylim, cex.axis=1.1, cex.lab=1.2, cex.main=1.25)
 
   axis(side=1, at=seq(1, 2, by=1), labels=names, font=2, cex.axis=1.2)
   axis(side=1, at=1, labels="n=1,904", line=1.3, col=NA, cex.axis=1.2)
   axis(side=1, at=2, labels="n=1,203", line=1.3, col=NA, cex.axis=1.2)
 
   text(2, 14, "*", col="black", cex=2.5)
   
   mtext(paste0("p-value = ", scientific(wilcox.test(expr ~ trait, exact=F)$p.value)), cex=1.25, line=0.3)
   dev.off()
}

###
##
file.name <- paste0("boxplot_nbl_tpm.gene.r5p47.rfd_IZ-NS-vs-IZ-S_*")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.nbl.ns, tpm.gene.log2.m.rfd.ctr.iz.nbl.s, main="NBL, SCLC-shared IZ (IZ-S)", names=c("IZ-NS", "IZ-S"), cols=c("red", "red"), ylim)

genes <- c("TERT")
for (g in 1:length(genes)) {
   id <- subset(ensGene, external_gene_name == genes[g])$ensembl_gene_id
   plotCYS(genes[g], as.numeric(tpm.gene.log2[id,]), samples$COR, 1, "black", "bottomright")
}

# -----------------------------------------------------------------------------
# 
# Last Modified: 01/21/20
# -----------------------------------------------------------------------------
tpm.gene.log2.m.rfd.ctr.iz.nbl <- tpm.gene.log2.m.rfd.ctr.iz
# > nrow(tpm.gene.log2.m.rfd.ctr.iz.nbl)
# [1] 3107
# > nrow(tpm.gene.log2.m.rfd.ctr.iz)
# [1] 2565

overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.nbl), rownames(tpm.gene.log2.m.rfd.ctr.iz.sclc))
diffs <- setdiff(rownames(tpm.gene.log2.m.rfd.ctr.iz.nbl), overlaps)
tpm.gene.log2.m.rfd.ctr.iz.nbl.s <- tpm.gene.log2.m.rfd.ctr.iz.nbl[diffs,]

tpm.gene.log2.m.rfd.ctr.iz.nbl.s <- tpm.gene.log2.m.rfd.ctr.iz.nbl.s[intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.nbl.s), rownames(de.tpm.gene)),]
tpm.gene.log2.m.rfd.ctr.iz.nbl.ns <- tpm.gene.log2.m.rfd.ctr.iz.nbl[overlaps,]
tpm.gene.log2.m.rfd.ctr.iz.nbl.ns <- tpm.gene.log2.m.rfd.ctr.iz.nbl.ns[intersect(rownames(tpm.gene.log2.m.rfd.ctr.iz.nbl.ns), rownames(de.tpm.gene)),]

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.nbl.s))
de.tpm.gene.rfd.iz.s <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.s$Q <- qvalue(de.tpm.gene.rfd.iz.s$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.s, file.path(wd.de.data, "de_nbl_tpm-gene-r5p47-rfd_src_q_iz_s_n39.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.iz.nbl.ns))
de.tpm.gene.rfd.iz.ns <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.iz.ns$Q <- qvalue(de.tpm.gene.rfd.iz.ns$P)$qvalue
writeTable(de.tpm.gene.rfd.iz.ns, file.path(wd.de.data, "de_nbl_tpm-gene-r5p47-rfd_src_q_iz_ns_n39.txt"), colnames=T, rownames=F, sep="\t")

##
xlab.text <- "NBL M2/M1 [log2FC]"

## IZ (S)
plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_r5p47_rfd_fdr0.09")
genes <- readTable(paste0(plot.de, "_iz_s.tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL and SCLC-shared initiated genes (IZ-S)", paste0("(n=1,203)"))
file.de <- paste0(plot.de, "_iz_s.pdf")
plotVolcano(de.tpm.gene.rfd.iz.s, 0.09, genes, file.de, file.main, xlab.text)

## IZ (S; testU)
de.tpm.gene.rfd.iz.s$RHO <- NA
de.tpm.gene.rfd.iz.s$P <- de.tpm.gene[rownames(de.tpm.gene.rfd.iz.s),]$P
de.tpm.gene.rfd.iz.s$Q <- qvalue(de.tpm.gene.rfd.iz.s$P)$qvalue

de.tpm.gene.rfd.iz.s <- de.tpm.gene.rfd.iz.s[order(de.tpm.gene.rfd.iz.s$P),]
writeTable(de.tpm.gene.rfd.iz.s, file.path(wd.de.data, "de_nbl_tpm-gene-r5p47-rfd_U_q_iz_s_n53.txt"), colnames=T, rownames=F, sep="\t")

plot.de <- file.path(wd.de.plots, "volcanoplot_nbl_r5p47_rfd_fdr0.09")
genes <- readTable(paste0(plot.de, "_iz_s.tab"), header=T, rownames=F, sep="\t")
file.main <- c("NBL and SCLC-shared initiated genes (IZ-S)", paste0("(n=1,203)"))
file.de <- paste0(plot.de, "_iz_s_RECQL4.pdf")
plotVolcano(de.tpm.gene.rfd.iz.s, 0.05, genes, file.de, file.main, xlab.text)






###
##
overlaps <- intersect(rownames(tpm.gene.log2.m.rfd.ctr.tz), rownames(tpm.gene.log2.m.rfd.ctr.tz.nl))
tpm.gene.log2.m.rfd.ctr.tz.sp <- tpm.gene.log2.m.rfd.ctr.tz[setdiff(rownames(tpm.gene.log2.m.rfd.ctr.tz), overlaps),]

tpm.gene.log2.m.rfd.ctr.tz.sp <- tpm.gene.log2.m.rfd.ctr.tz.sp[intersect(rownames(tpm.gene.log2.m.rfd.ctr.tz.sp), rownames(de.tpm.gene)),]

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.m.rfd.ctr.tz.sp))
de.tpm.gene.rfd.tz.sp <- de.tpm.gene[overlaps,]
de.tpm.gene.rfd.tz.sp$Q <- qvalue(de.tpm.gene.rfd.tz.sp$P)$qvalue
writeTable(de.tpm.gene.rfd.tz.sp, file.path(wd.de.data, "de_sclc_tpm-gene-r5p47-rfd_src_q_tz_sp_n70.txt"), colnames=T, rownames=F, sep="\t")


###
##
file.name <- paste0("boxplot_sclc_tpm.gene.r5p47.rfd_IZ-NSP-vs-IZ-SP_3")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.nsp, tpm.gene.log2.m.rfd.ctr.iz.sp, main="SCLC-specific IZ (IZ-S)", names=c("IZ-NS", "IZ-S"), cols=c("red", "red"), ylim)


###
## IZ vs TZ
# > median(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN)
# [1] 3.767014
# > median(tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 3.57522
# > testU(tpm.gene.log2.m.rfd.ctr.iz$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz$MEDIAN)
# [1] 0.02015283
file.name <- paste0("boxplot_nbl_tpm.gene.r5p47.rfd_IZ-vs-TZ")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz, tpm.gene.log2.m.rfd.ctr.tz, main="CTR", names=c("IZ", "TZ"), cols=c("red", "blue"), ylim)

## IZ (CD vs. HO)
# > median(tpm.gene.log2.m.rfd.ctr.iz.cd$MEDIAN)
# [1] 3.792572
# > median(tpm.gene.log2.m.rfd.ctr.iz.ho$MEDIAN)
# [1] 3.74394
# > testU(tpm.gene.log2.m.rfd.ctr.iz.cd$MEDIAN, tpm.gene.log2.m.rfd.ctr.iz.ho$MEDIAN)
# [1] 0.9836556
file.name <- paste0("boxplot_nbl_tpm.gene.r5p47.rfd_IZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.iz.ho, tpm.gene.log2.m.rfd.ctr.iz.cd, main="IZ", names=c("HO", "CD"), cols=c("red", "red"), ylim)

## TZ (CD vs. HO)
# > median(tpm.gene.log2.m.rfd.ctr.tz.cd$MEDIAN)
# [1] 3.502108
# > median(tpm.gene.log2.m.rfd.ctr.tz.ho$MEDIAN)
# [1] 3.661897
# > testU(tpm.gene.log2.m.rfd.ctr.tz.cd$MEDIAN, tpm.gene.log2.m.rfd.ctr.tz.ho$MEDIAN)
# [1] 0.5069909
file.name <- paste0("boxplot_nbl_tpm.gene.r5p47.rfd_TZ_HO-vs-CD")
plotBox(wd.de.plots, file.name, tpm.gene.log2.m.rfd.ctr.tz.ho, tpm.gene.log2.m.rfd.ctr.tz.cd, main="TZ", names=c("HO", "CD"), cols=c("blue", "blue"), ylim)

# -----------------------------------------------------------------------------
# Gene length vs. RFD slopes
# Last Modified: 10/01/20
# -----------------------------------------------------------------------------
plotTRC <- function(rfd, tpm, main.text, file.name, xlim, ylim, col, pos, ylab.text, isFlip=F) {
   xlab.text <- "Gene length [log10]"
 
   pdf(paste0(file.name, ".pdf"), height=6, width=6)
   plot(rfd ~ tpm, ylim=ylim, xlim=xlim, ylab="", xlab=xlab.text, main="", col=col[1], pch=1, cex=1.2, cex.axis=1.5, cex.lab=1.55, yaxt="n")
 
   lm.fit <- lm(rfd ~ tpm)
   abline(lm.fit, col=col[2], lwd=3)
 
   cor <- cor.test(rfd, tpm, method="spearman", exact=F)
   legend(pos, c(paste0("rho = ", round0(cor[[4]], digits=2)), paste0("p-value = ", scientific(cor[[3]], digits=2))), text.col=col[2], bty="n", cex=1.55)
 
   if (isFlip)
      axis(side=2, at=seq(0, 0.1, by=0.05), labels=c(0, -0.05, -0.1), cex.axis=1.5)
   else
      axis(side=2, at=seq(0, 0.1, by=0.05), labels=c(0, 0.05, 0.1), cex.axis=1.5)
   mtext(ylab.text, side=2, line=2.85, cex=1.5)
   mtext(main.text[1], line=1.8, cex=1.55, font=2)
   mtext(main.text[2], line=0.3, cex=1.55)
   dev.off()
}

###
##
tpm.gene.log2.m.rfd.ctr.iz <- cbind(tpm.gene.log2.m.rfd.ctr.iz, de.tpm.gene[rownames(tpm.gene.log2.m.rfd.ctr.iz),])
tpm.gene.log2.m.rfd.ctr.iz$length <- abs(tpm.gene.log2.m.rfd.ctr.iz$end_position - tpm.gene.log2.m.rfd.ctr.iz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.iz.cd <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.iz.ho <- subset(tpm.gene.log2.m.rfd.ctr.iz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.iz$length)))
ylim <- c(min(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD) * -1)
ylab.text <- "Initiation efficiency"

file.name <- file.path(wd.de.plots, "TRC3_LENGTH-vs-GENE-NRFD_IZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz$length), c("NBL initiation zone (IZ) genes", "(n=3,107)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "TRC3_LENGTH-vs-GENE-NRFD_IZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.cd$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.cd$length), c("NBL IZ (CD) genes", "(n=1,511)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

file.name <- file.path(wd.de.plots, "TRC3_LENGTH-vs-GENE-NRFD_IZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.iz.ho$GENE_NRFD, log10(tpm.gene.log2.m.rfd.ctr.iz.ho$length), c("NBL IZ (HO) genes", "(n=1,594)"), file.name, xlim, ylim, col=c("red", "black"), "topright", ylab.text)

##
tpm.gene.log2.m.rfd.ctr.tz <- cbind(tpm.gene.log2.m.rfd.ctr.tz, de.tpm.gene[rownames(tpm.gene.log2.m.rfd.ctr.tz),])
tpm.gene.log2.m.rfd.ctr.tz$length <- abs(tpm.gene.log2.m.rfd.ctr.tz$end_position - tpm.gene.log2.m.rfd.ctr.tz$start_position) / 1000
tpm.gene.log2.m.rfd.ctr.tz.cd <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC > 0)
tpm.gene.log2.m.rfd.ctr.tz.ho <- subset(tpm.gene.log2.m.rfd.ctr.tz, TRC < 0)

##
xlim <- c(min(log10(tpm.gene.log2.m.rfd.ctr.iz$length)), max(log10(tpm.gene.log2.m.rfd.ctr.tz$length)))
ylim <- c(max(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD), min(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD)) * -1
ylab.text <- "Termination efficiency"

file.name <- file.path(wd.de.plots, "TRC3_LENGTH-vs-GENE-NRFD_TZ")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz$length), c("NBL termination zone (TZ) genes", "(n=1,952)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "TRC3_LENGTH-vs-GENE-NRFD_TZ_CD")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.cd$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.cd$length), c("NBL TZ (CD) genes", "(n=963)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)

file.name <- file.path(wd.de.plots, "TRC3_LENGTH-vs-GENE-NRFD_TZ_HO")
plotTRC(tpm.gene.log2.m.rfd.ctr.tz.ho$GENE_NRFD * -1, log10(tpm.gene.log2.m.rfd.ctr.tz.ho$length), c("NBL TZ (HO) genes", "(n=989)"), file.name, xlim, ylim, col=c("blue", "black"), "topright", ylab.text, isFlip=T)
















# -----------------------------------------------------------------------------
# RFD vs. TPM
# Last Modified: 31/10/18
# -----------------------------------------------------------------------------
nrds.RT.NRFD <- nrds.RT.NRFD.nbl   ## MUY MUY IMPORTANTE!!

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.log2 <- getLog2andMedian(tpm.gene, 0.01)
# > nrow(tpm.gene.log2)
# [1] 18764

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
tpm.gene.log2 <- cbind(annot[rownames(tpm.gene.log2),], tpm.gene.log2[, "MEDIAN"], ensGene.bed[rownames(tpm.gene.log2),]) 
colnames(tpm.gene.log2)[8] <- "MEDIAN"

tpm.gene.log2$TSS_RFD <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$RFD
#tpm.gene.log2$TTS_RFD <- nrds.RT.NRFD[tpm.gene.log2$TTS,]$RFD
tpm.gene.log2$TSS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$SPLINE
#tpm.gene.log2$TTS_SPLINE <- nrds.RT.NRFD[tpm.gene.log2$TTS,]$SPLINE

tpm.gene.log2 <- tpm.gene.log2[!is.na(tpm.gene.log2$TSS_RFD),]
# > nrow(tpm.gene.log2)
# [1] 17773

##
tpm.gene.log2$TRC <- tpm.gene.log2$strand * tpm.gene.log2$TSS_RFD
tpm.gene.log2$TRC[which(tpm.gene.log2$TRC > 0)] <- 1
tpm.gene.log2$TRC[which(tpm.gene.log2$TRC < 0)] <- -1

## CTR (IZ)
tpm.gene.log2$TSS_NRFD <- nrds.RT.NRFD[tpm.gene.log2$TSS,]$NRFD
tpm.gene.log2.ctr <- subset(subset(tpm.gene.log2, TSS_RFD < 0.9), TSS_RFD > -0.9)
tpm.gene.log2.ctr.iz <- subset(tpm.gene.log2.ctr, TSS_NRFD > 0)
tpm.gene.log2.ctr.iz.cd <- subset(tpm.gene.log2.ctr.iz, TRC > 0)
tpm.gene.log2.ctr.iz.ho <- subset(tpm.gene.log2.ctr.iz, TRC < 0)

tpm.gene.log2.ctr.tz <- subset(tpm.gene.log2.ctr, TSS_NRFD < 0)
tpm.gene.log2.ctr.tz.cd <- subset(tpm.gene.log2.ctr.tz, TRC > 0)
tpm.gene.log2.ctr.tz.ho <- subset(tpm.gene.log2.ctr.tz, TRC < 0)

## TTR
tpm.gene.log2.ttr <- rbind(subset(tpm.gene.log2, TSS_RFD >= 0.9), subset(tpm.gene.log2, TSS_RFD <= -0.9))
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr$MEDIAN)
# [1] 4.887171e-05
testU(tpm.gene.log2.ctr.iz$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.02468772
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.iz$MEDIAN)
# [1] 4.337057e-06
testU(tpm.gene.log2.ttr$MEDIAN, tpm.gene.log2.ctr.tz$MEDIAN)
# [1] 0.2477697

## ALL
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2))
de.tpm.gene <- de.tpm.gene[overlaps,]
de.tpm.gene$Q <- qvalue(de.tpm.gene$P)$qvalue
#de.tpm.gene$FDR <- p.adjust(de.tpm.gene$P, method="BH", n=length(de.tpm.gene$P))
writeTable(de.tpm.gene, file.path(wd.de.data, "de_nbl_tpm.gene.r5p47-vs-wgs_src_q_n53.txt"), colnames=T, rownames=F, sep="\t")

## IZ
overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.iz))
de.tpm.gene.iz <- de.tpm.gene[overlaps,]
de.tpm.gene.iz$Q <- qvalue(de.tpm.gene.iz$P)$qvalue
writeTable(de.tpm.gene.iz, file.path(wd.de.data, "de_nbl_tpm.gene.r5p47-vs-wgs_src_q_iz_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.iz.ho))
de.tpm.gene.iz.ho <- de.tpm.gene[overlaps,]
de.tpm.gene.iz.ho$Q <- qvalue(de.tpm.gene.iz.ho$P)$qvalue
writeTable(de.tpm.gene.iz.ho, file.path(wd.de.data, "de_nbl_tpm.gene.r5p47-vs-wgs_src_q_iz_ho_n53.txt"), colnames=T, rownames=F, sep="\t")

overlaps <- intersect(rownames(de.tpm.gene), rownames(tpm.gene.log2.ctr.iz.cd))
de.tpm.gene.iz.cd <- de.tpm.gene[overlaps,]
de.tpm.gene.iz.cd$Q <- qvalue(de.tpm.gene.iz.cd$P)$qvalue
writeTable(de.tpm.gene.iz.cd, file.path(wd.de.data, "de_nbl_tpm.gene.r5p47-vs-wgs_src_q_iz_cd_n53_.txt"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# 
# Last Modified: 13/01/20
# -----------------------------------------------------------------------------
## PIF1
testU(tpm.gene.log2.sclc["ENSG00000140451",], tpm.gene.log2.nbl["ENSG00000140451",])
# [1] 0.4761733

## TOR1AIP1
testU(tpm.gene.log2.sclc["ENSG00000143337",], tpm.gene.log2.nbl["ENSG00000143337",])
# [1] 3.98806e-14

## BRCA2
testU(tpm.gene.log2.sclc["ENSG00000139618",], tpm.gene.log2.nbl["ENSG00000139618",])
# [1] 0.4270861









# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; 25 M2 vs 28 M1)
# Last Modified: 11/10/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : M2 as factor
argv      <- data.frame(predictor="M2", predictor.wt=0, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47_m2-m1_wilcox_q_n53")
file.main <- paste0("RT (n=25) vs WT (n=28) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; 14 Q4 vs 14 Q1)
# Last Modified: 11/10/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : Q4 as factor
samples <- subset(samples, Q4 %in% c(4,1))
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
# > dim(tpm.gene.log2)
# [1] 18764    28

argv      <- data.frame(predictor="Q4", predictor.wt=1, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47_q4-q1_wilcox_q_n25")
file.main <- paste0("Q4 (n=14) vs Q1 (n=14) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric; 11 Q3 vs 14 Q1)
# Last Modified: 11/10/19
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
##       Student's/t.test
## FDR : Q/BH
## DE  : Q4 as factor
samples <- subset(samples, Q4 %in% c(3,1))
tpm.gene.log2 <- tpm.gene.log2[, rownames(samples)]
# > dim(tpm.gene.log2)
# [1] 18764    25

argv      <- data.frame(predictor="Q4", predictor.wt=1, test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-r5p47_q3-q1_wilcox_q_n25")
file.main <- paste0("Q3 (n=11) vs Q1 (n=14) in ", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")












# -----------------------------------------------------------------------------
# Volcano plots of RB1-loss DE genes in LCNEC
# Figure(s)    : Figure 1 (A)
# Last Modified: 07/01/19
# -----------------------------------------------------------------------------
fdrToP <- function(fdr, de) {
   de.sig <- subset(de, FDR <= fdr)
   de.sig$log10P <- -log10(de.sig$P)
 
   return(max(de.sig$P))
}

plotVolcano <- function(de, fdr, genes, file.de, file.main) {
   pvalue <- fdrToP(fdr, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
 
   de$log10P <- -log10(de$P)
   xmax <- max(de$LOG2_FC)
   ymax <- max(de$log10P)
 
   pdf(file.de, height=7, width=7)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="CLL T20/T27 log2 fold change", ylab="-log10(p-value)", col="darkgray", main=file.main[1])

   abline(h=c(-log10(pvalue)), lty=5)
   text(xmax*-1 + 2*xmax/35, -log10(pvalue) + ymax/42, "FDR=0.02", cex=0.85)
   #abline(h=c(-log10(fdrToP(0.1, de))), lty=5, col="darkgray")
   #text(xmax*-1 + 2*xmax/50, -log10(fdrToP(0.1, de)) + ymax/42, "FDR=0.1", col="darkgray", cex=0.85)

   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="red")
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col="dodgerblue")
 
   for (g in 1:nrow(genes)) {
      gene <- subset(de, external_gene_name == genes[g,]$GENE)
      gene <- cbind(gene, genes[g,])
      
      if (nrow(gene) > 0) {
         points(gene$LOG2_FC, gene$log10P, pch=1, col="black")
         
         if (!is.na(gene$ADJ_1))
            if (is.na(gene$ADJ_2))
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=gene$ADJ_1, cex=0.75)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=0.75)
         else
            if (gene$LOG2_FC > 0)
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(0, -0.5), cex=0.75)
            else
               text(gene$LOG2_FC, gene$log10P, genes[g,]$GENE, col="black", adj=c(1, -0.5), cex=0.75)
      } else
         print(genes[g])
   }
   
   mtext(paste0("(", file.main[2], ")"), cex=1.2, line=0.3)
   legend("topright", legend=c("Upregulated", "Downregulated"), col=c("red", "dodgerblue"), pch=19)
   dev.off()
}

##
plot.main <- "Differential expression between CLL T29 and T33"
plot.de <- file.path(wd.de.plots, "volcanoplot-r50p100_cll_rt_q0.02")

## Natural killer cell receptor phenotypes
genes <- readTable(paste0(plot.de, "_nk.tab"), header=T, rownames=F, sep="\t")
file.main <- c(plot.main, "Natural killer cell receptors")
file.de <- paste0(plot.de, "_nk.pdf")
plotVolcano(de.tpm.gene, 0.02, genes, file.de, file.main)

# -----------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) on LCNEC RB1/WT ranked gene lists
# Figure(s)    : Figure S1 (A and B)
# Last Modified: 08/01/19
# -----------------------------------------------------------------------------
file.name <- paste0("de_nbl_tpm-gene-r5p47_rt_wilcox_q_n53")
writeRNKformat(de.tpm.gene, wd.de.gsea, file.name)

## Tirosh et al 2016 
load(file.path(wd.src.ref, "cycle.RData"))                           ## See guide-to-the/cycle.R
core.G1S <- intersect(core.G1S, rownames(tpm.gene.log2))             ## 41/43/43 are expressed in the dataset
core.G2M <- intersect(core.G2M, rownames(tpm.gene.log2))             ## 54/54/55
core.Stemness <- intersect(core.Stemness, rownames(tpm.gene.log2))   ## 54/58/63

writeGRPformat(core.G1S, wd.de.gsea, "core.G1-S")
writeGRPformat(core.G2M, wd.de.gsea, "core.G2-M")
writeGRPformat(core.Stemness, wd.de.gsea, "core.Stemness")

## Dominguez et al 2016
periodic.G1S <- intersect(periodic.G1S, rownames(tpm.gene.log2))   ## 262/304
periodic.G2M <- intersect(periodic.G2M, rownames(tpm.gene.log2))   ## 792/876

writeGRPformat(periodic.G1S, wd.de.gsea, "G1-S")
writeGRPformat(periodic.G2M, wd.de.gsea, "G2-M")










# -----------------------------------------------------------------------------
# Principal component analysis (PCA) of LCNEC samples on RB1-loss DE genes
# Figure(s)    : Figure 1 (B)
# Last Modified: 23/11/17
# -----------------------------------------------------------------------------
#load(file.path("/Users/tpyang/Work/uni-koeln/tyang2", "LCNEC", "analysis/expression/kallisto", paste0("lcnec", "-tpm-de/data/", "de_", "lcnec", "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))
genes.rt.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
# > length(genes.rb1.q0.05)
# [1] 145

##
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))

trait <- as.numeric(samples[,"RB1_MUT"])
trait[which(trait == 0)] <- "LCNEC (WT)"
trait[which(trait == 1)] <- "LCNEC (RB1)"
trait[which(is.na(trait))] <- "LCNEC (NA)"

file.main <- c("LCNEC samples on top 145 DE genes", "RB1-loss differential effect; FDR < 0.05")
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcnec_rb1_q0.05_145de", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NULL, flip.x=1, flip.y=-1)


# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
nrow(ensGene[intersect(genes.rb1.q0.05, core.G1S),])   ## 10/41 = 23.4%
nrow(ensGene[intersect(genes.rb1.q0.05, core.G2M),])   ##  4/54 = 7.4%
nrow(ensGene[intersect(genes.rb1.q0.05, core.SC),])    ##  0/54 = 0%

nrow(ensGene[intersect(genes.rb1.q0.05, genes.G1S),])   ## 15
nrow(ensGene[intersect(genes.rb1.q0.05, genes.G2M),])   ##  9

nrow(ensGene[intersect(genes.rb1.q0.05, unique(c(core.G1S, genes.G1S))),])   ## 19
nrow(ensGene[intersect(genes.rb1.q0.05, unique(c(core.G2M, genes.G2M))),])   ## 11

## Output
genes.rb1.n.s   <- rownames(subset(de.tpm.gene, FDR > 0.5))
genes.rb1.q0.5  <- rownames(subset(de.tpm.gene, FDR <= 0.5))
file.name <- paste0("de_", base, "_tpm-gene-r5p47_rb1_wilcox_q_n54_n.s.0.5")
#genes.rb1 <- genes.rb1.n.s
#genes.rb1 <- genes.rb1.q0.5
genes.rb1 <- genes.rb1.q0.1

output <- de.tpm.gene[intersect(genes.rb1, unique(c(core.G1S, genes.G1S))), c(1,2,8,9,12)]
output$Tirosh_2016 <- ""
output$Dominguez_2016 <- ""
output[intersect(genes.rb1, core.G1S),]$Tirosh_2016 <- "yes"
output[intersect(genes.rb1, genes.G1S),]$Dominguez_2016 <- "yes"
writeTable(output, file.path(wd.de.data, paste0(file.name, "_G1-S.txt")), colnames=T, rownames=F, sep="\t")

output <- de.tpm.gene[intersect(genes.rb1, unique(c(core.G2M, genes.G2M))), c(1,2,8,9,12)]
output$Tirosh_2016 <- ""
output$Dominguez_2016 <- ""
output[intersect(genes.rb1, core.G2M),]$Tirosh_2016 <- "yes"
output[intersect(genes.rb1, genes.G2M),]$Dominguez_2016 <- "yes"
writeTable(output, file.path(wd.de.data, paste0(file.name, "_G2-M.txt")), colnames=T, rownames=F, sep="\t")
