# =============================================================================
# Manuscript   : The dangeous case of DNA replication in SCLC
# Chapter II   : Correlating replicaiton-fork stress with transcriptional strand asymmetry
# Name         : manuscripts/asymmetries/sclc-asym-tx-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 30/01/18
# =============================================================================
#wd.source <- "/projects/cangen/tyang2/dev/R"              ## tyang2@cheops
#wd.source <- "/ngs/cangen/tyang2/dev/R"                   ## tyang2@gauss
wd.source <- "/Users/tpyang/Work/dev/R"                    ## tpyang@localhost

handbooks <- c("Common.R", "Asymmetry.R", "ReplicationTiming.R", "Rearrangement.R")   ## Required handbooks/libraries for the manuscript
invisible(sapply(handbooks, function(b) source(file.path(wd.source, "handbook-of", b))))

load(file.path(wd.source, "guide-to-the", "hg19.RData"))   ## The bioinformatician's guide to the reference genome
load(file.path(wd.source, "guide-to-the", "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd     <- "/ngs/cangen/tyang2/SCLC/analysis/"                ## tyang2@gauss
#wd.ngs <- "/ngs/cangen/tyang2/SCLC/ngs/WGS/" 
wd     <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/"   ## tpyang@localhost
wd.ngs <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/ngs/WGS/"

wd.rt       <- paste0(wd, "replication/sclc-wgs-rt/")
wd.rt.data  <- paste0(wd.rt, "data/")
wd.rt.plots <- paste0(wd.rt, "plots/")
wd.asym       <- paste0(wd, "asymmetries/sclc-wgs-asym/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 5: Define replicaiton timing direction for expressed genes (Following Step 4 in "asym-sclc-tx.R" and Step 5 from rt-sclc-wgs.R)
# Name: 
# Link(s): http://www.mun.ca/biology/scarr/2250_DNA_replication_&_transcription.html
#          https://stackoverflow.com/questions/43615469/how-to-calculate-the-slope-of-a-smoothed-curve-in-r
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
plotRT0 <- function(wd.rt.plots, title, chr, n, xmin, xmax, rpkms.chr.rt, bed.gc.chr, pair1, pair0, ext) {
   filename  <- paste0(wd.rt.plots, tolower(title), "_wgs_rt_", chr, "_", pair1, "-", pair0, "_n", n)
   main.text <- paste0("Read depth (CN-, GC-corrected RPKM) ratio (", pair1, "/", pair0, ") in ", title)
   xlab.text <- paste0("Chromosome ", gsub("chr", "", chr), " coordinate (Mb)")
   ylab.text <- "Replication time (log2 FC)"
 
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   ymin <- min(rpkms.chr.rt$MEDIAN)
   ymax <- max(rpkms.chr.rt$MEDIAN)
   if (!is.na(xmin) && !is.na(xmax)) filename <- paste0(filename, "_", xmin/1E6, "-", xmax/1E6, "Mb")
   if (is.na(xmin)) xmin <- 0
   if (is.na(xmax)) xmax <- subset(chromInfo, chrom == chr)$size

   if (ext == "pdf") {
      pdf(paste0(filename, ".pdf"), height=4, width=10)
   } else if (ext == "png")
      png(paste0(filename, ".png"), height=4, width=10, units="in", res=300)   ## ADD 16/05/17: res=300
   
   plot(NULL, ylim=c(ymin, ymax), xlim=c(xmin/1E6, xmax/1E6), xlab=xlab.text, ylab=ylab.text, main=main.text)
   points(bed.gc.chr$START/1E6, rpkms.chr.rt$MEDIAN, col="red", cex=0.3)
   abline(h=0, lwd=0.5, col="grey")
   lines(bed.gc.chr$START/1E6, smooth.spline(rpkms.chr.rt$MEDIAN)$y)
 
   #slopes <- diff(smooth.spline(rpkms.chr)$y)/diff((bed.gc.chr$START)/1E6)
   #slopes2 <- diff(smooth.spline(rpkms.chr)$y)/diff(smooth.spline(rpkms.chr)$x)
   
   #temp <- loess.smooth(bed.gc.chr$START, rpkms.chr)
   #slopes3 <- diff(temp$y)/diff(temp$x)
   
   for (c in 1:nrow(cytoBand.chr))
      abline(v=cytoBand.chr$chromEnd[c]/1E6, lty=5, lwd=0.4, col="lightgrey")
 
   dev.off()
}

getEnsGeneBED <- function(pos, bed.gc.chr) {
   bed.gc.chr.start <- subset(bed.gc.chr, pos >= START)
   bed.gc.chr.start.end <- subset(bed.gc.chr.start, pos <= END)
 
   return(rownames(bed.gc.chr.start.end))
}

BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
PAIR  <- paste0(PAIR1, "-", PAIR0)
CHR   <- 2
CUTOFF <- 0.15

###
##
bed.gc <- bed[which(bed$GC > 0),]   ## Only keep partitions (in the BED file) with a GC content
ensGene.tx <- ensGene[rownames(tpm.gene.sclc),]

ensGene.tx.rt <- ensGene.tx[1,]
ensGene.tx.rt$SLOPE_START <- 0
ensGene.tx.rt$SLOPE_END <- 0
ensGene.tx.rt <- ensGene.tx.rt[-1,]
for (c in 1:22) {
   #chr <- chrs[CHR]
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)

   ## Replication timing
   rpkms.chr <- readTable(paste0(wd.rt.data, tolower(BASE), "_rpkm.corr.gc.d.rt_", chr, "_", PAIR, "_n", length(samples), ".txt.gz"), header=T, rownames=T, sep="\t") 
   
   ##
   rpkms.chr.rt <- rpkms.chr[which(rpkms.chr$MEDIAN > -CUTOFF),]
   rpkms.chr.rt <- rpkms.chr.rt[which(rpkms.chr.rt$MEDIAN < CUTOFF),]
   bed.gc.chr <- bed.gc.chr[rownames(rpkms.chr.rt),]
   
   plotRT0(wd.rt.plots, BASE, chr, length(samples), NA, NA, rpkms.chr.rt, bed.gc.chr, PAIR1, PAIR0, "png")
   #plotRT0(wd.rt.plots, BASE, chr, length(samples), 50000000, 100000000, rpkms.chr.rt$MEDIAN, bed.gc.chr,PAIR1, PAIR0, "png")
   
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
save(ensGene.tx.rt, file=paste0(wd.asym.data, "sclc_asym_tx_rt.RData"))

# -----------------------------------------------------------------------------
# Step 6.0: Read in all the indels
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
tx.indel <- getMergedReport(samples[1], readTable(paste0(wd.asym.data, samples[1], "_mut.ens.tx.indel.vcf.gz"), header=T, rownames=F, sep=""))[1,]
tx.indel <- tx.indel[-1,]
tx.del   <- tx.indel
tx.in    <- tx.indel
for (s in 1:length(samples)) {
   sample <- samples[s]
   
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.tx.indel.vcf.gz"), header=T, rownames=F, sep="")
   vcf <- subset(vcf, CHROM %in% paste0("chr", c(1:22)))   ## Keep only autosomes
   indel <- getMergedReport(sample, vcf)
   tx.indel <- rbind(tx.indel, indel)
   
   indel$SIZE <- nchar(indel$ALT) - nchar(indel$REF)
   tx.del <- rbind(tx.del, subset(indel, SIZE < 0)[,-8])
   tx.in  <- rbind(tx.in,  subset(indel, SIZE > 0)[,-8]) 
}
save(ensGene.tx.rt, tx.indel, tx.del, tx.in, file=paste0(wd.asym.data, "sclc_asym_tx_rt_indel_del_in.RData"))

## Keep only indels on protein coding genes which has at least a SNV
#ens.tx.pcg.indel <- intersect(unique(tx.indel$ensembl_gene_id), rownames(tpm.gene.sclc.pcg))
# > length(ens.tx.pcg.indel)
# [1] 14136
ens.tx.pcg.snv.indel <- intersect(unique(tx.indel$ensembl_gene_id), ens.tx.pcg.snv)   ## From Step 3.0; ADD 02/02/18
# > length(ens.tx.pcg.snv.indel)
# [1] 14037

tx.pcg.snv.indel <- subset(tx.indel, ensembl_gene_id %in% ens.tx.pcg.snv.indel)
# > nrow(tx.indel)
# [1] 263649
# > nrow(tx.pcg.snv.indel)
# [1] 258738









# -----------------------------------------------------------------------------
# Step 6.0: Read in all the SVs
# Last Modified: 31/01/18
# -----------------------------------------------------------------------------
sv.bp.tpm.cna <- NA
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   ## INPUT: INV
   #sv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt")
   sv.file <- paste0(wd.ngs, "rearrangements/", sample, "_unsnarl.txt")
   sv.bp <- initBreakPoint(read.peiflyne.unsnarl.txt(sv.file, cutoff=10))   ## TO-DO

   sv.inter <- subset(sv.bp, Type1 == "Inter")   ## ADD: 20/04/17
   sv.intra <- subset(sv.bp, Type1 == "Intra")   ## ADD: 20/04/17
   sv.intra$Size <- abs(sv.intra$Exact_Pos_Pair_1 - sv.intra$Exact_Pos_Pair_2)
   sv.intra <- subset(sv.intra, Size >= 10000)   ## TO-DO
   sv.intra <- sv.intra[,-28]                    ## TO-DO
   sv.bp <- rbind(sv.intra, sv.inter)
   
   ## INPUT: Expression
   sv.bp.tpm <- setBreakPointTranscriptionEnsGene(sample, sv.bp, tpm.gene.sclc.log2, ensGene.tx.rt, pair=1)
   sv.bp.tpm <- setBreakPointTranscriptionEnsGene(sample, sv.bp.tpm, tpm.gene.sclc.log2, ensGene.tx.rt, pair=2)
   
   ## INPUT: CNA
   cna.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_iCN.seg")
   sv.bp.tpm <- setBreakPointCNA(sample, sv.bp.tpm, read.peiflyne.iCN.seg(cna.file), pair=1)
   sv.bp.tpm <- setBreakPointCNA(sample, sv.bp.tpm, read.peiflyne.iCN.seg(cna.file), pair=2)
   
   ##
   if (is.na(sv.bp.tpm.cna))
      sv.bp.tpm.cna <- sv.bp.tpm
   else
      sv.bp.tpm.cna <- rbind(sv.bp.tpm.cna, sv.bp.tpm)
}   

break1 <- sv.bp.tpm.cna[,c("BP_ENS_ID_1", "BP_iCN_1")]
break1 <- subset(break1, !is.na(BP_ENS_ID_1))
colnames(break1) <- c("ensembl_gene_id", "iCN")
break2 <- sv.bp.tpm.cna[,c("BP_ENS_ID_2", "BP_iCN_2")]
break2 <- subset(break2, !is.na(BP_ENS_ID_2))
colnames(break2) <- c("ensembl_gene_id", "iCN")
tx.bp <- rbind(break1, break2)

# -----------------------------------------------------------------------------
# Step 6.1: Read in all the SVs
# Last Modified: 31/01/18
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Step 7:
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
getHeadOnCollision <- function(ensGene.tx.rt, headon) {
   if (headon)
      return(rbind(subset(subset(ensGene.tx.rt, SLOPE_START > 0), strand == -1), subset(subset(ensGene.tx.rt, SLOPE_START < 0), strand == 1)))
   else
      return(rbind(subset(subset(ensGene.tx.rt, SLOPE_START > 0), strand == 1), subset(subset(ensGene.tx.rt, SLOPE_START < 0), strand == -1)))
}

getTxQ4RT <- function(tx.q4, ensGene.tx.rt, headon) {
   tx.q4.rt <- list()
   tx.q4.rt[[4]] <- intersect(tx.q4[[4]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
   tx.q4.rt[[3]] <- intersect(tx.q4[[3]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
   tx.q4.rt[[2]] <- intersect(tx.q4[[2]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
   tx.q4.rt[[1]] <- intersect(tx.q4[[1]], rownames(getHeadOnCollision(ensGene.tx.rt, headon)))
 
   return(tx.q4.rt)
}

###
##
ensGene.tx.rt <- subset(subset(ensGene.tx.rt, !is.na(SLOPE_START)), !is.na(SLOPE_END))
tx.input <- tx.indel   ## tx.indel   ## tx.del   ## tx.in
headon   <- T          ## Replication–transcription Head-on (T) / Co-directional (F) collisions
ensGene.tx.rt.input <- ensGene.tx.rt[unique(tx.input$ensembl_gene_id),]
tx.q4.rt.input <- getTxQ4RT(tx.q4, ensGene.tx.rt.input, headon=headon)

q4s <- list()
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_q4_bp_indel_ho.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("Q1", "Q2", "Q3", "Q4")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.rt.input[[q]]
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
      q4[1, q] <- getMutPerMbTxs(bp, unique(CntxGtx$ensembl_gene_id))
  
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      q4[2, q] <- getMutPerMbTxs(bp, unique(GntxCtx$ensembl_gene_id))
   }
   
   q4s[[i]] <- q4
   barplot(q4, ylab="Breaks/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3)
}
dev.off()

###
## Refining the plot
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_q4_bp_indel_ho.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                             RT-Tx (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                             RT-Tx (Co-directional)", "", "", "")
 
   barplot(q4, ylab="Breaks/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))   ##, ylim=c(0, 410))
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5: Log2 transcription-associated SNV asymmetry (Figure 3)
# Name: 
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_q4_bp_indel_ho_log2_ylim0.1.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                            RT-Tx (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                            RT-Tx (Co-directional)", "", "", "")
   
   cols <- c()   
   for (c in 1:ncol(q4))
      if (log2(q4[1,c]/q4[2,c]) > 0)
         cols <- c(cols, "#4D4D4D")
      else
         cols <- c(cols, "#E6E6E6")    
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.1, 0.1), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=cols)
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()





# -----------------------------------------------------------------------------
# Step: Test
# Name: 
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
getQ3ICNTx <- function(s1, txs) {
   return(quantile(subset(s1, ensembl_gene_id %in% txs)$iCN)[4])
}

getQ2ICNTx <- function(s1, txs) {
   return(quantile(subset(s1, ensembl_gene_id %in% txs)$iCN)[3])
}

getMaxICNTx <- function(s1, txs) {
   return(max(subset(s1, ensembl_gene_id %in% txs)$iCN))
}

getMedianICNTx <- function(s1, txs) {
   return(median(subset(s1, ensembl_gene_id %in% txs)$iCN))
}

q4s.icn <- list()
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_q4_icn_q3_ns.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4.icn <- as.matrix(toTable(0, 4, 2, c("Q1", "Q2", "Q3", "Q4")))
   rownames(q4.icn) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.rt.indel.ns[[q]]
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
      q4.icn[1, q] <- getQ3ICNTx(bp, unique(CntxGtx$ensembl_gene_id))
  
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      q4.icn[2, q] <- getQ3ICNTx(bp, unique(GntxCtx$ensembl_gene_id))
   }
 
   q4s.icn[[i]] <- q4.icn
   #colnames(q4.icn) <- c("                             Tx (Stressed)", "", "", "")
   colnames(q4.icn) <- c("                            Tx (Not stressed)", "", "", "")
   
   barplot(q4.icn, ylab="Q3 iCN", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"), ylim=c(0, 4))
   mtext(paste0("log2(", paste(rownames(q4.icn), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()






#q4s <- asym.q4s   ## TO-DO remove
save(tx.snv.s6, asyms, q4s, file=paste0(wd.asym.plots, "sclc_asym_tx_snv_s6_asyms_q4s.RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "sclc_asym_tx_snv_s6_asyms_q4s.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   #colnames(q4) <- c(paste0("                             ", paste(rownames(q4), collapse="/")), "", "", "")
   colnames(q4) <- c("", "         Tx", "", "")
   
   xlab <- ""
   if (i == 1 || i == 4)
      xlab <- paste(rownames(q4), collapse="/")
   barplot(q4, ylab="SNVs/Mb", xlab=xlab, main=paste(rownames(asyms[[i]]), collapse="/"), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))   ##, "white", "gray", "white", "gray", "white", "gray"))
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5: Log2 transcription-associated SNV asymmetry (Figure 3)
# Name: 
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "sclc_asym_tx_snv_s6_asyms_q4s_log2_ylim2.pdf"), height=4, width=6)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "        Tx", "", "")
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Log2 asymmetry", ylim=c(-2, 2), main=paste(rownames(asyms[[i]]), collapse="/"), beside=TRUE, width=.3, col=c("#E6E6E6", "#E6E6E6", "#E6E6E6", "#E6E6E6"))   ##, "white", "gray", "white", "gray", "white", "gray"))
}
dev.off()












for (s in 1:length(samples)) {
   sample <- samples[s]
   
   pdf(paste0(wd.asym.plots, "sclc_snv_s6_", sample, ".pdf"))
   par(mfrow=c(3,1))
   #par(mar=c(3.5, 3.5, 1,1), mgp=c(2, 0.65, 0), cex=0.9)
   for (i in 1:length(idxs)) {
      idx <- idxs[i]
  
      hist(wt)
      hist(mpg)
      hist(disp)
      
      tx.snv.s6[[i]][[1]][[1]]   ## E.g. C>A Tx(+)
      tx.snv.s6[[i]][[1]][[2]]   ##      C>A Tx(-)
  
      tx.snv.s6[[i]][[2]][[1]]   ##      G>T Tx(+)
      tx.snv.s6[[i]][[2]][[2]]   ##      G>T Tx(-)
   }
   dev.off()
}










# > getMutPerMb(tx.snv.s6[[i]][[1]][[1]])   ## E.g. C>A Tx(+)
# [1] 162.091
# >    getMutPerMb(tx.snv.s6[[i]][[1]][[2]])   ##      C>A Tx(-)
# [1] 293.3805
# >    getMutPerMb(tx.snv.s6[[i]][[2]][[1]])   ##      G>T Tx(+)
# [1] 281.7484
# >    getMutPerMb(tx.snv.s6[[i]][[2]][[2]])   ##      G>T Tx(-)
# [1] 175.6149
# >    
#  >    getMutPerMb0(tx.snv.s6[[i]][[1]][[1]])   ## E.g. C>A Tx(+)
# [1] 176.2969
# >    getMutPerMb0(tx.snv.s6[[i]][[1]][[2]])   ##      C>A Tx(-)
# [1] 264.0553
# >    getMutPerMb0(tx.snv.s6[[i]][[2]][[1]])   ##      G>T Tx(+)
# [1] 262.0042
# >    getMutPerMb0(tx.snv.s6[[i]][[2]][[2]])   ##      G>T Tx(-)
# [1] 194.1061








 
   snv.s6[[1]][[1]][[1]] <- vcf.CA.p   ## C>A/G>T
   snv.s6[[1]][[1]][[2]] <- vcf.CA.m
   snv.s6[[1]][[2]][[1]] <- vcf.GT.p
   snv.s6[[1]][[2]][[2]] <- vcf.GT.m
   
   snv.s6[[2]][[1]][[1]] <- vcf.CG.p   ## C>G/G>C
   snv.s6[[2]][[1]][[2]] <- vcf.CG.m
   snv.s6[[2]][[2]][[1]] <- vcf.GC.p
   snv.s6[[2]][[2]][[2]] <- vcf.GC.m
   
   snv.s6[[3]][[1]][[1]] <- vcf.CT.p   ## C>T/G>A
   snv.s6[[3]][[1]][[2]] <- vcf.CT.m
   snv.s6[[3]][[2]][[1]] <- vcf.GA.p
   snv.s6[[3]][[2]][[2]] <- vcf.GA.m
   
   snv.s6[[4]][[1]][[1]] <- vcf.CG.p
   snv.s6[[4]][[1]][[2]] <- vcf.CG.m
   snv.s6[[4]][[2]][[1]] <- vcf.GC.p
   snv.s6[[4]][[2]][[2]] <- vcf.GC.m
   
   
   snv.s6[[1]][[1]] <- vcf.CA
   snv.s6[[1]][[2]] <- vcf.GT
   snv.s6[[2]][[1]] <- vcf.CG
   snv.s6[[2]][[2]] <- vcf.GC
   
   
   ## C>A/G>T (Smoking)
   vcf.CA <- subset(subset(vcf, REF == "C"), ALT == "A")
   vcf.CA.p <- subset(vcf.CA, strand == 1)
   vcf.CA.m <- subset(vcf.CA, strand == -1)
 
   vcf.GT <- subset(subset(vcf, REF == "G"), ALT == "T")
   vcf.GT.p <- subset(vcf.GT, strand == 1)
   vcf.GT.m <- subset(vcf.GT, strand == -1)
}


##for (s in 1:length(samples)) {
for (s in 1:10) {
   #sample <- SAMPLE
   sample <- samples[s]
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep="")
  
   ## C>A/G>T (Smoking)
   vcf.CA <- subset(subset(vcf, REF == "C"), ALT == "A")
   vcf.CA.p <- subset(vcf.CA, strand == 1)
   vcf.CA.m <- subset(vcf.CA, strand == -1)
      
   vcf.GT <- subset(subset(vcf, REF == "G"), ALT == "T")
   vcf.GT.p <- subset(vcf.GT, strand == 1)
   vcf.GT.m <- subset(vcf.GT, strand == -1)

   ## C>G/G>C
   vcf.CG <- subset(subset(vcf, REF == "C"), ALT == "G")
   vcf.CG.p <- subset(vcf.CG, strand == 1)
   vcf.CG.m <- subset(vcf.CG, strand == -1)
   
   vcf.GC <- subset(subset(vcf, REF == "G"), ALT == "C")
   vcf.GC.p <- subset(vcf.GC, strand == 1)
   vcf.GC.m <- subset(vcf.GC, strand == -1)
   
   ## C>T/G>A (APOBEC/AID)
   vcf.CT <- subset(subset(vcf, REF == "C"), ALT == "T")
   vcf.CT.p <- subset(vcf.CT, strand == 1)
   vcf.CT.m <- subset(vcf.CT, strand == -1)
   
   vcf.GA <- subset(subset(vcf, REF == "G"), ALT == "A")
   vcf.GA.p <- subset(vcf.GA, strand == 1)
   vcf.GA.m <- subset(vcf.GA, strand == -1)

   ## T>A/A>T
   vcf.TA <- subset(subset(vcf, REF == "T"), ALT == "A")
   vcf.TA.p <- subset(vcf.TA, strand == 1)
   vcf.TA.m <- subset(vcf.TA, strand == -1)
   
   vcf.AT <- subset(subset(vcf, REF == "A"), ALT == "T")
   vcf.AT.p <- subset(vcf.AT, strand == 1)
   vcf.AT.m <- subset(vcf.AT, strand == -1)
   
   ## T>C/A>G (Sig. 5)
   vcf.TC <- subset(subset(vcf, REF == "T"), ALT == "C")
   vcf.TC.p <- subset(vcf.TC, strand == 1)
   vcf.TC.m <- subset(vcf.TC, strand == -1)
   
   vcf.AG <- subset(subset(vcf, REF == "A"), ALT == "G")
   vcf.AG.p <- subset(vcf.AG, strand == 1)
   vcf.AG.m <- subset(vcf.AG, strand == -1)
   
   ## T>G/A>C
   vcf.TG <- subset(subset(vcf, REF == "T"), ALT == "G")
   vcf.TG.p <- subset(vcf.TG, strand == 1)
   vcf.TG.m <- subset(vcf.TG, strand == -1)
   
   vcf.AC <- subset(subset(vcf, REF == "A"), ALT == "C")
   vcf.AC.p <- subset(vcf.AC, strand == 1)
   vcf.AC.m <- subset(vcf.AC, strand == -1)
}
