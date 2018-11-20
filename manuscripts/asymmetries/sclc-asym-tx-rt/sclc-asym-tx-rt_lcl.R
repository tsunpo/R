# =============================================================================
# Manuscript   : The dangeous case of DNA replication in SCLC
# Chapter II   : Correlating replicaiton-fork stress with transcriptional strand asymmetry
# Name         : manuscripts/asymmetries/sclc-asym-tx-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Asymmetry.R", "ReplicationTiming.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.data  <- file.path(wd.asym,  "data")
wd.asym.plots <- file.path(wd.asym,  "plots")

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 6.2: Divide genes into two groups (hand-on and co-directional)
# Last Modified: 20/02/18
# -----------------------------------------------------------------------------
#load(file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene <- tpm.gene[setdiff(rownames(tpm.gene), outliers1.5),]   ## ADD 26/06/18
#tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
#tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=F)
#tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

load(file.path(wd.asym.data, paste0(base, "_asym_tx_rt.RData")))
# > nrow(ensGene.tx.rt)   ## All genes
# [1] 18440
# > nrow(ensGene.tx.rt)   ## All protein coding
# [1] 16410
# > nrow(ensGene.tx.rt)   ## Nonredundant proten coding
# [1] 10604

ensGene.tx.rt.nona <- subset(subset(ensGene.tx.rt, !is.na(SLOPE_START)), !is.na(SLOPE_END))
ensGene.tx.rt.nona$SIGN <- sign(ensGene.tx.rt.nona$SLOPE_START / ensGene.tx.rt.nona$SLOPE_END)
ensGene.tx.rt.nona.sign <- subset(ensGene.tx.rt.nona, SIGN == 1)[,-11]
# > nrow(ensGene.tx.rt.nona)        ## All genes
# [1] 18440
# > nrow(ensGene.tx.rt.nona.sign)   ## All genes (LCL)
# [1] 17243

## TRICK: Assign replication has a "-" slope, therefore change/flip it's RT slope to "+" to be consistant with e.g. Tx(+)
ensGene.tx.rt.nona.sign$RT <- -1
ensGene.tx.rt.nona.sign[ensGene.tx.rt.nona.sign$SLOPE_START < 0,]$RT <- 1
#ensGene.tx.rt.nona.sign <- ensGene.tx.rt.nona.sign[setdiff(rownames(ensGene.tx.rt.nona.sign), outliers1.5),]   ## ADD 26/06/18
# > nrow(ensGene.tx.rt.nona.sign)
# [1] 9690

###
##
headon <- "ho"   ## Replication–transcription Head-on (ho) / Co-directional (cd) collisions

q4s <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_", headon, ".pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
   ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)
   GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
   
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   overlaps <-intersect(rownames(ensGene.tx.rt.nona.sign), ens.asym)
   ensGene.tx.rt.input <- ensGene.tx.rt.nona.sign[overlaps,]
   ens.asym.rt.q4 <- getTxQ4RT(ensGene.tx.rt.input, headon=headon, tpm.gene.input.log2)
   
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- ens.asym.rt.q4[[q]]
  
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)    ## SWAP 23/06/18
      CntxGtx.cg <- subset(CntxGtx, ensembl_gene_id %in% txs.cg)   ## SWAP 23/06/18
      q4[1, q] <- getMutPerMbTxs(CntxGtx.cg, txs.cg)
      
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      GntxCtx.gc <- subset(GntxCtx, ensembl_gene_id %in% txs.gc) 
      q4[2, q] <- getMutPerMbTxs(GntxCtx.gc, txs.gc)
        
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()
#save(q4s, file=file.path(wd.asym.data, paste0(base,"_asym_tx_rt_snv_q4s_rt_", prefix, ".RData")))

## ADD 22/02/18
if (headon == "cd") {
   q4s.rt[[2]] <- q4s
} else {
   q4s.rt[[3]] <- q4s
   save(q4s.rt, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt.RData")))
}

###
## Refining the plot
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_", headon, ".pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon == "ho")
      colnames(q4) <- c("                                     Lagging strand (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                     Leading strand (Co-directional)", "", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
 
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))   ##, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 6.3: Log2 transcription-associated SNV asymmetry (Figure S5)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_", headon, "_log2_ylim1.75.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon == "ho")
      colnames(q4) <- c("                                    Lagging strand (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                    Leading strand (Co-directional)", "", "", "")
 
   barplot(-log2(q4[1,]/q4[2,]), ylab="TCR efficiency", ylim=c(-0.5, 1.75), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=getLog2Colours(q4))
   mtext(paste0("-log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 7: Replication–transcription collision-induced SNV asymmetry (Figure 1)
# Last Modified: 22/02/18
# -----------------------------------------------------------------------------
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   ymax <- max(q4s.rt[[1]][[i]])
   for (j in 2:3)
      if (max(q4s.rt[[j]][[i]]) > ymax)
         ymax <- max(q4s.rt[[j]][[i]])
   
   pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], ".pdf")), height=4.5, width=7)
   par(mfrow=c(2, 3)) 
   for (j in 1:3) {
      q4 <- getRTTxQ4(j, i)
      
      barplot(q4, ylab="SNVs/Mb", main=getRTTxMain(i, asyms), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"), ylim=c(0, ymax))   ##, xlab="Number of genes")
      mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)      
   }
   
   for (j in 1:3) {
      q4 <- getRTTxQ4(j, i)
  
      barplot(-log2(q4[1,]/q4[2,]), ylab="TCR efficiency", ylim=c(-0.5, 1.75), main=getRTTxMain(i, asyms), beside=TRUE, width=.3, col=getLog2Colours(q4))
      mtext(paste0("-log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
   }
   dev.off()
}
