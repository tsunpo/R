# =============================================================================
# Manuscript   : The dangeous case of DNA replication in SCLC
# Chapter II   : Correlating replicaiton-fork stress with transcriptional strand asymmetry
# Name         : manuscripts/asymmetries/sclc-asym-tx-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 31/05/18
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
BASE <- "LUAD"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.data  <- file.path(wd.asym,  "data")
wd.asym.plots <- file.path(wd.asym,  "plots")

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "luad_wgs_n39-5.list"), header=F, rownames=F, sep="")

load(file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

# -----------------------------------------------------------------------------
# Step 6.2: Divide genes into two groups (hand-on and co-directional)
# Last Modified: 20/02/18
# -----------------------------------------------------------------------------
load(file.path(wd.asym.data, paste0(base, "_asym_tx_rt.RData")))
# > nrow(ensGene.tx.rt)
# [1] 10352   ## tpm.gene_r5_p47
# [1] 9394    ## tpm.gene_tpm0

ensGene.tx.rt.nona <- subset(subset(ensGene.tx.rt, !is.na(SLOPE_START)), !is.na(SLOPE_END))
ensGene.tx.rt.nona$SIGN <- sign(ensGene.tx.rt.nona$SLOPE_START / ensGene.tx.rt.nona$SLOPE_END)
ensGene.tx.rt.nona.sign <- subset(ensGene.tx.rt.nona, SIGN == 1)[,-10]
# > nrow(ensGene.tx.rt.nona)
# [1] 10255   ## tpm.gene_r5_p47
# [1] 9325
# > nrow(ensGene.tx.rt.nona.sign)
# [1] 9529   ## tpm.gene_r5_p47
# [1] 8646

## TRICK: Assign replication has a "-" slope, therefore change/flip it's RT slope to "+" to be consistant with e.g. Tx(+)
ensGene.tx.rt.nona.sign$RT <- -1
ensGene.tx.rt.nona.sign[ensGene.tx.rt.nona.sign$SLOPE_START < 0,]$RT <- 1

## ADD 13/05/18 from sclc-asym-tx-pcg-non.R
load(file.path(wd.asym.data, paste0(base, "_asym_tx_snv.RData")))
ens.tx.snv.input <- intersect(unique(tx.snv$ensembl_gene_id), rownames(tpm.gene.input.log2))   ## CHANGE 10/04/18 Needed in Step 5; ADD 02/02/18
tx.snv.input <- subset(tx.snv, ensembl_gene_id %in% ens.tx.snv.input)

###
##
headon <- T   ## Replication–transcription Head-on (T) / Co-directional (F) collisions
prefix <- "ho"
if (headon == F) prefix <- "cd"

ens.tx.rt.input <- intersect(unique(tx.snv.input$ensembl_gene_id), unique(ensGene.tx.rt.nona.sign$ensembl_gene_id))
ensGene.tx.rt.input <- ensGene.tx.rt.nona.sign[ens.tx.rt.input,]
tx.q4.rt.input <- getTxQ4RT(ensGene.tx.rt.input, headon=headon, tpm.gene.input.log2)
# > length(unique(tx.snv.input$ensembl_gene_id))
# [1] 9499   ## tpm.gene_r5_p47
# [1] 8675
# > nrow(ensGene.tx.rt.input)
# [1] 8725   ## tpm.gene_r5_p47
# [1] 7959

# for (q in 1:4)   ## CD
#    print(length(intersect(ens.tx.rt.input, tx.q4.rt.input[[q]])))
# [1] 1004   ## 1099
# [1] 1003
# [1] 1003
# [1] 1003

# for (q in 1:4)   ## HO
#    print(length(intersect(ens.tx.rt.input, tx.q4.rt.input[[q]])))
# [1] 987   ## 1083
# [1] 986
# [1] 986
# [1] 987

###
## (Figure S4)
#load(file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6_q4s.RData")))
#q4s.rt <- list()   ## ADD 31/05/18
#q4s.rt[[1]] <- q4s

q4s <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_", prefix, ".pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.rt.input[[q]]
  
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      #CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
      #q4[1, q] <- getMutPerMbTxs(bp, unique(CntxGtx$ensembl_gene_id))
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      #q4[1, q] <- getMutPerMbTxs(tx.input, txs.cg)
      q4[1, q] <- getMutPerMbTxs(CntxGtx, txs.cg)
      
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      #GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      #q4[2, q] <- getMutPerMbTxs(bp, unique(GntxCtx$ensembl_gene_id))
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      #q4[2, q] <- getMutPerMbTxs(tx.input, txs.gc)
      q4[2, q] <- getMutPerMbTxs(GntxCtx, txs.gc)
        
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()
save(q4s, file=file.path(wd.asym.data, paste0(base,"_asym_tx_rt_snv_q4s_rt_", prefix, ".RData")))

## ADD 22/02/18
if (headon == F) {
   q4s.rt[[2]] <- q4s
} else {
   q4s.rt[[3]] <- q4s
   save(q4s.rt, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt.RData")))
}

###
## Refining the plot
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_", prefix, ".pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
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
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_", prefix, "_log2_ylim1.5.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                                    Lagging strand (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                    Leading strand (Co-directional)", "", "", "")
 
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-1.5, 1.5), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=getLog2Colours(q4))
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 7: Replication–transcription collision-induced SNV asymmetry (Figure 1)
# Last Modified: 22/02/18
# -----------------------------------------------------------------------------
#for (i in 1:length(idxs)) {
for (i in 3:3) {
   idx <- idxs[i]
   ymax <- max(q4s.rt[[1]][[i]])
   for (j in 2:3)
      if (max(q4s.rt[[j]][[i]]) > ymax)
          ymax <- max(q4s.rt[[j]][[i]])
   
   pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_strands_", REFS[idx], ">", ALTS[idx], "_ylim0.4.pdf")), height=4.5, width=7)
   par(mfrow=c(2, 3)) 
   for (j in 1:3) {
      q4 <- getRTTxQ4(j, i)
      
      barplot(q4, ylab="SNVs/Mb", main=getRTTxMain(i, asyms), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"), ylim=c(0, ymax))   ##, xlab="Number of genes")
      mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)      
   }
   
   for (j in 1:3) {
      q4 <- getRTTxQ4(j, i)
  
      barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.4, 0.4), main=getRTTxMain(i, asyms), beside=TRUE, width=.3, col=getLog2Colours(q4))
      mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
   }
   dev.off()
}

# -----------------------------------------------------------------------------
# Step: Test
# Last Modified: 17/06/18
# -----------------------------------------------------------------------------
load(file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6_q4s.RData")))
q4s.rt.ho <- list()   ## ADD 31/05/18
q4s.rt.ho[[1]] <- q4s.rt[[3]]

q4s.rt.ho.tmp <- q4s.rt.ho
q4s.rt.ho.tmp[[2]] <- q4s.rt.ho[[3]]
q4s.rt.ho.tmp[[3]] <- q4s.rt.ho[[2]]

genes.ho.plus  <- rownames(subset(ensGene[genes.ho,], strand == 1))
genes.ho.minus <- rownames(subset(ensGene[genes.ho,], strand == -1))

q4s <- list()
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.rt.input[[q]]
      txs <- intersect(txs, genes.ho.plus)
  
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      q4[1, q] <- getMutPerMbTxs(CntxGtx, txs.cg)
  
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(GntxCtx, txs.gc)
  
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
}

##
for (i in 1:1) {
   idx <- idxs[i]
   ymax <- max(q4s.rt.ho[[1]][[i]])
   for (j in 2:3)
      if (max(q4s.rt.ho[[j]][[i]]) > ymax)
         ymax <- max(q4s.rt.ho[[j]][[i]])
 
   pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_strands_", REFS[idx], ">", ALTS[idx], "_ho_ylim1.5.pdf")), height=4.5, width=7)
   par(mfrow=c(2, 3)) 
   for (j in 1:3) {
      q4 <- getRTTxQ4_2(q4s.rt.ho, j, i)
  
      barplot(q4, ylab="SNVs/Mb", main=getRTTxMain_2(i, asyms), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"), ylim=c(0, ymax))   ##, xlab="Number of genes")
      mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)      
   }
 
   for (j in 1:3) {
      q4 <- getRTTxQ4_2(q4s.rt.ho, j, i)
  
      barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-1.5, 1.5), main=getRTTxMain_2(i, asyms), beside=TRUE, width=.3, col=getLog2Colours(q4))
      mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
   }
   dev.off()
}