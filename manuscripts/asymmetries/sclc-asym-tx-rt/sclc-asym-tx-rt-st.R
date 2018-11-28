# =============================================================================
# Manuscript   : The dangeous case of DNA replication in SCLC
# Chapter II   : Correlating replicaiton-fork stress with transcriptional strand asymmetry
# Name         : manuscripts/asymmetries/sclc-asym-tx-rt-st.R
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
BASE <- "SCLC"
base <- tolower(BASE)

#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd.anlys      <- file.path(wd, BASE, "analysis")
wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.data  <- file.path(wd.asym,  "data/bstrps")
wd.asym.plots <- file.path(wd.asym,  "plots/cut3/bstrps")

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 8: Fetching RT-Tx conflicts genes
# Last Modified: 17/06/18
# -----------------------------------------------------------------------------
load(file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt.RData")))
q4s.rt.cd <- list()   ## ADD 31/05/18
q4s.rt.cd[[1]] <- q4s.rt[[2]]
q4s.rt.ho <- list()   ## ADD 31/05/18
q4s.rt.ho[[1]] <- q4s.rt[[3]]

## Followed Step 6.2
headon <- "cd"   ## Replication–transcription Head-on (ho) / Co-directional (cd) collisions
strand <- "fw"   ## Forward transcription (fw) / Reverse transcription (re) strands

q4s <- list()
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   ## Cntx:Gtx  =  tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   CntxGtx <- tx.snv.s6[[i]][[1]][[1]]      ## ADD 23/06/18
   if (strand == "re")
      CntxGtx <- tx.snv.s6[[i]][[2]][[2]]   ## SWAP 23/06/18
    
   ## Gntx:Ctx  =  tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
   GntxCtx <- tx.snv.s6[[i]][[2]][[1]]
   if (strand == "re")             
      GntxCtx <- tx.snv.s6[[i]][[1]][[2]]
   
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ensGene.tx.rt.st <- ensGene.tx.rt.nona.sign[ens.asym,]
   ens.asym.rt.st.q4 <- getTxQ4RTStrand(ensGene.tx.rt.st, headon=headon, strand=strand, tpm.gene.input.log2)

   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- ens.asym.rt.st.q4[[q]]
  
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)    ## SWAP 23/06/18  
      CntxGtx.cg <- subset(CntxGtx, ensembl_gene_id %in% txs.cg)   ## SWAP 23/06/18
      q4[1, q] <- getMutPerMbTxs(CntxGtx.cg, txs.cg)
  
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      GntxCtx.gc <- subset(GntxCtx, ensembl_gene_id %in% txs.gc)
      q4[2, q] <- getMutPerMbTxs(GntxCtx.gc, txs.gc)
  
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
}
q4s.rt.cd[[2]] <- q4s

##
save(q4s.rt.cd, q4s.rt.ho, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt_st.RData")))

# -----------------------------------------------------------------------------
# Step 9: 
# Last Modified: 22/02/18
# -----------------------------------------------------------------------------
q4s.rt.st <- q4s.rt.ho

for (i in 1:6) {
   idx <- idxs[i]
   ymax <- max(q4s.rt.st[[1]][[i]])
   for (j in 2:3)
      if (max(q4s.rt.st[[j]][[i]]) > ymax)
         ymax <- max(q4s.rt.st[[j]][[i]])
 
   pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], "_", headon, ".pdf")), height=4.5, width=7)
   par(mfrow=c(2, 3)) 
   for (j in 1:3) {
      q4 <- getRTTxQ4_2(q4s.rt.st, j, i, headon)
  
      barplot(q4, ylab="SNVs/Mb", main=getRTTxMain_2(i, asyms, headon), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"), ylim=c(0, ymax))   ##, xlab="Number of genes")
      mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)      
   }
 
   for (j in 1:3) {
      q4 <- getRTTxQ4_2(q4s.rt.st, j, i, headon)
  
      barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-1.5, 1.5), main=getRTTxMain_2(i, asyms, headon), beside=TRUE, width=.3, col=getLog2Colours(q4))
      mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
   }
   dev.off()
}

# -----------------------------------------------------------------------------
# Step 10: Final reports
# Last Modified: 22/02/18
# -----------------------------------------------------------------------------
getRateQ4 <- function(mut.q4, len.q4) {
   rate.q4 <- list()
 
   for (q in 1:4) {
      muts <- as.numeric(mut.q4[[q]])
      lens <- as.numeric(len.q4[[q]])
      
      rate.q4[[q]] <- muts / lens * 1000000
   }
 
   return(rate.q4)
}

getLog2RatioQ4 <- function(cg.q4, gc.q4) {
   ratio.q4 <- list()
 
   for (q in 1:4) {
      cgs <- as.numeric(cg.q4[[q]])
      gcs <- as.numeric(gc.q4[[q]])

      ratio.q4[[q]] <- -log2(cgs / gcs)
   }
 
   return(ratio.q4)
}

getLog2RatioPerMbQ4 <- function(ratio.q4, len.q4) {
   ratio.permb.q4 <- list()
 
   for (q in 1:4) {
      ratios  <- as.numeric(ratio.q4[[q]])
      lengths <- as.numeric(len.q4[[q]])
  
      ratio.permb.q4[[q]] <- ratios / lengths * 1000000
   }
 
   return(ratio.permb.q4)
}


getLengthQ4 <- function(tx.q4) {
   length.q4 <- list()
   
   for (q in 1:4) {
      txs <- tx.q4[[q]]
  
      lengths <- list()
      for (g in 1:length(txs))
         lengths[[g]] <- getLengthTx(txs[g])
      
      length.q4[[q]] <- lengths  
   }
 
   return(length.q4)
}

getStrandQ4 <- function(tx.q4) {
   strand.q4 <- list()
 
   for (q in 1:4) {
      txs <- tx.q4[[q]]
  
      strands <- list()
      for (g in 1:length(txs))
         strands[[g]] <- ensGene[txs[g],]$strand
  
      strand.q4[[q]] <- strands  
   }
 
   return(strand.q4)
}

getMutationQ4 <- function(s1.all, txs.q4) {
   mut.q4 <- list()
   
   for (q in 1:4) {
      txs <- txs.q4[[q]]
      txs <- intersect(txs, unique(s1.all$ensembl_gene_id))   ## SWAP txs order BUG BUG BUG, SWAP 23/06/18
      s1 <- subset(s1.all, ensembl_gene_id %in% txs)      ## SWAP 23/06/18
      
      muts <- list()
      for (g in 1:length(txs))
         muts[[g]] <- nrow(subset(s1, ensembl_gene_id == txs[g]))   #getMutPerMbTx(s1, txs[g])

      mut.q4[[q]] <- muts
   }
 
   return(mut.q4)
}

getYlim <- function(q4s, start, end) {
   ymin <- 1000000000
   ymax <- 0
   
   for (i in start:end) {
      for (rt in 1:length(headons)) {
         for (st in 1:length(strands)) {
            for (q in 1:4) {
               y <- min(as.numeric(q4s[[i]][[rt]][[st]][[q]]))
               if (y < ymin)
                  ymin <- y
            }
         }
      }
   }
   
   for (i in start:end) {
      for (rt in 1:length(headons)) {
         for (st in 1:length(strands)) {
            for (q in 1:4) {
               y <- max(as.numeric(q4s[[i]][[rt]][[st]][[q]]))
               if (y != Inf)
                  if (y > ymax)
                     ymax <- y
            }
         }
      }
   }
   
   return(c(ymin, ymax))
}

plotQ4S <- function(q4s, file.name, file.main, mtext, ylim, ylab, isLog10) {
   counts <- c()
   quantiles <- c()
   for (q in 1:4) {
      counts    <- c(counts, as.numeric(q4s[[q]]))
      quantiles <- c(quantiles, mapply(x = 1:length(q4s[[q]]), function(x) paste0("q", q)))
   }
 
   pdf(file.name, height=6, width=4)
   names <- c("Q1", "Q2", "Q3", "Q4")
   if (isLog10)
      boxplot(log10(as.numeric(counts))~quantiles, ylab=ylab, xlab="Expression", names=names, main=file.main, ylim=log10(ylim))
   else
      boxplot(as.numeric(counts)~quantiles, ylab=ylab, xlab="Expression", names=names, main=file.main, ylim=ylim)
   mtext(mtext, cex=0.8, font=3, line=0.5)
   dev.off()
}

plotQ4SPlus <- function(q4s, file.name, file.main, mtext, ylim, ylab, isLog10, isPlus) {
   counts <- c()
   quantiles <- c()
   for (q in 1:4) {
      for (p in 1:length(q4s[[q]])) {  
         if (isPlus) {
            if (as.numeric(q4s[[q]][[p]]) > 0) {
               counts    <- c(counts, as.numeric(q4s[[q]][[p]]))
               quantiles <- c(quantiles, paste0("q", q))
            }
         } else {
            if (as.numeric(q4s[[q]][[p]]) < 0) {
               counts    <- c(counts, as.numeric(q4s[[q]][[p]]))
               quantiles <- c(quantiles, paste0("q", q))
            }
         }
      }
   }
 
   pdf(file.name, height=6, width=4)
   names <- c("Q1", "Q2", "Q3", "Q4")
   if (isLog10)
      boxplot(log10(as.numeric(counts))~quantiles, ylab=ylab, xlab="Expression", names=names, main=file.main, ylim=log10(ylim))
   else
      boxplot(as.numeric(counts)~quantiles, ylab=ylab, xlab="Expression", names=names, main=file.main, ylim=ylim)
   mtext(mtext, cex=0.8, font=3, line=0.5)
   dev.off()
}

##
headons <- c("cd", "ho")
strands <- c("fw", "re")
q4s.s6.rt.st.gen <- initTableS6RTSt()
q4s.s6.rt.st.len <- initTableS6RTSt()
#q4s.s6.rt.st.strand <- initTableS6RTSt()
q4s.s6.rt.st.mut <- initTableS6RTSt()
q4s.s6.rt.st.mut.cg <- initTableS6RTSt()
q4s.s6.rt.st.mut.gc <- initTableS6RTSt()
q4s.s6.rt.st.ratio  <- initTableS6RTSt()
#q4s.s6.rt.st.ratio.permb <- initTableS6RTSt()
q4s.s6.rt.st.rate <- initTableS6RTSt()
#q4s.s6.rt.st.g2g <- initTableS6RTSt()
for (i in 1:length(idxs)) {
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ensGene.tx.rt.st <- ensGene.tx.rt.nona.sign[ens.asym,]
   
   ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   #CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
   ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
   #GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
   
   for (rt in 1:length(headons)) {
      for (st in 1:length(strands)) {
         if (rt == 1 && st == 1) {
            CntxGtx <- tx.snv.s6[[i]][[1]][[1]]
            GntxCtx <- tx.snv.s6[[i]][[2]][[1]]         
         } else if (rt == 1 && st == 2) {
            CntxGtx <- tx.snv.s6[[i]][[2]][[2]]
            GntxCtx <- tx.snv.s6[[i]][[1]][[2]]
         } else if (rt == 2 && st == 1) {
            CntxGtx <- tx.snv.s6[[i]][[1]][[1]]
            GntxCtx <- tx.snv.s6[[i]][[2]][[1]]         
         } else if (rt == 2 && st == 2) {
            CntxGtx <- tx.snv.s6[[i]][[2]][[2]]
            GntxCtx <- tx.snv.s6[[i]][[1]][[2]]
         }
          
         txs <- getTxRTStrand(ensGene.tx.rt.st, headon=headons[rt], strand=strands[st])
         unique <- intersect(unique(CntxGtx$ensembl_gene_id), unique(GntxCtx$ensembl_gene_id))
         txs <- intersect(txs, unique)
         txs.q4 <- getTxQ4RTStrand(ensGene.tx.rt.st[txs,], headon=headons[rt], strand=strands[st], tpm.gene.input.log2)
         
         CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
         GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs)
         s1.all <- rbind(CntxGtx, GntxCtx)

         q4s.s6.rt.st.gen[[i]][[rt]][[st]] <- txs.q4
         q4s.s6.rt.st.len[[i]][[rt]][[st]] <- getLengthQ4(txs.q4)
         #q4s.s6.rt.st.strand[[i]][[rt]][[st]] <- getStrandQ4(txs.q4)
         q4s.s6.rt.st.mut[[i]][[rt]][[st]] <- getMutationQ4(s1.all, txs.q4)
         q4s.s6.rt.st.mut.cg[[i]][[rt]][[st]] <- getMutationQ4(CntxGtx, txs.q4)
         q4s.s6.rt.st.mut.gc[[i]][[rt]][[st]] <- getMutationQ4(GntxCtx, txs.q4)
         
         #tx.q4.all <- getTxQ4(NA, tpm.gene.input.log2)
         #g2g.q4.all <- getG2GQ4(tx.q4.all, NULL)
         
         #q4s.s6.rt.st.g2g[[i]][[rt]][[st]] <- getG2GQ4fromAll(txs.q4, tx.q4.all, g2g.q4.all)
         #q4s.s6.rt.st.g2g[[i]][[rt]][[st]] <- getG2GQ4(txs.q4, NULL)
      }
   }
}
save(q4s.s6.rt.st.gen, q4s.s6.rt.st.len, q4s.s6.rt.st.mut, q4s.s6.rt.st.mut.cg, q4s.s6.rt.st.mut.gc, q4s.s6.rt.st.ratio, q4s.s6.rt.st.rate, file=file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt_st_ratio.RData")))

getAll<- function(i, q4s.mut, as.numeric) {
   if (!is.null(i))  
      idxs <- idx[i]
 
   muts <- c()
   for (i in 1:length(idxs))
      for (rt in 1:length(headons))
         for (st in 1:length(strands))
            for (q in 1:4)
               if (as.numeric)
                  muts <- c(muts, as.numeric(q4s.mut[[i]][[rt]][[st]][[q]]))
               else
                  muts <- c(muts, q4s.mut[[i]][[rt]][[st]][[q]])
   return(muts)
}

###
##
gens <- getAll(1, q4s.s6.rt.st.gen, as.numeric=F)
muts <- getAll(1, q4s.s6.rt.st.mut, as.numeric=T)
muts.cg <- getAll(1, q4s.s6.rt.st.mut.cg, as.numeric=T)
muts.gc <- getAll(1, q4s.s6.rt.st.mut.gc, as.numeric=T)
lens <- getAll(1, q4s.s6.rt.st.len, as.numeric=T)
lens <- lens/1000000
exps <- as.numeric(tpm.gene.input.log2[gens,]$MEDIAN)

grep("ENSG00000174469", gens)
lens[1232]
muts[1232]
muts.cg[1232]
muts.gc[1232]

mtext <- "C>A/G>T"
file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_muts-lens.png"))
png(file.name, height=400, width=400)
plot(muts~lens, ylab="Mutation", xlab="Gene length (Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_muts-lens_log10.png"))
png(file.name, height=400, width=400)
plot(log10(muts)~lens, ylab="Mutation (log10)", xlab="Gene length", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_muts-lens_log10_log10.png"))
png(file.name, height=400, width=400)
plot(log10(muts)~log10(lens), ylab="Mutation (log10)", xlab="Gene length", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

lens <- lens/1000000
rates <- muts/lens

file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_rates-lens_log10_log10.png"))
png(file.name, height=400, width=400)
plot(log10(rates)~log10(lens), ylab="Mutation rate log10(SNV/Mb)", xlab="Gene length log10(Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_rates-lens_log10.png"))
png(file.name, height=400, width=400)
plot(log10(rates)~lens, ylab="Mutation rate log10(SNV/Mb)", xlab="Gene length (Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_rates-lens.png"))
png(file.name, height=400, width=400)
plot(rates~lens, ylab="Mutation rate (SNV/Mb)", xlab="Gene length (Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

###
## 09/10/18
ensGene.tx.rt.nona.sign.lcl.ca <- ensGene.tx.rt.nona.sign[gens,]
ensGene.tx.rt.nona.sign.lcl.ca$MUT <- muts
ensGene.tx.rt.nona.sign.lcl.ca$MUT_CG <- muts.cg
ensGene.tx.rt.nona.sign.lcl.ca$MUT_GC <- muts.gc
ensGene.tx.rt.nona.sign.lcl.ca$RATIO_LOG2 <- -log2(muts.cg / muts.gc)
ensGene.tx.rt.nona.sign.lcl.ca$MEDIAN <- exps
ensGene.tx.rt.nona.sign.lcl.ca$CD <- ensGene.tx.rt.nona.sign.lcl.ca$strand * ensGene.tx.rt.nona.sign.lcl.ca$RT

ratios <- -log2(muts.cg / muts.gc)
file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_ratios-lens_log10.png"))
png(file.name, height=400, width=400)
plot(ratios~log10(lens), ylab="TCR efficiency", xlab="Gene length log10(Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

## Whats the distribution of mutation numbers? 10/10/18
file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/asymmetries/sclc-asym-tx-rt/plots/cut3/lcl/sclc_tx_rt_freq_muts_cg_log10.pdf"
d <- getDensityCount(log10(ensGene.input$MUT_CG))
ymax <- max(d$y)
pdf(file.name, height=6, width=6)
plot(d, ylab="Frequency", xlab="Mutation (log10)", main="C>A", ylim=c(0, ymax))
dev.off()

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/asymmetries/sclc-asym-tx-rt/plots/cut3/lcl/sclc_tx_rt_freq_muts_gc_log10.pdf"
d <- getDensityCount(log10(ensGene.input$MUT_GC))
ymax <- max(d$y)
pdf(file.name, height=6, width=6)
plot(d, ylab="Frequency", xlab="Mutation (log10)", main="G>T", ylim=c(0, ymax))
dev.off()

file.name <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/asymmetries/sclc-asym-tx-rt/plots/cut3/lcl/sclc_tx_rt_freq_ratios.pdf"
d <- getDensityCount(ensGene.input$RATIO_LOG2)
ymax <- max(d$y)
pdf(file.name, height=6, width=6)
plot(d, ylab="Frequency", xlab="TCR efficiency", main="C>A / G>T", ylim=c(0, ymax))
dev.off()



## CUT2, TCR
ensGene.input <- ensGene.tx.rt.nona.sign.ca
ensGene.input <- subset(ensGene.input, MUT_CG >= 2)
ensGene.input <- subset(ensGene.input, MUT_GC >= 2)
#ensGene.input <- subset(ensGene.input, gene_biotype == "protein_coding")
txs.q4 <- getTxQ4(ensGene.input, NA)
q4ss <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

## 13/10/18 TEST LCL
overlaps <- intersect(rownames(ensGene.tx.rt.nona.sign.ca), rownames(ensGene.tx.rt.nona.sign))
overlaps4 <- setdiff(rownames(ensGene.tx.rt.nona.sign.ca), overlaps)
> length(overlaps)
[1] 8480     8480/9278=0.9139901
> length(overlaps4)
[1] 798
> 8480 + 798
[1] 9278

nrow(ensGene.tx.rt.nona.sign.ca)
# [1] 9278
# [1] 723   ## TCD; 798/9278 = 0.08600992
length(overlaps)
# [1] 8480
# [1] 659     ## TCD (LCL) 0.9114799
# 7821/8555   ## 0.9142022

tx.q4.fix <- getTxQ4(tpm.gene.input.log2[overlaps,], NA)
for (q in 1:4)
   print(length(tx.q4.fix[[q]]))

## 26/11/18
ensGene.input <- ensGene.tx.rt.nona.sign.ca
ensGene.input <- subset(ensGene.input, MUT_CG >= 3)
ensGene.input <- subset(ensGene.input, MUT_GC >= 3)

overlaps <- intersect(rownames(ensGene.rt.tx), rownames(ensGene.input))
ensGene.rt.tx.ca <- cbind(ensGene.rt.tx[overlaps,], ensGene.input[overlaps, 13:16])

overlaps <- intersect(sclc.tx.consist, overlaps)
ensGene.rt.tx.ca.consist <- ensGene.rt.tx.ca[overlaps,]
# > length(sclc.tx.consist.ca)
# [1] 4125
overlaps <- intersect(sclc.tx.inconsist, overlaps)
ensGene.rt.tx.ca.inconsist <- ensGene.rt.tx.ca[overlaps,]

wd.asym.plots <- file.path(wd.asym,  "plots/cut3/bstrps")
overlaps2 <- intersect(overlaps, rownames(subset(de.tpm.gene, LOG2_FC > 0)))
overlaps3 <- intersect(overlaps, rownames(subset(de.tpm.gene, LOG2_FC < 0)))
## TEST (TCR, no 0; 09/10/18); PCG 11/10/18

#ensGene.input <- ensGene.tx.rt.nona.sign.ca[overlaps2,]   ## TEST LCL, DIFF, RB1
#ensGene.input <- subset(ensGene.input, gene_biotype == "protein_coding")
#ensGene.input <- subset(ensGene.input, protein_coding_non_redundant == T)
#ensGene.input <- subset(ensGene.input, MUT_CG >= 3)
#ensGene.input <- subset(ensGene.input, MUT_GC >= 3)
#ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
ensGene.input <- ensGene.rt.tx.ca.consist
#ensGene.input <- ensGene.rt.tx.ca.consist[overlaps2,]
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
ensGene.input <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.ho.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.ho.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

#ensGene.input <- ensGene.tx.rt.nona.sign.ca
#ensGene.input <- ensGene.tx.rt.nona.sign.ca[overlaps33,]   ## TEST LCL
#ensGene.input <- subset(ensGene.input, gene_biotype == "protein_coding")
#ensGene.input <- subset(ensGene.input, protein_coding_non_redundant == T)
#ensGene.input <- subset(ensGene.input, MUT_CG >= 3)
#ensGene.input <- subset(ensGene.input, MUT_GC >= 3)
#ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
ensGene.input <- ensGene.rt.tx.ca.consist
#ensGene.input <- ensGene.rt.tx.ca.consist[overlaps2,]
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

####
###
## Consistent 26/11/18
overlaps1 <- intersect(rownames(ensGene.rt.tx.ca), sclc.tx.consist)
ensGene.rt.tx.ca.consist <- ensGene.rt.tx.ca[overlaps1,]
# > length(overlaps1)
# [1] 4125
ensGene.input <- ensGene.rt.tx.ca.consist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)

## Consistent > TCR
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input1 <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input1, tx.q4.fix)
q4ss.tcr.cd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.cd[[q]] <- ensGene.input1[txs.q4[[q]],]$RATIO_LOG2

ensGene.input2 <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input2, tx.q4.fix)
q4ss.tcr.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.ho[[q]] <- ensGene.input2[txs.q4[[q]],]$RATIO_LOG2

## Consistent > TCD
ensGene.input <- ensGene.rt.tx.ca.consist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input1 <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input1, tx.q4.fix)
q4ss.tcd.cd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd.cd[[q]] <- ensGene.input1[txs.q4[[q]],]$RATIO_LOG2

ensGene.input2 <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input2, tx.q4.fix)
q4ss.tcd.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd.ho[[q]] <- ensGene.input2[txs.q4[[q]],]$RATIO_LOG2

## Consistent > TCR > RB1
overlaps2 <- intersect(overlaps1, rownames(subset(de.tpm.gene, LOG2_FC > 0)))
ensGene.input <- ensGene.rt.tx.ca.consist[overlaps2,]
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)

txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcr.rb <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.rb[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input1 <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input1, tx.q4.fix)
q4ss.tcr.rb.cd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.rb.cd[[q]] <- ensGene.input1[txs.q4[[q]],]$RATIO_LOG2

ensGene.input2 <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input2, tx.q4.fix)
q4ss.tcr.rb.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.rb.ho[[q]] <- ensGene.input2[txs.q4[[q]],]$RATIO_LOG2

## Consistent > TCR > RB1
overlaps3 <- intersect(overlaps1, rownames(subset(de.tpm.gene, LOG2_FC < 0)))
ensGene.input <- ensGene.rt.tx.ca.consist[overlaps3,]
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)

txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcr.wt <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.wt[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input1 <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input1, tx.q4.fix)
q4ss.tcr.wt.cd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.wt.cd[[q]] <- ensGene.input1[txs.q4[[q]],]$RATIO_LOG2

ensGene.input2 <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input2, tx.q4.fix)
q4ss.tcr.wt.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.wt.ho[[q]] <- ensGene.input2[txs.q4[[q]],]$RATIO_LOG2

## Plots
idx <- 1
file.main <- getMain(rownames(asyms[[i]]))
q4 <- q4s[[i]]
mtext <- paste0("-log2(", paste(rownames(q4), collapse="/"), ")")

ensGene.input <- ensGene.rt.tx.ca.consist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
ylim.tcr <- c(min(ensGene.input$RATIO_LOG2), max(ensGene.input$RATIO_LOG2))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCR.pdf"))
plotQ4SS(q4ss.tcr, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCR_RB.pdf"))
plotQ4SS(q4ss.tcr.rb, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCR_WT.pdf"))
plotQ4SS(q4ss.tcr.wt, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)

## BTRAPS, TCR (with 0), CUT 3 (26/11/18)
testW(q4ss.tcr.cd[[1]], q4ss.tcr.ho[[1]])
testW(q4ss.tcr.cd[[2]], q4ss.tcr.ho[[2]])
# [1] 0.09873425
testW(q4ss.tcr.cd[[3]], q4ss.tcr.ho[[3]])
testW(q4ss.tcr.cd[[4]], q4ss.tcr.ho[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1
testW(q4ss.tcr.rb.cd[[1]], q4ss.tcr.rb.ho[[1]])
testW(q4ss.tcr.rb.cd[[2]], q4ss.tcr.rb.ho[[2]])
# [1] 0.05421148
testW(q4ss.tcr.rb.cd[[3]], q4ss.tcr.rb.ho[[3]])
testW(q4ss.tcr.rb.cd[[4]], q4ss.tcr.rb.ho[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), WT
testW(q4ss.tcr.wt.cd[[1]], q4ss.tcr.wt.ho[[1]])
testW(q4ss.tcr.wt.cd[[2]], q4ss.tcr.wt.ho[[2]])
testW(q4ss.tcr.wt.cd[[3]], q4ss.tcr.wt.ho[[3]])
testW(q4ss.tcr.wt.cd[[4]], q4ss.tcr.wt.ho[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1 vs WT
testW(q4ss.tcr.rb[[1]], q4ss.tcr.wt[[1]])
testW(q4ss.tcr.rb[[2]], q4ss.tcr.wt[[2]])
# [1] 1.2491e-06
testW(q4ss.tcr.rb[[3]], q4ss.tcr.wt[[3]])
testW(q4ss.tcr.rb[[4]], q4ss.tcr.wt[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1 vs WT (CD)
testW(q4ss.tcr.rb.cd[[1]], q4ss.tcr.wt.cd[[1]])
testW(q4ss.tcr.rb.cd[[2]], q4ss.tcr.wt.cd[[2]])
# [1] 1.939253e-05
testW(q4ss.tcr.rb.cd[[3]], q4ss.tcr.wt.cd[[3]])
testW(q4ss.tcr.rb.cd[[4]], q4ss.tcr.wt.cd[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1 vs WT (HO)
testW(q4ss.tcr.rb.ho[[1]], q4ss.tcr.wt.ho[[1]])
testW(q4ss.tcr.rb.ho[[2]], q4ss.tcr.wt.ho[[2]])
# [1] 0.009640875
testW(q4ss.tcr.rb.ho[[3]], q4ss.tcr.wt.ho[[3]])
testW(q4ss.tcr.rb.ho[[4]], q4ss.tcr.wt.ho[[4]])
####
###
##

####
###
## Inconsistent 26/11/18
overlaps1 <- intersect(rownames(ensGene.rt.tx.ca), sclc.tx.inconsist)
ensGene.rt.tx.ca.inconsist <- ensGene.rt.tx.ca[overlaps1,]
# > length(overlaps1)
# [1] 440

## Inconsistent > TCR
ensGene.input <- ensGene.rt.tx.ca.inconsist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input1 <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input1, tx.q4.fix)
q4ss.tcr.cd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.cd[[q]] <- ensGene.input1[txs.q4[[q]],]$RATIO_LOG2

ensGene.input2 <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input2, tx.q4.fix)
q4ss.tcr.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr.ho[[q]] <- ensGene.input2[txs.q4[[q]],]$RATIO_LOG2

## Inconsistent > TCD
ensGene.input <- ensGene.rt.tx.ca.inconsist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2
  
ensGene.input1 <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input1, tx.q4.fix)
q4ss.tcd.cd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd.cd[[q]] <- ensGene.input1[txs.q4[[q]],]$RATIO_LOG2

ensGene.input2 <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input2, tx.q4.fix)
q4ss.tcd.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd.ho[[q]] <- ensGene.input2[txs.q4[[q]],]$RATIO_LOG2

## Plots
idx <- 1
file.main <- getMain(rownames(asyms[[i]]))
q4 <- q4s[[i]]
mtext <- paste0("-log2(", paste(rownames(q4), collapse="/"), ")")

ensGene.input <- ensGene.rt.tx.ca.inconsist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 < 0)
ylim.tcr <- c(min(ensGene.input$RATIO_LOG2), max(ensGene.input$RATIO_LOG2))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_inconsist_TCD.pdf"))
plotQ4SS(q4ss.tcd, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_inconsist_TCD_CD.pdf"))
plotQ4SS(q4ss.tcd.cd, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_inconsist_TCD_HO.pdf"))
plotQ4SS(q4ss.tcd.ho, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)

## BTRAPS, TCR (with 0), CUT 3 (26/11/18)
testW(q4ss.tcr.cd[[1]], q4ss.tcr.ho[[1]])
testW(q4ss.tcr.cd[[2]], q4ss.tcr.ho[[2]])
# [1] 0.09873425
testW(q4ss.tcr.cd[[3]], q4ss.tcr.ho[[3]])
testW(q4ss.tcr.cd[[4]], q4ss.tcr.ho[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1
testW(q4ss.tcr.rb.cd[[1]], q4ss.tcr.rb.ho[[1]])
testW(q4ss.tcr.rb.cd[[2]], q4ss.tcr.rb.ho[[2]])
# [1] 0.05421148
testW(q4ss.tcr.rb.cd[[3]], q4ss.tcr.rb.ho[[3]])
testW(q4ss.tcr.rb.cd[[4]], q4ss.tcr.rb.ho[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), WT
testW(q4ss.tcr.wt.cd[[1]], q4ss.tcr.wt.ho[[1]])
testW(q4ss.tcr.wt.cd[[2]], q4ss.tcr.wt.ho[[2]])
testW(q4ss.tcr.wt.cd[[3]], q4ss.tcr.wt.ho[[3]])
testW(q4ss.tcr.wt.cd[[4]], q4ss.tcr.wt.ho[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1 vs WT
testW(q4ss.tcr.rb[[1]], q4ss.tcr.wt[[1]])
testW(q4ss.tcr.rb[[2]], q4ss.tcr.wt[[2]])
# [1] 
testW(q4ss.tcr.rb[[3]], q4ss.tcr.wt[[3]])
testW(q4ss.tcr.rb[[4]], q4ss.tcr.wt[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1 vs WT (CD)
testW(q4ss.tcr.rb.cd[[1]], q4ss.tcr.wt.cd[[1]])
testW(q4ss.tcr.rb.cd[[2]], q4ss.tcr.wt.cd[[2]])
# [1] 
testW(q4ss.tcr.rb.cd[[3]], q4ss.tcr.wt.cd[[3]])
testW(q4ss.tcr.rb.cd[[4]], q4ss.tcr.wt.cd[[4]])

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1 vs WT (HO)
testW(q4ss.tcr.rb.ho[[1]], q4ss.tcr.wt.ho[[1]])
testW(q4ss.tcr.rb.ho[[2]], q4ss.tcr.wt.ho[[2]])
# [1] 
testW(q4ss.tcr.rb.ho[[3]], q4ss.tcr.wt.ho[[3]])
testW(q4ss.tcr.rb.ho[[4]], q4ss.tcr.wt.ho[[4]])
####
###
##









txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

#ensGene.input <- ensGene.rt.tx.ca.consist[overlaps2,]
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
ensGene.input <- subset(ensGene.input, CD < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.ho.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
 q4ss.ho.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2


RB1
[1] 304
[1] 485
[1] 537
[1] 312

WT
[1] 573
[1] 537
[1] 352
[1] 214

wd.asym.plots <- file.path(wd.asym,  "plots/cut3/lcl/fixed/wt")
## TEST (TCD; 09/10/18)
#ensGene.input <- ensGene.tx.rt.nona.sign.ca
ensGene.input <- ensGene.tx.rt.nona.sign.ca[overlaps2,]   ## TEST LCL
#ensGene.input <- subset(ensGene.input, gene_biotype == "protein_coding")
#ensGene.input <- subset(ensGene.input, protein_coding_non_redundant == T)
ensGene.input <- subset(ensGene.input, MUT_CG >= 3)
ensGene.input <- subset(ensGene.input, MUT_GC >= 3)
ensGene.input <- subset(ensGene.input, RATIO_LOG2 < 0)
ensGene.input <- subset(ensGene.input, CD > 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.cd.tcd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.cd.tcd[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

#ensGene.input <- ensGene.tx.rt.nona.sign.ca
ensGene.input <- ensGene.tx.rt.nona.sign.ca[overlaps2,]   ## TEST LCL
#ensGene.input <- subset(ensGene.input, gene_biotype == "protein_coding")
#ensGene.input <- subset(ensGene.input, protein_coding_non_redundant == T)
ensGene.input <- subset(ensGene.input, MUT_CG >= 3)
ensGene.input <- subset(ensGene.input, MUT_GC >= 3)
ensGene.input <- subset(ensGene.input, RATIO_LOG2 < 0)
txs.q4 <- getTxQ4Fixed(ensGene.input, tx.q4.fix)
q4ss.tcd <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.tcd[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2


## TEST TO-DO (with same Q4 genes; 10/10/18)
ensGene.input <- ensGene.tx.rt.nona.sign.ca
txs.q4 <- getTxQ4(ensGene.input, NA)
#...

##
ensGene.input <- subset(ensGene.tx.rt.nona.sign.ca, CD > 0)
txs.q4 <- getTxQ4(ensGene.input, NA)
q4ss.ho <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.ho[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input <- subset(ensGene.tx.rt.nona.sign.ca, CD < 0)
ensGene.input <- subset(ensGene.input, strand < 0)
txs.q4 <- getTxQ4(ensGene.input, NA)
q4ss.ho.re <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.ho.re[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2

ensGene.input <- subset(ensGene.tx.rt.nona.sign.ca, CD < 0)
ensGene.input <- subset(ensGene.input, strand < 0)
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
txs.q4 <- getTxQ4(ensGene.input, NA)
q4ss.ho.re.tcr <- list(list(), list(), list(), list())
for (q in 1:4)
   q4ss.ho.re.tcr[[q]] <- ensGene.input[txs.q4[[q]],]$RATIO_LOG2



##
plotQ4SS <- function(q4ss, file.name, file.main, mtext, ylim, ylab, isLog10) {
   counts <- c()
   quantiles <- c()
   names <- c()
   for (q in 1:4) {
      counts    <- c(counts, q4ss[[q]])
      quantiles <- c(quantiles, mapply(x = 1:length(q4ss[[q]]), function(x) paste0("q", q)))
      names <- c(names, paste0("Q", q, "\nn=", length(q4ss[[q]])))
   }
 
   pdf(file.name, height=6, width=4)
   #par(mgp=c(3,2,0))
   if (isLog10)
      boxplot(log10(as.numeric(counts))~quantiles, ylab=ylab, xlab=c("", "Expression"), names=c("", "", "", ""), main=file.main, ylim=log10(ylim))
   else
      boxplot(as.numeric(counts)~quantiles, ylab=ylab, xlab=c("", "Expression"), names=c("", "", "", ""), main=file.main, ylim=ylim)
   mtext(mtext, cex=0.8, font=3, line=0.5)
   axis(1, at=1:4, names, line=1, lwd=0)
   dev.off()
}

file.main <- getMain(rownames(asyms[[i]]))
q4 <- q4s[[i]]
mtext <- paste0("-log2(", paste(rownames(q4), collapse="/"), ")")

ylim <- c(min(ensGene.input$RATIO_LOG2), max(ensGene.input$RATIO_LOG2))
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_CD.pdf"))
plotQ4SS(q4ss.cd, file.name, file.main, mtext, ylim=ylim, "TCR efficiency", isLog10=F)

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_HO_FW.pdf"))
plotQ4SS(q4ss.ho.fw, file.name, file.main, mtext, ylim=ylim, "TCR efficiency", isLog10=F)

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_HO_FW_TCR.pdf"))
plotQ4SS(q4ss.ho.fw.tcr, file.name, file.main, mtext, ylim=ylim, "TCR efficiency", isLog10=F)


## TEST (remove TCD AND 0; 09/10/18)
idx <- 1
file.main <- getMain(rownames(asyms[[i]]))
q4 <- q4s[[i]]
mtext <- paste0("-log2(", paste(rownames(q4), collapse="/"), ")")

ensGene.input <- ensGene.rt.tx.ca.consist
ensGene.input <- subset(ensGene.input, RATIO_LOG2 >= 0)
ylim.tcr <- c(min(ensGene.input$RATIO_LOG2), max(ensGene.input$RATIO_LOG2))

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCR_HO.pdf"))
plotQ4SS(q4ss.ho.tcr, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCR_CD.pdf"))
plotQ4SS(q4ss.cd.tcr, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCR.pdf"))
plotQ4SS(q4ss.tcr, file.name, file.main, mtext, ylim=ylim.tcr, "TCR efficiency", isLog10=F)

## BTRAPS, TCR (with 0), CUT 3 (26/11/18)
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
# [1] 0.2084318
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
# [1] 0.09873425
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
# [1] 0.4202235
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])
# [1] 0.3041756

## BTRAPS, TCR (with 0), CUT 3 (15/10/18), RB1
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
# [1] 0.2084318
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
# [1] 0.09873425
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
# [1] 0.4202235
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])
# [1] 0.3041756










## TEST (remove TCD AND 0; 09/10/18)
ylim.tcd <- c(min(ensGene.input$RATIO_LOG2), max(ensGene.input$RATIO_LOG2))
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCD.pdf"))
plotQ4SS(q4ss.tcd, file.name, file.main, mtext, ylim=ylim.tcd, "TCD", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCD_CD.pdf"))
plotQ4SS(q4ss.cd.tcd, file.name, file.main, mtext, ylim=ylim.tcd, "TCD", isLog10=F)
file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_", REFS[idx], ">", ALTS[idx], "_consist_TCD_HO.pdf"))
plotQ4SS(q4ss.ho.tcd, file.name, file.main, mtext, ylim=ylim.tcd, "TCD", isLog10=F)

test1 <- subset(tx.snv.input, ensembl_gene_id == "ENSG00000101773")
nrow(test1)
w <- subset(subset(test1, REF == REFS[idx]), ALT == ALTS[idx])
nrow(w)
c  <- subset(subset(test1, REF == REFS[idx+1]), ALT == ALTS[idx+1])
nrow(c)


## LCL, TCR (with 0), CUT 3, RB1 (15/10/18)
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
# [1] 0.5748316
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
# [1] 0.7136629
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
# [1] 0.4300532
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])
# [1] 0.5186372
 
testW(q4ss.cd.tcr[[1]], q4ss.cd.tcr[[2]])
# [1] 8.723108e-10
testW(q4ss.cd.tcr[[2]], q4ss.cd.tcr[[3]])
# [1] 0.1144443
testW(q4ss.cd.tcr[[3]], q4ss.cd.tcr[[4]])
# [1] 0.3392386

testW(q4ss.ho.tcr[[1]], q4ss.ho.tcr[[2]])
# [1] 1.304904e-11
testW(q4ss.ho.tcr[[2]], q4ss.ho.tcr[[3]])
# [1] 0.002632109
testW(q4ss.ho.tcr[[3]], q4ss.ho.tcr[[4]])
# [1] 0.01416655

median(q4ss.cd.tcr[[3]])
# [1] 1.874469
median(q4ss.ho.tcr[[3]])
# [1] 1.883777
median(q4ss.cd.tcr[[4]])
# [1] 1.678072
median(q4ss.ho.tcr[[4]])
# [1] 1.584963

nrow(ensGene.input)
# [1] 1638


testW(q4ss.tcr.wt[[1]], q4ss.tcr.rb[[1]])
# [1] 0.5205019
testW(q4ss.tcr.wt[[2]], q4ss.tcr.rb[[2]])
# [1] 2.650335e-07
testW(q4ss.tcr.wt[[3]], q4ss.tcr.rb[[3]])
# [1] 0.09743082
testW(q4ss.tcr.wt[[4]], q4ss.tcr.rb[[4]])
# [1] 0.5361507



## TCR (with 0), CUT 3, RB1 (15/10/18)
testW(q4ss.tcr[[1]], q4ss.tcr[[2]])
testW(q4ss.tcr[[2]], q4ss.tcr[[3]])
testW(q4ss.tcr[[3]], q4ss.tcr[[4]])

testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])

testW(q4ss.cd.tcr[[1]], q4ss.cd.tcr[[2]])
testW(q4ss.cd.tcr[[2]], q4ss.cd.tcr[[3]])
testW(q4ss.cd.tcr[[3]], q4ss.cd.tcr[[4]])

testW(q4ss.ho.tcr[[1]], q4ss.ho.tcr[[2]])
testW(q4ss.ho.tcr[[2]], q4ss.ho.tcr[[3]])
testW(q4ss.ho.tcr[[3]], q4ss.ho.tcr[[4]])

median(q4ss.cd.tcr[[3]])
median(q4ss.ho.tcr[[3]])
median(q4ss.cd.tcr[[4]])
median(q4ss.ho.tcr[[4]])

nrow(ensGene.input)

## TCD, LCL, CUT 3 (12/10/18)
testW(q4ss.cd.tcd[[1]], q4ss.ho.tcd[[1]])
testW(q4ss.cd.tcd[[2]], q4ss.ho.tcd[[2]])
testW(q4ss.cd.tcd[[3]], q4ss.ho.tcd[[3]])
testW(q4ss.cd.tcd[[4]], q4ss.ho.tcd[[4]])

testW(q4ss.tcd[[1]], q4ss.tcd[[2]])
testW(q4ss.tcd[[2]], q4ss.tcd[[3]])
testW(q4ss.tcd[[3]], q4ss.tcd[[4]])

testW(q4ss.cd.tcd[[3]], q4ss.cd.tcd[[4]])
testW(q4ss.ho.tcd[[3]], q4ss.ho.tcd[[4]])

nrow(ensGene.input)



## Remove TCD (with 0), CUT 3 (12/10/18)
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
[1] 0.07644623
> testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
[1] 0.0824414
> testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
[1] 0.8498209
> testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])
[1] 0.9416512
>  
 > testW(q4ss.tcr[[3]], q4ss.tcr[[4]])
[1] 0.6827462
> testW(q4ss.cd.tcr[[3]], q4ss.cd.tcr[[4]])
[1] 0.5119256
> testW(q4ss.ho.tcr[[3]], q4ss.ho.tcr[[4]])
[1] 0.7173703
>  
 > median(q4ss.cd.tcr[[3]])
[1] 1.736966
> median(q4ss.ho.tcr[[3]])
[1] 1.807355
> median(q4ss.cd.tcr[[4]])
[1] 1.736966
> median(q4ss.ho.tcr[[4]])
[1] 1.736966
> 
 > nrow(ensGene.input)
[1] 3660
> length(which(ensGene.input$RATIO_LOG2 == 0))
[1] 0

## Remove TCD (with 0), CUT 3 (11/10/18)
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])

testW(q4ss.cd.tcr[[3]], q4ss.cd.tcr[[4]])
testW(q4ss.ho.tcr[[3]], q4ss.ho.tcr[[4]])

median(q4ss.cd.tcr[[3]])
median(q4ss.ho.tcr[[3]])
median(q4ss.cd.tcr[[4]])
median(q4ss.ho.tcr[[4]])

nrow(ensGene.input)
length(which(ensGene.input$RATIO_LOG2 == 0))

## Remove TCD AND 0 (10/10/18)
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
# [1] 0.5699394
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
# [1] 0.2834457
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
# [1] 0.1916649
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])
# [1] 0.9355701

testW(q4ss.cd.tcr[[3]], q4ss.cd.tcr[[4]])
# [1] 0.004553327
testW(q4ss.ho.tcr[[3]], q4ss.ho.tcr[[4]])
# [1] 0.1420061

median(q4ss.cd.tcr[[3]])
# [1] 2.263034
median(q4ss.ho.tcr[[3]])
# [1] 2.169925
median(q4ss.cd.tcr[[4]])
# [1] 2
median(q4ss.ho.tcr[[4]])
# [1] 2

## Remove TCD but with 0
testW(q4ss.cd.tcr[[1]], q4ss.ho.tcr[[1]])
# [1] 0.05076076
testW(q4ss.cd.tcr[[2]], q4ss.ho.tcr[[2]])
# [1] 0.2485787
testW(q4ss.cd.tcr[[3]], q4ss.ho.tcr[[3]])
# [1] 0.1223452
testW(q4ss.cd.tcr[[4]], q4ss.ho.tcr[[4]])
# [1] 0.3311788

testW(q4ss.cd.tcr[[3]], q4ss.cd.tcr[[4]])
# [1] 0.0003360581
testW(q4ss.ho.tcr[[3]], q4ss.ho.tcr[[4]])
# [1] 0.003844611

median(q4ss.cd.tcr[[3]])
# [1] 2
median(q4ss.ho.tcr[[3]])
# [1] 2
median(q4ss.cd.tcr[[4]])
# [1] 1.925999
median(q4ss.ho.tcr[[4]])
# [1] 1.847997

## TCR + TCD
testW(q4ss.cd[[1]], q4ss.ho[[1]])
# [1] 0.2046874
testW(q4ss.cd[[2]], q4ss.ho[[2]])
# [1] 0.8209866
testW(q4ss.cd[[3]], q4ss.ho[[3]])
# [1] 0.2956054
testW(q4ss.cd[[4]], q4ss.ho[[4]])
# [1] 0.07701628

median(q4ss.cd[[4]])
# [1] 1.736966
median(q4ss.ho[[4]])
# [1] 1.678072

testW(q4ss.cd[[3]], q4ss.cd[[4]])
# [1] 0.01156121
testW(q4ss.ho[[3]], q4ss.ho[[4]])
# [1] 0.001568392





testW(q4ss.cd.fw[[1]], q4ss.cd.re[[1]])
# [1] 0.525491
testW(q4ss.cd.fw[[2]], q4ss.cd.re[[2]])
# [1] 0.003414948
testW(q4ss.cd.fw[[3]], q4ss.cd.re[[3]])
# [1] 0.8975287
testW(q4ss.cd.fw[[4]], q4ss.cd.re[[4]])
# [1] 0.4087991

testW(q4ss.ho.fw[[1]], q4ss.ho.re[[1]])
# [1] 0.7474507
testW(q4ss.ho.fw[[2]], q4ss.ho.re[[2]])
# [1] 0.1961129
testW(q4ss.ho.fw[[3]], q4ss.ho.re[[3]])
# [1] 0.1298423
testW(q4ss.ho.fw[[4]], q4ss.ho.re[[4]])
# [1] 0.2926868

testW(q4ss.cd.fw.tcd[[1]], q4ss.cd.re.tcd[[1]])
testW(q4ss.cd.fw.tcd[[2]], q4ss.cd.re.tcd[[2]])
testW(q4ss.cd.fw.tcd[[3]], q4ss.cd.re.tcd[[3]])
testW(q4ss.cd.fw.tcd[[4]], q4ss.cd.re.tcd[[4]])

testW(q4ss.ho.fw.tcd[[1]], q4ss.ho.re.tcd[[1]])
testW(q4ss.ho.fw.tcd[[2]], q4ss.ho.re.tcd[[2]])
testW(q4ss.ho.fw.tcd[[3]], q4ss.ho.re.tcd[[3]])
testW(q4ss.ho.fw.tcd[[4]], q4ss.ho.re.tcd[[4]])

testW(q4ss.ho.re.tcd[[1]], q4ss.cd.re.tcd[[1]])
testW(q4ss.ho.re.tcd[[2]], q4ss.cd.re.tcd[[2]])
testW(q4ss.ho.re.tcd[[3]], q4ss.cd.re.tcd[[3]])
testW(q4ss.ho.re.tcd[[4]], q4ss.cd.re.tcd[[4]])







###
##
outliers.idx <- which(lens >= 2000000)
outliers <- gens[outliers.idx]
# > ensGene[outliers,]
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name protein_coding_non_redundant
# ENSG00000174469 ENSG00000174469            chr7      1      145813453    148118090 protein_coding            CNTNAP2                         TRUE
# ENSG00000153707 ENSG00000153707            chr9     -1        8314246     10612723 protein_coding              PTPRD                         TRUE
# ENSG00000185565 ENSG00000185565            chr3     -1      115521235    117716095 protein_coding              LSAMP                         TRUE

###
##
getTitleRTTx <- function(headon, strand) {
   if (headon == "ho")
      if (strand == "fw")
         return("RT(L)-Tx(+)")
      else
         return("RT(R)-Tx(-)") 
   else
      if (strand == "fw")
         return("RT(R)-Tx(+)")
      else
         return("RT(L)-Tx(-)") 
}

plotSerenity <- function(q4s.mut, q4s.len, file.name, file.main, mtext, ylim, xlim, ylab, xlab) {
   for (q in 1:4) {
      file.name2 <- paste0(file.name, "_q", q, ".pdf")
      lens <- as.numeric(q4s.len[[q]]) / 1000000
      muts <- as.numeric(q4s.mut[[q]]) * lens
      
      pdf(file.name2, height=6, width=6)
      plot(muts~lens, ylab=ylab, xlab=xlab, main=file.main, ylim=NULL, xlim=NULL)
      mtext(mtext, cex=0.8, font=3, line=0.5)
      dev.off()  
   }
}

###
##
for (i in 1:length(idxs))
   for (rt in 1:length(headons))
      for (st in 1:length(strands))
         q4s.s6.rt.st.rate[[i]][[rt]][[st]] <- getRateQ4(q4s.s6.rt.st.mut[[i]][[rt]][[st]], q4s.s6.rt.st.len[[i]][[rt]][[st]])

for (i in 1:length(idxs))
   for (rt in 1:length(headons))
      for (st in 1:length(strands))
         q4s.s6.rt.st.ratio[[i]][[rt]][[st]] <- getLog2RatioQ4(q4s.s6.rt.st.mut.cg[[i]][[rt]][[st]], q4s.s6.rt.st.mut.gc[[i]][[rt]][[st]])
ratios <- getAll(1, q4s.s6.rt.st.ratio, as.numeric=T)

for (i in 1:length(idxs))
   for (rt in 1:length(headons))
      for (st in 1:length(strands))
         q4s.s6.rt.st.ratio.permb[[i]][[rt]][[st]] <- getLog2RatioPerMbQ4(q4s.s6.rt.st.ratio[[i]][[rt]][[st]], q4s.s6.rt.st.len[[i]][[rt]][[st]])

for (i in 1:1) {
  idx <- idxs[i] 
  
  for (rt in 1:length(headons)) {
      for (st in 1:length(strands)) {
         file.main <- getTitleRTTx(headons[rt], strands[st])
         mtext <- paste0(REFS[idx], ">", ALTS[idx], "/", REFS[idx+1], ">", ALTS[idx+1])
         
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_g2g_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.g2g[[s]][[rt]][[st]], file.name, file.main, getYlim(q4s.s6.rt.st.g2g), "Gene-to-gene min dist. (log10)")
       
         file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_len_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         plotQ4S(q4s.s6.rt.st.len[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.len, 1, 1), "Gene length (log10)", isLog10=T)
         
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_mut_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.mut[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.mut, 1, 1), "Mutation (log10)", isLog10=T)
         
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_rate_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.rat[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.rat, 1, 1), "Mutation rate (log10)", isLog10=T)
         
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_ratio_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.ratio[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.ratio, 1, 1), "TCR efficiency (-log2 Cntx:Gtx/Gntx:Ctx)", isLog10=F)

         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_ratio_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], "_TCR-plus.pdf"))
         #plotQ4SPlus(q4s.s6.rt.st.ratio[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.ratio, 1, 1), "TCR (-log2 Cntx:Gtx/Gntx:Ctx)", isLog10=F, isPlus=T)

         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_ratio_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], "_TCP-minus.pdf"))
         #plotQ4SPlus(q4s.s6.rt.st.ratio[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.ratio, 1, 1), "TCD (-log2 Cntx:Gtx/Gntx:Ctx)", isLog10=F, isPlus=F)
        
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_g2g_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.g2g[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.g2g, 1, 1), "G2G min distance (log10)", isLog10=T)
      }
   }
}

###
## TCR + TCD
testW(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[1]][[4]])
# [1] 0.04540133
testW(q4s.s6.rt.st.ratio[[1]][[1]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[2]][[4]])
# [1] 0.0839085
testW(q4s.s6.rt.st.ratio[[1]][[2]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[1]][[4]])
# [1] 0.07596404
testW(q4s.s6.rt.st.ratio[[1]][[2]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[2]][[4]])
# [1] 0.05387942

###
## TCR and TCD
testWPlus <- function(a, b, isPlus) {
   if (isPlus)
      return(testW(a[which(a > 0)], b[which(b > 0)]))
   else
      return(testW(a[which(a < 0)], b[which(b < 0)]))
}

###
## Without ratio=0
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[1]][[4]], isPlus=T)
# [1] 0.004155493
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[2]][[4]], isPlus=T)
# [1] 0.5974946
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[1]][[4]], isPlus=T)
# [1] 0.2342088
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[2]][[4]], isPlus=T)
# [1] 0.7891488

##
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[1]][[4]], isPlus=F)
# [1] 0.6681504
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[2]][[4]], isPlus=F)
# [1] 0.087542
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[1]][[4]], isPlus=F)
# [1] 0.960271
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[2]][[4]], isPlus=F)
# [1] 0.001752346

# -----------------------------------------------------------------------------
# Step ?: Verify 
# Last Modified: 09/10/18
# -----------------------------------------------------------------------------
for (i in 1:length(idxs))
   for (rt in 1:length(headons))
      for (st in 1:length(strands))
         q4s.s6.rt.st.ratio[[i]][[rt]][[st]] <- getLog2RatioQ4(q4s.s6.rt.st.mut.cg[[i]][[rt]][[st]], q4s.s6.rt.st.mut.gc[[i]][[rt]][[st]])












###
##
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[1]][[4]], isPlus=T)
# [1] 0.5267008
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[2]][[4]], isPlus=T)
# [1] 0.08196797
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[1]][[4]], isPlus=T)
# [1] 0.377691
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[2]][[4]], isPlus=T)
# [1] 0.4985916

##
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[1]][[4]], isPlus=F)
# [1] 0.1955011
testWPlus(q4s.s6.rt.st.ratio[[1]][[1]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[2]][[4]], isPlus=F)
# [1] 0.5707632
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[1]][[4]], isPlus=F)
# [1] 0.003357378
testWPlus(q4s.s6.rt.st.ratio[[1]][[2]][[2]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[2]][[4]], isPlus=F)
# [1] 0.7491214


testW(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[1]][[1]][[4]])

testW(q4s.s6.rt.st.ratio[[1]][[1]][[1]][[3]], q4s.s6.rt.st.ratio[[1]][[2]][[1]][[3]])



for (i in 1:1) {
   idx <- idxs[i] 
 
   for (rt in 1:length(headons)) {
      for (st in 1:length(strands)) {
         file.main <- getTitleRTTx(headons[rt], strands[st])
         mtext <- paste0(REFS[idx], ">", ALTS[idx], "/", REFS[idx+1], ">", ALTS[idx+1])
      }
   }
}

## 12/07/18 SRC
src.rt.st.s1 <- toTable(0, 10, 8, c("Q", "RT", "ST", "TCD", "LM_SLOPE", "LM_P", "LM_Q", "SRC_RHO", "SRC_P", "SRC_Q"))
for (i in 1:1) {
   idx <- idxs[i] 
 
   row <- 1
   for (q in 3:3)
      for (rt in 1:length(headons))
         for (st in 1:length(strands)) {
            for (tcd in 1:length(tcds)) {
               src.rt.st.s1$Q[row] <- paste0("Q", q)
               src.rt.st.s1$RT[row] <- headons[rt]
               src.rt.st.s1$ST[row] <- strands[st]
               src.rt.st.s1$TCD[row] <- tcds[tcd]
            
               tcd.idx <- c()
               if (tcd == 1)
                  tcd.idx <- which(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]] > 0)
               else
                  tcd.idx <- which(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]] < 0)
             
               src.rt.st.s1$SRC_RHO[row] <- cor.test(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]][tcd.idx]), log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]][tcd.idx])), method="spearman", exact=F)[[4]]
               src.rt.st.s1$SRC_P[row]   <- cor.test(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]][tcd.idx]), log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]][tcd.idx])), method="spearman", exact=F)[[3]]
         
               fit <- lm(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]][tcd.idx])~log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]][tcd.idx])))
               src.rt.st.s1$LM_P[row]     <- summary(fit)[[4]][2,4]
               src.rt.st.s1$LM_SLOPE[row] <- as.numeric(fit$coefficients[2])
               
               file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_ratio_permb-len_", REFS[idx], ">", ALTS[idx], "_q", q, "_", headons[rt], "_", strands[st], "_", tcds[tcd], ".pdf"))
               pdf(file.name, height=6, width=6)
               plot(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]][tcd.idx])~log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]][tcd.idx])))
               abline(fit)
               dev.off()  
            
               row <- row+1
            }
         }
}

## 12/07/18 SRC
src.rt.st.s1 <- toTable(0, 10, 4, c("Q", "RT", "ST", "TCD", "LM_SLOPE", "LM_P", "LM_Q", "SRC_RHO", "SRC_P", "SRC_Q"))
for (i in 1:1) {
 idx <- idxs[i] 
 
 row <- 1
 for (q in 3:3)
  for (rt in 1:length(headons))
   for (st in 1:length(strands)) {
     src.rt.st.s1$Q[row] <- paste0("Q", q)
     src.rt.st.s1$RT[row] <- headons[rt]
     src.rt.st.s1$ST[row] <- strands[st]
     src.rt.st.s1$TCD[row] <- tcds[tcd]
     
     src.rt.st.s1$SRC_RHO[row] <- cor.test(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]]), log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]])), method="spearman", exact=F)[[4]]
     src.rt.st.s1$SRC_P[row]   <- cor.test(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]]), log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]])), method="spearman", exact=F)[[3]]
     
     fit <- lm(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]])~log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]])))
     src.rt.st.s1$LM_P[row]     <- summary(fit)[[4]][2,4]
     src.rt.st.s1$LM_SLOPE[row] <- as.numeric(fit$coefficients[2])
     
     file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_ratio_permb-len_", REFS[idx], ">", ALTS[idx], "_q", q, "_", headons[rt], "_", strands[st], ".pdf"))
     pdf(file.name, height=6, width=6)
     plot(as.numeric(q4s.s6.rt.st.ratio[[i]][[rt]][[st]][[q]])~log10(as.numeric(q4s.s6.rt.st.g2g[[i]][[rt]][[st]][[q]])))
     abline(fit)
     dev.off()  
     
     row <- row+1
   }
}



##
file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_ratios-lens.png"))
png(file.name, height=400, width=400)
plot(ratios~lens, ylab="TCR ratio (log2)", xlab="Gene length (Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

##
file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_ratios-exps.png"))
png(file.name, height=400, width=400)
plot(ratios~exps, ylab="TCR ratio (log2)", xlab="Gene expression log2(TPM + 0.01)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()

##
file.name <- file.path(wd.asym.plots, paste0(base, "_tx_rt_ratios-tpm_CNTNAP2.png"))
png(file.name, height=400, width=400)
plot(ratios~lens, ylab="TCR ratio (log2)", xlab="Gene length (Mb)", main="SCLC")
mtext(mtext, cex=0.8, font=3, line=0.5)
dev.off()





##
for (i in 1:1) {
   idx <- idxs[i] 
 
   for (rt in 1:1) {
      for (st in 1:2) {
         file.main <- getTitleRTTx(headons[rt], strands[st])
         mtext <- paste0(REFS[idx], ">", ALTS[idx], "/", REFS[idx+1], ">", ALTS[idx+1])
         
         file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_mut-len6_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st]))
         plotSerenity(q4s.s6.rt.st.mut[[i]][[rt]][[st]], q4s.s6.rt.st.len[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.mut), xlim=getYlim(q4s.s6.rt.st.len), ylab="Mutation", xlab="Gene length (Mb)")
      }
   }
}





# -----------------------------------------------------------------------------
# Step: Remove outliers
# Last Modified: 22/06/18
# -----------------------------------------------------------------------------
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
   ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
   GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
 
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes

   ##
   CntxGtx.cg <- subset(CntxGtx, ensembl_gene_id %in% ens.asym)
      q4[1, q] <- getMutPerMbTxs(CntxGtx.cg, txs.cg)   ## ADD 22/06/18: Also implemented txs.cg in getMutPerMbTxs() now
  
      GntxCtx.gc <- subset(GntxCtx, ensembl_gene_id %in% txs)
      txs.gc <- intersect(unique(GntxCtx.gc$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(GntxCtx.gc, txs.gc)   ## ADD 22/06/18: Also implemented txs.gc in getMutPerMbTxs() now
  
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(q4)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()




 
   #mut.q4 <- getMutQ4(i, tx.snv.s6, ens.asym.rt.st.q4)
   
 
   #file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], "_", "ho", "_", "re", ".pdf"))
   #file.main <- "C>A/G>T - RT(R)-Tx(-)"
   #plotMutQ4(mut.q4, file.name, file.main)
 
   file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], "_", "ho", "_", "re_size", ".pdf"))
   file.main <- "C>A/G>T - RT(R)-Tx(-)"
   plotSizeQ4(size.q4, file.name, file.main)
}





##
headon <- T
strand <- F   ## Forward transcription (T) / Reverse transcription (F) strands

s6.rt.st.g2g <- initTableS6RTSt()
s6.rt.st.len <- initTableS6RTSt()
s6.rt.st.mut <- initTableS6RTSt()
#for (i in 1:length(idxs)) {
for (i in 1:1) {
   idx <- idxs[i]
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ensGene.tx.rt.st <- ensGene.tx.rt.nona.sign[ens.asym,]
   ens.asym.rt.st.q4 <- getTxQ4RTStrand(ensGene.tx.rt.st, headon=headon, strand=strand, tpm.gene.input.log2)

   #mut.q4 <- getMutQ4(i, tx.snv.s6, ens.asym.rt.st.q4)
   size.q4 <- getSizeQ4(ens.asym.rt.st.q4)
   
   #file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], "_", "ho", "_", "re", ".pdf"))
   #file.main <- "C>A/G>T - RT(R)-Tx(-)"
   #plotMutQ4(mut.q4, file.name, file.main)
   
   file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], "_", "ho", "_", "re_size", ".pdf"))
   file.main <- "C>A/G>T - RT(R)-Tx(-)"
   plotSizeQ4(size.q4, file.name, file.main)
}

file.name <- file.path(wd.asym.plots, paste0(base, "_asym_tx_rt_snv_q4s_collision_", REFS[idx], ">", ALTS[idx], "_", "ho", "_", "fw", ".pdf"))
file.main <- "C>A/G>T - RT(L)-Tx(+)"
plotMutQ4(muts.q4, file.name, file.main)










ens.tx.rt.input <- intersect(unique(tx.snv.input$ensembl_gene_id), unique(ensGene.tx.rt.nona.sign$ensembl_gene_id))
ensGene.tx.rt.input <- ensGene.tx.rt.nona.sign[ens.tx.rt.input,]
tx.q4.rt.ho.plus  <- getTxQ4RTStrand(ensGene.tx.rt.input, headon=T, plus=T, tpm.gene.input.log2)
tx.q4.rt.ho.minus <- getTxQ4RTStrand(ensGene.tx.rt.input, headon=T, plus=F, tpm.gene.input.log2)
tx.q4.rt.cd.plus  <- getTxQ4RTStrand(ensGene.tx.rt.input, headon=F, plus=T, tpm.gene.input.log2)
tx.q4.rt.cd.minus <- getTxQ4RTStrand(ensGene.tx.rt.input, headon=F, plus=F, tpm.gene.input.log2)

#tx.q4.rt.ho.plus.q4 <- tx.q4.rt.ho.plus[[4]]
#writeTable(ensGene.tx.rt.input[tx.q4.rt.ho.plus.q4,], file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt_ho_plus_q4.txt")), colnames=T, rownames=F, sep="\t")
#tx.q4.rt.ho.minus.q4 <- tx.q4.rt.ho.minus[[4]]
#writeTable(ensGene.tx.rt.input[tx.q4.rt.ho.minus.q4,], file.path(wd.asym.data, paste0(base, "_asym_tx_rt_snv_q4s_rt_ho_minus_q4.txt")), colnames=T, rownames=F, sep="\t")






# =============================================================================
# Inner Class  : Test gene to gene distance 
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/05/18
# =============================================================================











# =============================================================================
# Inner Class  : Test gene to gene distance 
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/05/18
# =============================================================================
tx.q4.rt.ho <- list()
for (q in 1:4)
   tx.q4.rt.ho[[q]] <- intersect(ens.tx.rt.input, tx.q4.rt.input[[q]])

tx.q4.rt.cd <- list()
for (q in 1:4)
   tx.q4.rt.cd[[q]] <- intersect(ens.tx.rt.input, tx.q4.rt.input[[q]])

tx.q4.rt <- getTxQ4(ens.tx.rt.input, tpm.gene.pcg.non.sclc.log2)
# > for (q in 1:4)
#  + print(length(intersect(ens.tx.rt.input, tx.q4.rt[[q]])))
# [1] 2373
# [1] 2373
# [1] 2372
# [1] 2373

###
##
getTSS <- function(ensGene.genes) {
   ensGene.genes$TSS <- 0
 
   for (g in 1:nrow(ensGene.genes)) {
      if (ensGene.genes$strand[g] > 0)
         ensGene.genes$TSS[g] <- ensGene.genes$start_position[g]
      else
         ensGene.genes$TSS[g] <- ensGene.genes$end_position[g]
   }
 
   return(ensGene.genes)
}

g2g.rt.collision <- list()

## All genes
g2g.rt <- list()
for (q in 1:4) {
   genes <- tx.q4.rt[[q]]
   ensGene.genes <- ensGene[genes,]
   ensGene.genes <- getTSS(ensGene.genes)
 
   g2g <- list()
   for (g in 1:length(genes)) {
      ensGene.gene <- ensGene.genes[genes[g],]
      ensGene.genes.chr <- subset(ensGene.genes, chromosome_name == ensGene.gene$chromosome_name)
      ensGene.genes.chr <- subset(ensGene.genes.chr, ensembl_gene_id != ensGene.gene$ensembl_gene_id)
    
      g2g[[g]] <- min(abs(ensGene.gene$TSS - ensGene.genes.chr$TSS))
   }
   g2g.rt[[q]] <- g2g
}
g2g.rt.collision[[1]] <- g2g.rt

## CD genes
g2g.rt <- list()
for (q in 1:4) {
   genes <- tx.q4.rt.cd[[q]]
   ensGene.genes <- ensGene[genes,]
   ensGene.genes <- getTSS(ensGene.genes)
 
   g2g <- list()
   for (g in 1:length(genes)) {
      ensGene.gene <- ensGene.genes[genes[g],]
      ensGene.genes.chr <- subset(ensGene.genes, chromosome_name == ensGene.gene$chromosome_name)
      ensGene.genes.chr <- subset(ensGene.genes.chr, ensembl_gene_id != ensGene.gene$ensembl_gene_id)
      
      g2g[[g]] <- min(abs(ensGene.gene$TSS - ensGene.genes.chr$TSS))
   }
   g2g.rt[[q]] <- g2g
}
g2g.rt.collision[[2]] <- g2g.rt

##
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_gene-to-gene.pdf"), height=2.25, width=7)
par(mfrow=c(1, 3))

for (i in 1:3) {
 g2g.rt <- g2g.rt.collision[[i]]
 
 boxplot(as.numeric(g2g.rt), ylab="Distance", main="", beside=TRUE, width=.3)  #, ylim=c(ymin, ymax))
}

dev.off()

## Moment of truth
median(as.numeric(g2g.rt.collision[[1]][[1]]))
median(as.numeric(g2g.rt.collision[[1]][[2]]))
median(as.numeric(g2g.rt.collision[[1]][[3]]))
median(as.numeric(g2g.rt.collision[[1]][[4]]))

median(as.numeric(g2g.rt.collision[[2]][[1]]))
median(as.numeric(g2g.rt.collision[[2]][[2]]))
median(as.numeric(g2g.rt.collision[[2]][[3]]))
median(as.numeric(g2g.rt.collision[[2]][[4]]))

median(as.numeric(g2g.rt.collision[[3]][[1]]))
median(as.numeric(g2g.rt.collision[[3]][[2]]))
median(as.numeric(g2g.rt.collision[[3]][[3]]))
median(as.numeric(g2g.rt.collision[[3]][[4]]))

boxplot(log10(as.numeric(g2g.rt.collision[[1]][[1]])), ylab="Distance", main="")
boxplot(log10(as.numeric(g2g.rt.collision[[1]][[2]])), ylab="Distance", main="")
boxplot(log10(as.numeric(g2g.rt.collision[[1]][[3]])), ylab="Distance", main="")
boxplot(log10(as.numeric(g2g.rt.collision[[1]][[4]])), ylab="Distance", main="")

plot(density(as.numeric(g2g.rt.collision[[1]][[1]])/1000000), ylab="Distance", main="", beside=TRUE, width=.3) 
plot(density(as.numeric(g2g.rt.collision[[1]][[2]])/1000000), ylab="Distance", main="", beside=TRUE, width=.3) 
plot(density(as.numeric(g2g.rt.collision[[1]][[3]])/1000000), ylab="Distance", main="", beside=TRUE, width=.3) 
plot(density(as.numeric(g2g.rt.collision[[1]][[4]])/1000000), ylab="Distance", main="", beside=TRUE, width=.3)

plot(density(log10(as.numeric(g2g.rt.collision[[1]][[1]]))), ylab="Distance", main="", beside=TRUE, width=.3) 
plot(density(log10(as.numeric(g2g.rt.collision[[1]][[2]]))), ylab="Distance", main="", beside=TRUE, width=.3) 
plot(density(log10(as.numeric(g2g.rt.collision[[1]][[3]]))), ylab="Distance", main="", beside=TRUE, width=.3) 
plot(density(log10(as.numeric(g2g.rt.collision[[1]][[4]]))), ylab="Distance", main="", beside=TRUE, width=.3)

###
##
xmax <- 0
ymax <- 0
for (q in 1:4) {
   d <- density(log10(as.numeric(g2g.rt.collision[[1]][[q]])))
 
   if (max(d$x) > xmax) {
      xmax <- max(d$x)
   } else if (max(d$y) > ymax) {
      ymax <- max(d$y)
   }
}

plot(density(log10(as.numeric(g2g.rt.collision[[1]][[1]]))), ylab="Distance", ylim=c(0, ymax+0.05), main="") 
lines(density(log10(as.numeric(g2g.rt.collision[[1]][[2]])))) 
lines(density(log10(as.numeric(g2g.rt.collision[[1]][[3]])))) 
lines(density(log10(as.numeric(g2g.rt.collision[[1]][[4]]))))

###
##
distances <- c()
quantiles <- c()
for (q in 1:4) {
   distances <- c(distances, as.numeric(g2g.rt.collision[[1]][[q]]))
   quantiles <- c(quantiles, mapply(x = 1:length(g2g.rt.collision[[1]][[q]]), function(x) paste0("q", q)))
}

pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_gene-to-gene.pdf"), height=5, width=5)
boxplot(log10(distances)~quantiles, ylab="Minumum gene-to-gene distance (log10)", xlab="Expression", main="SCLC")
dev.off()





boxplot(log10(as.numeric(g2g.rt.collision[[1]][[1]])), ylab="Distance", main="")
boxplot(log10(as.numeric(g2g.rt.collision[[1]][[2]])), ylab="Distance", main="")
boxplot(log10(as.numeric(g2g.rt.collision[[1]][[3]])), ylab="Distance", main="")
boxplot(log10(as.numeric(g2g.rt.collision[[1]][[4]])), ylab="Distance", main="")



# -----------------------------------------------------------------------------
# Examples: RB1
# Last Modified: 13/04/18
# -----------------------------------------------------------------------------
# > ensGene.chr.start.end
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name overlapping
# ENSG00000238086 ENSG00000238086           chr13     -1       48890999     48894611     pseudogene          PPP1R26P1          NA
# ENSG00000177197 ENSG00000177197           chr13     -1       48902220     48902700     pseudogene             PCNPP5          NA
# ENSG00000139687 ENSG00000139687           chr13      1       48877887     49056122 protein_coding                RB1        TRUE
# ENSG00000139679 ENSG00000139679           chr13     -1       48963707     49018840 protein_coding              LPAR6        TRUE

# > ensGene.tx.rt.nona.sign[c("ENSG00000139687", "ENSG00000139679"),]
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name SLOPE_START SLOPE_END RT
# ENSG00000139687 ENSG00000139687           chr13      1       48877887     49056122 protein_coding                RB1   0.1551074 0.2126124 -1
# ENSG00000139679 ENSG00000139679           chr13     -1       48963707     49018840 protein_coding              LPAR6   0.1981649 0.2111415 -1

# > nrow(subset(tx.snv, ensembl_gene_id == "ENSG00000139687"))
# [1] 102
# > nrow(subset(tx.snv, ensembl_gene_id == "ENSG00000139679"))
# [1] 17

# -----------------------------------------------------------------------------
# Examples: TP53
# Last Modified: 13/04/18
# -----------------------------------------------------------------------------
# > ensGene.chr.start.end
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name overlapping
# ENSG00000262251 ENSG00000262251           chr17     -1        7588578      7589689 sense_intronic      RP11-199F11.2          NA
# ENSG00000141510 ENSG00000141510           chr17     -1        7565097      7590856 protein_coding               TP53        TRUE
# ENSG00000141499 ENSG00000141499           chr17      1        7589389      7606820 protein_coding             WRAP53        TRUE

# > ensGene.tx.rt.nona.sign[c("ENSG00000141510", "ENSG00000141499", "ENSG00000262251"),]
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name SLOPE_START    SLOPE_END RT
# ENSG00000141510 ENSG00000141510           chr17     -1        7565097      7590856 protein_coding               TP53 -0.05183592 -0.019433001  1
# ENSG00000141499 ENSG00000141499           chr17      1        7589389      7606820 protein_coding             WRAP53 -0.02064783 -0.002023487  1
# NA                         <NA>            <NA>     NA             NA           NA           <NA>               <NA>          NA           NA NA

# > nrow(subset(tx.snv, ensembl_gene_id == "ENSG00000141510"))
# [1] 88
# > nrow(subset(tx.snv, ensembl_gene_id == "ENSG00000141499"))
# [1] 16

# -----------------------------------------------------------------------------
# Examples: MYCN and NCYM
# Last Modified: 13/04/18
# -----------------------------------------------------------------------------
# > ensGene.chr.start.end
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name overlapping
# ENSG00000233718 ENSG00000233718            chr2     -1       16061159     16082371      antisense             MYCNOS          NA
# ENSG00000134323 ENSG00000134323            chr2      1       16080686     16087129 protein_coding               MYCN       FALSE

# > ensGene.tx.rt.nona.sign[c("ENSG00000134323", "ENSG00000233718"),]
# ensembl_gene_id chromosome_name strand start_position end_position   gene_biotype external_gene_name SLOPE_START  SLOPE_END RT
# ENSG00000134323 ENSG00000134323            chr2      1       16080686     16087129 protein_coding               MYCN  -0.2278129 -0.2309268  1
# NA                         <NA>            <NA>     NA             NA           NA           <NA>               <NA>          NA         NA NA






# =============================================================================
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
#getHeadOnCollision <- function(ensGene.tx.rt, headon) {
#   if (headon)
#      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), SLOPE_START > 0), subset(subset(ensGene.tx.rt, strand == 1), SLOPE_END < 0)))
#   else
#      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), SLOPE_END < 0), subset(subset(ensGene.tx.rt, strand == 1), SLOPE_START > 0)))
#}










# -----------------------------------------------------------------------------
# Step 6: Read in all the SVs
# Last Modified: 31/01/18
# -----------------------------------------------------------------------------
sv.bp.tpm.cna <- NA
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   ## INPUT: INV
   #sv.file <- paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_unsnarl.txt")
   sv.file <- paste0(wd.ngs, "rearrangements/", sample, "_unsnarl.txt")
   sv.bp <- initBreakPoint(read.peiflyne.unsnarl.txt(sv.file, cutoff=10))   ## Filter 1: TO-DO
   #sv.bp <- subset(sv.bp, Length_BP >= 10000)                              ## Filter 2: TO-DO
 
   ## INPUT: Expression
   sv.bp.tpm <- setBreakPointTranscriptionEnsGene(sample, sv.bp,     tpm.gene.sclc.log2, ensGene.tx.rt, pair=1)
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

break1 <- sv.bp.tpm.cna[,c("BP_ENS_ID_1", "BP_iCN_1")]   ## TO-DO 06/02/18: add sample IDs!!!
break1 <- subset(break1, !is.na(BP_ENS_ID_1))
colnames(break1) <- c("ensembl_gene_id", "iCN")

break2 <- sv.bp.tpm.cna[,c("BP_ENS_ID_2", "BP_iCN_2")]
break2 <- subset(break2, !is.na(BP_ENS_ID_2))
colnames(break2) <- c("ensembl_gene_id", "iCN")

tx.bp <- rbind(break1, break2)

###
##
ens.tx.pcg.snv.bp <- intersect(unique(tx.bp$ensembl_gene_id), unique(tx.pcg.snv$ensembl_gene_id))
tx.pcg.snv.bp <- subset(tx.bp, ensembl_gene_id %in% ens.tx.pcg.snv.bp)
# > nrow(tx.pcg.snv.bp)
# [1] 10639

ens.tx.pcg.snv.indel.bp <- intersect(unique(tx.bp$ensembl_gene_id), unique(tx.pcg.snv.indel$ensembl_gene_id))
tx.pcg.snv.indel.bp <- subset(tx.bp, ensembl_gene_id %in% ens.tx.pcg.snv.indel.bp)
# > nrow(tx.pcg.snv.indel.bp)
# [1] 10303









###
##
ens.tx.pcg.snv.indel.bp <- intersect(unique(tx.bp$ensembl_gene_id), unique(tx.pcg.snv.indel$ensembl_gene_id))
tx.pcg.snv.indel.bp <- subset(tx.bp, ensembl_gene_id %in% ens.tx.pcg.snv.indel.bp)
# > nrow(tx.pcg.snv.indel.bp)
# [1] 10303

nhej <- setdiff(tx.pcg.snv.bp[,1], tx.pcg.snv.indel.bp[,1])
for (q in 1:4)
   print(length(intersect(nhej, tx.q4[[q]])))
# [1] 41
# [1] 47
# [1] 59
# [1] 60
# > length(nhej)
# [1] 207

### 
##
test <- setdiff(ens.tx.pcg.snv, ens.tx.pcg.snv.indel)
for (q in 1:4)
   print(length(intersect(test, tx.q4[[q]])))
# [1] 959
# [1] 706
# [1] 666
# [1] 946

for (q in 1:4)
   print(length(intersect(ens.tx.pcg.snv.indel.bp, tx.q4[[q]])))






# -----------------------------------------------------------------------------
# Step 7:
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
getHeadOnCollision <- function(ensGene.tx.rt, headon) {
   if (headon)
      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), SLOPE_START > 0), subset(subset(ensGene.tx.rt, strand == 1), SLOPE_END < 0)))
   else
      return(rbind(subset(subset(ensGene.tx.rt, strand == -1), SLOPE_END < 0), subset(subset(ensGene.tx.rt, strand == 1), SLOPE_START > 0)))
}

getHeadOnCollision0 <- function(ensGene.tx.rt, headon) {
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

#getLaggingStrandColour <- function(headon) {
#   if (headon)
#      return(c("salmon", "#E6E6E6"))
#   else
#      return(c("#E6E6E6", "salmon"))
#}

###
##
ensGene.tx.rt.nona <- subset(subset(ensGene.tx.rt, !is.na(SLOPE_START)), !is.na(SLOPE_END))
ensGene.tx.rt.nona$SIGN <- sign(ensGene.tx.rt.nona$SLOPE_START / ensGene.tx.rt.nona$SLOPE_END)
ensGene.tx.rt.nona.sign <- subset(ensGene.tx.rt.nona, SIGN == 1)[,-10]
# > nrow(ensGene.tx.rt.nona)
# [1] 18093
# > nrow(ensGene.tx.rt.nona.sign)
# [1] 16896

tx.input <- tx.pcg.snv.indel.bp
headon <- F   ## Replication–transcription Head-on (T) / Co-directional (F) collisions
prefix <- "cd"
ens.tx.pcg.snv.indel.bp.rt <- intersect(unique(tx.input$ensembl_gene_id), unique(ensGene.tx.rt.nona.sign$ensembl_gene_id))
ensGene.tx.rt.input <- ensGene.tx.rt.nona.sign[ens.tx.pcg.snv.indel.bp.rt,]
tx.q4.rt.input <- getTxQ4RT(tx.q4, ensGene.tx.rt.input, headon=headon)
# > length(unique(tx.input$ensembl_gene_id))
# [1] 2885
# > nrow(ensGene.tx.rt.input)
# [1] 2483

q4s <- list()
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, ".pdf"), height=4.5, width=7)
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
      q4[1, q] <- getMutPerMbTxs(tx.input, txs.cg)

      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      #GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      #q4[2, q] <- getMutPerMbTxs(bp, unique(GntxCtx$ensembl_gene_id))
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(tx.input, txs.gc)
      
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
   
   q4s[[i]] <- q4
   barplot(q4, ylab="Breaks/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))
}
dev.off()
save(q4s, file=paste0(wd.asym.data, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, ".RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                                     RT-Tx (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                     RT-Tx (Co-directional)", "", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
 
   barplot(q4, ylab="Breaks/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("#E6E6E6", "#4D4D4D"))   ##, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 8: Log2 transcription-associated SNV asymmetry (Figure 3)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "sclc_asym_tx_rt_snv_indel_bp_q4s_", prefix, "_log2_ylim0.2.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (headon)
      colnames(q4) <- c("                                     RT-Tx (Head-on)", "", "", "")
   else
      colnames(q4) <- c("                                     RT-Tx (Co-directional)", "", "", "")
   
   cols <- c()   
   for (c in 1:ncol(q4))
      if (log2(q4[1,c]/q4[2,c]) > 0)
         cols <- c(cols, "#4D4D4D")
      else
         cols <- c(cols, "#E6E6E6")    
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.2, 0.2), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=cols)
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()







