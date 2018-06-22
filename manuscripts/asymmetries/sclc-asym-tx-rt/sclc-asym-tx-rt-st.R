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
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ensGene.tx.rt.st <- ensGene.tx.rt.nona.sign[ens.asym,]
   ens.asym.rt.st.q4 <- getTxQ4RTStrand(ensGene.tx.rt.st, headon=headon, strand=strand, tpm.gene.input.log2)

   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- ens.asym.rt.st.q4[[q]]
  
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
      CntxGtx <- subset(CntxGtx, ensembl_gene_id %in% txs)
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      q4[1, q] <- getMutPerMbTxs(CntxGtx, txs.cg)
  
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
      GntxCtx <- subset(GntxCtx, ensembl_gene_id %in% txs) 
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(GntxCtx, txs.gc)
  
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
q4s.rt.st <- q4s.rt.cd

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
getLengthQ4 <- function(tx.q4) {
   length.q4 <- list()
   
   for (q in 1:4) {
      txs <- tx.q4[[q]]
  
      lengths <- list()
      for (g in 1:length(txs))
         lengths[[g]] <- getMbTx(txs[g])
      
      length.q4[[q]] <- lengths  
   }
 
   return(length.q4)
}

getMutQ4 <- function(s1.all, txs.q4) {
   mut.q4 <- list()
   
   for (q in 1:4) {
      txs <- txs.q4[[q]]
      s1 <- subset(s1.all, ensembl_gene_id %in% txs)
      txs <- intersect(unique(s1$ensembl_gene_id), txs)

      muts <- list()
      for (g in 1:length(txs))
         muts[[g]] <- getMutPerMbTx(s1, txs[g])

      mut.q4[[q]] <- muts
   }
 
   return(mut.q4)
}

getYlim <- function(q4s) {
   ymin <- 1000000
   ymax <- 0
   
   for (s in 1:length(idxs)) {
      for (rt in 1:length(headons)) {
         for (st in 1:length(strands)) {
            for (q in 1:4) {
               y <- log10(min(as.numeric(q4s[[s]][[rt]][[st]][[q]])))
               if (y < ymin)
                  ymin <- y
            }
         }
      }
   }
   
   for (s in 1:length(idxs)) {
      for (rt in 1:length(headons)) {
         for (st in 1:length(strands)) {
            for (q in 1:4) {
               y <- log10(max(as.numeric(q4s[[s]][[rt]][[st]][[q]])))
               if (y > ymax)
                  ymax <- y
            }
         }
      }
   }
   
   return(c(ymin, ymax))
}

plotQ4S <- function(q4s, file.name, file.main, mtext, ylim, ylab) {
   counts <- c()
   quantiles <- c()
   for (q in 1:4) {
      counts    <- c(counts, as.numeric(q4s[[q]]))
      quantiles <- c(quantiles, mapply(x = 1:length(q4s[[q]]), function(x) paste0("q", q)))
   }
 
   pdf(file.name, height=6, width=4)
   names <- c("25%", "50%", "75%", "100%")
   boxplot(log10(as.numeric(counts))~quantiles, ylab=ylab, xlab="Expression", names=names, main=file.main, ylim=ylim)
   mtext(mtext, cex=0.8, font=3, line=0.5)
   dev.off()
}

##
headons <- c("ho", "cd")
strands <- c("fw", "re")
#q4s.s6.rt.st.g2g <- initTableS6RTSt()
#q4s.s6.rt.st.len <- initTableS6RTSt()
q4s.s6.rt.st.mut <- initTableS6RTSt()

for (i in 1:length(idxs)) {
   ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
   ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
   GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
   s1.all <- rbind(CntxGtx, GntxCtx)

   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ensGene.tx.rt.st <- ensGene.tx.rt.nona.sign[ens.asym,]
   
   for (rt in 1:length(headons)) {
      for (st in 1:length(strands)) {
         txs.q4 <- getTxQ4RTStrand(ensGene.tx.rt.st, headon=headons[rt], strand=strands[st], tpm.gene.input.log2)

         #q4s.s6.rt.st.g2g[[i]][[rt]][[st]] <- getG2GQ4(s1.all, txs.q4)
         #q4s.s6.rt.st.len[[i]][[rt]][[st]] <- getLengthQ4(txs.q4)
         q4s.s6.rt.st.mut[[i]][[rt]][[st]] <- getMutQ4(s1.all, txs.q4)
      }
   }
}

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

###
##
for (i in 1:6) {
  idx <- idxs[i] 
  
  for (rt in 1:length(headons)) {
      for (st in 1:length(strands)) {
         file.main <- getTitleRTTx(headons[rt], strands[st])
         mtext <- paste0(REFS[idx], ">", ALTS[idx], "/", REFS[idx+1], ">", ALTS[idx+1])
         
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_g2g_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.g2g[[s]][[rt]][[st]], file.name, file.main, getYlim(q4s.s6.rt.st.g2g), "Gene-to-gene min dist. (log10)")
       
         #file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_len_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         #plotQ4S(q4s.s6.rt.st.len[[s]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.len), "Gene length (log10)")
         
         file.name <- file.path(wd.asym.plots, paste0(base, "_s6_tx_rt_mut_", REFS[idx], ">", ALTS[idx], "_", headons[rt], "_", strands[st], ".pdf"))
         plotQ4S(q4s.s6.rt.st.mut[[i]][[rt]][[st]], file.name, file.main, mtext, ylim=getYlim(q4s.s6.rt.st.mut), "Mutation rate (log10)")
      }
   }
}
   
 
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







