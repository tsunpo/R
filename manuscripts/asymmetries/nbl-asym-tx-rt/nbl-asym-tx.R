# =============================================================================
# Manuscript   : The dangeous case of transcription in NBL
# Chapter I    : Transcriptional strand asymmetry
# Name         : manuscripts/asymmetries/sclc-asym-tx.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 21/04/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "Asymmetry.R", "Mutation.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
BASE <- "NBL"
#wd     <- paste0("/ngs/cangen/tyang2/", BASE, "/analysis/")                   ## tyang2@gauss
#wd.ngs <- paste0("/ngs/cangen/tyang2/", BASE, "/ngs/WGS/")
wd     <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/")   ## tpyang@localhost
wd.ngs <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/ngs/WGS/")

wd.asym       <- paste0(wd, "asymmetries/", tolower(BASE), "-asym-tx/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "nbl_wgs_n57.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 1: Finding mutations locate within Ensembl genes
# Last Modified: 23/01/18
# -----------------------------------------------------------------------------
for (s in 1:length(samples)) {
   #sample <- SAMPLE 
   sample <- samples[s]
   vcf <- read.peiflyne.mutcall.filtered.vcf(paste0(wd.ngs, sample, "/", sample, "_ANALYSIS/", sample, "_mutcall_filtered.vcf"), pass=T, rs=F)

   colnames <- c(colnames(ensGene), colnames(vcf))
   vcf.gene <- toTable("", length(colnames), 0, colnames)
   #vcf.gene.10kb <- NA
   for (c in 1:24) {
      chr <- chrs[c]
      vcf.chr <- subset(vcf, CHROM == chr)
      ensGene.chr <- subset(ensGene, chromosome_name == chr)  
      
      vcf.gene.chr <- toTable("", length(colnames), 0, colnames)
      for (g in 1:nrow(ensGene.chr)) {
         ensGene.chr.gene <- ensGene.chr[g,]
         
         vcf.gene.chr0 <- getEnsGeneSNVs(vcf.chr, ensGene.chr.gene)
         if (nrow(vcf.gene.chr0) != 0) {
            if (nrow(vcf.gene.chr) == 0)
               vcf.gene.chr <- getMergedTable(ensGene.chr.gene, vcf.gene.chr0)
            else
               vcf.gene.chr <- rbind(vcf.gene.chr, getMergedTable(ensGene.chr.gene, vcf.gene.chr0))
            
            #vcf.chr.gene.10kb <- getEnsGeneTSS(vcf.chr, ensGene.chr.gene, 10000)
         }
      }
      
      if (nrow(vcf.gene.chr) != 0) {
         vcf.gene.chr <- vcf.gene.chr[order(vcf.gene.chr$POS),]
      
         if (nrow(vcf.gene) == 0)
            vcf.gene <- vcf.gene.chr
         else
            vcf.gene <- rbind(vcf.gene, vcf.gene.chr)
      }
   }
   
   writeTable(vcf.gene, gzfile(paste0(wd.asym.data, sample, "_mut.ens.vcf.gz")), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 2: Keep only genes that are transcribed in NBL 
# Last Modified: 24/01/18
# -----------------------------------------------------------------------------   
load("/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/expression/kallisto/nbl-tpm-de/data/nbl_kallisto_0.43.1_tpm.gene_r5_p47.RData")   ## tpyang@localhost
tpm.gene.nbl <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F, nonOverlappingOnly=F)   ## CHANGE 02/04/18

lookups <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/metadata/Peifer 2015/SEQC_RNAseq_ID_table.txt", header=T, rownames=2, sep="")
lookups$Patient_ID_WGS <- paste0("P", lookups$Patient_ID_WGS)
colnames(tpm.gene.nbl) <- lookups[colnames(tpm.gene.nbl),]$Patient_ID_WGS

tpm.gene.input <- tpm.gene.nbl   ## ADD 21/04/18
for (s in 1:length(samples)) {
   #sample <- SAMPLE 
   sample <- samples[s]
   
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.vcf.gz"), header=T, rownames=F, sep="")
   ens.tx <- intersect(unique(vcf$ensembl_gene_id), rownames(tpm.gene.input))   ## All SNVs on "expressed" genes (regardless protein coding or not)
   vcf.tx <- subset(vcf, ensembl_gene_id %in% ens.tx)
   writeTable(vcf.tx, gzfile(paste0(wd.asym.data, sample, "_mut.ens.tx.vcf.gz")), colnames=T, rownames=F, sep="\t")
   
   ## Seperate SNVs
   vcf.tx.snv <- subset(vcf.tx,     REF %in% c("A", "T", "C", "G"))
   vcf.tx.snv <- subset(vcf.tx.snv, ALT %in% c("A", "T", "C", "G"))   
   writeTable(vcf.tx.snv, gzfile(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz")), colnames=T, rownames=F, sep="\t")
   
   ## Seperate Indels
   vcf.tx$TMP <- paste0(vcf.tx$REF, vcf.tx$ALT)
   vcf.tx.indel <- vcf.tx[which(nchar(vcf.tx$TMP) != 2),]
   writeTable(vcf.tx.indel[,-16], gzfile(paste0(wd.asym.data, sample, "_mut.ens.tx.indel.vcf.gz")), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 3.1: Read in all the SNVs
# Last Modified: 02/02/18
# -----------------------------------------------------------------------------
tx.snv <- initMergedReport(F)
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   vcf <- readTable(paste0(wd.asym.data, sample, "_mut.ens.tx.snv.vcf.gz"), header=T, rownames=F, sep="")
   vcf <- subset(vcf, CHROM %in% paste0("chr", c(1:22)))   ## Keep only autosomes
   snv <- getMergedReport(sample, vcf)
   tx.snv <- rbind(tx.snv, snv)
}
save(tx.snv, file=paste0(wd.asym.data, "nbl_asym_tx_snv.RData"))   ## All SNVs on "expressed" genes (regardless protein coding or not)
# > nrow(tx.snv)   ## BWA-MEM (n=57)
# [1] 52468
# > nrow(tx.snv)   ## BWA (n=56)
# [1] 47207

###
## Build up S6 table
tx.snv.s6 <- getTableS6(tx.snv, isExon=F)
save(tx.snv.s6, file=paste0(wd.asym.data, "nbl_asym_tx_snv_s6.RData"))   ## All SNVs on "expressed" genes (regardless protein coding or not)

# -----------------------------------------------------------------------------
# Step 4.1: Keep only SNVs on expressed genes
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
ens.tx.snv.input <- intersect(unique(tx.snv$ensembl_gene_id), rownames(tpm.gene.input))   ## Needed in Step 5; ADD 02/02/18
tx.snv.input <- subset(tx.snv, ensembl_gene_id %in% ens.tx.snv.input)
save(ens.tx.snv.input, tx.snv.input, file=paste0(wd.asym.data, "nbl_asym_tx_snv_input.RData"))   ## All SNVs on "expressed" genes (regardless protein coding or not)
# > nrow(tx.snv)
# [1] 52468
# > nrow(tx.snv.input)   ## BWA-MEM
# [1] 52468
# > length(ens.tx.snv.input)
# [1] 10857
# > nrow(tpm.gene.input)
# [1] 18081

# > nrow(tx.snv.input)   ## BWA
# [1] 47207
# > length(ens.tx.snv.input)
# [1] 10516
# > nrow(tpm.gene.input)
# [1] 18081

###
## Compare to SCLC
# > nrow(tx.snv.sclc)
# [1] 1441303
# > length(ens.tx.snv.sclc)
# [1] 17332
# > nrow(tpm.gene.sclc)
# [1] 18440

# > length(intersect(ens.tx.snv.nbl, ens.tx.snv.sclc))   ## TO-DO
# [1] 10240

###
## Keep only genes with at least one SNV
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

tx.q4 <- getTxQ4(ens.tx.snv.input, tpm.gene.input.log2)   ## ADD 02/04/18; Divide Q4 based on genes with at least one SNV (not all expressed genes)
save(tx.q4, file=paste0(wd.asym.data, "nbl_asym_tx_q4.RData"))
# for (q in 1:4)
#    print(length(intersect(ens.tx.snv.input, tx.q4[[q]])))
# [1] 2715
# [1] 2714
# [1] 2714
# [1] 2714

# -----------------------------------------------------------------------------
# Step 5.1: SNV asymetrey (Figure S1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
asyms <- list()
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_s6.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]

   asym <- as.matrix(toTable(0, 2, 2, c("Tx(+)", "Tx(-)")))
   rownames(asym) <- c(paste0(REFS[idx], ">", ALTS[idx]), paste0(REFS[idx+1], ">", ALTS[idx+1]))
   for (j in 1:2)
      for (k in 1:2)
         asym[j, k] <- getMutPerMb(tx.snv.s6[[i]][[j]][[k]], tpm.gene.input)   ## Only look at SNVs on expressed genes
               ## E.g. getMutPerMb(tx.snv.s6[[1]][[1]][[1]])   ## E.g. C>A Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[1]][[2]])   ##      C>A Tx(-)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[1]])   ##      G>T Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[2]])   ##      G>T Tx(-)
   
   asyms[[i]] <- asym
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()
save(tx.snv.s6, asyms, file=paste0(wd.asym.data, "nbl_asym_tx_snv_s6.RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_s6_ylim10.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   asym <- asyms[[i]]
   idx <- idxs[i]
   
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"), ylim=c(0, 10))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5.2: Transcription-associated SNV asymmetry (Figure S2)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
q4s <- list()
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_s6_q4s.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4[[q]]
      
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
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(q4)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()
save(tx.snv.s6, asyms, q4s, file=paste0(wd.asym.data, "nbl_asym_tx_snv_s6_q4s.RData"))

## ADD 16/03/18
q4s.exon <- list()
q4s.exon[[1]] <- q4s

###
## Refining the plot
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_s6_q4s.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "            Tx", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
   
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))   #, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5.3: Log2 transcription-associated SNV asymmetry (Figure S3)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_s6_q4s_log2_ylim0.25.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "           Tx", "", "")
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-0.25, 0.25), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=getLog2Colours(q4))   ## ADD 08/03/18
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()
