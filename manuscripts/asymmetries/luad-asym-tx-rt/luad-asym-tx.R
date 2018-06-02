# =============================================================================
# Manuscript   : The dangeous case of DNA replication
# Chapter      : Transcriptional strand asymmetry in LUAD
# Name         : manuscripts/asymmetries/cll-asym-tx.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 27/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "Asymmetry.R", "Mutation.R", "DifferentialExpression.R", "Mutation.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "LUAD"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.files <- file.path(wd.asym, "files")
wd.asym.data  <- file.path(wd.asym, "data")
wd.asym.plots <- file.path(wd.asym, "plots")
setwd(wd.asym)

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "luad_wgs_n39-5.list"), header=F, rownames=F, sep="")

# -----------------------------------------------------------------------------
# Step 1: Finding mutations locate within Ensembl genes
# Last Modified: 23/01/18
# -----------------------------------------------------------------------------
for (s in 1:length(samples)) {
   sample <- samples[s]

   vcf <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_mutcall_filtered.vcf")), pass=T, rs=F)
   vcf.gene <- getSNVinEnsGene(sample, vcf, ensGene)
   
   writeTable(vcf.gene, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.vcf.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 2: Keep only genes that are transcribed in this cancer type 
# Last Modified: 24/01/18
# -----------------------------------------------------------------------------   
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)   ## CHANGE 02/04/18

for (s in 1:length(samples)) {
   #sample <- SAMPLE 
   sample <- samples[s]
   
   vcf <- readTable(file.path(wd.asym.files, paste0(sample, "_mut.ens.vcf.gz")), header=T, rownames=F, sep="")
   ens.tx <- intersect(unique(vcf$ensembl_gene_id), rownames(tpm.gene.input))   ## All SNVs on "expressed" genes (regardless protein coding or not)
   vcf.tx <- subset(vcf, ensembl_gene_id %in% ens.tx)
   writeTable(vcf.tx, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.vcf.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Seperate SNVs
   vcf.tx.snv <- subset(vcf.tx,     REF %in% c("A", "T", "C", "G"))
   vcf.tx.snv <- subset(vcf.tx.snv, ALT %in% c("A", "T", "C", "G"))   
   writeTable(vcf.tx.snv, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.snv.vcf.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Seperate Indels
   vcf.tx$TMP <- paste0(vcf.tx$REF, vcf.tx$ALT)
   vcf.tx.indel <- vcf.tx[which(nchar(vcf.tx$TMP) != 2),]
   writeTable(vcf.tx.indel[,-16], gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.indel.vcf.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 3.1: Read in all the SNVs
# Last Modified: 02/02/18
# -----------------------------------------------------------------------------
tx.snv <- initMergedReport(F)
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   vcf <- readTable(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.snv.vcf.gz")), header=T, rownames=F, sep="")
   vcf <- subset(vcf, CHROM %in% paste0("chr", c(1:22)))   ## Keep only autosomes
   snv <- getMergedReport(sample, vcf)
   tx.snv <- rbind(tx.snv, snv)
}
save(tx.snv, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv.RData")))   ## All SNVs on "expressed" genes (regardless protein coding or not)
# > nrow(tx.snv)
# [1] 395454

###
## Build up S6 table
tx.snv.s6 <- getTableS6(tx.snv, isExon=F)
save(tx.snv.s6, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6.RData")))   ## All SNVs on "expressed" genes (regardless protein coding or not)

# -----------------------------------------------------------------------------
# Step 4.1: Keep only SNVs on expressed, "non-redundant" protein coding genes
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)   ## CHANGE 12/04/18

ens.tx.snv.input <- intersect(unique(tx.snv$ensembl_gene_id), rownames(tpm.gene.input))   ## Needed in Step 5; ADD 02/02/18
tx.snv.input <- subset(tx.snv, ensembl_gene_id %in% ens.tx.snv.input)
save(ens.tx.snv.input, tx.snv.input, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_input.RData")))
# > nrow(tx.snv.input)   ## All SNVs on expressed, "non-redundant" protein coding genes
# [1] 272763

###
## Keep only genes with at least one SNV
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

tx.q4 <- getTxQ4(ens.tx.snv.input, tpm.gene.input.log2)   ## ADD 02/04/18; Divide Q4 based on genes with at least one SNV (not all expressed genes)
save(tx.q4, file=file.path(wd.asym.data, paste0(base, "_asym_tx_q4.RData")))
# for (q in 1:4)
#    print(length(intersect(ens.tx.snv.input, tx.q4[[q]])))
# [1] 2375
# [1] 2375
# [1] 2374
# [1] 2375

# -----------------------------------------------------------------------------
# Step 5.1: SNV asymetrey (Figure S1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
asyms <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6.pdf")), height=4.5, width=7)
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
save(tx.snv.s6, asyms, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6.RData")))

###
## Refining the plot
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_ylim85.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   asym <- asyms[[i]]
   idx <- idxs[i]
   
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"), ylim=c(0, 85))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5.2: Transcription-associated SNV asymmetry (Figure S2)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
q4s <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s.pdf")), height=4.5, width=7)
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
save(tx.snv.s6, asyms, q4s, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6_q4s.RData")))

## ADD 16/03/18
#q4s.exon <- list()
#q4s.exon[[1]] <- q4s

###
## Refining the plot
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s.pdf")), height=4.5, width=7)
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
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_log2_ylim1.5.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "           Tx", "", "")
   
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-1.5, 1.5), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=getLog2Colours(q4))   ## ADD 08/03/18
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()
