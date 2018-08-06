# =============================================================================
# Manuscript   : The dangeous case of DNA replication
# Chapter      : Transcriptional strand asymmetry in LUAD
# Name         : manuscripts/asymmetries/cll-asym-tx.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/18
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
BASE <- "SCLC"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.asym  <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-tx-rt"))
wd.asym.files <- file.path(wd.asym, "files")
wd.asym.data  <- file.path(wd.asym, "data")
wd.asym.plots <- file.path(wd.asym, "plots")

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "sclc_wgs_n101.list"), header=F, rownames=F, sep="")

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
# [1] 1441303

###
## Build up S6 table
tx.snv.s6 <- getTableS6(tx.snv, isExon=F)
save(tx.snv.s6, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6.RData")))   ## All SNVs on "expressed" genes (regardless protein coding or not)

# -----------------------------------------------------------------------------
# Step 4.1: Keep only SNVs on expressed, "non-redundant" protein coding genes
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
#tpm.gene <- tpm.gene[setdiff(rownames(tpm.gene), outliers1.5),]   ## ADD 26/06/18
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)   ## CHANGE 12/04/18
#tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=F)   ## TEST 29/06/18
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

ens.tx.snv.input <- intersect(rownames(tpm.gene.input), unique(tx.snv$ensembl_gene_id))   ## Needed in Step 5; ADD 02/02/18
tx.snv.input <- subset(tx.snv, ensembl_gene_id %in% ens.tx.snv.input)
save(ens.tx.snv.input, tx.snv.input, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_input.RData")))
# > nrow(tx.snv.input)   ## All SNVs on expressed, "all" protein coding genes
# [1] 1414777
# > nrow(tx.snv.input)   ## All SNVs on expressed, "non-redundant" protein coding genes
# [1] 1021827
# > nrow(tx.snv.input)   ## After removing CNTNAP2, PTPRD, LSAMP
# [1] 992324
# > nrow(tx.snv.input)   ## After removing CNTNAP2, PTPRD, LSAMP and RBFOX1
# [1] 988794

###
## Keep only genes with at least one SNV                   ## REMOVED 20/06/18
#tx.q4 <- getTxQ4(ens.tx.snv.input, tpm.gene.input.log2)   ## ADD 02/04/18; Divide Q4 based on genes with at least one SNV (not all expressed genes)
#save(tx.q4, file=file.path(wd.asym.data, paste0(base, "_asym_tx_q4.RData")))
# for (q in 1:4)
#    print(length(intersect(ens.tx.snv.input, tx.q4[[q]])))
# [1] 2587
# [1] 2587
# [1] 2586
# [1] 2587

# -----------------------------------------------------------------------------
# Step 5.1: SNV asymetrey (Figure S1)
# Last Modified: 20/06/18
# -----------------------------------------------------------------------------
ens.asyms <- list()   ## ADD 20/06/18
asyms     <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   
   txs.onboth.fw <- intersect(unique(tx.snv.s6[[i]][[1]][[1]]$ensembl_gene_id), unique(tx.snv.s6[[i]][[2]][[1]]$ensembl_gene_id))   ## ADD 27/06/18
   txs.onboth.re <- intersect(unique(tx.snv.s6[[i]][[1]][[2]]$ensembl_gene_id), unique(tx.snv.s6[[i]][[2]][[2]]$ensembl_gene_id))   
   ens.asyms[[i]] <- intersect(rownames(tpm.gene.input), c(txs.onboth.fw, txs.onboth.re))   ## CHANGE 27/06/18; ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   
   asym <- as.matrix(toTable(0, 2, 2, c("Tx(+)", "Tx(-)")))
   rownames(asym) <- c(paste0(REFS[idx], ">", ALTS[idx]), paste0(REFS[idx+1], ">", ALTS[idx+1]))
   for (j in 1:2)
      for (k in 1:2) {
         #ens.asym   <- c(ens.asym, getMutGenes(tx.snv.s6[[i]][[j]][[k]], tpm.gene.input))   ## ADD 20/06/18
         asym[j, k] <- getMutPerMb(tx.snv.s6[[i]][[j]][[k]], tpm.gene.input[ens.asyms[[i]],])   ## Only look at SNVs on expressed genes, which is with at least one SNV on each strand
               ## E.g. getMutPerMb(tx.snv.s6[[1]][[1]][[1]])   ## E.g. C>A Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[1]][[2]])   ##      C>A Tx(-)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[1]])   ##      G>T Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[2]])   ##      G>T Tx(-)
      }

   #ens.asyms[[i]] <- unique(ens.asym)   ## BUG 27/06/18
   asyms[[i]] <- asym
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()
save(tx.snv.s6, ens.asyms, asyms, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6.RData")))

###
## Refining the plot
max(asyms[[1]])
# [1] 313.1942
# [1] 326.1996   ## After removing genes with only SNVs on one strand
# [1] 318.4453   ## After removing genes with only SNVs on one strand, and genes longer than 1.5MB (CNTNAP2, PTPRD, LSAMP and RBFOX1)
# [1] 305.9242   ## After removing genes longer than 2MB (CNTNAP2, PTPRD and LSAMP)
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_ylim318.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   asym <- asyms[[i]]
   idx <- idxs[i]
   
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"), ylim=c(0, 318.4453))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5.2: Transcription-associated SNV asymmetry (Figure S2)
# Last Modified: 22/06/18
# -----------------------------------------------------------------------------
q4s <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
   ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
   GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
   
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ens.asym.q4 <- getTxQ4(ens.asym, tpm.gene.input.log2) 
   
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- ens.asym.q4[[q]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
      
      CntxGtx.cg <- subset(CntxGtx, ensembl_gene_id %in% txs)
      txs.cg <- intersect(txs, unique(CntxGtx.cg$ensembl_gene_id))   ## SWAP txs order BUG BUG BUG 29/06/18
      q4[1, q] <- getMutPerMbTxs(CntxGtx.cg, txs.cg)   ## ADD 22/06/18: Also implemented txs.cg in getMutPerMbTxs() now

      GntxCtx.gc <- subset(GntxCtx, ensembl_gene_id %in% txs)
      txs.gc <- intersect(txs, unique(GntxCtx.gc$ensembl_gene_id))   ## SWAP txs order BUG BUG BUG 29/06/18
      q4[2, q] <- getMutPerMbTxs(GntxCtx.gc, txs.gc)   ## ADD 22/06/18: Also implemented txs.gc in getMutPerMbTxs() now
      
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
   
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(q4)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()
save(tx.snv.s6, ens.asyms, asyms, q4s, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6_q4s.RData")))

## ADD 16/03/18
q4s.rt <- list()
q4s.rt[[1]] <- q4s

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






# -----------------------------------------------------------------------------
# Test: EZH2 ~ CNTNAP2 ?
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2-EZH2.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2~EZH2, ylab="CNTNAP2", xlab="EZH2", main="Gene expression")
dev.off()

CNTNAP2.muts <- table[overlaps,]$Freq

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2-CNTNAP2.muts.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2~CNTNAP2.muts, ylab="CNTNAP2", xlab="CNTNAP2 mutation", main="Gene expression")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_EZH2-CNTNAP2.muts.pdf"))
pdf(file.name, height=6, width=6)
plot(EZH2~CNTNAP2.muts, ylab="EZH2", xlab="CNTNAP2 mutation", main="Gene expression")
dev.off()

##
CNTNAP2 <- as.numeric(tpm.gene.input.log2["ENSG00000174469", overlaps])

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2-CNTNAP2_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2.ratios~CNTNAP2, ylab="CNTNAP2 TCR ratio (log2)", xlab="CNTNAP2 log2(TPM + 0.01)", main="SCLC (n=70)")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2_mut-CNTNAP2_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2.ratios~CNTNAP2.muts, ylab="CNTNAP2 TCR ratio (log2)", xlab="CNTNAP2 mutation", main="SCLC (n=70)")
dev.off()

##
RB1 <- as.numeric(tpm.gene.input.log2["ENSG00000139687", overlaps])

file.name <- file.path(wd.asym.plots, paste0(base, "_RB1-RB1_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(RB1.ratios~RB1, ylab="RB1 TCR ratio (log2)", xlab="RB1 log2(TPM + 0.01)", main="SCLC (n=70)")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2_mut-CNTNAP2_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2.ratios~CNTNAP2.muts, ylab="CNTNAP2 TCR ratio (log2)", xlab="CNTNAP2 mutation", main="SCLC (n=70)")
dev.off()

###
## 34 R(L)-Tx(+)/Q3/TCD genes
i <- 1
ratios <- log2(as.numeric(q4s.s6.rt.st.mut.cg[[1]][[1]][[1]][[3]])/as.numeric(q4s.s6.rt.st.mut.gc[[1]][[1]][[1]][[3]]))
idx <- which(ratios > 0)
genes <- q4s.s6.rt.st.gen[[1]][[1]][[1]][[3]][idx]

CntxGtx      <- tx.snv.s6[[i]][[1]][[1]]      ## ADD 23/06/18
GntxCtx      <- tx.snv.s6[[i]][[2]][[1]]
n <- c()
for (g in 1:length(genes)) {
   gene <- genes[g]

   CntxGtx.gene <- subset(CntxGtx, ensembl_gene_id == gene)
   GntxCtx.gene <- subset(GntxCtx, ensembl_gene_id == gene)
   s1 <- rbind(CntxGtx.gene, GntxCtx.gene)
   samples <- unique(s1$SAMPLE)
   n <- c(n, length(samples))
}

de$SRC_RHO <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[4]])
de$SRC_P   <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[3]])

