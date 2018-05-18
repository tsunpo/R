# =============================================================================
# Manuscript   : The dangeous case of transcription in NBL
# Chapter III  : 
# Name         : manuscripts/asymmetries/nbl-asym-tx-exon.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 21/04/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "Asymmetry.R", "DifferentialExpression.R")
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
# Step 3.2: Determin SNVs are on exons or introns (following Step 3.1 in nbl-asym-tx.R)
# Last Modified: 15/03/18
# -----------------------------------------------------------------------------
load(file=paste0(wd.asym.data, "nbl_asym_tx_snv.RData"))   ## All SNVs on "expressed" genes, regardless protein coding or not
# > nrow(tx.snv)
# [1] 52468
#ensGene.transcript.exon <- subset(ensGene.transcript.exon, exonFrame != -1)   ## exonFrames value of -1 means that the exon is UTR
# > nrow(ensGene.transcript.exon)
# [1] 724012

tx.snv.exon <- tx.snv
tx.snv.exon$exon <- F
tx.snv.exon$exon <- mapply(x = 1:nrow(tx.snv), function(x) isSNVonExon(tx.snv$CHROM[x], tx.snv$POS[x], ensGene.transcript.exon))
save(tx.snv.exon, file=paste0(wd.asym.data, "nbl_asym_tx_snv_exon.RData"))   ## All SNVs on expressed genes

utr <- setdiff(rownames(subset(tx.snv.exon, exon == F)), rownames(subset(tx.snv.exon.utr, exon == F)))
tx.snv.exon[utr,]$exon <- -1
save(tx.snv.exon, file=paste0(wd.asym.data, "nbl_asym_tx_snv_exon.RData"))   ## All SNVs on expressed genes

# > nrow(subset(tx.snv.exon, exon == T))   ## But still with SNVs on UTRs (exonFrame == -1)
# [1] 2631
# > nrow(subset(tx.snv.exon, exon == F))
# [1] 48241
# nrow(subset(tx.snv.exon, exon == -1))
# [1] 1596

# > length(unique(tx.snv.exon$ensembl_gene_id))
# [1] 10516
# > length(unique(subset(tx.snv.exon, exon == T)$ensembl_gene_id))
# [1] 2242
# > length(unique(subset(tx.snv.exon, exon == F)$ensembl_gene_id))
# [1] 9827
# > length(unique(subset(tx.snv.exon, exon == -1)$ensembl_gene_id))
# [1] 1331
# > length(unique(subset(tx.snv.exon, exon == T)$SAMPLE))
# [1] 57

# > genes.utr <- unique(subset(tx.snv.exon, exon == -1)$ensembl_gene_id)
# > genes.exon <- unique(subset(tx.snv.exon, exon == T)$ensembl_gene_id)
# > genes.intron <- unique(subset(tx.snv.exon, exon == F)$ensembl_gene_id)
# > length(intersect(genes.utr, genes.exon))
# [1] 210
# > length(intersect(genes.utr, genes.intron))
# [1] 924
# > length(intersect(intersect(genes.utr, genes.intron), genes.exon))
# [1] 167
# > genes.all <- intersect(intersect(genes.utr, genes.intron), genes.exon)
# > for (q in 1:4)
#  + print(length(intersect(genes.all, tx.q4[[q]])))
# [1] 32
# [1] 48
# [1] 38
# [1] 49
# > for (q in 1:4)
#  + print(length(intersect(genes.utr, tx.q4[[q]])))
# [1] 303
# [1] 317
# [1] 323
# [1] 388
# > for (q in 1:4)
#  + print(length(intersect(genes.intron, tx.q4[[q]])))
# [1] 2432
# [1] 2511
# [1] 2469
# [1] 2415
# > for (q in 1:4)
#  + print(length(intersect(genes.exon, tx.q4[[q]])))
# [1] 587
# [1] 583
# [1] 512
# [1] 560
# > cor.test(c(303, 317, 323, 388), c(25, 50, 75, 100), method="spearman", exact=F)
# 
# Spearman's rank correlation rho
# 
# data:  c(303, 317, 323, 388) and c(25, 50, 75, 100)
# S = 0, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 1 

# > cor.test(c(2432, 2511, 2469, 2415), c(25, 50, 75, 100), method="spearman", exact=F)
# 
# Spearman's rank correlation rho
# 
# data:  c(2432, 2511, 2469, 2415) and c(25, 50, 75, 100)
# S = 14, p-value = 0.6
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#  rho 
# -0.4 









###
## Build up S6 table (with exon/intron information)
tx.snv.exon.s6 <- getTableS6(tx.snv.exon, T)
save(tx.snv.exon.s6, file=paste0(wd.asym.data, "nbl_asym_tx_snv_exon_s6.RData"))   ## All SNVs on "expressed" genes, regardless protein coding or not

# -----------------------------------------------------------------------------
# Step 4.2: Divide SNVs into two groups (exons and introns)
# Last Modified: 16/03/18
# -----------------------------------------------------------------------------
e <- F
prefix <- "intron"
tx.snv.input <- subset(tx.snv.exon, exon == e)

ens.tx.snv.input <- unique(tx.snv.input$ensembl_gene_id)
# > length(ens.tx.snv.input)   ## Exon
# [1] 2115
# > length(ens.tx.snv.input)   ## Intron
# [1] 9899

tx.q4.input <- getTxQ4(ens.tx.snv.input, tpm.gene.nbl.log2)
# > for (q in 1:4)   ## Exons
#  +  print(length(intersect(ens.tx.snv.input, tx.q4.input[[q]])))
# [1] 529
# [1] 529
# [1] 528
# [1] 529
# 
# > for (q in 1:4)   ## Introns
#  +    print(length(intersect(ens.snv.input, tx.q4.input[[q]])))
# [1] 2475
# [1] 2475
# [1] 2474
# [1] 2475

# -----------------------------------------------------------------------------
# Step 8.1: SNV asymetrey (Exon/Intron)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
q4s <- list()
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
 
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- tx.q4.input[[q]]
  
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.exon.s6[[i]][[1]][[1]],          tx.snv.exon.s6[[i]][[2]][[2]])
      CntxGtx.e <- subset(CntxGtx, exon == e)
      txs.cg <- intersect(unique(CntxGtx.e$ensembl_gene_id), txs)
      q4[1, q] <- getMutPerMbTxs(CntxGtx.e, txs.cg)
      
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.exon.s6[[i]][[2]][[1]],          tx.snv.exon.s6[[i]][[1]][[2]])
      GntxCtx.e <- subset(GntxCtx, exon == e)
      txs.gc <- intersect(unique(GntxCtx.e$ensembl_gene_id), txs)
      q4[2, q] <- getMutPerMbTxs(GntxCtx.e, txs.gc)
        
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
 
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=paste(rownames(q4), collapse="/"), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()
save(q4s, file=paste0(wd.asym.data, "nbl_asym_tx_snv_q4s_", prefix, ".RData"))

## ADD 16/03/18
#q4s.exon[[2]] <- q4s
q4s.exon[[3]] <- q4s
save(q4s.exon, file=paste0(wd.asym.data, "nbl_asym_tx_snv_q4s_exon.RData"))

###
## Refining the plot
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_q4s_", prefix, ".pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (e)
      colnames(q4) <- c("                                     Tx (Exon)", "", "", "")
   else
      colnames(q4) <- c("                                     Tx (Intron)", "", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
 
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))   ##, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 8.2: Log2 transcription-associated SNV asymmetry (Figure S5)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(paste0(wd.asym.plots, "nbl_asym_tx_snv_q4s_", prefix, "_log2_ylim1.25.pdf"), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   if (e)
      colnames(q4) <- c("                                    Tx (Exon)", "", "", "")
   else
      colnames(q4) <- c("                                    Tx (Intron)", "", "", "")
 
   barplot(log2(q4[1,]/q4[2,]), ylab="Asymmetry", ylim=c(-1.25, 1.25), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=getLog2Colours(q4))
   mtext(paste0("log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 9: Find genes with intronic T-to-G transversions, which might result in splice defect
# Last Modified: 21/03/18
# -----------------------------------------------------------------------------
i <- 6
genes <- c()
snvs <- NA
#for (i in 1:length(idxs)) {
   #idx <- idxs[i]
   for (q in 4:4) {
      txs <- tx.q4.input[[q]]
 
      ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
      CntxGtx <- rbind(tx.snv.exon.s6[[i]][[1]][[1]],          tx.snv.exon.s6[[i]][[2]][[2]])
      txs.cg <- intersect(unique(CntxGtx$ensembl_gene_id), txs)
      CntxGtx.q <- subset(CntxGtx, ensembl_gene_id %in% txs.cg)
 
      ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
      GntxCtx <- rbind(tx.snv.exon.s6[[i]][[2]][[1]],          tx.snv.exon.s6[[i]][[1]][[2]])
      txs.gc <- intersect(unique(GntxCtx$ensembl_gene_id), txs)
      GntxCtx.q <- subset(GntxCtx, ensembl_gene_id %in% txs.gc)
 
      genes <- c(genes, txs.cg, txs.gc)
   
      if (is.na(snvs)[1])
         snvs <- rbind(CntxGtx.q, GntxCtx.q)
      else
         snvs <- rbind(snvs, CntxGtx.q, GntxCtx.q)
   }
#}
genes <- unique(genes)
writeTable(genes, paste0(wd.asym.data, "genes_T>G_Q4.txt"), colnames=F, rownames=F, sep="\t")

genes.table <- data.frame(table(snvs$ensembl_gene_id))
genes.table <- genes.table[order(genes.table$Freq, decreasing=T),]
genes.table$Var1 <- as.vector(genes.table$Var1)

genes.f2 <- subset(genes.table, Freq >= 2)$Var1
writeTable(genes.f2, paste0(wd.asym.data, "genes_T>G_Q4_F2.txt"), colnames=F, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Replace Ensembl Gene IDs to gene name in Reactome results
# -----------------------------------------------------------------------------
wd.reactome <- "/Users/tpyang/Work/uni-koeln/tyang2/NBL/analysis/asymmetries/nbl-wgs-asym/pathway/genes_T>G_Q4/"
list <- readTable(paste0(wd.reactome, "genes_T>G_Q4.txt"), header=F, rownames=F, sep="\t")
list <- ensGene[list, c("ensembl_gene_id",	"external_gene_name")]

reactome <- read.csv(paste0(wd.reactome, "result.csv"))
colnames(reactome) <- gsub("X.", "", colnames(reactome))
reactome$Submitted.entities.found <- as.vector(reactome$Submitted.entities.found)
for (r in 1:nrow(reactome)) {
   ids <- as.vector(reactome$Submitted.entities.found[r])
   ids <- unlist(strsplit(ids, ";"))
 
   for (i in 1:length(ids))
      if (nrow(list[ids[i],]) != 0)
         ids[i] <- list[ids[i],]$external_gene_name
 
   reactome$Submitted.entities.found[r] <- paste(ids, collapse=";")
}
writeTable(reactome, paste0(wd.reactome, "result.tsv"), colnames=T, rownames=F, sep="\t")

##
file.de <- paste0(wd.reactome, "genes_T>G_Q4.pdf")
pdf(file.de, height=7, width=7.5)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,27,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()
















##
file.de <- paste0(wd.reactome, "genes_intron_T>G_Tntx_Q4.pdf")
pdf(file.de, height=7, width=7.5)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,27,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

##
file.de <- paste0(wd.reactome, "genes_exon_T>G_Q4.pdf")
pdf(file.de, height=7, width=8.7)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,33,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

##
file.de <- paste0(wd.reactome, "genes_exon_T>G_Ttx_Q4.pdf")
pdf(file.de, height=7, width=8.9)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,34,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()

##
file.de <- paste0(wd.reactome, "genes_intron_Q4.pdf")
pdf(file.de, height=7, width=7.5)
par(las=1)             ## make label text perpendicular to axis
par(mar=c(5,27,4,2))   ## increase y-axis margin
barplot(-log10(reactome[20:1,]$Entities.pValue), names.arg=reactome[20:1,]$Pathway.name, xlab="Significance (-log10 P-value)", horiz=T)
dev.off()



# -----------------------------------------------------------------------------
# Recruitment of NuMA to mitotic centrosomes (n=9)
# -----------------------------------------------------------------------------
genes <- c("ENSG00000155093")   ## PTPRN2
test <- unique(as.vector(subset(tx.pcg.snv, ensembl_gene_id %in% genes)$SAMPLE))




genes <- c("ENSG00000234745")   ##
genes <- c("ENSG00000155366")   ## RHOC

subset(tx.input, ensembl_gene_id %in% genes)

## 16 genes (9 samples)
genes <- c("ENSG00000072864", "ENSG00000104833", "ENSG00000197102", "ENSG00000175216", "ENSG00000136861", "ENSG00000137822", "ENSG00000136811", "ENSG00000170027", "ENSG00000141551", "ENSG00000054282", "ENSG00000077380", "ENSG00000088986", "ENSG00000078674", "ENSG00000258947", "ENSG00000127914", "ENSG00000160299")

#genes <- c("ENSG00000088986", "ENSG00000104833", "ENSG00000258947", "ENSG00000197102", "ENSG00000127914", "ENSG00000160299", "ENSG00000136861", "ENSG00000137822", "ENSG00000136811")
#genes <- c("ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811")
samples <- unique(subset(tx.input, ensembl_gene_id %in% genes)$SAMPLE)

## 19 genes (13 samples)
genes <- c("ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000150093", "ENSG00000182871", "ENSG00000142798", "ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000063660", "ENSG00000142798", "ENSG00000063660", "ENSG00000142798", "ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000149182", "ENSG00000104833", "ENSG00000151150", "ENSG00000197102", "ENSG00000173898", "ENSG00000063660", "ENSG00000142798", "ENSG00000063660", "ENSG00000142798", "ENSG00000063660", "ENSG00000142798", "ENSG00000104833", "ENSG00000197102", "ENSG00000127914", "ENSG00000136861", "ENSG00000136811", "ENSG00000176884", "ENSG00000127914", "ENSG00000176884", "ENSG00000127914", "ENSG00000107643", "ENSG00000177889", "ENSG00000068796", "ENSG00000150093", "ENSG00000104833", "ENSG00000197102", "ENSG00000018236", "ENSG00000110492", "ENSG00000106367", "ENSG00000068796", "ENSG00000104833", "ENSG00000197102", "ENSG00000173898")
samples <- unique(subset(tx.input, ensembl_gene_id %in% genes)$SAMPLE)
samples <- samples[-13]

# > samples <- c("P15239", "P15403", "P23122")
# > subset(subset(tx.input, ensembl_gene_id %in% genes), SAMPLE %in% samples)
# SAMPLE CHROM       POS REF ALT ensembl_gene_id strand exon
# 3297  P15239 chr21  47851564   G   A ENSG00000160299      1 TRUE
# 5450  P15403 chr19   6502292   T   G ENSG00000104833     -1 TRUE
# 34569 P23122  chr9 131219552   A   C ENSG00000136811      1 TRUE
# 35173 P23122 chr12 120934118   T   G ENSG00000088986      1 TRUE

# -----------------------------------------------------------------------------
# Pre-D.E.
# -----------------------------------------------------------------------------
pheno <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/metadata/Peifer 2015/nature14980-s1.txt", header=T, rownames=T, sep="\t")
lookups <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/metadata/Peifer 2015/Samples_NEW_NUMBERS.txt", header=T, rownames=2, sep="")
rownames(pheno) <- lookups[rownames(pheno),]$OLD

pheno$Sample_Old <- ""
pheno$Sample_Old <- rownames(pheno)
pheno <- pheno[c(1,6,2:5)]

## T>G SNVs (Recruitment of NuMA to mitotic centrosomes; n=9)
pheno$T_C_Q4 <- "no"
pheno[match(samples, rownames(pheno)),]$T_C_Q4 <- "yes"

## PTPRN2
subset(ensGene, external_gene_name == "PTPRN2")
quantile(tpm.gene.nbl.log2["ENSG00000155093",])
ptprn2 <- colnames(tpm.gene.nbl.log2["ENSG00000155093", which(tpm.gene.nbl.log2["ENSG00000155093",] >= 5.68718),])

pheno$PTPRN2 <- "low"
pheno[match(ptprn2, rownames(pheno)),]$PTPRN2 <- "high"

## PTPRN2_Freq5
pheno$PTPRN2_mut <- "wt"
pheno[match(samples, rownames(pheno)),]$PTPRN2_mut <- "mut-tg"

pheno$PTPRN2_combined <- pheno$PTPRN2
n7 <- rownames(subset(subset(pheno, PTPRN2 == "high"), PTPRN2_mut == "mut-tg"))
pheno[match(n7, rownames(pheno)),]$PTPRN2_combined <- "low"

## MYCN expression
subset(ensGene, external_gene_name == "MYCN")
quantile(tpm.gene.nbl.log2["ENSG00000134323",])
mycn <- colnames(tpm.gene.nbl.log2["ENSG00000134323", which(tpm.gene.nbl.log2["ENSG00000134323",] >= 4.890977),])

pheno$MYCN <- "low"
pheno[match(mycn, rownames(pheno)),]$MYCN <- "high"

## ACD
acd <- unique(c(rownames(subset(pheno, MYCN.amp == "yes")), rownames(subset(pheno, T_C_Q4 == "yes"))))

pheno$ACD <- "yes"
pheno[match(acd, rownames(pheno)),]$ACD <- "no"

## PTPRN2
PTPRN2 <- unique(subset(rbind(CntxGtx.q, GntxCtx.q), ensembl_gene_id == "ENSG00000155093")$SAMPLE)

pheno$PTPRN2 <- "wt"
pheno[match(PTPRN2, rownames(pheno)),]$PTPRN2 <- "mut-tg"

## PTPRN2_Freq5
pheno$PTPRN2_Freq5 <- "wt"
pheno[match(samples, rownames(pheno)),]$PTPRN2_Freq5 <- "mut-tg"

## PTPRN2_Freq3
pheno$PTPRN2_Freq3 <- "wt"
samples <- intersect(samples, rownames(pheno))
pheno[match(samples, rownames(pheno)),]$PTPRN2_Freq3 <- "mut-tg"

## SCD
scd <- unique(c(rownames(subset(pheno, MYCN.amp == "yes")), rownames(subset(pheno, T_C_Q4 == "yes")), rownames(subset(pheno, PTPRN2 == "mut-tg"))))

pheno$SCD <- "no"
pheno[match(scd, rownames(pheno)),]$SCD <- "yes"

## SCD2
scd2 <- unique(c(rownames(subset(pheno, MYCN.amp == "yes")), rownames(subset(pheno, PTPRN2 == "mut-tg"))))

pheno$SCD2 <- "no"
pheno[match(scd2, rownames(pheno)),]$SCD2 <- "yes"

## SCD3
scd3 <- unique(c(rownames(subset(pheno, MYCN.amp == "yes")), rownames(subset(pheno, PTPRN2_Freq5 == "mut-tg"))))

pheno$SCD3 <- "no"
pheno[match(scd3, rownames(pheno)),]$SCD3 <- "yes"






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







