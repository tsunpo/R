# =============================================================================
# Library      : Transcriptional Strand Asymmetry
# Name         : handbook-of/Asymmetry.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/01/18
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Finding mutations locate within Ensembl genes (Step 1)
# Last Modified: 29/01/18
# -----------------------------------------------------------------------------
getEnsGeneSNV <- function(pos, ensGene.chr) {
   ensGene.chr.start <- subset(ensGene.chr, pos >= start_position)
   ensGene.chr.start.end <- subset(ensGene.chr.start, pos <= end_position)
 
   return(ensGene.chr.start.end)
}

getEnsGeneSNVs <- function(vcf.chr, ensGene.chr.gene) {
   vcf.chr.start <- subset(vcf.chr, POS >= ensGene.chr.gene$start_position)
   vcf.chr.start.end <- subset(vcf.chr.start, POS <= ensGene.chr.gene$end_position)
 
   return(vcf.chr.start.end)
}

getEnsGeneTSS <- function(vcf.chr, ensGene.chr.gene, distance) {
   if (ensGene.chr.gene$strand == 1) {
      vcf.chr.start <- subset(vcf.chr, POS >= ensGene.chr.gene$start_position - distance)
      vcf.chr.start.end <- subset(vcf.chr.start, POS <= ensGene.chr.gene$start_position)
   } else {
      vcf.chr.start <- subset(vcf.chr, POS >= ensGene.chr.gene$end_position)
      vcf.chr.start.end <- subset(vcf.chr.start, POS <= ensGene.chr.gene$end_position + distance)
   }
 
   return(vcf.chr.start.end)
}

getMergedTable <- function(ensGene.chr.gene, vcf.gene.chr0) {
   ensGene.chr.genes <- ensGene.chr.gene
   if (nrow(vcf.gene.chr0) >= 2)
      for (v in 2:nrow(vcf.gene.chr0))
      ensGene.chr.genes <- rbind(ensGene.chr.genes, ensGene.chr.gene)
  
   return(cbind(ensGene.chr.genes, vcf.gene.chr0))
}

# -----------------------------------------------------------------------------
# Method: (Step 3.1)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
getMergedReport <- function(sample, vcf) {
   vcf <- cbind(toTable(sample, 1, nrow(vcf), "SAMPLE"), vcf)
   
   return(vcf[,c("SAMPLE", "CHROM", "POS", "REF", "ALT", "ensembl_gene_id", "strand")])
}

# -----------------------------------------------------------------------------
# Method: SNV asymetrey (Step 3.2)
# Last Modified: 25/01/18
# -----------------------------------------------------------------------------
getMutPerMb <- function(s1, tpm.gene) {
   txs <- intersect(unique(s1$ensembl_gene_id), rownames(tpm.gene))
   s1 <- subset(s1, ensembl_gene_id %in% txs)
   
   mb <- 0
   for (t in 1:length(txs))
      mb <- mb + abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)
 
   return(nrow(s1) / mb * 1000000)
}

getMain <- function(rownames) {
   return(paste0(rownames[1], " (", rownames[2], ")"))
}

# -----------------------------------------------------------------------------
# Method: Transcription-associated SNV asymmetry (Step 4.1)
# Last Modified: 01/02/18
# -----------------------------------------------------------------------------
getMutPerMbTx <- function(s1, tx) {
   mb <- abs(ensGene[tx,]$end_position - ensGene[tx,]$start_position)
 
   return(nrow(subset(s1, ensembl_gene_id == tx)) / mb * 1000000)
}

getMutPerMbTxs <- function(s1, txs) {
   mb <- 0
   for (t in 1:length(txs))
      mb <- mb + abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)
 
   return(nrow(subset(s1, ensembl_gene_id %in% txs)) / mb * 1000000)
}

# =============================================================================
# Inner Class  : Replicative Strand Asymmetry
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/01/18
# =============================================================================

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/02/18
# =============================================================================
#getMutPerMb0 <- function(s1) {
#   txs <- unique(s1$ensembl_gene_id)
#   mpm <- 0
#   for (t in 1:length(txs))
#      mpm <- mpm + ((nrow(subset(s1, ensembl_gene_id == txs[t])) / abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)) * 1000000)
# 
#   return(mpm / length(txs))
#}

#getMutPerMbTx0 <- function(s1, txs) {
#   mpm <- 0
#   for (t in 1:length(txs))
#      mpm <- mpm + ((nrow(subset(s1, ensembl_gene_id == txs[t])) / abs(ensGene[txs[t],]$end_position - ensGene[txs[t],]$start_position)) * 1000000)
# 
#   return(mpm / length(txs))
#}
