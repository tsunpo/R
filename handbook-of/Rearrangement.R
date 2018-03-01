# =============================================================================
# Library      : Structural Genomic Rearrangement
# Name         : handbook-of/Rearrangement.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Re-annotate breakpoints
# Last Modified: 22/04/17
# -----------------------------------------------------------------------------
initBreakPoint <- function(sv) {
   colnames <- c("Length_Pair_1", "Length_Pair_2", "Length_BP")
   bp <- toTable(NA, length(colnames), nrow(sv), colnames)
   bp <- cbind(sv, bp)
   
   bp$Length_Pair_1 <- bp$Max_Pos_Pair1 - bp$Min_Pos_Pair1
   bp$Length_Pair_2 <- bp$Max_Pos_Pair2 - bp$Min_Pos_Pair2
   bp.intra.idx <- which(bp$Type1 == "Intra")
   bp[bp.intra.idx,]$Length_BP <- bp[bp.intra.idx,]$Exact_Pos_Pair_2 - bp[bp.intra.idx,]$Exact_Pos_Pair_1
 
   return(bp)
}

# -----------------------------------------------------------------------------
# Methods: Re-annotate breakpoints and add TPM vaules
# Last Modified: 24/04/17
# -----------------------------------------------------------------------------
getBreakPointEnsGene <- function(chr, pos, ensGene.gene.tpm) {
   ensGene.gene.chr <- subset(ensGene.gene.tpm, chromosome_name == chr)
   ensGene.gene.chr.start <- subset(ensGene.gene.chr, pos >= start_position)
   ensGene.gene.chr.start.end <- subset(ensGene.gene.chr.start, pos <= end_position)
   
   return(ensGene.gene.chr.start.end)
}

getBreakPointICN <- function(chr, pos, cna) {
   cna.chr <- subset(cna, Chr == chr)
   cna.chr.start <- subset(cna.chr, pos >= Start)
   cna.chr.start.end <- subset(cna.chr.start, pos <= End)
 
   return(cna.chr.start.end)
}

getBreakPointCytoBand <- function(chr, pos) {
   cytoBand.chr <- subset(cytoBand, chrom == chr)
   cytoBand.chr.start <- subset(cytoBand.chr, pos >= chromStart)
   cytoBand.chr.start.end <- subset(cytoBand.chr.start, pos <= chromEnd)
 
   return(cytoBand.chr.start.end)
}

setBreakPointTranscriptionEnsGene <- function(sample, bp, tpm.gene, ensGene.gene.tpm, pair) {
   colnames.expr <- c("BP_N", "BP_BAND", "BP_ENS", "BP_ENS_ID", "BP_ENS_BIO", "BP_TPM", "BP_TPM_ALL")
   colnames <- paste0(colnames.expr, "_", pair)
   bp.tpm <- cbind(bp, toTable(NA, length(colnames), nrow(bp), colnames))
   
   for (b in 1:nrow(bp.tpm)) {
      cytoBand.bp <- getBreakPointCytoBand(bp.tpm[b, paste0("_Pair", pair)], bp.tpm[b, paste0("Exact_Pos_Pair_", pair)])
      if (nrow(cytoBand.bp) != 0)
         if (nrow(cytoBand.bp) == 1)
         bp.tpm[b, paste0("BP_BAND_", pair)] <- cytoBand.bp$name2
      else
         bp.tpm[b, paste0("BP_BAND_", pair)] <- "--"
    
      ##
      ensGene.gene.tpm.bp <- getBreakPointEnsGene(bp.tpm[b, paste0("Chr_Pair", pair)], bp.tpm[b, paste0("Exact_Pos_Pair_", pair)], ensGene.gene.tpm)
      if (nrow(ensGene.gene.tpm.bp) != 0) {
         bp.tpm[b, paste0("BP_N_", pair)]       <- nrow(ensGene.gene.tpm.bp)
         bp.tpm[b, paste0("BP_ENS_", pair)]     <- paste(ensGene.gene.tpm.bp$external_gene_name, collapse=",")
         bp.tpm[b, paste0("BP_ENS_ID_", pair)]  <- paste(ensGene.gene.tpm.bp$ensembl_gene_id, collapse=",")
         bp.tpm[b, paste0("BP_ENS_BIO_", pair)] <- paste(ensGene.gene.tpm.bp$gene_biotype, collapse=",")
         
         if (nrow(ensGene.gene.tpm.bp) == 1) {
            if (is.null(tpm.gene[ensGene.gene.tpm.bp$ensembl_gene_id, sample]))
               bp.tpm[b, paste0("BP_TPM_", pair)] <- "--"
            else
               bp.tpm[b, paste0("BP_TPM_", pair)] <- as.numeric(tpm.gene[ensGene.gene.tpm.bp$ensembl_gene_id, sample])
            bp.tpm[b, paste0("BP_TPM_ALL_", pair)] <- median(as.numeric(tpm.gene[ensGene.gene.tpm.bp$ensembl_gene_id,]))
         } ##else TO-DO?
      }
   }
    
   return(bp.tpm)
}

setBreakPointCNA <- function(sample, bp.tpm, cna, pair) {
   colnames.cna <- paste0("BP_iCN_", pair)
   bp.tpm.cna <- cbind(bp.tpm, toTable(NA, length(colnames.cna), nrow(bp.tpm), colnames.cna))
 
   for (b in 1:nrow(bp.tpm.cna)) {
      cna.bp <- getBreakPointICN(bp.tpm.cna[b, paste0("Chr_Pair", pair)], bp.tpm.cna[b, paste0("Exact_Pos_Pair_", pair)], cna)
      if (nrow(cna.bp) != 0) {
         if (nrow(cna.bp) == 1)
            bp.tpm.cna[b, paste0("BP_iCN_", pair)] <- cna.bp$iCN
      } else {
         bp.tpm.cna[b, paste0("BP_iCN_", pair)] <- "--"
      }
   }
 
   return(bp.tpm.cna)
}

setAnalysisReady <- function(bp.tpm.cna, pair) {
   colnames.pos <- c("Chr_Pair", "Exact_Pos_Pair_")
   colnames.expr <- c("BP_N", "BP_BAND", "BP_ENS", "BP_ENS_ID", "BP_ENS_BIO", "BP_TPM", "BP_TPM_ALL", "BP_iCN")
   bp.tpm.cna.pair <- bp.tpm.cna[,c("Sample", paste0(colnames.pos, pair), paste0(colnames.expr, "_", pair))]
   colnames(bp.tpm.cna.pair) <- c("Sample", "CHR", "POS", colnames.expr)
   
   return(bp.tpm.cna.pair)
}

getAnalysisReady <- function(sv.bp.tpm.cna) {
   ## Check if two breakpoints are all NA
   sv.bp.tpm.cna$BP_ENS_NA <- T
   sv.bp.tpm.cna$BP_ENS_NA[intersect(which(is.na(sv.bp.tpm.cna$BP_ENS_1)), which(is.na(sv.bp.tpm.cna$BP_ENS_2)))] <- F
   sv.bp.tpm.cna.cln <- subset(sv.bp.tpm.cna, BP_ENS_NA == T)
   
   ## Check if two breakpoints are on the same gene
   #sv.bp.tpm.cna$BP_ENS_12 <- T
   #sv.bp.tpm.cna$BP_ENS_12[which(sv.bp.tpm.cna$BP_ENS_1 == sv.bp.tpm.cna$BP_ENS_2)] <- F
   #sv.bp.tpm.cna.cln <- subset(sv.bp.tpm.cna.cln, BP_ENS_12 == T)
   
   ##
   sv.bp.tpm.cna.cln.all <- rbind(setAnalysisReady(sv.bp.tpm.cna.cln, pair=1), setAnalysisReady(sv.bp.tpm.cna.cln, pair=2))
   
   ## Check if breakpoint is on two genes
   sv.bp.tpm.cna.cln.all <- subset(sv.bp.tpm.cna.cln.all, BP_N == 1)
   
   return(sv.bp.tpm.cna.cln.all[,-4])
}

getBreakPointFreqTable <- function(sv.bp.tpm.cln.pcg, cutoff) {
   sv.bp.tpm.cln.pcg.freq <- as.data.frame(table(sv.bp.tpm.cln.pcg$BP_ENS))
   sv.bp.tpm.cln.pcg.freq <- sv.bp.tpm.cln.pcg.freq[order(sv.bp.tpm.cln.pcg.freq$Freq, decreasing=T),]
   colnames(sv.bp.tpm.cln.pcg.freq) <- c("BP_ENS", "BP_N")
   sv.bp.tpm.cln.pcg.freq[,1] <- as.vector(sv.bp.tpm.cln.pcg.freq[,1])
 
   colnames.freq <- c("BP_SN", "BP_TPM", "BP_TPM_ALL", "BP_TPM_EFF", "BP_iCN", "BP_iCN_SAMPLES", "BP_SAMPLES", "BP_TPM_SAMPLES", "BP_BAND", "BP_ENS_ID", "CHR", "POS_SAMPLES")
   freq <- toTable(0, length(colnames.freq), nrow(sv.bp.tpm.cln.pcg.freq), colnames.freq)
   sv.bp.tpm.cln.pcg.freq <- cbind(sv.bp.tpm.cln.pcg.freq, freq)
   for (g in 1:nrow(sv.bp.tpm.cln.pcg.freq)) {
      gene <- sv.bp.tpm.cln.pcg.freq[g, "BP_ENS"]
      sv.bp.tpm.cln.pcg.gene <- subset(sv.bp.tpm.cln.pcg, BP_ENS == gene)
  
      sv.bp.tpm.cln.pcg.freq$BP_SN[g] <- length(unique(sv.bp.tpm.cln.pcg.gene$Sample))
      sv.bp.tpm.cln.pcg.freq$BP_TPM[g] <- median(as.numeric(sv.bp.tpm.cln.pcg.gene$BP_TPM[which(sv.bp.tpm.cln.pcg.gene$BP_TPM != "--")]))
      sv.bp.tpm.cln.pcg.freq$BP_TPM_ALL[g] <- sv.bp.tpm.cln.pcg.gene$BP_TPM_ALL[1]
      
      sv.bp.tpm.cln.pcg.freq$BP_iCN[g] <- median(as.numeric(sv.bp.tpm.cln.pcg.gene$BP_iCN))
      sv.bp.tpm.cln.pcg.freq$BP_iCN_SAMPLES[g] <- paste(sv.bp.tpm.cln.pcg.gene$BP_iCN, collapse=";")
      
      sv.bp.tpm.cln.pcg.freq$BP_BAND[g] <- sv.bp.tpm.cln.pcg.gene$BP_BAND[1]  
      sv.bp.tpm.cln.pcg.freq$BP_ENS_ID[g] <- sv.bp.tpm.cln.pcg.gene$BP_ENS_ID[1]
      sv.bp.tpm.cln.pcg.freq$CHR[g] <- sv.bp.tpm.cln.pcg.gene$CHR[1]
      sv.bp.tpm.cln.pcg.freq$POS_SAMPLES[g] <- paste(sv.bp.tpm.cln.pcg.gene$POS, collapse=";")
      
      sv.bp.tpm.cln.pcg.freq$BP_SAMPLES[g] <- paste(sv.bp.tpm.cln.pcg.gene$Sample, collapse=";")
      sv.bp.tpm.cln.pcg.freq$BP_TPM_SAMPLES[g] <- paste(sv.bp.tpm.cln.pcg.gene$BP_TPM, collapse=";")
   }
   sv.bp.tpm.cln.pcg.freq <- sv.bp.tpm.cln.pcg.freq[order(sv.bp.tpm.cln.pcg.freq$BP_SN, decreasing=T),]
   sv.bp.tpm.cln.pcg.freq <- sv.bp.tpm.cln.pcg.freq[which(!is.na(sv.bp.tpm.cln.pcg.freq$BP_TPM)),]
   
   ## Only report BP_SN >= 2
   sv.bp.tpm.cln.pcg.freq <- subset(sv.bp.tpm.cln.pcg.freq, BP_SN >= cutoff)
   sv.bp.tpm.cln.pcg.freq$BP_TPM_EFF <- sv.bp.tpm.cln.pcg.freq$BP_TPM - sv.bp.tpm.cln.pcg.freq$BP_TPM_ALL
   
   return(sv.bp.tpm.cln.pcg.freq)
}

# =============================================================================
# Inner Class  : PeifLyne FileReader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 24/04/17
# =============================================================================
read.peiflyne.iCN.seg <- function(cna.file) {
   cna <- readTable(cna.file, header=F, rownames=F, sep="")
   colnames(cna) <- c("Sample", "Chr", "Start", "End", "V5", "iCN")
 
   return(cna)
}

read.peiflyne.unsnarl.txt <- function(sv.file, cutoff) {
   sv <- readTable(sv.file, header=T, rownames=F, sep="")[,-c(25,26)]
   sv.qc <- subset(sv, N_Clipped_Normal <= cutoff)
   
   return(sv.qc)
}
