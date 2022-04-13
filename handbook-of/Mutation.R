# =============================================================================
# Library      : Somatic Mutation
# Name         : handbook-of/Mutation.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/04/17
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Finding mutations locate within Ensembl genes
# Last Modified: 23/01/18
# -----------------------------------------------------------------------------
getSNVinEnsGene <- function(vcf, ensGene) {
   colnames <- colnames(vcf)
   vcf.gene <- toTable("", length(colnames), 0, colnames)
   
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
   
   return(vcf.gene)
}

getSNVEnsGeneTSS <- function(vcf, ensGene, window=1000000) {
   colnames <- colnames(vcf)
   vcf.gene <- toTable("", length(colnames), 0, colnames)
 
   for (c in 1:24) {
      chr <- chrs[c]
      vcf.chr <- subset(vcf, CHROM == chr)
      ensGene.chr <- subset(ensGene, chromosome_name == chr)  
  
      vcf.gene.chr <- toTable("", length(colnames), 0, colnames)
      for (g in 1:nrow(ensGene.chr)) {
         ensGene.chr.gene <- ensGene.chr[g,]
   
         vcf.gene.chr0 <- getEnsGeneTSS(vcf.chr, ensGene.chr.gene, window)
         if (nrow(vcf.gene.chr0) != 0) {
            if (nrow(vcf.gene.chr) == 0)
               vcf.gene.chr <- getMergedTable(ensGene.chr.gene, vcf.gene.chr0)
            else
               vcf.gene.chr <- rbind(vcf.gene.chr, getMergedTable(ensGene.chr.gene, vcf.gene.chr0))
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
 
   return(vcf.gene)
}

filtered <- function(df, colnames, cutoff) {
   for (c in 1:length(colnames))
      df <- df[which(df[,colnames[c]] >= cutoff),]
   return(df)
}

# =============================================================================
# Inner Class  : PeifLyne File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/03/19; 01/05/17
# =============================================================================
read.peiflyne.mutcall.filtered.vcf <- function(vcf.file, pass, rs) {
   vcf <- readTable(vcf.file, header=F, rownames=F, sep="")
   colnames(vcf) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")
 
   if (pass)
      vcf <- subset(vcf, FILTER == "PASS")
   if (!rs)
      vcf <- subset(vcf, ID == ".")
   return(vcf)
}

read.peiflyne.muts.txt <- function(txt.file, type) {
   txt <- readTable(txt.file, header=T, rownames=F, sep="\t")

   txt <- subset(txt, Type_2 == type)
   return(txt)
}

# =============================================================================
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
