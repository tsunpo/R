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

# -----------------------------------------------------------------------------
# Method: CNV drivers
# Last Modified: 03/07/22
# -----------------------------------------------------------------------------
getGeneCNV <- function(wd.driver.data, samples) {
   colnames <- rownames(samples)
   rownames <- rownames(ensGene)
   cna.gene <- toTable(NA, length(colnames), length(rownames), colnames)
   rownames(cna.gene) <- rownames
 
   for (s in 1:nrow(samples)) {
      sample <- rownames(samples)[s]
  
      cna.dupl <- readTable(file.path(wd.driver.data, paste0(sample, "_ens.iCN.seg.gz")), header=T, rownames=F, sep="")
      cna <- as.data.frame(sort(table(cna.dupl$ensembl_gene_id), decreasing=T))
      rownames(cna) <- cna$Var1
  
      cna.1 <- subset(cna, Freq == 1)
      cna.1$CN <- mapply(x = 1:nrow(cna.1), function(x) subset(cna.dupl, ensembl_gene_id == rownames(cna.1)[x])$CN)
  
      cna.2 <- subset(cna, Freq > 1)
      cna.2$CN <- mapply(x = 1:nrow(cna.2), function(x) median(subset(cna.dupl, ensembl_gene_id == rownames(cna.2)[x])$CN))
  
      cna <- rbind(cna.1[,c(1,3)], cna.2[,c(1,3)])
      colnames(cna) <- c("ensembl_gene_id", "CN")
  
      overlaps <- intersect(rownames(cna.gene), cna$ensembl_gene_id)
      cna.gene[overlaps, s] <- cna[overlaps,]$CN
   }
   
   return(cna.gene) 
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
