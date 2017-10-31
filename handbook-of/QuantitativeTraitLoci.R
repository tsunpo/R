# =============================================================================
# Library: Handbook of Quantitative Trait Loci
# Name: handbook-of/QuantitativeTraitLoci.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/17
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Calculate ratio of raw sequencing reads between tumour and matched normal
# Last Modified: 20/06/17
# -----------------------------------------------------------------------------
initRTEQTL <- function(bed.gc.chr, ensGene.gene.chr) {
   colnames <- c("GxE", "GENOTYPE", "EXPRESSION", "P", "FDR", "RHO")
   qtl <- toTable(NA, length(colnames), nrow(bed.gc.chr), colnames)
 
   qtl$GENOTYPE <- rownames(bed.gc.chr)
   qtl$EXPRESSION <- ensGene.gene.chr$ensembl_gene_id
   qtl$GxE <- paste0(qtl$GENOTYPE, "-", qtl$EXPRESSION)
   qtl <- cbind(qtl, bed.gc.chr)
   rownames(qtl) <- qtl$GxE
 
   qtl$chromosome_name <- y.ens$chromosome_name
   qtl$strand <- y.ens$strand
   qtl$start_position <- y.ens$start_position
   qtl$end_position   <- y.ens$end_position
   qtl$gene_biotype <- y.ens$gene_biotype
   qtl$external_gene_name <- y.ens$external_gene_name
 
   return(qtl)
}

getGenotypes <- function(ensGene.gene.chr, bed.gc.chr, window) {
   chr <- ensGene.gene.chr$chromosome_name
   strand <- ensGene.gene.chr$strand
   if (strand > 0)
      tss <- ensGene.gene.chr$start_position
   else
      tss <- ensGene.gene.chr$end_position
 
   start <- tss - window
   end   <- tss + window
   bed.gc.chr <- subset(bed.gc.chr, CHR == chr)
   bed.gc.chr.end <- subset(bed.gc.chr, END >= start)
   bed.gc.chr.end.start <- subset(bed.gc.chr.end, START <= end)   
 
   return(bed.gc.chr.end.start)
}


# =============================================================================
# Inner Class: In-house FileReader
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/17
# =============================================================================
read.rt.eqtl.txt.gz <- function(qtl.file) {
   return(readTable(qtl.file, header=T, rownames=T, sep=""))
}
 

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
