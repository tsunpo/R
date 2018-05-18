# =============================================================================
# Library      : Somatic Mutation
# Name         : handbook-of/Mutation.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 13/04/17
# =============================================================================

# =============================================================================
# Inner Class  : PeifLyne File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 01/05/17
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

# =============================================================================
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
