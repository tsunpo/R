# ============================================================================= #
# Library: Handbook of Cancer Evolution
# Name: handbook-of/Clonality.R (TO-DO)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 17/02/17
# ============================================================================= #

# -----------------------------------------------------------------------------
# Methods: Generate *_muts_expAF.c1 & c2
# -----------------------------------------------------------------------------
initExpAF <- function(nrow) {
   return(toTable(0, 12, nrow, c("Mut_ID", "Chr", "Position", "Wt", "Mut", "AF_obs", "Coverage", "AF_exp", "Mut_Copies", "Mut_Copies_Raw", "Is_Subclonal_CN", "iCN")));
}

initExpAFfromVCF <- function(sample, vcf.segment) { 
   expAF <- initExpAF(nrow(vcf.segment))
   
   expAF$Chr <- vcf.segment$CHROM[1]
   expAF$Position <- vcf.segment$POS
   expAF$Wt  <- vcf.segment$REF
   expAF$Mut <- vcf.segment$ALT
   expAF$Mut_ID <- paste(sample, paste(expAF$Chr, expAF$Position, sep=":"), "SNM", sep="_")
   
   ## VAFobs
   expAF$AF_obs <- mapply(v = 1:nrow(vcf.segment), function(v) obsVAF(vcf.segment$INFO[v]))
   expAF$Coverage <- mapply(v = 1:nrow(vcf.segment), function(v) coverage(vcf.segment$INFO[v]))

   return(expAF);
}

vcfGetSegment <- function(vcf, chromosome, start, end) {
   vcf.segment <- vcf[vcf$CHROM == paste("chr", chromosome, sep=""),]
   vcf.segment <- vcf.segment[vcf.segment $POS >= start,]
   vcf.segment <- vcf.segment[vcf.segment $POS <= end,]
   
   return(vcf.segment);
}

obsVAF <- function(format) {
   format <- unlist(strsplit(format, ";"))
   
   for (f in 1:length(format)) {
   	  value <- unlist(strsplit(format[f], "="))
      if (value[1] == "AF")
         return(value[2])
   }
}

coverage <- function(format) {
   format <- unlist(strsplit(format, ";"))
   
   for (f in 1:length(format)) {
   	  value <- unlist(strsplit(format[f], "="))
      if (value[1] == "DP")
         return(value[2])
   }
}

multiplicity <- function(p, CN, obsVAF) {
   return( round((( 2*(1-p) + CN*p ) * as.numeric(obsVAF)) / p, 7 ));
}

multiplicityHat <- function(majorCN, m) {
   mHat <- floor(m+0.5)
   
   if (mHat == 0) {
      return(1);
   } else if (mHat > majorCN) {
      return(majorCN);
   } else
      return(mHat);
}

multiplicityHatC2 <- function(theta1, m) {
   mHat <- 1
   mHatPlusTheta1 <- mHat + theta1$cellular_prevalence
   diff <- min(abs(m - mHat), abs(m - mHatPlusTheta1))
   
   for (mHat.tmp in 2:theta1$major_cn) {
      mHatPlusTheta1.tmp <- mHat.tmp + theta1$cellular_prevalence
   	  diff.tmp <- min(abs(m - mHat.tmp), abs(m - mHatPlusTheta1.tmp))
   	  
   	  if (diff.tmp < diff) {
   	  	  mHat <- mHat.tmp
   	  	  mHatPlusTheta1 <- mHatPlusTheta1.tmp
   	  	  diff <- diff.tmp
   	  }
   }
   
   return(mHatPlusTheta1)
}

expVAF <- function(p, CN, mHat) {
   return( round(mHat * p / ( 2*(1-p) + CN*p ), 7) );
}

expVAFC2 <- function(p, CN, mHatPlusTheta1, m, theta1) {
   mHat <- mHatPlusTheta1 - theta1$cellular_prevalence
   r <- mHat
   
   if (abs(m - mHatPlusTheta1) < abs(m - mHat))
      r <- mHatPlusTheta1
   	
   return( round(p * r / ( 2*(1-p) + CN*p ), 7) );
}

# =============================================================================
# Inner Class: ICGC FileReader
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================

# =============================================================================
# Inner Class: Sclust FileReader
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
