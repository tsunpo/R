# =============================================================================
# Library: Handbook of Mutational Signatures
# Name: handbook-of/Mutation.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 17/02/17
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Mutated sample numbers per gene table
# Last Modified: 11/02/17
# -----------------------------------------------------------------------------
getMutatedSampleNumbersPerGeneTable <- function(muts.raw, samples, muts) {
   muts.num <- toTable(NA, length(muts), length(samples), muts)   ## UPDATE 11/02/17: Change to NA (Not avaialbe) from 0 (WT)
   rownames(muts.num) <- samples
 
   for (s in 1:length(samples)) {
      sample <- samples[s]
      mutations <- subset(muts.raw, PAT_ID == sample)
  
      if (nrow(mutations) != 0)   ## BUG FIX 11/02/17
        for (m in 1:length(muts))
      muts.num[s, m] <- nrow(subset(mutations, Gene_Hugo == muts[m]))
   }
   return(muts.num)
}

# -----------------------------------------------------------------------------
# Method: Mutated sample numbers per gene table
# Last Modified: 11/02/17
# -----------------------------------------------------------------------------
getMutatedSampleNumbersPerGeneTable <- function(muts.raw, samples, muts) {
   muts.num <- toTable(NA, length(muts), length(samples), muts)   ## UPDATE 11/02/17: Change to NA (Not avaialbe) from 0 (WT)
   rownames(muts.num) <- samples
 
   for (s in 1:length(samples)) {
      sample <- samples[s]
      mutations <- subset(muts.raw, PAT_ID == sample)
  
      if (nrow(mutations) != 0)   ## BUG FIX 11/02/17
         for (m in 1:length(muts))
            muts.num[s, m] <- nrow(subset(mutations, Gene_Hugo == muts[m]))
   }
   return(muts.num)
}

# -----------------------------------------------------------------------------
# Methods: Fisher's exact test
# Last Modified: 17/02/17
# -----------------------------------------------------------------------------
getFishersTable <- function(pheno, trait1, trait2) {
   pheno.nona <- removeNA(pheno, trait1)
   pheno.nona <- removeNA(pheno.nona, trait2)   
   trait1.v <- sort(unique(pheno.nona[,trait1]))
   trait2.v <- sort(unique(pheno.nona[,trait2]))
                    
   table <- toTable(0, length(trait1.v), length(trait2.v), paste0(trait1, "_", trait1.v))
   rownames(table) <- paste0(trait2, "_", trait2.v)
   for (x in 1:length(trait1.v))
      for (y in 1:length(trait2.v))
         table[y, x] <- nrow(subset0(subset0(pheno.nona, trait1, trait1.v[x]), trait2, trait2.v[y]))
   
   return(table)
}

testFishers <- function(table) {
   return(fisher.test(table)$p.value)
}

testChisquared <- function(table) {
   return(chisq.test(table)$p.value)
}

# -----------------------------------------------------------------------------
# Method: Mutational signatures
# Last Modified: 02/03/17
# -----------------------------------------------------------------------------

# =============================================================================
# Inner Class: PeifLyne FileReader
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================

# =============================================================================
# Inner Class: Collections of test/obsolete/deprecated methods
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
