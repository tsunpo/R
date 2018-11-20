# =============================================================================
# Library      : Common Daily R Meal
# Name         : handbook-of/Analytics.R 
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 16/02/17
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Read data from a file
# Tricks: 1) Assign a specified column as row names if a column number in 'rownames' is given
#         2) Return a vector if only one column is in the file.
# Usage: table <- readTable("file.txt", header=T, rownames=T, sep="\t")
#        table <- readTable("file.txt", header=T, rownames=2, sep="\t")
#        list  <- readTable("list.txt", header=F, rownames=F, sep="")
# Last Modified: 22/10/16
# -----------------------------------------------------------------------------
readTable <- function(file, header, rownames, sep) {
   table <- read.table(file, header=header, sep=sep, fill=T, as.is=T, comment.char="#", na.strings="NA") 

   if (rownames == T)
   	  rownames(table) <- table[,1]
   else if (is.numeric(rownames))         ## Assign a specified column as row names if a column number in 'rownames' is given
   	  rownames(table) <- table[,rownames]

   if (header == F && ncol(table) == 1)   ## Return a vector if only one column is in the file
      table <- as.vector(table[,1])

   return(table)
}

# -----------------------------------------------------------------------------
# Method: Output table to a file
# Tricks: Automatically add an additional "ID_REF" for the first column name if rownames=T
# Usage: writeTable(table, "file.txt", colnames=T, rownames=T, sep="\t")
# Last Modified: 18/01/13
# -----------------------------------------------------------------------------
writeTable <- function(table, file, colnames, rownames, sep) {
   if (rownames == T) {
   	  names <- c("ID_REF", colnames(table))
   	  table$ID_REF <- rownames(table)
   	  
   	  write.table(table[,names], file=file, row.names=F, col.names=colnames, sep=sep, quote=F, na="NA")
   } else
      write.table(table, file=file, row.names=rownames, col.names=colnames, sep=sep, quote=F, na="NA")
}

# -----------------------------------------------------------------------------
# Method: Initialise a data-frame matrix
# Usage: toTable(NA, 8, 0, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
#        toTable("", length(genes), length(samples), genes)
#        toTable(0, length(genes), length(samples), genes)
# Last Modified: 02/02/17
# -----------------------------------------------------------------------------
toTable <- function(value, ncol, nrow, colnames) {
   table <- data.frame(matrix(value, nrow, ncol))
   names(table) <- colnames
   
   if (value == "" && !is.na(value))   ## ADD "" 02/02/17; ADD !is.na() 11/02/17 
      for (c in 1:ncol)
         table[,c] <- as.vector(table[,c])

   return(table)
}

# -----------------------------------------------------------------------------
# Method: Modification of subset(refGene, name2 == "MYCL")
# Usage: subset0(refGene, "name2", "MYCL")
#        subset0(refGene, "name2", c("MYCL", "MYCN", "MYC"))   ## BUG FIX 16/02/17
# Last Modified: 16/02/17
# -----------------------------------------------------------------------------
subset0 <- function(df, filter, value) {
   if (length(value) == 1) {
      return(df[which(value == df[,filter]),])   ## BUG FIX 16/02/17: Dose not work on length(value) != 1
   } else {
      subsets <- df[which(value[1] == df[,filter]),]
      for (v in 2:length(value))
         subsets <- rbind(subsets, df[which(value[v] == df[,filter]),])
      
      return(subsets)
   }
}

separator <- function(number) {
   return(formatC(number, format="f", big.mark=",", digits=0))
}
