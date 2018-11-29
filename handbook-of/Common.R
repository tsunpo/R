# =============================================================================
# Library      : Daily R Meal
# Name         : handbook-of/Common.R 
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 29/11/18
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Read data from a file
# Tricks: 1) Assign a specified column as row names if a column number in 'rownames' is given
#         2) Return a vector if only one column is in the file.
# Usages: table <- readTable("file.txt", header=T, rownames=T, sep="\t")
#         table <- readTable("file.txt", header=T, rownames=2, sep="\t")
#         array <- readTable("list.txt", header=F, rownames=F, sep="")
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
# Usages: writeTable(table, "file.txt", colnames=T, rownames=T, sep="\t")
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
# Usages: toTable(NA, 8, 0, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
#         toTable("", length(genes), length(samples), genes)
#         toTable(0, length(genes), length(samples), genes)
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
# Method: Modification of subset(ensGene, external_gene_name %in% "MYCL")
# Usages: subset0(ensGene, "external_gene_name", "MYCL")
#         subset0(ensGene, "external_gene_name", c("MYC", "MYCL", "MYCN"))
# Last Modified: 29/11/18
# -----------------------------------------------------------------------------
subset0 <- function(df, colname, values) {
   return(df[which(df[, colname] %in% values),])   ## NEW VERSION: 29/11/18
}

separator <- function(number) {
   return(formatC(number, format="f", big.mark=",", digits=0))
}

# -----------------------------------------------------------------------------
# Method: Shortcuts for fetching an Ensembl gene
# Usages: getGene("ENSG00000157168")
#         getGene("NRG")
#         getGene("NRG1")
#         getGene("NRG1", beginWith=T)
#         getGene("-AS2", endWith=F)
#         getGene("chr8", 32500000, 33500000, protein_coding=T)
# Last Modified: 29/11/18
# -----------------------------------------------------------------------------
getGene <- function(gene, start=NA, end=NA, beginWith=F, endWith=F, protein_coding=F) {
   if (grepl("^ENSG", gene)) {
      return(ensGene[gene,])
   } else if (grepl("^chr", gene, ignore.case=F)) {
      if (is.na(start) || is.na(end))
         return(ensGene[1,][-1,])
      else {
         ensGene.chr     <- subset(ensGene, chromosome_name == gene)
         ensGene.chr.s   <- subset(ensGene.chr,   end_position >= start)   ## Find genes that overlapping to the searched region;
         ensGene.chr.s.e <- subset(ensGene.chr.s, start_position <= end)
         
         if (protein_coding)
            ensGene.chr.s.e <- subset(ensGene.chr.s.e, gene_biotype == "protein_coding")
         return(ensGene.chr.s.e[order(ensGene.chr.s.e$start_position),])   ## and sort by start_position
      }
   } else {
      genes <- subset(ensGene, external_gene_name == gene)
      if (nrow(genes) == 0)                                           ## If there is no direct match of the input gene name;
         genes <- ensGene[grepl(gene, ensGene$external_gene_name),]   ## report all gene names begin with the searched gene;
 
      if (beginWith) {
         if (nrow(ensGene[grepl(gene, ensGene$external_gene_name),]) > 1)
            genes <- ensGene[grepl(gene, ensGene$external_gene_name),]
         genes <- genes[grepl(paste0("^", gene), genes$external_gene_name),]
      } else if (endWith)
         genes <- genes[grepl(paste0(gene, "$"), genes$external_gene_name),] 
 
      if (protein_coding)
         genes <- subset(genes, gene_biotype == "protein_coding")
      return(genes[order(genes$external_gene_name),])                 ## and sort by gene names
   }
}
