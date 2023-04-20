# =============================================================================
# Library      : Reusable common methods
# Name         : handbook-of/Commons.R 
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
   table <- read.table(file, header=header, sep=sep, fill=T, as.is=T, comment.char="#", na.strings="NA", quote="") 

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
# Usages: subset00(ensGene, "external_gene_name", "MYCL")
#         subset00(ensGene, "external_gene_name", c("MYC", "MYCL", "MYCN"))
# Last Modified: 209/03/22; 9/11/18
# -----------------------------------------------------------------------------
grep0 <- function(patterns, vectors) {
   grep(paste(patterns, collapse="|"), vectors)
}

subset0 <- function(table, colname, values) {
   return(table[which(table[, colname] == values),])   ## NEW VERSION: 09/03/22
}

subset00 <- function(table, colname, values) {
   return(table[which(table[, colname] %in% values),])   ## NEW VERSION: 29/11/18
}

round0 <- function(number, digits) {
   return(sprintf(paste0("%.", digits, "f"), round(number, digits=digits)))
}

separator <- function(number) {
   return(formatC(number, format="f", big.mark=",", digits=0))
}

scientific <- function(number, digits=2) {
   return(formatC(number, format="E", digits=digits))
}

strsplit0 <- function(list, sep, pos) {
   return(mapply(x = 1:length(list), function(x) unlist(strsplit(list[x], sep))[pos]))
}

which0 <- function(df, value) {
   return(which(rownames(df) == value))
}

firstup <- function(x) {
	  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
	  return(x)
}

toString0 <- function(x, sep) {
	  a <- unlist(strsplit(x, sep))
	  n <- length(a)
	  return(paste0(c(firstup(a[1]), a[2:n]), collapse=" "))
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

getGenes <- function(genes) {
   outs <- ensGene[1,][-1,]
   
   for (g in 1:length(genes)) {
      out <- getGene(genes[g])
      if (nrow(out) > 0)
         outs <- rbind(outs, out)
   }
   
   return(outs)
}

# -----------------------------------------------------------------------------
# Collections: Texts for plot
# Last Modified: 01/04/22
# -----------------------------------------------------------------------------
text.In.silico  <- expression(italic("In silico") ~ "sorting")

text.Log10.P    <- expression("-log" * ""[10] * "(" * italic("P") * ")")
text.Log2.TPM.1 <- expression("log" * ""[2] * "(TPM + 1)")
text.Log2       <- expression("log" * ""[2])

text.SCLC <- "Lung-SCLC"

firstup <- function(x) {
   x <- tolower(x)
   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
   
   return(x)
}

# -----------------------------------------------------------------------------
# Collections: Colours for plot
# Links: https://blog.datawrapper.de/beautifulcolors/
#        https://medium.com/nightingale/how-to-choose-the-colors-for-your-data-visualizations-50b2557fa335
#        https://cran.r-project.org/web/packages/unikn/vignettes/colors.html
#        https://standards.aviva.com/global-experience-principles/style-guide/colour/ (AVIVA)
#        https://brandportal.basf.com/global/en/create/corporate-design/colors.html (BASF)
#        https://nasa.github.io/nasawds-site/components/colors/ (NASA)
#        https://www.dpdhl-brands.com/dhl/en/guides/design-basics/colors-materials/colors.html (DHL)
#        https://www.imperial.ac.uk/brand-style-guide/visual-identity/brand-colours/ (Imperial)
#        https://colorswall.com/palette/38/ (Google)
#        https://colorswall.com/palette/72/ (Microsoft)
#        https://en.wikipedia.org/wiki/Shades_of_green#Kelly_green
#        https://www.rgbtohex.net/
#        https://www.hexcolortool.com/
# Last Modified: 22/11/20
# -----------------------------------------------------------------------------
## AVIVA
aviva.red     <- "#B10101"
aviva.orange  <- "#FFA000"

aviva.lightyellow <- "#FFE340"
aviva.yellow      <- "#FFD900"
aviva.darkyellow  <- "#FFC100"

aviva.green   <- "#3E812C"
aviva.skyblue <- "#44C0FF"

aviva.lightblue <- "#407BC9"
aviva.blue      <- "#004FB6"
aviva.darkblue  <- "#00359D"

aviva.lightplum <- "#7B5182"
aviva.plum      <- "#4E1758"
aviva.darkplum  <- "#340D3D"

## BASF
basf.red        <- "#C50022"
basf.orange     <- "#F39500"
basf.lightgreen <- "#65AC1E"
basf.darkgreen  <- "#00793A"
basf.lightblue  <- "#21A0D2"
basf.darkblue   <- "#004A96"

## Google
google.red    <- "#EA4335"
google.yellow <- "#FBBC05"
google.green  <- "#34A853"
google.blue   <- "#4285F4"

## DHL
dhl.red    <- "#D40511"
dhl.yellow <- "#FFCC00"

## NASA
nasa.blue <- "#105bd8"

## FlowJo
flowjo.blue <- "#989aff"
flowjo.red  <- "#ff9899"
flowjo.grey <- "#b7b7b7"

## Color wheel
colorwheel.orange     <- "#FFA700"
colorwheel.lightblue  <- "#7DC0F7"
colorwheel.mediumblue <- "#005FCF"

# -----------------------------------------------------------------------------
# Method: My colours
# Last Modified: 22/11/20
# -----------------------------------------------------------------------------
## RT/RFD
red       <- google.red
orange    <- colorwheel.orange
yellow    <- dhl.yellow
green     <- basf.lightgreen
lightblue <- colorwheel.lightblue
blue      <- nasa.blue
purple    <- aviva.lightplum
grey      <- "darkgray"

## PCA
colorwheel.orange.lighter20 <- "#FFDA33"
basf.orange.lighter20     <- "#FFC833"
basf.darkgreen.lighter20  <- "#33AC6D"
basf.darkgreen.lighter30  <- "#4DC687"
basf.lightgreen.lighter10 <- "#7FC638"
basf.lightgreen.lighter20 <- "#98DF51"

#google.red.lighter60 <- "#FFDCCE"
google.red.lighter50 <- "#FFC2B4"
#google.red.lighter40 <- "#FFA99B"
#google.red.lighter35 <- "#FF9C8E"
#google.red.lighter30 <- "#FF9082"
google.red.lighter25 <- "#FF8375"
nasa.blue.lighter25  <- "#509BFF"
nasa.blue.lighter50  <- "#8FDAFF"
#nasa.blue.lighter60  <- "#A9F4FF"

orange.lighter <- basf.orange.lighter20
green.lighter  <- basf.lightgreen.lighter10

#red.lightest <- google.red.lighter60
red.lighter  <- google.red.lighter50
red.light    <- google.red.lighter25
blue.light   <- nasa.blue.lighter25
blue.lighter <- nasa.blue.lighter50
#blue.lightest <- nasa.blue.lighter60

# -----------------------------------------------------------------------------
# TEMP
# -----------------------------------------------------------------------------
#basf.red.lighter50 <- "#FF7FA1"
#dhl.red.lighter50  <- "#FF8490"

#green      <- "#01DF01"
#lightgreen <- "#55C8A9"
#darkgreen  <- "#00A64C"
#imperial.green.lime   <- "#BBCEOO"
#wiki.kelly.green      <- "#4CBB17"
#imperial.green.kermit <- "#66A40A"
#microsoft.green <- "#7FBA00"
#twitter.blue <- "#00ACEE"

#pinky1   <- "#F3BFCB"
#pinky2   <- "#ECA0B2"
#seablau2 <- "#A6E1F4"
#seablau3 <- "#59C7EB"

#bbc.red    <- "#BE002C"
#bbc.sports <- "#ffd230"
#bbc.science.blue <- "#0066c2"

#aviva.lightwarmgrey <- "#858392"
#aviva.warmgrey      <- "#5C596D"
#aviva.darkwarmgrey  <- "#4B485B"

#aviva.lightdarkgrey <- "#737373"
#aviva.darkgrey      <- "#444444"
#aviva.darkdarkgrey  <- "#2D2D2D"

#nasa.gray      <- "#5b616b"
#nasa.lightgray <- "#aeb0b5"
