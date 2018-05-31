# =============================================================================
# Library      : DNA Replication Origin
# Name         : handbook-of/ReplicationOrigin.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 18/05/18
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Gene-to-gene minimum distance
# Last Modified: 18/05/18
# -----------------------------------------------------------------------------
getTSS <- function(ensGene.genes) {
   ensGene.genes$TSS <- 0
 
   for (g in 1:nrow(ensGene.genes)) {
      if (ensGene.genes$strand[g] > 0)
         ensGene.genes$TSS[g] <- ensGene.genes$start_position[g]
      else
         ensGene.genes$TSS[g] <- ensGene.genes$end_position[g]
   }
 
   return(ensGene.genes)
}

getG2GQ4 <- function(tx.q4) {
   g2g.q4 <- list()
 
   for (q in 1:4) {
      genes <- tx.q4[[q]]
      ensGene.genes <- ensGene[genes,]
      ensGene.genes <- getTSS(ensGene.genes)
  
      g2g <- list()
      for (g in 1:length(genes)) {
         ensGene.gene <- ensGene.genes[genes[g],]
         ensGene.genes.chr <- subset(ensGene.genes, chromosome_name == ensGene.gene$chromosome_name)
         ensGene.genes.chr <- subset(ensGene.genes.chr, ensembl_gene_id != ensGene.gene$ensembl_gene_id)
   
         g2g[[g]] <- min(abs(ensGene.gene$TSS - ensGene.genes.chr$TSS))
      }
      g2g.q4[[q]] <- g2g
   }
   
   return(g2g.q4)
}

plotG2GQ4 <- function(g2g.q4, file.name, file.main, ylim) {
   distances <- c()
   quantiles <- c()
   for (q in 1:4) {
      distances <- c(distances, as.numeric(g2g.q4[[q]]))
      quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
   }
 
   pdf(file.name, height=6, width=4)
   names <- c("25%", "50%", "75%", "100%")
   boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=names, ylim=ylim, main=file.main)
   dev.off()
}

testT <- function(q3, q4) {
   return(t.test(as.numeric(q3), as.numeric(q4))$p.value)
}

testT <- function(q3, q4) {
 return(t.test(as.numeric(qX), as.numeric(q4))$p.value)
}

# -----------------------------------------------------------------------------
# Pipelines: Gene-to-gene minimum distance
# Last Modified: 18/05/18
# -----------------------------------------------------------------------------
pipeTPM <- function(wd, BASE) {
   base <- tolower(BASE)
 
   load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5_p47.RData")))
   tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)

   return(tpm.gene.input)
}

pipeRO <- function(tpm.gene.input.log2) {
   tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
   g2g.q4 <- getG2GQ4(tx.q4)
   
   return(g2g.q4)
}

# -----------------------------------------------------------------------------
# Methods: Density/Count
# http://www.sthda.com/english/articles/32-r-graphics-essentials/133-plot-one-variable-frequency-graph-density-distribution-and-more/
# Last Modified: 26/05/18
# -----------------------------------------------------------------------------
getDensity <- function(distances, count) {
   d <- density(log10(as.numeric(distances)))
   if (count)
      d$y <- d$y * d$n
    
   return(d)
}

plotDensity <- function(g2g.q4, file.name, file.main, count) {
   ymax <- max(getDensity(g2g.q4[[1]], count)$y)
   for (q in 2:4)
      if (max(getDensity(g2g.q4[[q]], count)$y) > ymax)
         ymax <- max(getDensity(g2g.q4[[q]], count)$y)
      
   xmax <- max(getDensity(g2g.q4[[1]], count)$x)
   for (q in 2:4)
      if (max(getDensity(g2g.q4[[q]], count)$x) > xmax)
         xmax <- max(getDensity(g2g.q4[[q]], count)$x)   
      
   cols <- c("blue", "deepskyblue", "salmon", "red")
   pdf(file.name, height=6, width=6)
   plot(getDensity(g2g.q4[[1]], count), col=cols[1], ylab="Frequency", xlab="Gene-to-gene min dist. (log10)", ylim=c(0, ymax), xlim=c(0, xmax), main=file.main)
   for (q in 2:4)
      lines(getDensity(g2g.q4[[q]], count), col=cols[q])

   legend("topleft", legend=c("25%", "50%", "75%", "100%"), levels(cols), fill=cols) 
   dev.off()
}

# =============================================================================
# Inner Class  : In-house File Reader
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 
# =============================================================================

# =============================================================================
# Inner Class  : Collections of test/obsolete/deprecated methods
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified:
# =============================================================================
