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

plotG2GQ4 <- function(wd.asym.plots, BASE, sample.size, g2g.q4) {
   distances <- c()
   quantiles <- c()
   for (q in 1:4) {
      distances <- c(distances, as.numeric(g2g.q4[[q]]))
      quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
   }
 
   pdf(paste0(wd.asym.plots, tolower(BASE), "_asym_tx_g2g.pdf"), height=6, width=4)
   boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=c("25%", "50%", "75%", "100%"), main=paste0(BASE, " (n=", sample.size, ")"))
   dev.off()
}

testT <- function(qX, q4) {
   return(t.test(as.numeric(qX), as.numeric(q4))$p.value)
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
