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

getG2G <- function(genes) {
   genes <- ensGene[genes,]  
   genes <- getTSS(genes)
   
   g2g <- list()
   for (g in 1:nrow(genes)) {
      gene <- genes[g,]
    
      genes.chr <- subset(genes, chromosome_name == gene$chromosome_name)
      if (nrow(genes.chr) != 1)   ## ADD 09/09/18
         genes.chr <- subset(genes.chr, ensembl_gene_id != gene$ensembl_gene_id)
    
      g2g[[g]] <- min(abs(gene$TSS - genes.chr$TSS))
   }
   
   return(g2g)
}

getTxQ4G2G <- function(tx.q4) {
   tx.q4.g2g <- list()

   for (q in 1:4)
      tx.q4.g2g[[q]] <- getG2G(tx.q4[[q]])
   
   return(tx.q4.g2g)
}

boxplotTxQ4G2G <- function(wd.de.plots, base, BASE, tx.q4.g2g, samples) {
   distances <- c()
   quantiles <- c()
   for (q in 1:4) {
      distances <- c(distances, as.numeric(tx.q4.g2g[[q]]))
      quantiles <- c(quantiles, mapply(x = 1:length(tx.q4.g2g[[q]]), function(x) paste0("q", q)))
   }
 
   file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_g2g.pdf"))
   file.main <- paste0(BASE, " (n=", nrow(samples), ")")
   pdf(file.name, height=6, width=4)
   names <- c("Q1", "Q2", "Q3", "Q4")
   boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=names, main=file.main)
   dev.off()
}

testT <- function(q3, q4) {
   return(t.test(as.numeric(q3), as.numeric(q4))$p.value)
}

testW <- function(q3, q4) {
   trait <- rep(0, length(q3))
   trait <- c(trait, rep(1, length(q4)))
   trait <- as.factor(trait)
 
   expr <- as.numeric(c(q3, q4))
   return(wilcox.test(expr ~ trait, exact=F)$p.value)
}

testOnewayANOVA <- function(g2g.q4) {
   distances <- c()
   quantiles <- c()
   for (q in 1:4) {
      distances <- c(distances, as.numeric(g2g.q4[[q]]))
      quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
   }
   
   res.aov <- aov(distances ~ quantiles)
   return(summary(res.aov)[[1]]$Pr[1])
}

# -----------------------------------------------------------------------------
# Methods: Density/Count
# http://www.sthda.com/english/articles/32-r-graphics-essentials/133-plot-one-variable-frequency-graph-density-distribution-and-more/
# Last Modified: 26/05/18
# -----------------------------------------------------------------------------
getTxQ4G2GDensity <- function(distances, count) {
   d <- density(log10(as.numeric(distances)))
   if (count)
      d$y <- d$y * d$n
    
   return(d)
}

plotTxQ4G2GDensity <- function(wd.de.plots, base, BASE, tx.q4.g2g, samples, count, ymax) {
   if (is.null(ymax)) {
      ymax <- max(getTxQ4G2GDensity(tx.q4.g2g[[1]], count)$y)
      for (q in 2:4)
         if (max(getTxQ4G2GDensity(tx.q4.g2g[[q]], count)$y) > ymax)
            ymax <- max(getTxQ4G2GDensity(tx.q4.g2g[[q]], count)$y)
   }

   xmax <- max(getTxQ4G2GDensity(tx.q4.g2g[[1]], count)$x)
   for (q in 2:4)
      if (max(getTxQ4G2GDensity(tx.q4.g2g[[q]], count)$x) > xmax)
         xmax <- max(getTxQ4G2GDensity(tx.q4.g2g[[q]], count)$x)   
   
   file.name <- file.path(wd.de.plots, paste0("plot_", base, "_genes_tx_q4_g2g_density.pdf"))
   file.main <- paste0(BASE, " (n=", nrow(samples), ")")   
   cols <- c("blue", "deepskyblue", "salmon", "red")
   pdf(file.name, height=6, width=6)
   plot(getTxQ4G2GDensity(tx.q4.g2g[[1]], count), col=cols[1], ylab="Frequency", xlab="Gene-to-gene min dist. (log10)", ylim=c(0, ymax), xlim=c(0, xmax), main=file.main)
   for (q in 2:4)
      lines(getTxQ4G2GDensity(tx.q4.g2g[[q]], count), col=cols[q])

   legend("topleft", legend=c("Q1", "Q2", "Q3", "Q4"), levels(cols), fill=cols) 
   dev.off()
}

# -----------------------------------------------------------------------------
# Pipelines: Gene-to-gene minimum distance
# Last Modified: 18/05/18
# -----------------------------------------------------------------------------
pipeTPM <- function(wd, BASE) {
   base <- tolower(BASE)
 
   load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_tpm0.RData")))
   tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)

   return(tpm.gene.input)
}

pipeTxQ4G2G <- function(wd, BASE, genes.list) {
   tpm.gene.input <- pipeTPM(wd, BASE)
   tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input[genes.list,])
   tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
   tx.q4.g2g <- getTxQ4G2G(tx.q4)

   p3 <- testW(tx.q4.g2g[[3]], tx.q4.g2g[[4]])
   p2 <- testW(tx.q4.g2g[[2]], tx.q4.g2g[[4]])
   p1 <- testW(tx.q4.g2g[[1]], tx.q4.g2g[[4]])
   c(p1, p2, p3)
   p.adjust(c(p1, p2, p3), method="bonferroni")
   
   return(tx.q4.g2g)
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
# getQ1G2GQ4 <- function(tx.q4) {
#  g2g.q4 <- list()
#  
#  for (q in 1:4) {
#   genes <- tx.q4[[q]]
#   ensGene.genes <- ensGene[genes,]
#   ensGene.genes <- getTSS(ensGene.genes)
#   
#   g2g <- list()
#   for (g in 1:length(genes)) {
#    ensGene.gene <- ensGene.genes[genes[g],]
#    ensGene.genes.chr <- subset(ensGene.genes, chromosome_name == ensGene.gene$chromosome_name)
#    ensGene.genes.chr <- subset(ensGene.genes.chr, ensembl_gene_id != ensGene.gene$ensembl_gene_id)
#    
#    g2g[[g]] <- min(abs(ensGene.gene$TSS - ensGene.genes.chr$TSS))
#   }
#   
  ## So far same as in getG2GQ4, but only keep the least G2G distance gene
#   q1 <- quantile(as.numeric(g2g))[2]
#   g2g.q4[[q]] <- g2g[as.numeric(g2g) <= q1]
#  }
#  
#  return(g2g.q4)
# }

# get1KBG2GQ4 <- function(tx.q4, min) {
#  g2g.q4 <- list()
#  
#  for (q in 1:4) {
#   genes <- tx.q4[[q]]
#   ensGene.genes <- ensGene[genes,]
#   ensGene.genes <- getTSS(ensGene.genes)
#   
#   g2g <- list()
#   for (g in 1:length(genes)) {
#    ensGene.gene <- ensGene.genes[genes[g],]
#    ensGene.genes.chr <- subset(ensGene.genes, chromosome_name == ensGene.gene$chromosome_name)
#    ensGene.genes.chr <- subset(ensGene.genes.chr, ensembl_gene_id != ensGene.gene$ensembl_gene_id)
#    
#    g2g[[g]] <- min(abs(ensGene.gene$TSS - ensGene.genes.chr$TSS))
#   }
#   
  ## So far same as in getG2GQ4, but only keep the least G2G distance gene
#   g2g.q4[[q]] <- g2g[as.numeric(g2g) <= min]
#  }
#  
#  return(g2g.q4)
# }
