# =============================================================================
# Library      : Cell Cycle
# Name         : handbook-of/CellCycle.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 07/09/18
# =============================================================================

# -----------------------------------------------------------------------------
# Distribution of cell cycle genes in four quantiles
# Last Modified: 07/09/18
# -----------------------------------------------------------------------------
getTxQ4Cycle <- function(tx.q4, core.G1S, core.G2M, periodic.G1S, periodic.G2M) {
   tx.q4.cycle   <- list(list(), list(), list(), list())
   for (q in 1:4) {
      genes <- tx.q4[[q]]

      genes.G1S <- intersect(genes, unique(c(core.G1S, periodic.G1S)))
      genes.G2M <- intersect(genes, unique(c(core.G2M, periodic.G2M)))

      tx.q4.cycle[[q]][[1]] <- genes.G1S
      tx.q4.cycle[[q]][[2]] <- genes.G2M
      tx.q4.cycle[[q]][[3]] <- setdiff(genes, c(genes.G1S, genes.G2M))
   }

   return(tx.q4.cycle)
}

barplotTxQ4Cycle <- function(wd.de.plots, base, BASE, tx.q4.cycle, beside) {
   colnames <- c("Q1", "Q2", "Q3", "Q4")
   rownames <- c("G1/S", "G2/M", "Others")
   
   percent <- toTable(0, length(colnames), 3, colnames)
   counts  <- toTable(0, length(colnames), 3, colnames)
   for (q in 1:4) {
      sum <- length(tx.q4.cycle[[q]][[1]])
      for (r in 2:3)
         sum <- sum + length(tx.q4.cycle[[q]][[r]])
  
      for (r in 1:3) {
         percent[r, q] <- length(tx.q4.cycle[[q]][[r]]) / sum * 100
         counts[r, q]  <- length(tx.q4.cycle[[q]][[r]]) 
      }
   }
   writeTable(counts,  file.path(wd.de.plots, paste0(base, "_genes_tx_q4_cycle_counts.txt")), colnames=T, rownames=T, sep="\t")
   writeTable(percent, file.path(wd.de.plots, paste0(base, "_genes_tx_q4_cycle_percent.txt")), colnames=T, rownames=T, sep="\t")
   percent <- as.matrix(percent)
   
   ##
   file.name <- file.path(wd.de.plots, paste0("barplot_", base, "_genes_tx_q4_cycle_percent.pdf"))
   cols <- c("cornsilk", "darkgray")
   pdf(file.name, height=6, width=4)
   barplot(percent[-3,], ylab="Cell cycle genes (%)", xlab="Expression", col=cols, main=BASE, beside=beside)
   legend("topleft", legend=rownames[-3], fill=cols)
   dev.off()
}

# -----------------------------------------------------------------------------
# Gene length analysis
# Last Modified: 07/09/18
# -----------------------------------------------------------------------------
initLength0 <- function() {
   colnames <- c(colnames(ensGene), "Length", "Group")
   return(toTable("", length(colnames), 0, colnames))
}

initLength <- function(genes, group) {
   ens.genes <- ensGene[genes,]
   ens.genes$Length <- mapply(x = 1:nrow(ens.genes), function(x) getLengthTx(genes[x]))
   ens.genes$Group  <- group
 
   return(ens.genes)
}

getTxQ4Length <- function(tx.q4) {
   tx.q4.length <- initLength0()
   for (q in 1:4) {
      genes <- tx.q4[[q]]
      tx.q4.length <- rbind(tx.q4.length, initLength(genes, q))
   }
   tx.q4.length$Group <- as.factor(tx.q4.length$Group)
   tx.q4.length$Length <- log10(tx.q4.length$Length)
   
   return(tx.q4.length)
}

boxplotTxQ4Length <- function(wd.de.plots, base, BASE, tx.q4.length) {
   colnames <- c("Q1", "Q2", "Q3", "Q4")

   file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_length.pdf"))
   pdf(file.name, height=6, width=4)
   ymin <- min(tx.q4.length$Length)
   ymax <- max(tx.q4.length$Length)
   boxplot(Length ~ Group, data=tx.q4.length, outline=T, names=colnames, ylim=c(ymin, ymax), ylab="Gene length (log10)", xlab="Expression", main=BASE)
   dev.off()
}
