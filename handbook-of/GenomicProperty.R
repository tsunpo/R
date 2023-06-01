# =============================================================================
# Library      : Genomic Properties
# Name         : handbook-of/GenomicProperty.R
# Author       : Tsun-Po Yang (ty2@sanger.ac.uk)
# Last Modified: 30/05/23
# =============================================================================

# =============================================================================
# Methods: Density plot
# Last Modified: 30/05/23
# =============================================================================
plotSizeDensity <- function(sizes, main.text, file.name, title) {
   xlab.text <- "Size [kb]"
   ylab.text <- "Density"
   d <- density(sizes)

   pdf(file.name, height=6, width=6)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(d, xlab=xlab.text, ylab=ylab.text, main=main.text, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)

   rug(jitter(sizes))
   dev.off()
}
