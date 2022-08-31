## http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
plotPCA <- function(x, y, pca, trait, wd.de.data, file.name, size, file.main, legend, legend.title, cols, flip.x, flip.y) {
   scores <- pcaScores(pca)
   trait[is.na(trait)] <- "NA"
   trait.v <- legend.title
   
   #if (isNA(cols))
   #   cols <- c("red", "deepskyblue", "forestgreen", "purple3", "blue", "gold", "lightsalmon", "turquoise1", "limegreen")   #, "salmon", "tomato", "steelblue2", "cyan")
   #if (length(trait.v) > length(cols))
   #   cols <- rainbow(length(trait.v))
   #else
   #   cols <- cols[1:length(trait.v)]
   trait.col <- mapply(x = 1:length(trait), function(x) return(cols[which(trait[x] == trait.v)]))   ## Assign colours to each subtypes
   
   cols[which(trait.v == "NA")] <- "lightgray"
   trait.col[which(trait == "NA")] <- "lightgray"
   xlab.txt <- paste0("PC", x, " (", pcaProportionofVariance(pca, x), "%)")
   ylab.txt <- paste0("PC", y, " (", pcaProportionofVariance(pca, y), "%)")
   
   pdf(file.path(wd.de.data, paste0(file.name, "_", names(scores)[x], "-", names(scores)[y], ".pdf")), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   #plot(scores[, x]*flip.x, scores[, y]*flip.y, col=trait.col, pch=19, cex=2, main=file.main[1], xlab=xlab.txt, ylab=ylab.txt, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   plot(scores[, x]*flip.x, scores[, y]*flip.y, col=NA, pch=19, cex=2, main=file.main[1], xlab=xlab.txt, ylab=ylab.txt, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   idx <- which(trait == "non-WGD")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   idx <- which(trait == "WGD")
   points(scores[idx, x]*flip.x, scores[idx, y]*flip.y, col=trait.col[idx], pch=19, cex=2)
   
   #mtext(file.main[2], cex=1.2, line=0.3)
   for (l in 1:length(trait.v))
      trait.v[l] <- paste0(trait.v[l], " (n=", length(which(trait == trait.v[l])), ")")
      #trait.v[l] <- trait.v[l]
   
   #trait.v[5] <- "Others"   ## For PCA ALL
   #if (BASE != "") {
   #   trait.v <- trait.v[1:4]
   #   trait.v <- paste0(BASE, " ", trait.v)
   #   cols <- cols[1:4]
   #}
   #if (is.na(legend.title))
   #   legend(legend, trait.v, col=cols, pch=19, pt.cex=2.5, cex=1.8)   ##bty="n")
   #else
      legend(legend, trait.v, col=cols, pch=19, pt.cex=2.5, cex=1.6)
   #mtext(ylab.txt, side=2, line=2.75, cex=1.9)
   dev.off()
}
