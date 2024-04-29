# =============================================================================
# Library      : Survival Probability
# Name         : handbook-of/Survival.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 09/08/21; 03/05/19
# =============================================================================
#library(survival)
#library(survminer)
#library(tidyverse)
#library(broom)
#library(grid)
#library(dplyr)
#library(ggalt)
#library(gridExtra)

# -----------------------------------------------------------------------------
# Methods: Kaplan-Meier survival analysis
# Link(s): https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_survival_analysis.pdf
#          http://www.stat.columbia.edu/~madigan/W2025/notes/survival.pdf
#          https://support.minitab.com/en-us/minitab/18/help-and-how-to/modeling-statistics/reliability/supporting-topics/reliability-metrics/what-is-the-survival-probability/
# Last Modified: 03/05/19
# -----------------------------------------------------------------------------
plotOS <- function(file.name, main.text, xlab.text, ylab.text, x, y, p, rho, lwd=4) {
   ymax <- -log10(p)
   ymax <- ymax + ymax/5.5
   xmin <- min(x)
    
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(y ~ x, ylim=c(0, ymax), xaxt="n", xlab=xlab.text, ylab=ylab.text, col="white", main=main.text, pch=1, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   idx.down   <- which(x < rho)
   points(x[idx.down], y[idx.down], pch=16, col=blue, cex=2)
   idx.up   <- which(x >= rho)
   points(x[idx.up],   y[idx.up],   pch=16, col=red, cex=2)   
   
   points(rho, -log10(p), pch=1, col=red, cex=2, lwd=lwd)
   abline(v=rho, lty=5, lwd=5, col=red)
   text(xmin, ymax - ymax/10, expression(italic("    P")~"="), col=red, pos=3, cex=1.7)
   text(xmin, ymax - ymax/10, paste0("                           ", scientific(p, digits=2)), col=red, pos=3, cex=1.7)
   
   axis(side=1, at=seq(-0.4, 0.4, by=0.4), labels=c(-0.4, 0, 0.4), cex.axis=1.7)
   axis(side=1, at=rho, labels=round(rho, 2), cex.axis=1.7, col.axis=red, font.axis=2)
 
   dev.off()
}

plotSurvfit <- function(fit, file.name, main.text, legend.labs, name, strata, cols, size) {
   legends <- c()
   strata.cols <- cols
   for (l in 1:length(legend.labs)) {
      strata.idx <- grep(strata[l], gsub(name, "", names(fit$strata)))
    
      legends <- c(legends, paste0(legend.labs[l], " (n=", fit[[1]][strata.idx], ")"))
      strata.cols[strata.idx] <- cols[l]
   }
 
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="Probability [%]", col=strata.cols, main=main.text[1], mark.time=T, lwd=3, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   legend("topright", legend=legends, lwd=5, col=cols, cex=1.65)

   max <- max(fit[[2]])
   #text(max/5.5, 0.18, "                HR = 1.90 (1.61 - 2.25)", cex=1.7)
   text(max/4, 0.07, expression(italic('P')~"="~"                "), cex=1.7)
   text(max/4, 0.07, paste0("      ", scientific(surv_pvalue(fit)$pval, digits=2)), cex=1.7)
 
   dev.off()
}


plotSurvfit55 <- function(fit, file.name, main.text, legends.original, legends.text, cols) {
   legends <- c()
   for (l in 1:length(legends.text))
      legends <- c(legends, paste0(legends.text[l], " (n=", fit[[1]][l], ")"))
 
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="", col=cols, main=main.text[1], mark.time=T, lwd=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   legend("topright", legend=legends, lwd=4, col=cols, cex=1.7)
   #mtext(main.text[2], cex=1.25, line=0.2)
 
   mtext("% OS", side=2, line=2.85, cex=1.7)
 
   max <- max(fit[[2]])
   #text(max/5.5, 0.18, "                HR = 1.90 (1.61 - 2.25)", cex=1.7)
   text(max/4, 0.07, expression(italic('P')~"="~"                "), cex=1.7)
   text(max/4, 0.07, paste0("      ", scientific(surv_pvalue(fit)$pval, digits=2)), cex=1.7)

   dev.off()
}

plotStripchartV <- function(wd.rt.plots, file.name, main.text, samples.v, hists.surv, samples.surv) {
   adjustcolor.red  <- adjustcolor(red, alpha.f=0.5)
   adjustcolor.blue <- adjustcolor(blue, alpha.f=0.5)
   cols <- c(adjustcolor.blue, adjustcolor.red)

   pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=12, width=8)
   par(mar = c(5, 14.5, 4, 2))
   boxplot(COR ~ CANCER, data=subset(samples.v, histology_abbreviation %in% hists.surv), horizontal=T, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2)
   #text(labels=labels, x=1:26, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, side=1, cex=1.8)

   #abline(v=sample$COR, lty=5, lwd=5, col=green)

   #stripchart(COR ~ CANCER, data=samples.h,                           method="jitter", cex=1.5, pch=19, col="lightgray", vertical=F, add=T, at=c(1:length(hists)))
   stripchart(COR ~ CANCER, data=subset(samples.surv, SORTING == "G1"), method="jitter", cex=1.5, pch=19, col=cols[1], vertical=F, add=T, at=c(1:length(hists.surv)))
   stripchart(COR ~ CANCER, data=subset(samples.surv, SORTING == "S"),  method="jitter", cex=1.5, pch=19, col=cols[2], vertical=F, add=T, at=c(1:length(hists.surv)))

   axis(side=1, at=seq(-0.4, 0.4, by=0.4), labels=c(-0.4, 0, 0.4), cex.axis=1.8)
   #axis(side=1, at=sample$COR, labels=round(sample$COR, 2), cex.axis=1.8, col.axis=green, font.axis=2)

   mtext(expression(italic('In silico')~"sorting"), side=1, line=3.5, cex=1.9)
   #mtext("", cex=1.2, line=0.3)
   mtext(rev(hists.surv), side=2, line=0.5, cex=1.9, at=1:length(hists.surv), las=1)

   legend("bottomright", legend=c("Non-proliferative", "              "), pch=19, pt.cex=2.5, col=c(blue, "white"), text.col=c("black", "white"), cex=1.7, horiz=T, text.width=0.59)
   legend("bottomright", legend="Proliferative", pch=19, pt.cex=2.5, col=red, text.col="black", bty="n", cex=1.7, horiz=T, text.width=0.49)
   dev.off()
}

# -----------------------------------------------------------------------------
# Driver PMR (Proliferation Mutational Ratio)
# Last Modified: 03/04/22
# -----------------------------------------------------------------------------
getMutGene <- function(annotmuts.impact) {
	  colnames <- unique(annotmuts.impact$sampleID)
	  rownames <- rownames(sort(table(annotmuts.impact$gene), decreasing=T))
	  mut.gene <- toTable(0, length(colnames), length(rownames), colnames)
	  rownames(mut.gene) <- rownames
	  for (c in 1:length(colnames)) {
		    sample.mut <- subset(annotmuts.impact, sampleID == colnames[c])
	    	mut <- as.data.frame(sort(table(sample.mut$gene), decreasing=T))
		    if (nrow(mut) == 1) {
			      colnames(mut) <- "Freq"
		    } else
			      rownames(mut) <- mut$Var1
		
		    overlaps <- intersect(rownames, rownames(mut))
		    if (nrow(mut) == 1) {
			      mut.gene[overlaps, c] <- mut[overlaps,]
		    } else
			      mut.gene[overlaps, c] <- mut[overlaps,]$Freq
	  }
	  
	  return(mut.gene)
}

getPdNdS <- function(mut.gene.n, mut.gene.s, samples.mut.hist, s=0, n=1, t="Fisher") {
	  pmr.mut.gene <- toTable(0, 14, nrow(mut.gene.n), c("MUT", "FREQ_N", "FREQ_S", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "P", "FDR", "G1MR", "SMR", "PMR", "TYPE_N", "TYPE_S"))
   pmr.mut.gene$MUT <- rownames(mut.gene.n)
   rownames(pmr.mut.gene) <- pmr.mut.gene$MUT
   
   pmr.mut.gene$FREQ_N <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene.n[x,])))
   pmr.mut.gene$FREQ_S <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene.s[x,])))
   #pmr.mut.gene <- subset(pmr.mut.gene, FREQ_SAMPLE >= s)
 
   for (g in 1:nrow(pmr.mut.gene)) {
      gene <- pmr.mut.gene$MUT[g]
      ids.n  <- colnames(mut.gene.n)[which(mut.gene.n[gene,] != 0)]
      ids.s  <- colnames(mut.gene.s)[which(mut.gene.s[gene,] != 0)]
      pmr.mut.gene$TYPE_N[g] <- paste(unique(samples.mut.hist[ids.n,]$histology_abbreviation), collapse=",")
      pmr.mut.gene$TYPE_S[g] <- paste(unique(samples.mut.hist[ids.s,]$histology_abbreviation), collapse=",")
      
      samples.mut.mut <- subset(samples.mut.hist, icgc_specimen_id %in% ids.n)
      samples.mut.wt  <- subset(samples.mut.hist, icgc_specimen_id %in% ids.s)
      
      mut.g1 <- rownames(subset(samples.mut.mut, M2 == 1))
      mut.s  <- rownames(subset(samples.mut.mut, M2 == 2))
      wt.g1  <- rownames(subset(samples.mut.wt,  M2 == 1))
      wt.s   <- rownames(subset(samples.mut.wt,  M2 == 2))
      
      test <- toTable(0, 2, 2, c("G1", "S"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- sum(mut.gene.n[g, mut.g1])
      test[1, 2] <- sum(mut.gene.n[g, mut.s])
      test[2, 1] <- sum(mut.gene.s[g, wt.g1])
      test[2, 2] <- sum(mut.gene.s[g, wt.s])

      if ((test[1, 1] >= n) && (test[1, 2] >= n) && (test[2, 1] >= n) && (test[2, 2] >= n)) {
         pmr.mut.gene$MUT_G1[g] <- test[1, 1]
         pmr.mut.gene$MUT_S[g]  <- test[1, 2]
         pmr.mut.gene$WT_G1[g]  <- test[2, 1]
         pmr.mut.gene$WT_S[g]   <- test[2, 2]
         if (t == "USP")    pmr.mut.gene$P[g] <- USP.test(test)[[1]]
         if (t == "Fisher") pmr.mut.gene$P[g] <- fisher.test(test)[[1]]
         if (t == "Chisq")  pmr.mut.gene$P[g] <- chisq.test(test)[[3]]
   
         pmr.mut.gene$G1MR[g] <- pmr.mut.gene$MUT_G1[g] / pmr.mut.gene$WT_G1[g]
         pmr.mut.gene$SMR[g]  <- pmr.mut.gene$MUT_S[g]  / pmr.mut.gene$WT_S[g]
         pmr.mut.gene$PMR[g]  <- pmr.mut.gene$SMR[g]    / pmr.mut.gene$G1MR[g]
      }
   }
   pmr.mut.gene <- subset(pmr.mut.gene, P != 0)
   pmr.mut.gene <- pmr.mut.gene[order(pmr.mut.gene$P, decreasing=F),]
   if (nrow(pmr.mut.gene) >= 4) {
      pmr.mut.gene$FDR <- testFDR(pmr.mut.gene$P, "BH")
   }
   #pmr <- pmr.mut.gene
   #annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
   #pmr.mut.gene <- cbind(annot[rownames(pmr),], pmr.mut.gene[, -1])
   
   return(pmr.mut.gene)
}

getPMR <- function(mut.gene, samples.mut, s=0, n=1, t="Fisher") {
	  pmr.mut.gene <- toTable(0, 13, nrow(mut.gene), c("MUT", "FREQ_MUT", "FREQ_SAMPLE", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "P", "FDR", "G1MR", "SMR", "PMR", "TYPE"))
	  pmr.mut.gene$MUT <- rownames(mut.gene)
	  rownames(pmr.mut.gene) <- pmr.mut.gene$MUT
	
	  pmr.mut.gene$FREQ_MUT <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene[x,])))
	  pmr.mut.gene$FREQ_SAMPLE <- mapply(x = 1:nrow(pmr.mut.gene), function(x) length(which(mut.gene[x,] != 0)))
	  #pmr.mut.gene <- subset(pmr.mut.gene, FREQ_SAMPLE >= s)
	
	  for (g in 1:nrow(pmr.mut.gene)) {
		    gene <- pmr.mut.gene$MUT[g]
		    ids  <- colnames(mut.gene)[which(mut.gene[gene,] != 0)]
		    pmr.mut.gene$TYPE[g] <- paste(unique(samples.mut[ids,]$histology_abbreviation), collapse=",")
		
		    samples.mut.mut <- subset(samples.mut, icgc_specimen_id %in% ids)
		    samples.mut.wt  <- samples.mut[setdiff(rownames(samples.mut), ids),]
		
		    mut.g1 <- rownames(subset(samples.mut.mut, M2 == 1))
		    mut.s  <- rownames(subset(samples.mut.mut, M2 == 2))
		    wt.g1  <- rownames(subset(samples.mut.wt,  M2 == 1))
		    wt.s   <- rownames(subset(samples.mut.wt,  M2 == 2))
		
		    if ((length(mut.g1) >= n) && (length(mut.s) >= n) && (length(wt.g1) >= n) && (length(wt.s) >= n)) {
			      test <- toTable(0, 2, 2, c("G1", "S"))
			      rownames(test) <- c("MUT", "WT")
			      test[1, 1] <- length(mut.g1)
			      test[1, 2] <- length(mut.s)
			      test[2, 1] <- length(wt.g1)
			      test[2, 2] <- length(wt.s)
			
			      pmr.mut.gene$MUT_G1[g] <- length(mut.g1)
			      pmr.mut.gene$MUT_S[g]  <- length(mut.s)
			      pmr.mut.gene$WT_G1[g]  <- length(wt.g1)
			      pmr.mut.gene$WT_S[g]   <- length(wt.s)
			      if (t == "USP")    pmr.mut.gene$P[g] <- USP.test(test)[[1]]
			      if (t == "Fisher") pmr.mut.gene$P[g] <- fisher.test(test)[[1]]
			      if (t == "Chisq")  pmr.mut.gene$P[g] <- chisq.test(test)[[3]]

			      pmr.mut.gene$G1MR[g] <- pmr.mut.gene$MUT_G1[g] / pmr.mut.gene$WT_G1[g]
			      pmr.mut.gene$SMR[g]  <- pmr.mut.gene$MUT_S[g]  / pmr.mut.gene$WT_S[g]
			      pmr.mut.gene$PMR[g]  <- pmr.mut.gene$SMR[g]    / pmr.mut.gene$G1MR[g]
		    }
	  }
	  pmr.mut.gene <- subset(pmr.mut.gene, P != 0)
	  pmr.mut.gene <- pmr.mut.gene[order(pmr.mut.gene$P, decreasing=F),]
	  if (nrow(pmr.mut.gene) >= 4) {
		    pmr.mut.gene$FDR <- testFDR(pmr.mut.gene$P, "BH")
	  }
	  #pmr <- pmr.mut.gene
	  #annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
	  #pmr.mut.gene <- cbind(annot[rownames(pmr),], pmr.mut.gene[, -1])
	
	  return(pmr.mut.gene)
}

plotPMR <- function(file.name, main.text, xlab.text, ylab.text, pmr.mut.gene.hist, xmax=NA, ymax=NA, size=6, genes=c("TP53"), cols=c(adjustcolor.blue, adjustcolor.red), h=3) {
	  x <- pmr.mut.gene.hist$PMR
	  y <- -log10(pmr.mut.gene.hist$P)
	  #genes <- pmr.mut.gene.hist$MUT
	
	  xlim <- c(0, max(x))
	  if (!is.na(xmax)) xlim <- c(0, xmax)
	  ylim <- c(0, max(y))
	  if (!is.na(ymax)) ylim <- c(0, ymax)
	
	  pdf(paste0(file.name, ".pdf"), height=size, width=size)
	  par(mar=c(5.1, 4.8, 4.1, 1.3))
	  plot(y ~ x, ylab=ylab.text, xlim=xlim, ylim=ylim, xlab=xlab.text, main=main.text, pch=16, cex=2, col="white", cex.axis=1.9, cex.lab=2, cex.main=2.1)
	
	  idx.down   <- which(x < 1)
	  points(x[idx.down], y[idx.down], pch=16, col=cols[1], cex=2)
	  idx.up   <- which(x >= 1)
	  points(x[idx.up],   y[idx.up],   pch=16, col=cols[2], cex=2)   
	
	  abline(v=1, lty=5, col=red, lwd=2)
	  abline(h=h, lty=5, col="black", lwd=2)
	
	  #genes <- c("TP53")
	  par(xpd=T)
	  if (length(genes) != 0) {
		    for (g in 1:length(genes)) {
			      if (!is.na(genes[g])) {
				        idx <- which(pmr.mut.gene.hist$MUT == genes[g])
				
				        if (length(idx) > 0) {
					          if (x[idx] > 1) {
					            	points(x[idx], y[idx], pch=1, col="black", cex=2)
					            	#text(x[idx], y[idx], genes[g], col="black", adj=c(0, -0.5), cex=1.8)
							           text(x[idx], y[idx], genes[g], col="black", pos=3, cex=2)
					          } else {
						            points(x[idx], y[idx], pch=1, col="black", cex=2)
					          	  #text(x[idx], y[idx], genes[g], col="black", adj=c(1, -0.5), cex=1.8)
							           text(x[idx], y[idx], genes[g], col="black", pos=3, cex=2)
					          }
            }
			      }
		    }
	  }
	
	  dev.off()
}

plotdNdSCOR <- function(file.name, main.text, xlab.text, ylab.text, pmr.mut.gene.hist, xmax=NA, ymax=NA, size=6, cols=c(adjustcolor.blue, adjustcolor.red), h=3) {
	  x <- pmr.mut.gene.hist$dNdS_RHO
	  y <- -log10(pmr.mut.gene.hist$dNdS_P)
	  genes <- pmr.mut.gene.hist$MUT
	
	  xlim <- c(min(x), max(x))
	  if (!is.na(xmax)) xlim <- c(0, xmax)
	  ylim <- c(0, max(y))
	  if (!is.na(ymax)) ylim <- c(0, ymax)
	
	  pdf(paste0(file.name, ".pdf"), height=size, width=size)
	  par(mar=c(5.1, 4.7, 4.1, 1.4))
	  plot(y ~ x, ylab=ylab.text, xlim=xlim, ylim=ylim, xlab=xlab.text, main=main.text, pch=16, cex=2, col="white", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	  idx.down   <- which(x < 0)
	  points(x[idx.down], y[idx.down], pch=16, col=cols[1], cex=2)
	  idx.up   <- which(x >= 0)
	  points(x[idx.up],   y[idx.up],   pch=16, col=cols[2], cex=2)   
	
	  #abline(v=1, lty=5, col=red, lwd=2)
	  abline(h=h, lty=5, col="black", lwd=2)
	
	  par(xpd=T)
	  if (length(genes) != 0) {
		    for (g in 1:5) {
			      if (!is.na(genes[g])) {
				        idx <- which(pmr.mut.gene.hist$MUT == genes[g])
				
				        if (length(idx) > 0) {
					          if (x[idx] > 1) {
						            points(x[idx], y[idx], pch=1, col="black", cex=2)
						            #text(x[idx], y[idx], genes[g], col="black", adj=c(0, -0.5), cex=1.8)
						            text(x[idx], y[idx], genes[g], col="black", pos=3, cex=1.8)
					          } else {
						            points(x[idx], y[idx], pch=1, col="black", cex=2)
						            #text(x[idx], y[idx], genes[g], col="black", adj=c(1, -0.5), cex=1.8)
						            text(x[idx], y[idx], genes[g], col="black", pos=3, cex=1.8)
					          }
				        }
			      }
		    }
	  }
	
	  dev.off()
}

plotPMRBCL2 <- function(file.name, main.text, xlab.text, ylab.text, pmr.mut.gene.hist, xmax=NA, ymax=NA, size=6, cols=c(adjustcolor.blue, adjustcolor.red), h=3) {
	  x <- pmr.mut.gene.hist$PMR
	  y <- -log10(pmr.mut.gene.hist$P)
	  genes <- pmr.mut.gene.hist$MUT
 
   xlim <- c(0, max(x))
   if (!is.na(xmax)) xlim <- c(0, xmax)
   ylim <- c(0, max(y))
   if (!is.na(ymax)) ylim <- c(0, ymax)
 
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(y ~ x, ylab=ylab.text, xlim=xlim, ylim=ylim, xlab=xlab.text, main=main.text, pch=16, cex=2, col="white", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   idx.down   <- which(x < 1)
   points(x[idx.down], y[idx.down], pch=16, col=cols[1], cex=2)
   idx.up   <- which(x >= 1)
   points(x[idx.up],   y[idx.up],   pch=16, col=cols[2], cex=2)   
 
   abline(v=1, lty=5, col=red, lwd=2)
   abline(h=h, lty=5, col="black", lwd=2)
   
   par(xpd=T)
   if (length(genes) != 0) {
   	  for (g in 1:5) {
   	  	  if (!is.na(genes[g])) {
   	        if ((genes[g] == "BCL2") || (genes[g] == "CREBBP" || (genes[g] == "TP53") || (genes[g] == "MYC"))) {
               idx <- which(pmr.mut.gene.hist$MUT == genes[g])

               if (length(idx) > 0) {
                  if (x[idx] > 1) {
                     points(x[idx], y[idx], pch=1, col="black", cex=2)
                     if ((genes[g] == "MYC") || (genes[g] == "TP53"))
                        text(x[idx], y[idx], genes[g], col="black", pos=3, cex=1.8)
                     else
                        text(x[idx], y[idx], genes[g], col="black", adj=c(0, -0.5), cex=1.8)
                  } else {
                     points(x[idx], y[idx], pch=1, col="black", cex=2)
                  	  if ((genes[g] == "CREBBP") || (genes[g] == "BCL2"))
                  		    text(x[idx], y[idx], genes[g], col="black", pos=3, cex=1.8)
                  	  else
                        text(x[idx], y[idx], genes[g], col="black", adj=c(1, -0.5), cex=1.8)
                  }
               }
   	        }
   	     }
   	  }
   }
 
   dev.off()
}

plotPMRHIST <- function(file.name, main.text, xlab.text, ylab.text, pmr.mut.gene.hist, xmax=NA, ymax=NA, size=6, cols=c(adjustcolor.blue, adjustcolor.red), h=3) {
	x <- pmr.mut.gene.hist$PMR
	y <- -log10(pmr.mut.gene.hist$P)
	genes <- pmr.mut.gene.hist$MUT
	
	xlim <- c(0, max(x))
	if (!is.na(xmax)) xlim <- c(0, xmax)
	ylim <- c(0, max(y))
	if (!is.na(ymax)) ylim <- c(0, ymax)
	
	pdf(paste0(file.name, ".pdf"), height=size, width=size)
	par(mar=c(5.1, 4.7, 4.1, 1.4))
	plot(y ~ x, ylab=ylab.text, xlim=xlim, ylim=ylim, xlab=xlab.text, main=main.text, pch=16, cex=2, col="white", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
	
	idx.down   <- which(x < 1)
	points(x[idx.down], y[idx.down], pch=16, col=cols[1], cex=2)
	idx.up   <- which(x >= 1)
	points(x[idx.up],   y[idx.up],   pch=16, col=cols[2], cex=2)   
	
	abline(v=1, lty=5, col=red, lwd=2)
	abline(h=h, lty=5, col="black", lwd=2)
	
	par(xpd=T)
	if (length(genes) != 0) {
		for (g in 1:7) {
			if (!is.na(genes[g])) {
				#if ((genes[g] == "BCL2") || (genes[g] == "CREBBP" || (genes[g] == "TP53") || (genes[g] == "MYC"))) {
					idx <- which(pmr.mut.gene.hist$MUT == genes[g])
					
					if (length(idx) > 0) {
						if (x[idx] > 1) {
							points(x[idx], y[idx], pch=1, col="black", cex=2)
							if ((genes[g] == "Uterus-AdenoCA") || (genes[g] == "") || (genes[g] == "") || (genes[g] == ""))
								text(x[idx], y[idx], genes[g], col="black", pos=3, cex=1.8)
							else
								text(x[idx], y[idx], genes[g], col="black", adj=c(0, -0.5), cex=1.8)
						} else {
							points(x[idx], y[idx], pch=1, col="black", cex=2)
							if ((genes[g] == "") || (genes[g] == ""))
								text(x[idx], y[idx], genes[g], col="black", pos=3, cex=1.8)
							else
								text(x[idx], y[idx], genes[g], col="black", adj=c(1, -0.5), cex=1.8)
						}
				#}
				}
			}
		}
	}
	
	dev.off()
}




getEPMR <- function(mut.gene.tpm.hist, samples.mut.tpm.hist) {
   de.mut.tpm <- toTable(0, 21, nrow(mut.gene.tpm.hist), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "G1MR", "SMR", "PMR", "P", "FDR", "MUT_G1_2", "MUT_S_2", "WT_G1_2", "WT_S_2", "G1MR_2", "SMR_2", "PMR_2", "P_2", "FDR_2"))
   de.mut.tpm$MUT <- rownames(mut.gene.tpm.hist)
   rownames(de.mut.tpm) <- de.mut.tpm$MUT
   
   de.mut.tpm$MUT_FREQ <- mapply(x = 1:nrow(de.mut.tpm), function(x) sum(as.numeric(mut.gene.tpm.hist[x,])))
   de.mut.tpm$SAMPLE_FREQ <- mapply(x = 1:nrow(de.mut.tpm), function(x) length(which(mut.gene.tpm.hist[x,] != 0)))
   de.mut.tpm <- subset(de.mut.tpm, SAMPLE_FREQ >= 20)
   
   for (g in 1:nrow(de.mut.tpm)) {
      gene <- de.mut.tpm$MUT[g]
      ids  <- colnames(mut.gene.tpm.hist)[which(mut.gene.tpm.hist[gene,] != 0)]
  
      samples.mut.mut <- subset(samples.mut.tpm.hist, icgc_specimen_id %in% ids)
      samples.mut.wt  <- samples.mut.tpm.hist[setdiff(rownames(samples.mut.tpm.hist), rownames(samples.mut.mut)),]
  
      mut.g1 <- rownames(subset(samples.mut.mut, SORTING == "G1"))
      mut.s  <- rownames(subset(samples.mut.mut, SORTING == "S"))
      wt.g1  <- rownames(subset(samples.mut.wt,  SORTING == "G1"))
      wt.s   <- rownames(subset(samples.mut.wt,  SORTING == "S"))
  
      if ((length(mut.g1) > 0) && (length(mut.s) > 0) && (length(wt.g1) > 0) && (length(wt.s) > 0)) {
         test <- toTable(0, 2, 2, c("G1", "S"))
         rownames(test) <- c("MUT", "WT")
         test[1, 1] <- length(mut.g1)
         test[1, 2] <- length(mut.s)
         test[2, 1] <- length(wt.g1)
         test[2, 2] <- length(wt.s)
   
         de.mut.tpm$MUT_G1[g] <- length(mut.g1)
         de.mut.tpm$MUT_S[g]  <- length(mut.s)
         de.mut.tpm$WT_G1[g]  <- length(wt.g1)
         de.mut.tpm$WT_S[g]   <- length(wt.s)
         #de.mut.tpm$P[g]     <- fisher.test(test)[[1]]
         de.mut.tpm$P[g]      <- chisq.test(test)[[3]]
   
         de.mut.tpm$G1MR[g] <- de.mut.tpm$MUT_G1[g] / de.mut.tpm$WT_G1[g]
         de.mut.tpm$SMR[g]  <- de.mut.tpm$MUT_S[g]  / de.mut.tpm$WT_S[g]
         de.mut.tpm$PMR[g]  <- de.mut.tpm$SMR[g]    / de.mut.tpm$G1MR[g]
   
         ##
         test2 <- toTable(0, 2, 2, c("G1", "S"))
         rownames(test2) <- c("MUT", "WT")
         test2[1, 1] <- median(as.numeric(tpm.gene.mut[gene, mut.g1]))
         test2[1, 2] <- median(as.numeric(tpm.gene.mut[gene, mut.s]))
         test2[2, 1] <- median(as.numeric(tpm.gene.mut[gene, wt.g1]))
         test2[2, 2] <- median(as.numeric(tpm.gene.mut[gene, wt.s]))
   
         if ((test2[1, 1] > 0) && (test2[1, 2] > 0) && (test2[2, 1] > 0) && (test2[2, 2] > 0)) {
            de.mut.tpm$MUT_G1_2[g] <- median(as.numeric(tpm.gene.mut[gene, mut.g1]))
            de.mut.tpm$MUT_S_2[g]  <- median(as.numeric(tpm.gene.mut[gene, mut.s]))
            de.mut.tpm$WT_G1_2[g]  <- median(as.numeric(tpm.gene.mut[gene, wt.g1]))
            de.mut.tpm$WT_S_2[g]   <- median(as.numeric(tpm.gene.mut[gene, wt.s]))
            #de.mut.tpm$P[g]      <- fisher.test(test)[[1]]
            de.mut.tpm$P_2[g]      <- chisq.test(test2)[[3]]
          
            de.mut.tpm$G1MR_2[g] <- de.mut.tpm$MUT_G1_2[g] / de.mut.tpm$WT_G1_2[g]
            de.mut.tpm$SMR_2[g]  <- de.mut.tpm$MUT_S_2[g]  / de.mut.tpm$WT_S_2[g]
            de.mut.tpm$PMR_2[g]  <- de.mut.tpm$SMR_2[g]    / de.mut.tpm$G1MR_2[g]
         }
      }
   }
   de.mut.tpm <- subset(de.mut.tpm, P   != 0)
   de.mut.tpm <- subset(de.mut.tpm, P_2 != 0)
   de.mut.tpm <- de.mut.tpm[order(de.mut.tpm$P, decreasing=F),]
   
   de <- de.mut.tpm
   annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
   de.mut.tpm <- cbind(annot[rownames(de),], de.mut.tpm[, -1])
   
   return(de.mut.tpm)
}

plotPMRDriver <- function(file.name, main.text, xlab.text, ylab.text, de.mut.tpm, genes, size=6, xmax="", ymax="", cols=c(adjustcolor.blue, adjustcolor.red)) {
   x <- de.mut.tpm$PMR
   y <- -log10(de.mut.tpm$P)
   
   xlim <- c(0, max(x))
   if (xmax != "") xlim <- c(0, xmax)
   ylim <- c(0, max(y))
   if (ymax != "") ylim <- c(0, ymax)
   
   pdf(paste0(file.name, ".pdf"), height=size, width=size)
   par(mar=c(5.1, 4.7, 4.1, 1.4))
   plot(y ~ x, ylab=ylab.text, xlim=xlim, ylim=ylim, xlab=xlab.text, main=main.text, pch=1, cex=2, col="white", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)

   idx.down   <- which(de.mut.tpm$PMR < 1)
   points(x[idx.down], y[idx.down], pch=1, col=cols[1], cex=2)
   idx.up   <- which(de.mut.tpm$PMR >= 1)
   points(x[idx.up],   y[idx.up],   pch=1, col=cols[2], cex=2)   
   
   abline(v=1, lty=5, col="black", lwd=2)

   if (length(genes) != 0) {
      for (g in 1:length(genes)) {
         idx <- which(de.mut.tpm$MUT == genes[g])
         #idx <- which(de.mut.tpm$external_gene_name == genes[g])
         
         if (length(idx) > 0) {
            if (x[idx] > 1) {
               points(x[idx], y[idx], pch=1, col=cols[2], cex=2)
               if (genes[g] == "TP53")
                  text(x[idx], y[idx], genes[g], col="black", adj=-0.12, cex=1.8)
               else
                  text(x[idx], y[idx], genes[g], col="black", adj=c(0, -0.5), cex=1.8)
            } else {
               points(x[idx], y[idx], pch=1, col=cols[1], cex=2)
               #text(x[idx], y[idx], genes[g], col="black", adj=c(0, -0.5), cex=1.8)
               text(x[idx], y[idx], genes[g], col="black", adj=c(1, -0.5), cex=1.8)
            }
         }
      }
   }
   
   dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified:
# -----------------------------------------------------------------------------
survESAD <- function(purities2) {
   #purities2 <- purities2[!is.na(purities2$os_days),]
 
   purities2$Tumor_content_pathology <- as.numeric(purities2$Tumor_content_pathology)
 
   ## https://docs.icgc.org/submission/guide/donor-clinical-data-guidelines/#donor-clinical-data-guidelines
   purities2$OS_month <- purities2$os_days / 30.44
 
   purities2$OS_censor <- purities2$death
   purities2$OS_censor <- as.numeric(purities2$OS_censor)
 
   purities2$Survival <- "Short"
   idx <- which(purities2$os_days >= 730)
   purities2[idx,]$Survival <- "Long"
   purities2$Survival <- as.factor(purities2$Survival)
 
   purities2$Response <- gsub("Complete Responder", "Complete", purities2$Response)
   purities2$Response <- gsub("Majorresponder",     "Major",    purities2$Response)
   purities2$Response <- gsub("Minorresponder",     "Minor",    purities2$Response)
   purities2$Response <- as.factor(purities2$Response)
 
   rownames(purities2) <- paste0(purities2$Patient_ID, "_B")
   return(purities2)
}

## https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/Surv
## https://media.nature.com/original/nature-assets/nature/journal/v524/n7563/extref/nature14664-s1.xlsx
survSCLC <- function(phenos, samples, isCensored) {
   overlaps <- intersect(rownames(samples), rownames(phenos))   ## ADD 08/05/19
   phenos <- cbind(samples[overlaps,], phenos[overlaps,])

   ## Stage UICC
   phenos$STAGE <- NA
   phenos[which(phenos$stage_UICC == "I"), ]$STAGE <- "I-II"
   phenos[which(phenos$stage_UICC == "Ia"),]$STAGE <- "I-II"
   phenos[which(phenos$stage_UICC == "Ib"),]$STAGE <- "I-II"
   #phenos[which(phenos$stage_UICC == "IB"),]$STAGE <- "I-II"
   
   phenos[which(phenos$stage_UICC == "II"), ]$STAGE <- "I-II"
   phenos[which(phenos$stage_UICC == "IIa"),]$STAGE <- "I-II"
   #phenos[which(phenos$stage_UICC == "2b"),]$STAGE  <- "I-II"
   phenos[which(phenos$stage_UICC == "IIb"),]$STAGE <- "I-II"
   
   #phenos[which(phenos$stage_UICC == "III"), ]$STAGE <- "III-IV"
   phenos[which(phenos$stage_UICC == "IIIa"),]$STAGE <- "III-IV"
   phenos[which(phenos$stage_UICC == "IIIb"),]$STAGE <- "III-IV"
   
   phenos[which(phenos$stage_UICC == "IV"),]$STAGE <- "III-IV"
   phenos$STAGE <- as.factor(phenos$STAGE)
   
   ## Q4
   #phenos$Q4 <- as.factor(phenos$Q4)
   #phenos$RT <- as.factor(phenos$RT)
   phenos$M2 <- phenos$M2 + 1
    
   ## For Cox regression model
   #phenos$Surgery <- "no"
   #phenos$Surgery[which(phenos$tissue.sampling == "surgical resection")] <- "yes"
   #phenos$Chemotherapy <- phenos$chemotherapy..yes.no.
   #phenos$Radiation    <- phenos$radiation..yes.no.
   
   ## OS_censor
   phenos.surv <- phenos[!is.na(phenos$overall_survival..months.),]
   phenos.surv$OS_month <- phenos.surv$overall_survival..months.
   #phenos.surv$Groups   <- phenos.surv$Stage
 
   phenos.surv$OS_censor <- phenos.surv$Status..at.time.of.last.follow.up.
   phenos.surv$OS_censor <- gsub("dead",  2, phenos.surv$OS_censor)   ## BUG FIX 07/05/19: 0=alive, 1=dead
   phenos.surv$OS_censor <- gsub("alive", 1, phenos.surv$OS_censor)
   phenos.surv$OS_censor <- as.numeric(phenos.surv$OS_censor)
   
   #phenos.surv$SG1 <- "G1"
   #phenos.surv[which(phenos.surv$COR > -0.6503083),]$SG1 <- "S"
   #phenos$SEX   <- phenos$sex
   #phenos$AGE   <- phenos$age
   phenos$SORTING <- NA
    
   if (isCensored) {
      return(phenos.surv)
   } else
      return(phenos)
}

survNBL <- function(samples.nbl) {
   ## OS_censor
   samples.nbl$OS_month <- samples.nbl$OS_d / 30.44

   samples.nbl$OS_censor <- samples.nbl$OS_bin
   samples.nbl$OS_censor <- as.numeric(samples.nbl$OS_censor)
 
   samples.nbl$AGE <- as.numeric(samples.nbl$pat_age)
   samples.nbl$SORTING <- NA
 
   return(samples.nbl)
}

## https://docs.icgc.org/submission/guide/donor-clinical-data-guidelines/#donor-clinical-data-guidelines
## https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/Surv
## https://media.nature.com/original/nature-assets/nature/journal/v524/n7563/extref/nature14664-s1.xlsx
survICGC <- function(samples.hist) {
   samples.hist <- samples.hist[!is.na(samples.hist$donor_vital_status),]
   samples.hist <- subset(samples.hist, donor_vital_status != "")

   samples.hist <- removeNA(samples.hist, "donor_survival_time")
   samples.hist <- subset(samples.hist, donor_survival_time != 0)
 
   samples.hist <- removeNA(samples.hist, "donor_age_at_diagnosis")
   samples.hist <- subset(samples.hist, donor_age_at_diagnosis != 0)
   
   if (nrow(samples.hist) != 0) {
      samples.hist$OS_censor <- samples.hist$donor_vital_status
      samples.hist$OS_censor <- gsub("alive",    0, samples.hist$OS_censor)   ## BUG FIX 07/05/19: 0=alive, 1=dead
      samples.hist$OS_censor <- gsub("deceased", 1, samples.hist$OS_censor)
      samples.hist$OS_censor <- as.numeric(samples.hist$OS_censor)
 
      ## https://docs.icgc.org/submission/guide/donor-clinical-data-guidelines/#donor-clinical-data-guidelines
      samples.hist$OS_month <- samples.hist$donor_survival_time / 30.44
      
      samples.hist$donor_sex <- as.factor(samples.hist$donor_sex)
      
      ## Variables for processing coxph() and ggforest()
      samples.hist$SEX <- samples.hist$donor_sex
      samples.hist$SEX <- gsub("female", "F", samples.hist$SEX)
      samples.hist$SEX <- gsub("male",   "M", samples.hist$SEX)
      samples.hist$AGE <- samples.hist$donor_age_at_diagnosis
      samples.hist$SORTING <- NA
      
      return(samples.hist)
   } else
      return(samples.hist[1,][-1,])
}

getSurvfitPvals <- function(pvals, samples.surv.sclc) {
   pvals <- c()
   
   for (s in 1:nrow(samples.surv.sclc)) {
      samples.surv.sclc$SORTING <- "G1"
      idx <- which(samples.surv.sclc$COR >= samples.surv.sclc$COR[s])
      if (length(idx) != 0)
         samples.surv.sclc[idx,]$SORTING <- "S"
      samples.surv.sclc$SORTING <- as.factor(samples.surv.sclc$SORTING)
  
      if (length(unique(samples.surv.sclc$SORTING)) != 1){
         fit <- survfit(Surv(OS_month, OS_censor) ~ SORTING, data=samples.surv.sclc)
         pvals <- c(pvals, surv_pvalue(fit)$pval)
      } else {
         pvals <- c(pvals, NA)
      }
   }
   
   return(pvals)
}

setProliferation <- function(samples.sg1, cor) {
   samples.sg1$SG1 <- "G1"
   idx <- which(samples.sg1$COR >= cor)
   if (length(idx) != 0)
      samples.sg1[idx,]$SG1 <- "S"
   samples.sg1$SG1 <- as.factor(samples.sg1$SG1)
   
   return(samples.sg1)
}

# -----------------------------------------------------------------------------
# Methods: 
# Last Modified: 28/03/22
# -----------------------------------------------------------------------------
## https://www.rdocumentation.org/packages/survminer/versions/0.4.3/topics/surv_pvalue
get_surv_pvalue <- function(fit) {
   #surv_pvalue(fit, method="log-rank")$pval.txt
   pvalue <- surv_pvalue(fit, method="log-rank")$pval
 
   if (surv_pvalue(fit, method="log-rank")$pval < 0.001) {
      return("log-rank p-value < 0.001")
   } else
      return(paste0("log-rank p-value = ", round0(pvalue, digits=3)))
}
