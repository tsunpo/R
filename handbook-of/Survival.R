# =============================================================================
# Library      : Survival Probability
# Name         : handbook-of/Survival.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 09/08/21; 03/05/19
# =============================================================================

# -----------------------------------------------------------------------------
# Methods: Kaplan-Meier survival analysis
# Link(s): https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_survival_analysis.pdf
#          http://www.stat.columbia.edu/~madigan/W2025/notes/survival.pdf
#          https://support.minitab.com/en-us/minitab/18/help-and-how-to/modeling-statistics/reliability/supporting-topics/reliability-metrics/what-is-the-survival-probability/
# Last Modified: 03/05/19
# -----------------------------------------------------------------------------
plotOS <- function(file.name, main.text, xlab.text, ylab.text, x, y, p, rho, line=2.75, lwd=4) {
   ymax <- -log10(p)
   ymax <- ymax + ymax/5.5
   xmin <- min(x)
    
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(y ~ x, ylim=c(0, ymax), ylab="", xaxt="n", xlab=xlab.text, col="dimgray", main=main.text, pch=1, cex=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   points(rho, -log10(p), pch=1, col=green, cex=2, lwd=lwd)
   abline(v=rho, lty=5, lwd=5, col=green)
   text(xmin, ymax - ymax/10, expression(italic("    P")~"="), pos=3, cex=1.7)
   text(xmin, ymax - ymax/10, paste0("                           ", scientific(surv_pvalue(fit)$pval, digits=2)), pos=3, cex=1.7)
   
   axis(side=1, at=seq(-0.4, 0.4, by=0.4), labels=c(-0.4, 0, 0.4), cex.axis=1.7)
   axis(side=1, at=rho, labels=round(rho, 2), cex.axis=1.7, col.axis=green, font.axis=2)
 
   mtext(ylab.text, side=2, line=2.4, cex=1.8)
   dev.off()
}

plotSurvfit <- function(fit, file.name, main.text, legends.text, cols) {
   legends <- c()
   for (l in 1:length(legends.text))
      legends <- c(legends, paste0(legends.text[l], " (n=", fit[[1]][l], ")"))
 
   pdf(paste0(file.name, ".pdf"), height=5, width=6)
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="", col=cols, main=main.text[1], mark.time=T, lwd=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   legend("topright", legend=legends, lwd=4, col=cols, cex=1.7)
   #mtext(main.text[2], cex=1.25, line=0.2)
   
   mtext("% OS", side=2, line=2.85, cex=1.7)
   
   max <- max(fit[[2]])
   text(max/5.5, 0.18, "                HR = 1.90 (1.61 - 2.25)", cex=1.7)
   text(max/5.5, 0.07, expression(italic('P')~"="~"                "), cex=1.7)
   text(max/5.5, 0.07, paste0("      ", scientific(surv_pvalue(fit)$pval, digits=2)), cex=1.7)
   #legend("bottomleft", expression(italic('P')~"="~scientific(surv_pvalue(fit)$pval)), bty="n", cex=1.6)
   #legend("bottomleft", expression(italic('P')~"="~"2.62e-18"), bty="n", cex=1.6)
   
   dev.off()
}

plotSurvfit55 <- function(fit, file.name, main.text, legends.text, cols) {
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
getPMR <- function(mut.gene.test, samples.test) {
   de.mut.tpm <- toTable(0, 12, nrow(mut.gene.test), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "P", "FDR", "G1MR", "SMR", "PMR"))
   de.mut.tpm$MUT <- rownames(mut.gene.test)
   rownames(de.mut.tpm) <- de.mut.tpm$MUT
 
   de.mut.tpm$MUT_FREQ <- mapply(x = 1:nrow(de.mut.tpm), function(x) sum(as.numeric(mut.gene.test[x,])))
   de.mut.tpm$SAMPLE_FREQ <- mapply(x = 1:nrow(de.mut.tpm), function(x) length(which(mut.gene.test[x,] != 0)))
   de.mut.tpm <- subset(de.mut.tpm, SAMPLE_FREQ >= 20)
 
   for (g in 1:nrow(de.mut.tpm)) {
      gene <- de.mut.tpm$MUT[g]
      ids  <- colnames(mut.gene.test)[which(mut.gene.test[gene,] != 0)]
  
      samples.mut.mut <- subset(samples.test, icgc_specimen_id %in% ids)
      samples.mut.wt  <- samples.test[setdiff(rownames(samples.test), rownames(samples.mut.mut)),]
  
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
      }
   }
   de.mut.tpm <- subset(de.mut.tpm, P != 0)
   de.mut.tpm <- de.mut.tpm[order(de.mut.tpm$P, decreasing=F),]
 
   de <- de.mut.tpm
   annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
   de.mut.tpm <- cbind(annot[rownames(de),], de.mut.tpm[, -1])
 
   return(de.mut.tpm)
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

plotPMR <- function(file.name, main.text, xlab.text, ylab.text, x, y, line=2.4) {
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(y ~ x, ylab="", xlab=xlab.text, main=main.text, pch=1, cex=2, col="dimgray", cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
 
   abline(v=1, lty=5, col="black", lwd=2)

   if (nrow(genes) != 0) {
      for (g in 1:nrow(genes)) {
         gene <- subset(de, external_gene_name == genes[g,]$GENE)
         gene <- cbind(gene, genes[g,])
     
         if (nrow(gene) > 0) {
            points(gene$LOG2_FC, gene$log10P, pch=1, col="black", cex=1.4)
      
            if (!is.na(gene$ADJ_1))
               if (is.na(gene$ADJ_2))
                  text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=gene$ADJ_1, cex=1.2)
               else
                  text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(gene$ADJ_1, gene$ADJ_2), cex=1.2)
            else
               if (gene$LOG2_FC > 0)
                  text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(0, -0.5), cex=1.2)
               else
                  text(gene$LOG2_FC, gene$log10P, paste(genes[g,]$GENE, genes[g,]$ADJ_3), col="black", adj=c(1, -0.5), cex=1.2)
         } else
            print(genes[g,])
      }
   }
   
   mtext(ylab.text, side=2, line=line, cex=1.8)
   dev.off()
}

# -----------------------------------------------------------------------------
# 
# Last Modified:
# -----------------------------------------------------------------------------
## https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/Surv
## https://media.nature.com/original/nature-assets/nature/journal/v524/n7563/extref/nature14664-s1.xlsx
survSCLC <- function(phenos, samples, isCensored) {
   overlaps <- intersect(rownames(samples), rownames(phenos))   ## ADD 08/05/19
   phenos <- cbind(phenos[overlaps,], samples[overlaps,])

   ## Stage UICC
   phenos$Stage <- NA
   phenos[which(phenos$stage_UICC == "I"), ]$Stage <- 1
   phenos[which(phenos$stage_UICC == "Ia"),]$Stage <- 1
   phenos[which(phenos$stage_UICC == "Ib"),]$Stage <- 1
   #phenos[which(phenos$stage_UICC == "IB"),]$Stage <- 1
   
   phenos[which(phenos$stage_UICC == "II"), ]$Stage <- 2
   phenos[which(phenos$stage_UICC == "IIa"),]$Stage <- 2
   phenos[which(phenos$stage_UICC == "IIb"),]$Stage <- 2
   
   #phenos[which(phenos$stage_UICC == "III"), ]$Stage <- 3
   phenos[which(phenos$stage_UICC == "IIIa"),]$Stage <- 3
   phenos[which(phenos$stage_UICC == "IIIb"),]$Stage <- 3
   
   phenos[which(phenos$stage_UICC == "IV"),]$Stage <- 4
   #phenos$Stage <- as.factor(phenos$Stage)
   
   ## Q4
   #phenos$Q4 <- as.factor(phenos$Q4)
   #phenos$RT <- as.factor(phenos$RT)
   
   ## For Cox regression model
   phenos$Surgery <- "no"
   phenos$Surgery[which(phenos$tissue.sampling == "surgical resection")] <- "yes"
   phenos$Sex          <- phenos$sex
   phenos$Chemotherapy <- phenos$chemotherapy..yes.no.
   phenos$Radiation    <- phenos$radiation..yes.no.
   
   ## OS_censor
   phenos.surv <- phenos[!is.na(phenos$overall_survival..months.),]
   phenos.surv$OS_month <- phenos.surv$overall_survival..months.
   phenos.surv$Groups   <- phenos.surv$Stage
 
   phenos.surv$OS_censor <- phenos.surv$Status..at.time.of.last.follow.up.
   phenos.surv$OS_censor <- gsub("dead",  1, phenos.surv$OS_censor)   ## BUG FIX 07/05/19: 0=alive, 1=dead
   phenos.surv$OS_censor <- gsub("alive", 0, phenos.surv$OS_censor)
   phenos.surv$OS_censor <- as.numeric(phenos.surv$OS_censor)
   
   phenos.surv$SG1 <- "G1"
   phenos.surv[which(phenos.surv$COR > -0.6503083),]$SG1 <- "S"
   
   if (isCensored) {
      return(phenos.surv)
   } else
      return(phenos)
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

getSurvfitPvals <- function(samples.surv) {
   pvals <- c()
   
   for (s in 1:nrow(samples.surv)) {
      samples.surv$SORTING <- "G1"
      idx <- which(samples.surv$COR >= samples.surv$COR[s])
      if (length(idx) != 0)
         samples.surv[idx,]$SORTING <- "S"
      samples.surv$SORTING <- as.factor(samples.surv$SORTING)
  
      fit <- survfit(Surv(OS_month, OS_censor) ~ SORTING, data=samples.surv)
      if (length(unique(samples.surv$SORTING)) != 1){
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
