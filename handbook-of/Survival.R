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
## https://www.rdocumentation.org/packages/survminer/versions/0.4.3/topics/surv_pvalue
get_surv_pvalue <- function(fit) {
   #surv_pvalue(fit, method="log-rank")$pval.txt
   pvalue <- surv_pvalue(fit, method="log-rank")$pval
   
   if (surv_pvalue(fit, method="log-rank")$pval < 0.001) {
      return("log-rank p-value < 0.001")
   } else
      return(paste0("log-rank p-value = ", round0(pvalue, digits=3)))
}

plotSurvfit <- function(fit, file.name, main.text, legends.text, cols) {
   legends <- c()
   for (l in 1:length(legends.text))
      legends <- c(legends, paste0(legends.text[l], " (n=", fit[[1]][l], ")"))
 
   pdf(paste0(file.name, ".pdf"), height=5, width=6)
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="", col=cols, main=main.text[1], mark.time=T, lwd=2, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
   legend("topright", legend=legends, lwd=4, col=cols, cex=1.7)
   #mtext(main.text[2], cex=1.25, line=0.2)
   
   mtext("OS [%]", side=2, line=2.85, cex=1.7)
   
   max <- max(fit[[2]])
   text(max/5.5, 0.18, "                HR = 1.90 (1.61 - 2.25)", cex=1.7)
   text(max/5.5, 0.07, expression(italic('P')~"="~"                "), cex=1.7)
   text(max/5.5, 0.07, paste0("      ", scientific(surv_pvalue(fit)$pval, digits=2)), cex=1.7)
   #legend("bottomleft", expression(italic('P')~"="~scientific(surv_pvalue(fit)$pval)), bty="n", cex=1.6)
   #legend("bottomleft", expression(italic('P')~"="~"2.62e-18"), bty="n", cex=1.6)
   
   dev.off()
}

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
