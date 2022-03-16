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
library(survival)
library(survminer)
library(broom)
library(grid)

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
 
   pdf(paste0(file.name, ".pdf"), height=5, width=5)
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="Overall survival (%)", col=cols, main=main.text[1], mark.time=T, lwd=1.5, cex.axis=1.3, cex.lab=1.4, cex.main=1.5)
   legend("topright", legend=legends, lwd=2, col=cols, cex=1.25)
   mtext(main.text[2], cex=1.25, line=0.2)
   text(0, 0, get_surv_pvalue(fit), adj= c(-0.05, 0.1), col="black", cex=1.25)
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

## https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/Surv
## https://media.nature.com/original/nature-assets/nature/journal/v524/n7563/extref/nature14664-s1.xlsx
survICGC <- function(samples.hist) {
   samples.hist <- removeNA(samples.hist, "donor_survival_time")
   samples.hist <- subset(samples.hist, donor_survival_time != 0)
   samples.hist <- removeNA(samples.hist, "donor_age_at_diagnosis")
   samples.hist <- subset(samples.hist, donor_age_at_diagnosis != 0)
   
   if (nrow(samples.hist) != 0) {
      samples.hist$donor_survival_time <- samples.hist$donor_survival_time / 30.44
      #samples.hist <- samples.hist[!is.na(phenos.surv$donor_vital_status),]
      
      samples.hist$OS_month <- samples.hist$donor_survival_time
      
      samples.hist$OS_censor <- samples.hist$donor_vital_status
      samples.hist$OS_censor <- gsub("deceased", 1, samples.hist$OS_censor)   ## BUG FIX 07/05/19: 0=alive, 1=dead
      samples.hist$OS_censor <- gsub("alive",    0, samples.hist$OS_censor)
      samples.hist$OS_censor <- as.numeric(samples.hist$OS_censor)
 
      samples.hist$donor_sex <- as.factor(samples.hist$donor_sex)
      #samples.hist$M2 <- as.factor(samples.hist$M2)
      #samples.hist$Q4 <- as.factor(samples.hist$Q4)
      
      #samples.hist$In_silico <- samples.hist$COR
      #samples.hist$G1_vs_S <- samples.hist$SG1
      samples.hist$SEX     <- samples.hist$donor_sex
      samples.hist$AGE     <- samples.hist$donor_age_at_diagnosis
   
      return(samples.hist)
   } else
      return(samples.hist[1,][-1,])
}

# -----------------------------------------------------------------------------
# Methods:  analysis
# Last Modified: 09/08/21
# -----------------------------------------------------------------------------
