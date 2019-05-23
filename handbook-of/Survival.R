# =============================================================================
# Library      : Survival Probability
# Name         : handbook-of/Survival.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 03/05/19
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
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="Overall survival (%)", col=cols, main=main.text[1], mark.time=T)
   legend("topright", legend=legends, lwd=1, col=cols)
   mtext(main.text[2], cex=1.2, line=0.3)
   text(0, 0, get_surv_pvalue(fit), adj= c(-0.05, 0.1), col="black")
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
   
   phenos[which(phenos$stage_UICC == "II"), ]$Stage <- 1
   phenos[which(phenos$stage_UICC == "IIa"),]$Stage <- 1
   phenos[which(phenos$stage_UICC == "IIb"),]$Stage <- 1
   
   #phenos[which(phenos$stage_UICC == "III"), ]$Stage <- 2
   phenos[which(phenos$stage_UICC == "IIIa"),]$Stage <- 2
   phenos[which(phenos$stage_UICC == "IIIb"),]$Stage <- 2
   
   phenos[which(phenos$stage_UICC == "IV"),]$Stage <- 2
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
   
   if (isCensored) {
      return(phenos.surv)
   } else
      return(phenos)
}
