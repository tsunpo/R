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
   
   return(paste0("log-rank p-value = ", round0(pvalue, digits=2)))
}

plotSurvfitRT <- function(fit, file.name, main.text, isRT) {
   cols    <- c("blue", "skyblue3", "lightcoral", "red")
   extension <- "Q4"
   legends <- c()
   for (q in 1:4)
      legends <- c(legends, paste0("Q", q, " (n=", fit[[1]][q], ")"))
   if (isRT) {
      cols    <- c("blue", "red")
      extension <- "RT"
      legends <- c()
      for (q in 1:2)
         legends <- c(legends, paste0("M", q, " (n=", fit[[1]][q], ")"))
   }
 
   pdf(paste0(file.name, "_", extension, ".pdf"), height=6, width=6)
   plot(fit, ylim=c(0, 1), xlab="Months", ylab="Overall survival (OS)", col=cols, main=main.text)
   legend("topright", legend=legends, lwd=1, col=cols)
   mtext(get_surv_pvalue(fit), cex=1.2, line=0.3)
   dev.off()
}
