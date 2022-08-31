# =============================================================================
# Manuscript   : 
# Chapter      : Reconstruction of replication timing profile in tumour cells
# Name         : manuscripts/replication/cll-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 26/02/19
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Survival.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "SCLC"
base <- tolower(BASE)

wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "George 2015")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
phenos  <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# 
# Last Modified: 15/04/22
# -----------------------------------------------------------------------------
samples.surv.sclc <- survSCLC(phenos, samples, isCensored=T)
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

idx <- which(is.na(pvals))
s <- which(pvals == min(pvals[-idx]))
rho.sclc <- samples.surv.sclc$COR[s]
samples.surv.sclc$SORTING <- "G1"
idx <- which(samples.surv.sclc$COR >= rho.sclc)
if (length(idx) != 0)
   samples.surv.sclc[idx,]$SORTING <- "S"
samples.surv.sclc$SORTING <- as.factor(samples.surv.sclc$SORTING)

fit <- survfit(Surv(OS_month, OS_censor) ~ SORTING, data=samples.surv.sclc)
pval <- surv_pvalue(fit)$pval
file.name <- file.path(paste0(wd.rt.plots, "/hists/survfit_in-silico_samples.surv.hist_", text.SCLC, "G1-S"))
plotSurvfit(fit, file.name, text.SCLC, legend.labs=c("Resting", "Proliferative"), name="SORTING", strata=c("G1", "S"), cols=c(blue, red), size=5)

file.name <- file.path(paste0(wd.rt.plots, "/hists/survfit_in-silico_samples.surv.hist_", text.SCLC, "S-G1"))
plotSurvfit(fit, file.name, text.SCLC, legend.labs=c("Proliferative", "Resting"), name="SORTING", strata=c("S", "G1"), cols=c(red, blue), size=5)

idx <- which(is.na(pvals))
x <- samples.surv.sclc$COR[-idx]
y <- -log10(pvals[-idx])
file.name <- paste0(wd.rt.plots, "/hists/correlation_in-silico_P_KM_OS_", text.SCLC)
plotOS(file.name, text.SCLC, text.In.silico, text.Log10.P, x, y, pvals[s], rho.sclc, lwd=3)

res.cox <- coxph(Surv(OS_month, OS_censor) ~ SORTING + STAGE + SEX + AGE, data=samples.surv.sclc)
pdf(paste0(wd.rt.plots, "/hists/hazard_in-silico_samples.surv.hist_", text.SCLC, "_.pdf"), height=2.7, width=5)
ggforest(res.cox, data=samples.surv.sclc, main=text.SCLC, cpositions=c(0.02, 0.15, 0.35), fontsize=2.5)
dev.off()

###
##
size <- nrow(samples.surv.sclc)
sigs <- 0
for (p in 1:10000) {
   idx <- sample(1:size, size, replace=F)
   fit <- survfit(Surv(OS_month[idx], OS_censor[idx]) ~ SORTING, data=samples.surv.sclc)
 
   if (surv_pvalue(fit)$pval <= pval)
      sigs <- sigs + 1
}
sigs/10000

###
##
test <- toTable(0, 2, 2, c("I-II", "III-IV"))
rownames(test) <- c("S", "G1")

test[1, 1] <- nrow(subset(subset(samples.surv.sclc, SORTING == "S"), STAGE == "I-II"))
test[1, 2] <- nrow(subset(subset(samples.surv.sclc, SORTING == "S"), STAGE == "III-IV"))
test[2, 1] <- nrow(subset(subset(samples.surv.sclc, SORTING == "G1"), STAGE == "I-II"))
test[2, 2] <- nrow(subset(subset(samples.surv.sclc, SORTING == "G1"), STAGE == "III-IV"))
fisher.test(test)[[1]]
# [1] 0.6590286

file.name <- paste0("boxplot_STAGAE_SORTING")
ylab.text <- expression(italic('in silico')~"sorting [rho]")
wd <- "/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/expression/kallisto/sclc-tpm-de/plots"
plotBox02(wd, file.name, subset(samples.surv.sclc, STAGE == "I-II")$COR, subset(samples.surv.sclc, STAGE == "III-IV")$COR, "Lung-SCLC", names=c("I-II", "III-IV"), cols=c("lightgray", "dimgray"), ylab.text)



###
##
fit <- survfit(Surv(OS_month, OS_censor) ~ STAGE, data=samples.surv.sclc)
pval <- surv_pvalue(fit)$pval
file.name <- file.path(paste0(wd, "/survfit_STAGE_samples.surv.hist_", text.SCLC))
plotSurvfit55(fit, file.name, text.SCLC, c("I-II", "III-IV"), c("lightgray", "dimgray"))

idx <- which(is.na(pvals))
x <- samples.surv.sclc$COR[-idx]
y <- -log10(pvals[-idx])
file.name <- paste0(wd.rt.plots, "/hists/correlation_in-silico_P_KM_OS_", text.SCLC)
plotOS(file.name, text.SCLC, text.In.silico, text.Log10.P, x, y, pvals[s], rho, lwd=3)

res.cox <- coxph(Surv(OS_month, OS_censor) ~ STAGE + SORTING + SEX + AGE, data=samples.surv.sclc)
pdf(paste0(wd.rt.plots, "/hazard_in-silico_samples.surv.hist_", text.SCLC, "_.pdf"), height=2.7, width=5)
ggforest(res.cox, data=samples.surv.sclc, main=text.SCLC, cpositions=c(0.02, 0.15, 0.35), fontsize=2.5)
dev.off()

# -----------------------------------------------------------------------------
# 
# Last Modified: 27/09/19
# -----------------------------------------------------------------------------
phenos.surv <- survSCLC(phenos, samples, isCensored=T)

test <- toTable(0, 4, 4, c("Q4", "Q3", "Q2", "Q1"))
rownames(test) <- c("IV", "III", "II", "I")

test[1, 1] <- nrow(subset(subset(phenos.surv, Q4 == 4), Stage == 4))
test[1, 2] <- nrow(subset(subset(phenos.surv, Q4 == 3), Stage == 4))
test[1, 3] <- nrow(subset(subset(phenos.surv, Q4 == 2), Stage == 4))
test[1, 4] <- nrow(subset(subset(phenos.surv, Q4 == 1), Stage == 4))
test[2, 1] <- nrow(subset(subset(phenos.surv, Q4 == 4), Stage == 3))
test[2, 2] <- nrow(subset(subset(phenos.surv, Q4 == 3), Stage == 3))
test[2, 3] <- nrow(subset(subset(phenos.surv, Q4 == 2), Stage == 3))
test[2, 4] <- nrow(subset(subset(phenos.surv, Q4 == 1), Stage == 3))
test[3, 1] <- nrow(subset(subset(phenos.surv, Q4 == 4), Stage == 2))
test[3, 2] <- nrow(subset(subset(phenos.surv, Q4 == 3), Stage == 2))
test[3, 3] <- nrow(subset(subset(phenos.surv, Q4 == 2), Stage == 2))
test[3, 4] <- nrow(subset(subset(phenos.surv, Q4 == 1), Stage == 2))
test[4, 1] <- nrow(subset(subset(phenos.surv, Q4 == 4), Stage == 1))
test[4, 2] <- nrow(subset(subset(phenos.surv, Q4 == 3), Stage == 1))
test[4, 3] <- nrow(subset(subset(phenos.surv, Q4 == 2), Stage == 1))
test[4, 4] <- nrow(subset(subset(phenos.surv, Q4 == 1), Stage == 1))
fisher.test(test)[[1]]
# [1] 0.7301972

## Surgery
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("Yes", "No")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Surgery == "yes"))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Surgery == "yes"))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Surgery == "no"))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Surgery == "no"))
fisher.test(test)[[1]]
# [1] 0.06256722
# > test
#     M2 M1
# Yes 35 45
# No   9  3

## Surgery (Yes) -> UICC stage
phenos.surv <- subset(phenos.surv, Surgery == "yes")

test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("II", "I")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 4))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 4))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 2))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 2))
fisher.test(test)[[1]]
# [1] 1
# > test
#    M2 M1
# II 14 17
# I  21 27

## Surgery (No) -> UICC stage
phenos.surv <- subset(phenos.surv, Surgery == "no")

test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("II", "I")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 4))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 4))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 2))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 2))
fisher.test(test)[[1]]
# [1] 1
# > test
#    M2 M1
# II  9  3
# I   0  0


## UICC stage
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("II", "I")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 4))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 4))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 2))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 2))
fisher.test(test)[[1]]
# [1] 0.4043607
# > test
#    M2 M1
# II 23 20
# I  21 27

## UICC stage (Female)
phenos.surv <- subset(phenos.surv, sex == "female")
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("II", "I")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 4))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 4))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 2))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 2))
fisher.test(test)[[1]]
# [1] 0.181141
# > test
#    M2 M1
# II 10  6
# I   7 12

## UICC stage (Male)
phenos.surv <- subset(phenos.surv, sex == "male")
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("II", "I")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 4))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 3)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 4))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 1), Stage == 2))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Stage == 1)) + nrow(subset(subset(phenos.surv, M2 == 0), Stage == 2))
fisher.test(test)[[1]]
# [1] 1
# > test
#    M2 M1
# II 13 14
# I  14 15
x
# -----------------------------------------------------------------------------
# Cox regression model
# Link(s): http://www.sthda.com/english/wiki/cox-proportional-hazards-model
#          http://www.sthda.com/english/wiki/survminer-0-3-0
#          http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/BS704_Survival6.html
# Last Modified: 08/05/19
# -----------------------------------------------------------------------------
phenos.surv <- survSCLC(phenos, samples, isCensored=T)

phenos.surv$Stage[which(phenos.surv$Stage == 1)] <- "II"
phenos.surv$Stage[which(phenos.surv$Stage == 2)] <- "I"
phenos.surv$RT[which(phenos.surv$RT == 0)] <- "M1"
phenos.surv$RT[which(phenos.surv$RT == 1)] <- "M2"
phenos.surv$RT <- phenos.surv$Q4
#phenos.surv$RT <- as.factor(phenos.surv$RT)
phenos.surv$RT[which(phenos.surv$Q4 == 1)] <- "Q1"
phenos.surv$RT[which(phenos.surv$Q4 == 2)] <- "Q2"
phenos.surv$RT[which(phenos.surv$Q4 == 3)] <- "Q3"
phenos.surv$RT[which(phenos.surv$Q4 == 4)] <- "Q4"
phenos.surv$RT <- as.factor(phenos.surv$RT)
phenos.surv$Sex[which(phenos.surv$Sex == "male")] <- "boy"
phenos.surv$Sex[which(phenos.surv$Sex == "female")] <- "girl"
#phenos.surv$Sex <- as.factor(phenos.surv$Sex)

#phenos.surv <- subset(phenos.surv, Surgery == "yes")
#phenos.surv <- subset(phenos.surv, sex == "male")
#phenos.surv <- subset(phenos.surv, Chemotherapy == "yes")
#phenos.surv <- subset(phenos.surv, Stage == "I")
#phenos.surv <- subset(phenos.surv, RT == "Q2+Q3")

res.cox <- coxph(Surv(OS_month, OS_censor) ~ Surgery + Stage + RT + Chemotherapy + Radiation + Sex, data=phenos.surv)
#res.cox <- coxph(Surv(OS_month, OS_censor) ~ Stage + RT + Radiation + Chemotherapy + Sex, data=phenos.surv)
#res.cox <- coxph(Surv(OS_month, OS_censor) ~ RT + Stage + Chemotherapy + Radiation, data=phenos.surv)
ggforest(res.cox, data=phenos.surv, main = "Hazard ratio in SCLC patients", cpositions = c(0.02, 0.22, 0.4))
# > summary(res.cox)
# Call:
#  coxph(formula = Surv(OS_month, OS_censor) ~ Surgical + Chemotherapy + 
#         Radiation + Q4 + Sex + Stage, data = phenos.surv)
# 
# n= 73, number of events= 46 
# (19 observations deleted due to missingness)
# 
# coef exp(coef) se(coef)      z Pr(>|z|)    
# Surgicalyes     -2.56216   0.07714  0.60743 -4.218 2.46e-05 ***
#  Chemotherapyyes -0.60082   0.54836  0.43068 -1.395  0.16301    
# Radiationyes    -0.34751   0.70644  0.39709 -0.875  0.38149    
# Q4Q2            -0.73596   0.47904  0.43892 -1.677  0.09359 .  
# Q4Q3            -0.90180   0.40584  0.51654 -1.746  0.08084 .  
# Q4Q4            -0.18829   0.82838  0.46183 -0.408  0.68350    
# Sexmale          0.21805   1.24364  0.33206  0.657  0.51142    
# StageIII-IV      0.99061   2.69289  0.35669  2.777  0.00548 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Surgicalyes       0.07714    12.9638   0.02345    0.2537
# Chemotherapyyes   0.54836     1.8236   0.23576    1.2755
# Radiationyes      0.70644     1.4155   0.32439    1.5384
# Q4Q2              0.47904     2.0875   0.20266    1.1324
# Q4Q3              0.40584     2.4640   0.14746    1.1170
# Q4Q4              0.82838     1.2072   0.33506    2.0480
# Sexmale           1.24364     0.8041   0.64870    2.3842
# StageIII-IV       2.69289     0.3713   1.33845    5.4179
# 
# Concordance= 0.75  (se = 0.048 )
# Rsquare= 0.364   (max possible= 0.99 )
# Likelihood ratio test= 33.05  on 8 df,   p=6.023e-05
# Wald test            = 34.16  on 8 df,   p=3.797e-05
# Score (logrank) test = 48.33  on 8 df,   p=8.559e-08

# -----------------------------------------------------------------------------
# Kaplan-Meier survival analysis (OS~kmeans): RT and Stage
# Last Modified: 07/05/19
# -----------------------------------------------------------------------------
phenos.surv <- survSCLC(phenos, samples)

phenos.surv <- subset(phenos.surv, tissue.sampling == "surgical resection")
phenos.surv <- subset(phenos.surv, sex == "male")
#phenos.surv <- subset(phenos.surv, Q2 == "Q1+Q4")
phenos.surv <- subset(phenos.surv, chemotherapy..yes.no. == "yes")
#phenos.surv <- subset(phenos.surv, radiation..yes.no. == "yes")
phenos.surv <- subset(phenos.surv, Stage == 2)
#phenos.surv <- subset(phenos.surv, RT == 0)
patient.text <- paste0(BASE, ", surgical, male, chemo, III-IV patients (n=", nrow(phenos.surv), ")")
#patient.text <- paste0(BASE, ", surgical, female, chemo patients (n=", nrow(phenos.surv), ")")

##
phenos.surv$Groups <- phenos.surv$Stage
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Stage"))
main.text <- c(patient.text, "Stage UICC")
plotSurvfit(fit, file.name, main.text, c("  I-II", "III-IV"), c("blue", "red"))

##
phenos.surv$Groups <- phenos.surv$sex
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Chemo_Sex"))
main.text <- c(patient.text, "Sex")
plotSurvfit(fit, file.name, main.text, c("F", "M"), c("lightcoral", "skyblue3"))


##
phenos.surv$Groups <- phenos.surv$Q4
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Male_Chemo_Stage-III-IV_Q4"))
main.text <- c(patient.text, "Replication timing")
plotSurvfit(fit, file.name, main.text, c("Q1", "Q2", "Q3", "Q4"), c("blue", "skyblue3", "lightcoral", "red"))

##
phenos.surv$Groups[which(phenos.surv$Q4 == 1)] <- "Q1+Q4"
phenos.surv$Groups[which(phenos.surv$Q4 == 2)] <- "Q2+Q3"
phenos.surv$Groups[which(phenos.surv$Q4 == 3)] <- "Q2+Q3"
phenos.surv$Groups[which(phenos.surv$Q4 == 4)] <- "Q1+Q4"
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Male_Chemo_Stage-III-IV_Q2"))
main.text <- c(patient.text, "Replication timing")
plotSurvfit(fit, file.name, main.text, c("Q1+Q4", "Q2+Q3"), c("red", "blue"))

##
phenos.surv$Groups <- "no"
phenos.surv$Groups[which(phenos.surv$tissue.sampling == "surgical resection")] <- "yes"
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery"))
main.text <- c(patient.text, "Surgery")
plotSurvfit(fit, file.name, main.text, c("No", "Yes"), c("red", "blue"))

##
phenos.surv <- phenos.surv[!is.na(phenos.surv$radiation..yes.no.),]
phenos.surv$Groups <- phenos.surv$radiation..yes.no.
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Male_I-II_Radio"))
main.text <- c(patient.text, "Radiotherapy")
plotSurvfit(fit, file.name, main.text, c("No", "Yes"), c("red", "blue"))

##
phenos.surv <- phenos.surv[!is.na(phenos.surv$chemotherapy..yes.no.),]
phenos.surv$Groups <- phenos.surv$chemotherapy..yes.no.
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Male_I-II_Chemo"))
main.text <- c(patient.text, "Chemotherapy")
plotSurvfit(fit, file.name, main.text, c("No", "Yes"), c("red", "blue"))

##
phenos.surv$Groups <- phenos.surv$RT
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Surgery_Female_Chemo_RT"))
main.text <- c(patient.text, "Replication timing")
plotSurvfit(fit, file.name, main.text, c("M1", "M2"), c("blue", "red"))

# -----------------------------------------------------------------------------
# Kaplan-Meier survival analysis (OS~kmeans): MUT genes
# Last Modified: 02/05/19
# -----------------------------------------------------------------------------
results.gene$Survival <- NA
results.gene$Survival_BH <- NA
results.gene$Survival_SC <- NA
results.gene$Survival_SC_BH <- NA
for (g in 1:nrow(results.gene)) {
   gene <- rownames(results.gene)[g]
   #gene <- "COL22A1"
   samples.mut <- unique(subset(txt, Gene_Hugo == gene)$PAT_ID)
   phenos.surv.mut <- phenos.surv
   
   phenos.surv.mut$MUT <- 0
   phenos.surv.mut[samples.mut,]$MUT <- 1
   #phenos <- subset(phenos, sex == "male")
   phenos.surv.mut <- subset(phenos.surv.mut, tissue.sampling == "surgical resection")
   phenos.surv.mut <- subset(phenos.surv.mut, chemotherapy..yes.no. == "yes")
   phenos.surv.mut$Groups <- phenos.surv.mut$MUT
   extension <- paste0("_", gene, "_S+C")
   
   ##
   fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv.mut)
   file.name <- file.path(wd.rt.plots, "survfit", "MUT", "S+C", paste0("survfit_", base, "_OS~Groups", extension))
   main.text <- c(paste0(BASE, " overall survival (OS)"), paste0(gene, " status (in sugery and chemo)"))
   cols    <- c("blue", "red")
   legends <- c(paste0("WT", " (n=", fit[[1]][1], ")"), paste0("MUT", " (n=", fit[[1]][2], ")"))

   if (!is.na(fit[[1]][2])) {
      pdf(paste0(file.name, ".pdf"), height=5, width=5)
      plot(fit, ylim=c(0, 1), xlab="Months", ylab="Survival probability", col=cols, main=main.text[1], mark.time=T)
      legend("topright", legend=legends, lwd=1, col=cols)
      mtext(main.text[2], cex=1.2, line=0.3)
      text(0, 0, get_surv_pvalue(fit), adj= c(-0.05, 0.1), col="black")
      dev.off()
   
      #results.gene.sub$Survival[g] <- surv_pvalue(fit, method="log-rank")$pval
   }
}
results.gene.sub$Survival_BH <- p.adjust(results.gene.sub$Survival, "BH")
writeTable(results.gene.sub, file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6_M2-M1_Survival_S+C.txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Kaplan-Meier survival analysis (OS~kmeans)
# Last Modified: 02/05/19
# -----------------------------------------------------------------------------
results.gene$Survival <- NA
#results.gene$Survival_SC <- 0
for (g in 1:nrow(results.gene.sub)) {
   gene <- rownames(results.gene.sub)[g]
 
   samples.mut <- unique(subset(txt, Gene_Hugo == gene)$PAT_ID)     ## ADD 06/05/19
   phenos.surv.mut <- phenos.surv[intersect(rownames(phenos.surv), samples.mut),]  ## BUG 07/05/19; ADD 06/05/19
   #phenos.surv.mut <- phenos.surv[setdiff(phenos.surv$Sample.ID, samples.mut),]   ## ADD 05/05/19
   #phenos.sub <- subset(phenos.sub, sex == "female")
   #phenos.surv.mut <- subset(phenos.surv.mut, tissue.sampling == "surgical resection")
   #phenos.surv.mut <- subset(phenos.surv.mut, chemotherapy..yes.no. == "yes")

   ##
   phenos.surv.mut$Groups <- phenos.surv.mut$Stage
   extension <- paste0("_Stage_", gene, "-S+C")
   fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv.mut)
   file.name <- file.path(wd.rt.plots, "survfit", "MUT", "MUT", paste0("survfit_", base, "_OS~Groups", extension))
   main.text <- c(paste0(BASE, " overall survival (OS)"), paste0(gene, " sugery and chemo patients (n=", nrow(phenos.surv.mut), ")"))
   if (!is.na(fit[[1]][1]) && !is.na(fit[[1]][2]))
      plotSurvfitStage(fit, file.name, main.text)
   
   ##
   phenos.surv.mut$Groups <- phenos.surv.mut$RT
   extension <- paste0("_RT_", gene, "-S+C")
   fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv.mut)
   file.name <- file.path(wd.rt.plots, "survfit", "MUT", "MUT", paste0("survfit_", base, "_OS~Groups", extension))
   main.text <- c(paste0(BASE, " overall survival (OS)"), paste0(gene, " sugery and chemo patients (n=", nrow(phenos.surv.mut), ")"))
   if (!is.na(fit[[1]][1]) && !is.na(fit[[1]][2]))
      plotSurvfitRT(fit, file.name, main.text, isRT)
   
   ##
   phenos.surv.mut$Groups <- phenos.surv.mut$sex
   extension <- paste0("_Sex_", gene, "-S+C")
   fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.surv.mut)
   file.name <- file.path(wd.rt.plots, "survfit", "MUT", "MUT", paste0("survfit_", base, "_OS~Groups", extension))
   main.text <- c(paste0(BASE, " overall survival (OS)"), paste0(gene, " sugery and chemo patients (n=", nrow(phenos.surv.mut), ")"))
   if (!is.na(fit[[1]][1]) && !is.na(fit[[1]][2]))
      plotSurvfitSex(fit, file.name, main.text)
}

# -----------------------------------------------------------------------------
# Kaplan-Meier survival analysis (OS~kmeans)
# Last Modified: 02/05/19
# -----------------------------------------------------------------------------
phenos <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
phenos <- cbind(phenos[samples$SAMPLE_ID,], samples)
#extension <- ""
#phenos <- subset(phenos, sex == "female")
#extension <- "_sex_male"
#phenos <- subset(phenos, tissue.sampling == "surgical resection")
#extension <- "_sex_male_surgery"
#phenos <- subset(phenos, chemotherapy..yes.no. == "yes")
#extension <- "_sex_male_surgery_chemotherapy_yes"
#phenos <- subset(phenos, Stage == 2)
#extension <- "_Stage3+4_Female"
isRT <- T

## OS_censor
phenos.sub <- phenos[!is.na(phenos$overall_survival..months.),]
phenos.sub$OS_month <- phenos.sub$overall_survival..months.
phenos.sub$Groups <- phenos.sub$Q4
if (isRT)
   phenos.sub$Groups <- phenos.sub$RT

phenos.sub$OS_censor <- phenos.sub$Status..at.time.of.last.follow.up.
phenos.sub$OS_censor <- gsub("dead", 0, phenos.sub$OS_censor)
phenos.sub$OS_censor <- gsub("alive", 1, phenos.sub$OS_censor)
phenos.sub$OS_censor <- as.numeric(phenos.sub$OS_censor)

##
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", extension))
main.text <- c(paste0(BASE, " overall survival (OS)"), paste0("PDE4DIP patients (n=", nrow(phenos.sub), ")"))
plotSurvfitRT(fit, file.name, main.text, isRT)

# -----------------------------------------------------------------------------
# Stage UICC
# Last Modified: 05/05/19
# -----------------------------------------------------------------------------
phenos <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
phenos <- cbind(phenos[samples$SAMPLE_ID,], samples)
phenos <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")

##
phenos$Stage <- NA
phenos[which(phenos$stage_UICC == "I"),]$Stage <- 1
phenos[which(phenos$stage_UICC == "Ia"),]$Stage <- 1
phenos[which(phenos$stage_UICC == "Ib"),]$Stage <- 1
phenos[which(phenos$stage_UICC == "IB"),]$Stage <- 1

phenos[which(phenos$stage_UICC == "II"),]$Stage <- 1
phenos[which(phenos$stage_UICC == "IIa"),]$Stage <- 1
phenos[which(phenos$stage_UICC == "IIb"),]$Stage <- 1

phenos[which(phenos$stage_UICC == "III"),]$Stage <- 2
phenos[which(phenos$stage_UICC == "IIIa"),]$Stage <- 2
phenos[which(phenos$stage_UICC == "IIIb"),]$Stage <- 2

phenos[which(phenos$stage_UICC == "IV"),]$Stage <- 2
phenos <- phenos[!is.na(phenos$Stage),]

## OS_censor
phenos.sub <- phenos[!is.na(phenos$overall_survival..months.),]
phenos.sub$OS_month <- phenos.sub$overall_survival..months.
phenos.sub$Groups <- phenos.sub$Stage

phenos.sub$OS_censor <- phenos.sub$Status..at.time.of.last.follow.up.
phenos.sub$OS_censor <- gsub("dead",  1, phenos.sub$OS_censor)
phenos.sub$OS_censor <- gsub("alive", 0, phenos.sub$OS_censor)
phenos.sub$OS_censor <- as.numeric(phenos.sub$OS_censor)

##
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Stage_UICC_Mark_Switch"))
main.text <- c(paste0(BASE, " overall survival (OS)"), paste0("All patients (n=", nrow(phenos.sub), ")"))
cols    <- c("blue", "red")
legends <- c(paste0("  I-II (n=", fit[[1]][1], ")"), paste0("III-IV (n=", fit[[1]][2], ")"))

pdf(paste0(file.name, ".pdf"), height=5, width=5)
plot(fit, ylim=c(0, 1), xlab="Months", ylab="Survival probability", col=cols, main=main.text[1], mark.time=T)
legend("topright", title="Stage UICC", legend=legends, lwd=1, col=cols, title.adj=0.19)
mtext(main.text[2], cex=1.2, line=0.3)
#legend("bottomleft", legend=get_surv_pvalue(fit), bty="n", text.col="black")
text(0, 0, get_surv_pvalue(fit), adj= c(-0.05, 0.1), col="black")
dev.off()

# -----------------------------------------------------------------------------
# Stage T
# Last Modified: 05/05/19
# -----------------------------------------------------------------------------
phenos <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
phenos <- cbind(phenos[samples$SAMPLE_ID,], samples)

phenos$Stage <- NA
phenos[which(phenos$stage_T == "1"),]$Stage <- 1
phenos[which(phenos$stage_T == "1a"),]$Stage <- 1
phenos[which(phenos$stage_T == "1b"),]$Stage <- 1

phenos[which(phenos$stage_T == "2"),]$Stage <- 1
phenos[which(phenos$stage_T == "2a"),]$Stage <- 1
phenos[which(phenos$stage_T == "2b"),]$Stage <- 1

phenos[which(phenos$stage_T == "3"),]$Stage <- 2

phenos[which(phenos$stage_T == "4"),]$Stage <- 2
phenos <- phenos[!is.na(phenos$Stage),]

## OS_censor
phenos.sub <- phenos[!is.na(phenos$overall_survival..months.),]
phenos.sub$OS_month <- phenos.sub$overall_survival..months.
phenos.sub$Groups <- phenos.sub$Stage

phenos.sub$OS_censor <- phenos.sub$Status..at.time.of.last.follow.up.
phenos.sub$OS_censor <- gsub("dead", 0, phenos.sub$OS_censor)
phenos.sub$OS_censor <- gsub("alive", 1, phenos.sub$OS_censor)
phenos.sub$OS_censor <- as.numeric(phenos.sub$OS_censor)

##
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub)
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", "_Stage_T"))
main.text <- c(paste0(BASE, " overall survival (OS)"), paste0("All patients (n=", nrow(phenos.sub), ")"))
cols    <- c("blue", "red")
legends <- c(paste0("I-II (n=", fit[[1]][1], ")"), paste0("III-IV (n=", fit[[1]][2], ")"))

pdf(paste0(file.name, ".pdf"), height=6, width=6)
plot(fit, ylim=c(0, 1), xlab="Months", ylab="Survival probability", col=cols, main=main.text[1])
legend("topright", title="Stage T", legend=legends, lwd=1, col=cols, title.adj=0.13)
mtext(main.text[2], cex=1.2, line=0.3)
#legend("bottomleft", legend=get_surv_pvalue(fit), bty="n", text.col="black")
text(0, 0, get_surv_pvalue(fit), adj= c(-0.05, 0.1), col="black")
dev.off()










##
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("III-IV", "I-II")
 
test[1, 2] <- nrow(subset(subset(phenos, RT == 0), Stage2 %in% c(3,4)))
test[1, 1] <- nrow(subset(subset(phenos, RT == 1), Stage2 %in% c(3,4)))
test[2, 2] <- nrow(subset(subset(phenos, RT == 0), Stage2 %in% c(1,2)))
test[2, 1] <- nrow(subset(subset(phenos, RT == 1), Stage2 %in% c(1,2)))
 
fisher.test(test)[[1]]
writeTable(test, file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6_Q4-Q1.txt"), colnames=T, rownames=T, sep="\t")











pdf()
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub) 
plot(fit, ylim=c(0, 1), xlab="Months", ylab="Survival probability", col=cols, main=)
legend("topright", legend=c("Q1", "Q2", "Q3", "Q4"), lwd=1, col=cols)
mtext(paste0("log-rank P = ", surv_pvalue(fit, method="log-rank")$pval), cex=1.2, line=0.3)
dev.off()

## OS~kmeans
cols <- c("blue", "red")
pdf(file.path(wd.rt.plots, "survfit_sclc_OS~Groups_sex_male_RT.pdf"), height=6, width=6)
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub) 
plot(fit, ylim=c(0, 1), xlab="Months", ylab="Survival probability", col=cols, main="SCLC overall survival (n=57; male)")
legend("topright", legend=c("M1", "M2"), lwd=1, col=cols)
mtext(paste0("log-rank p-value = ", surv_pvalue(fit)$pval), cex=1.2, line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# Kaplan-Meier survival analysis
# Last Modified: 02/05/19
# -----------------------------------------------------------------------------
phenos <- readTable(file.path(wd.meta, "nature14664-s1_ST2.txt"), header=T, rownames=T, sep="\t")
phenos <- phenos[samples$SAMPLE_ID,]
phenos$Groups <- samples[rownames(phenos),]$Q4

phenos <- phenos[phenos$CCF_per_read >= 0.086,]

median(subset(phenos, Groups == 1)$Het_Score_Muts)
median(subset(phenos, Groups == 4)$Het_Score_Muts)
testU(subset(phenos, Groups == 1)$Het_Score_Muts, subset(phenos, Groups == 4)$Het_Score_Muts)

median(subset(phenos, Groups == 1)$CCF_per_read)
median(subset(phenos, Groups == 4)$CCF_per_read)
testU(subset(phenos, Groups == 1)$CCF_per_read, subset(phenos, Groups == 4)$CCF_per_read)







## OS_censor
phenos.sub <- phenos[!is.na(phenos$overall_survival..months.),]
phenos.sub$OS_month <- phenos.sub$overall_survival..months.
#phenos.sub$Groups <- samples[rownames(phenos.sub),]$Q4






###
## PFS_censor
phenos.sub <- phenos[!is.na(phenos$progression.free_survival..months.),]
#phenos.sub$Groups <- samples[rownames(phenos.sub),]$Q4
phenos.sub$Groups <- samples[rownames(phenos.sub),]$RT

phenos.sub$OS_censor <- phenos.sub$Status..at.time.of.last.follow.up.
phenos.sub$OS_censor <- gsub("dead", 0, phenos.sub$OS_censor)
phenos.sub$OS_censor <- gsub("alive", 1, phenos.sub$OS_censor)
phenos.sub$OS_censor <- as.numeric(phenos.sub$OS_censor)

phenos.sub$OS_month <- phenos.sub$progression.free_survival..months.





## PFS~kmeans
cols <- c("blue", "skyblue3", "lightcoral", "red")
pdf(file.path(wd.rt.plots, "survfit_sclc_PFS~Groups_Q4.pdf"))
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub) 
plot(fit, ylim=c(0, 1), xlab="Months", ylab="PFS", col=cols, main="SCLC progression free survival (n=37)")
legend("topright", legend=c("Q1", "Q2", "Q3", "Q4"), lwd=1, col=cols)
mtext(surv_pvalue(fit)$pval.txt, cex=1.2, line=0.3)
dev.off()

##
cols <- c("blue", "red")
pdf(file.path(wd.rt.plots, "survfit_sclc_PFS~Groups_RT.pdf"))
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub) 
plot(fit, ylim=c(0, 1), xlab="Months", ylab="PFS", col=cols, main="SCLC progression free survival (n=37)")
legend("topright", legend=c("M1", "M2"), lwd=1, col=cols)
mtext(surv_pvalue(fit)$pval.txt, cex=1.2, line=0.3)
dev.off()
