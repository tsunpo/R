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
BASE <- "NBL"
base <- tolower(BASE)

wd.wgs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "Peifer 2015")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples <- readTable(file.path(wd.wgs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")
phenos  <- readTable(file.path(wd.meta, "NB_WGS_Mutations_Overview.txt"), header=T, rownames=T, sep="")

# -----------------------------------------------------------------------------
# 
# Last Modified: 16/08/18
# -----------------------------------------------------------------------------
## Risk group
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("High", "Low")

test[1, 1] <- nrow(subset(phenos[c(q4,q3),], risk == "high"))
test[1, 2] <- nrow(subset(phenos[c(q2,q1),], risk == "high"))
test[2, 1] <- nrow(subset(phenos[c(q4,q3),], risk == "low"))
test[2, 2] <- nrow(subset(phenos[c(q2,q1),], risk == "low"))
fisher.test(test)[[1]]
# [1] 4.382596e-07

## MYCN amp
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("MYCN_amp", "MYCN")

test[1, 1] <- nrow(subset(phenos[c(q4,q3),], MYCN_amp == 1))
test[1, 2] <- nrow(subset(phenos[c(q2,q1),], MYCN_amp == 1))
test[2, 1] <- nrow(subset(phenos[c(q4,q3),], MYCN_amp == 0))
test[2, 2] <- nrow(subset(phenos[c(q2,q1),], MYCN_amp == 0))
# > test
# M2 M1
# MYCN_amp  9  1
# MYCN     19 27
# > fisher.test(test)[[1]]
# [1] 0.01159974

## MYCN amp
test <- toTable(0, 2, 2, c("Q4", "Q3"))
rownames(test) <- c("MYCN_amp", "MYCN")

test[1, 1] <- nrow(subset(phenos[q4,], MYCN_amp == 1))
test[1, 2] <- nrow(subset(phenos[q3,], MYCN_amp == 1))
test[2, 1] <- nrow(subset(phenos[q4,], MYCN_amp == 0))
test[2, 2] <- nrow(subset(phenos[q3,], MYCN_amp == 0))
# > test
# Q4 Q3
# MYCN_amp  1  8
# MYCN     13  6
# > fisher.test(test)[[1]]
# [1] 0.01275362

## TERT_Rearrangement
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("TERT_re", "TERT")

test[1, 1] <- nrow(subset(phenos[c(q4,q3),], TERT_Rearrangement == 1))
test[1, 2] <- nrow(subset(phenos[c(q2,q1),], TERT_Rearrangement == 1))
test[2, 1] <- nrow(subset(phenos[c(q4,q3),], TERT_Rearrangement == 0))
test[2, 2] <- nrow(subset(phenos[c(q2,q1),], TERT_Rearrangement == 0))
# > test
# M2 M1
# TERT_re  9  4
# TERT    19 24
# > fisher.test(test)[[1]]
# [1] 0.2046831




test[1, 2] <- nrow(subset(phenos, risk == "low"))

test[2, 2] <- nrow(subset(subset(phenos, RT == 0), Stage2 %in% c(1,2)))
test[2, 1] <- nrow(subset(subset(phenos, RT == 1), Stage2 %in% c(1,2)))

fisher.test(test)[[1]]





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
