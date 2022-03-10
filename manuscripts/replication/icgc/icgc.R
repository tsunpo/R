# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/icgc-wgs-rt/icgc-wgs.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 19/02/22
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "ReplicationTiming.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.icgc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 18/02/22
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "ICGC"
base  <- tolower(BASE)

wd.ngs <- file.path(wd, BASE, "ngs/WGS")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.rt.plots <- "/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/analysis/replication/icgc-wgs-rt/data/"

# -----------------------------------------------------------------------------
# 
# Last Modified: 18/02/22
# -----------------------------------------------------------------------------
raw <- readTable(file.path(wd.meta, "PCAWG.raw.data.copy_number.converted_data.list"), header=F, rownames=F, sep="\t")
# > length(raw)
# [1] 2951
seg <- readTable(file.path(wd.meta, "PCAWG.raw.data.copy_number.sclust_final_copy_number_analysis_files.list"), header=F, rownames=F, sep="\t")
# > length(seg)
# [1] 2778
overlaps <- intersect(seg, gsub("_CONVERTED", "", raw))
# > length(overlaps)
# [1] 2777

table <- readTable(file.path(wd.meta, "PCAWG.wgs_id_histology_appreviations.txt"), header=T, rownames=T, sep="\t")
dim(table)
overlaps.2 <- intersect(rownames(table), overlaps)
table <- table[overlaps.2,]
table$project_code <- ""

mapping <- readTable(file.path(wd.meta, "PCAWG.wgs_id_histology_table_mached.txt"), header=T, rownames=F, sep="\t")
# > dim(mapping)
# [1] 3911   29

idx <- which(table$wgs_id %in% mapping$pcawg_wgs_id)
table <- table[idx,]
rownames(table) <- table$specimen_id

for (s in 1:nrow(table)) {
   codes <- mapping[which(mapping$pcawg_wgs_id == table$wgs_id[s]),]$project_code
   if (length(codes) > 1) {
      if (codes[1] != codes[2])
         print(codes)
      table$project_code[s] <- codes[1]
   } else {
      table$project_code[s] <- codes
   }
}
sort(table(table$project_code), decreasing=T)
sort(table(table$histology_abbreviation), decreasing=T)

###
##
clinicals <- readTable(file.path("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/metadata/DCC_DATA_RELEASE", "pcawg_donor_clinical_August2016_v9.txt"), header=T, rownames=F, sep="\t")
# > dim(clinicals)
# [1] 1462   18

idx2 <- which(clinicals$icgc_donor_id %in% mapping$icgc_donor_id)
clinicals <- clinicals[idx2,]
clinicals$icgc_specimen_id <- NA
for (s in 1:nrow(clinicals)) {
   ids <- mapping[which(mapping$icgc_donor_id == clinicals$icgc_donor_id[s]),]$icgc_specimen_id
   if (length(ids) > 1) {
      clinicals$icgc_specimen_id[s] <- paste0(paste(unique(ids), collapse=","), ",")
   } else {
      clinicals$icgc_specimen_id[s] <- paste0(ids, ",")
   }
}
rownames(clinicals) <- clinicals$icgc_specimen_id
writeTable(clinicals, "/Users/tpyang/Work/uni-koeln/tyang2/ICGC/metadata/DCC_DATA_RELEASE/pcawg_donor_clinical_August2016_v9_tyang2.txt", colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
projects <- unique(table$histology_abbreviation)
removes  <- c("SoftTissue-Liposarc", "Cervix-SCC", "CNS-Oligo", "Bone-Benign", "Myeloid-AML", "SoftTissue-Leiomyo", "Breast-LobularCA", "Bone-Epith", "Breast-DCIS", "Myeloid-MDS", "Cervix-AdenoCA")
diffs    <- setdiff(projects, removes)
table.histology <- subset(table, histology_abbreviation %in% diffs)

histologies <- unique(table.histology$histology_abbreviation)
icgc <- toTable(0, 2, length(histologies), c("MEDIAN", "N"))
rownames(icgc) <- histologies
for (h in 1:length(histologies)) {
   samples   <- subset(table.histology, histology_abbreviation == histologies[h])
   icgc$N[h] <- nrow(samples)
   
   sum <- c()
   for (s in 1:nrow(samples)) {
      sample <- rownames(samples)[s]
      load(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/analysis/replication/icgc-wgs-rt/data/samples/rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData"))
      
      sum <- c(sum, as.numeric(cor))
   }
   icgc$MEDIAN[h] <- median(sum)
}
icgc <- icgc[order(icgc$MEDIAN, decreasing=T),]
writeTable(icgc, "/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/analysis/replication/icgc-wgs-rt/data/icgc_histology_abbreviation>20.txt", colnames=T, rownames=T, sep="\t")

###
##
wd.rt.data <- "/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/analysis/replication/icgc-wgs-rt/data/"
donors <- c("icgc_donor_id", "donor_sex", "donor_vital_status", "donor_age_at_diagnosis", "donor_survival_time")
samples <- toTable(0, 11, sum(icgc$N), c("CANCER", "COR", "Q4", "M2", "icgc_specimen_id", "histology_abbreviation", donors))
idx.cor    <- 0
idx.sample <- 0
for (h in nrow(icgc):1) {
   samples.list <- subset(table.histology, histology_abbreviation == rownames(icgc)[h])
   samples.q4   <- setSamplesQ4(wd.rt.data, rownames(samples.list))
   
   for (s in 1:nrow(samples.list)) {
      sample <- rownames(samples.list)[s]
      load(paste0(wd.rt.data, "/samples/rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData"))
  
      samples$COR[idx.sample + s] <- as.numeric(cor)
      samples$icgc_specimen_id[idx.sample + s] <- sample
      samples$histology_abbreviation[idx.sample + s] <- rownames(icgc)[h]
      
      samples[idx.sample + s, c(donors)] <- clinicals[grep(paste0(sample, ","), clinicals$icgc_specimen_id), c(donors)]
   }
   samples$CANCER[(idx.sample+1):(idx.sample+nrow(samples.list))] <- idx.cor
   samples$Q4[(idx.sample+1):(idx.sample+nrow(samples.list))] <- samples.q4$Q4
   samples$M2[(idx.sample+1):(idx.sample+nrow(samples.list))] <- samples.q4$M2
   
   idx.cor <- idx.cor + 1
   idx.sample <- idx.sample + nrow(samples.list)
}

###
##
file.name <- "stripchart_ICGC_H"
main.text <- expression(bold(~bolditalic('in silico')~"sorting of PCAWG samples"))
cols <- c(blue, blue.lighter, red.lighter, red)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=18, width=9)
par(mar = c(5, 14.5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h, horizontal=T, yaxt="n", xaxt="n", ylab="", main=main.text, col="white", outline=F, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
#text(labels=labels, x=1:26, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, side=1, cex=1.8)

stripchart(COR ~ CANCER, data=subset(samples.h, Q4 == 1), method="jitter", cex=1.5, pch=19, col=cols[1], vertical=F, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples.h, Q4 == 2), method="jitter", cex=1.5, pch=19, col=cols[2], vertical=F, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples.h, Q4 == 3), method="jitter", cex=1.5, pch=19, col=cols[3], vertical=F, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples.h, Q4 == 4), method="jitter", cex=1.5, pch=19, col=cols[4], vertical=F, add=T, at=c(1:26))

axis(side=1, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.7)
mtext("Spearman's rho", side=1, line=3.5, cex=1.8)
#mtext("", cex=1.2, line=0.3)
mtext(rev(rownames(icgc.tmp)), side=2, line=0.5, cex=1.8, at=1:26, las=1)

legend("bottomright", legend=c("Q1", "Q2", "Q3", "Q4"), pch=19, pt.cex=2.5, col=cols, cex=1.8, horiz=T)
dev.off()

###
##
file.name <- "stripchart_ICGC_V"
main.text <- expression(bold(~bolditalic('in silico')~"sorting of 2,612 ICGC PCAWG samples"))
#labels <- paste0(labels=rownames(icgc), " (n=", icgc$N, ")")
labels <- rownames(icgc)
cols <- c(blue, blue.lighter, red.lighter, red)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
par(mar = c(11, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.v, yaxt="n", xaxt="n", ylab="", main=main.text, col="white", outline=F, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
text(labels=labels, x=1:26, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, side=1, cex=1.8)

stripchart(COR ~ CANCER, data=subset(samples.v, Q4 == 1), method="jitter", cex=1.5, pch=19, col=cols[1], vertical=T, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples.v, Q4 == 2), method="jitter", cex=1.5, pch=19, col=cols[2], vertical=T, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples.v, Q4 == 3), method="jitter", cex=1.5, pch=19, col=cols[3], vertical=T, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples.v, Q4 == 4), method="jitter", cex=1.5, pch=19, col=cols[4], vertical=T, add=T, at=c(1:26))

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.7)
mtext("Spearman's rho", side=2, line=3.5, cex=1.8)
#mtext("", cex=1.2, line=0.3)
#mtext(text=rownames(icgc), side=1, cex=1.9, line=1.3, at=par("usr")[3], srt=35, adj=0, xpd=T)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.8)
dev.off()
save(table, table.histology, icgc, samples.h, samples.v, file=file.path(wd.rt.plots, paste0("stripchart_ICGC.RData")))

###
## 
colname <- "donor_age_at_diagnosis"
icgc$age_N <- NA
icgc$age_rho <- NA
icgc$age_P   <- NA
for (h in nrow(icgc):1) {
   hist <- rownames(icgc)[h]
   samples.hist <- subset(samples, histology_abbreviation == hist)
   samples.hist <- removeNA(samples.hist, colname)
   samples.hist <- subset(samples.hist, donor_age_at_diagnosis != 0)
   
   if (nrow(samples.hist) > 1) {
      x <- samples.hist[, colname]
      y <- samples.hist$COR
      file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/correlation_in-silico_vs_", colname, "_", hist)
      plotCorrelation(file.name, hist, "Age", expression(bolditalic('in silico')~"sorting"), x, y, "topright", line=2.4)
      
      cor <- cor.test(y, x, method="spearman", exact=F)
      icgc$age_rho[h] <- round0(cor[[4]], digits=2)
      icgc$age_P[h] <- scientific(cor[[3]], digits=2)
      icgc$age_N[h] <- nrow(samples.hist)
   }
}

###
##
colname <- "donor_survival_time"
icgc$survival_N <- NA
icgc$survival_rho <- NA
icgc$survival_P   <- NA
for (h in 17:nrow(icgc)) {
   hist <- rownames(icgc)[h]
   samples.hist <- subset(samples, histology_abbreviation == hist)
   samples.hist <- removeNA(samples.hist, colname)
   samples.hist <- subset(samples.hist, donor_survival_time != 0)

   if (nrow(samples.hist) > 10) {
      x <- samples.hist[, colname]
      y <- samples.hist$COR
      file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/correlation_in-silico_vs_", colname, "_", hist)
      plotCorrelation(file.name, hist, "Survival time", expression(bolditalic('in silico')~"sorting"), x, y, "topright", line=2.4)
  
      cor <- cor.test(y, x, method="spearman", exact=F)
      icgc$survival_rho[h] <- round0(cor[[4]], digits=2)
      icgc$survival_P[h] <- scientific(cor[[3]], digits=2)
      icgc$survival_N[h] <- nrow(samples.hist)
   }
}
save(table, table.histology, icgc, samples, samples.h, samples.v, clinicals, file=file.path(wd.rt.plots, paste0("stripchart_ICGC_surv.RData")))

###
##
colname <- "donor_survival_time"
library(survival)
library(survminer)
for (h in 1:nrow(icgc)) {
   hist <- rownames(icgc)[h]
   samples.hist <- subset(samples, histology_abbreviation == hist)
   samples.hist.surv <- survICGC(samples.hist, colname)
   
   if (nrow(samples.hist.surv) > 10) {
      ## Cox regression model
   
      print(hist)
      #res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR + donor_sex + donor_age_at_diagnosis, data=samples.hist.surv)
      res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR, data=samples.hist.surv)
      print(res.cox)
      print("----------------------------------------------------------------")
      #pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/hazard_", hist, ".pdf"), height=3, width=5)
      #ggforest(res.cox, data=samples.hist.surv, main=paste0("Hazard ratio in ", hist), cpositions = c(0.02, 0.22, 0.4), fontsize=0.7)
      #dev.off()
      #ggforest(res.cox)
      
      ##
      fit <- survfit(Surv(OS_month, OS_censor) ~ M2, data=samples.hist.surv)
      file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/M2/survfit_", hist, "_M2"))
      main.text <- c(hist, expression(bolditalic('in silico')~"sorting"))
      plotSurvfit(fit, file.name, main.text, c("M1", "M2"), c(blue, red))
      
      ##
      fit <- survfit(Surv(OS_month, OS_censor) ~ Q4, data=samples.hist.surv)
      file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/Q4/survfit_", hist, "_Q4"))
      main.text <- c(hist, expression(bolditalic('in silico')~"sorting"))
      plotSurvfit(fit, file.name, main.text, c("Q1", "Q2", "Q3", "Q4"), c(blue, blue.lighter, red.lighter, red))
   }
}

###
## Wow
samples.surv <- survICGC(samples)
res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR + donor_sex + donor_age_at_diagnosis, data=samples.surv)
print(res.cox)

#samples.surv <- samples.surv.neg
colname <- "donor_survival_time"
x <- samples.surv[, colname]
y <- samples.surv$COR
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/correlation_in-silico_vs_", colname, "_", "ALL")
plotCorrelation(file.name, "1,683 PCAWG samples", "Survival time", expression(italic('in silico')~"sorting"), x, y, "topright", line=2.4)

colname <- "donor_age_at_diagnosis"
x <- samples.surv[, colname]
y <- samples.surv$COR
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/correlation_in-silico_vs_", colname, "_", "ALL")
plotCorrelation(file.name, "1,683 PCAWG samples", "Age at diagnosis", expression(italic('in silico')~"sorting"), x, y, "topright", line=2.4)

fit <- survfit(Surv(OS_month, OS_censor) ~ donor_sex, data=samples.surv)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/survfit_Sex_ALL"))
main.text <- c("Sex", "")
plotSurvfit(fit, file.name, main.text, c("F", "M"), c(red.lighter, blue.lighter))

## Wow
file.name <- paste0("boxplot_PCAWG_F-vs-M_COR")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/", file.name, subset(samples.surv, donor_sex == "male")$COR, subset(samples.surv, donor_sex == "female")$COR, expression(italic('in silico')~"sorting"), names=c("Male", "Female"), cols=c(blue.lighter, red.lighter))

file.name <- paste0("boxplot_PCAWG_F-vs-M_Age")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/", file.name, subset(samples.surv, donor_sex == "male")$donor_age_at_diagnosis, subset(samples.surv, donor_sex == "female")$donor_age_at_diagnosis, "All", names=c("Male", "Female"), cols=c(blue.lighter, red.lighter))

file.name <- paste0("boxplot_PCAWG_G1-vs-S_Age")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/", file.name, subset(samples.surv, SG1 == "G1")$donor_age_at_diagnosis, subset(samples.surv, SG1 == "S")$donor_age_at_diagnosis, "All", names=c("G1", "S"), cols=c(blue, red))


###
## ???
res.cox <- coxph(Surv(OS_month, OS_censor) ~ SG1 + COR + donor_sex + donor_age_at_diagnosis, data=samples.surv)

fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/survfit_", "ALL", "_S-vs-G1"))
main.text <- c("Cell cycle statue", "")
plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

icgc.pos <- subset(icgc, survival_rho > 0)
icgc.neg <- subset(icgc, survival_rho < 0)
samples.surv.pos <- subset(samples.surv, histology_abbreviation %in% rownames(icgc.pos))
samples.surv.neg <- subset(samples.surv, histology_abbreviation %in% rownames(icgc.neg))

res.cox <- coxph(Surv(OS_month, OS_censor) ~ SG1 + COR + donor_sex + donor_age_at_diagnosis, data=samples.surv.pos)
fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv.pos)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/survfit_", "ALL", "_S-vs-G1_RHO>0"))
main.text <- c("Cell cycle statue", "RHO > 0")
plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

res.cox <- coxph(Surv(OS_month, OS_censor) ~ SG1 + COR + donor_sex + donor_age_at_diagnosis, data=samples.surv.neg)
fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv.neg)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/survfit_", "ALL", "_S-vs-G1_RHO<0"))
main.text <- c("Cell cycle statue", "RHO < 0")
plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

## ???
test <- toTable(0, 2, 2, c("G1", "S"))
rownames(test) <- c("M", "F")
test[1, 1] <- nrow(subset(subset(samples.surv.pos, donor_sex == "male"), SG1 == "G1"))
test[1, 2] <- nrow(subset(subset(samples.surv.pos, donor_sex == "male"), SG1 == "S"))
test[2, 1] <- nrow(subset(subset(samples.surv.pos, donor_sex == "female"), SG1 == "G1"))
test[2, 2] <- nrow(subset(subset(samples.surv.pos, donor_sex == "female"), SG1 == "S"))
fisher.test(test)[[1]]

test <- toTable(0, 2, 2, c("G1", "S"))
rownames(test) <- c("M", "F")
test[1, 1] <- nrow(subset(subset(samples.surv.neg, donor_sex == "male"), SG1 == "G1"))
test[1, 2] <- nrow(subset(subset(samples.surv.neg, donor_sex == "male"), SG1 == "S"))
test[2, 1] <- nrow(subset(subset(samples.surv.neg, donor_sex == "female"), SG1 == "G1"))
test[2, 2] <- nrow(subset(subset(samples.surv.neg, donor_sex == "female"), SG1 == "S"))
fisher.test(test)[[1]]

# > hist
# [1] "Kidney-RCC"
test <- toTable(0, 2, 2, c("G1", "S"))
rownames(test) <- c("M", "F")
test[1, 1] <- nrow(subset(subset(samples.hist, donor_sex == "male"), SG1 == "G1"))
test[1, 2] <- nrow(subset(subset(samples.hist, donor_sex == "male"), SG1 == "S"))
test[2, 1] <- nrow(subset(subset(samples.hist, donor_sex == "female"), SG1 == "G1"))
test[2, 2] <- nrow(subset(subset(samples.hist, donor_sex == "female"), SG1 == "S"))
fisher.test(test)[[1]]




## donor_sex
## Surgery
test <- toTable(0, 2, 2, c("M2", "M1"))
rownames(test) <- c("Yes", "No")

test[1, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Surgery == "yes"))
test[1, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Surgery == "yes"))
test[2, 1] <- nrow(subset(subset(phenos.surv, M2 == 1), Surgery == "no"))
test[2, 2] <- nrow(subset(subset(phenos.surv, M2 == 0), Surgery == "no"))
fisher.test(test)[[1]]
## donor_vital_status








###
##
file.name <- "beeswarm_ICGC_H"
main.text <- expression(bold(~bolditalic('in silico')~"sorting of 2,612 ICGC PCAWG samples"))
#labels <- paste0(labels=rownames(icgc), " (n=", icgc$N, ")")
labels <- rownames(icgc)
cols <- c(blue, blue.lighter, red.lighter, red)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
par(mar = c(11, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h, yaxt="n", xaxt="n", ylab="", main=main.text, col="white", outline=F, cex.axis=1.7, cex.lab=1.8, cex.main=1.9)
text(labels=rev(labels), x=1:26, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, side=1, cex=1.8)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col=blue, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col=blue.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col=red.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col=red, pch=19, cex=1, add=T)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.7)
mtext("Spearman's rho", side=2, line=3.5, cex=1.8)
#mtext("", cex=1.2, line=0.3)
#mtext(text=rownames(icgc), side=1, cex=1.9, line=1.3, at=par("usr")[3], srt=35, adj=0, xpd=T)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.8)
dev.off()






#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_ICGC_Spearman's.pdf"), height=9, width=45)
ymax <- 0.8   #max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=NA, ylim=c(ymin, ymax), ylab="", main=expression(bold(~bolditalic('in silico')~"sorting of ICGC PCAWG samples")), yaxt="n", xaxt="n", cex.axis=1.8, cex.lab=1.9, cex.main=2.1)
abline(h=0, lty=5, lwd=2)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.9)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col=blue, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col=blue.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col=red.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col=red, pch=19, cex=1, add=T)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Spearman's rho", side=2, line=2.7, cex=1.9)
#mtext("", cex=1.2, line=0.3)
mtext(text=rownames(icgc), side=1, cex=1.9, line=1.3, at=c(1:26), srt=35, adj=1.1, xpd=T)
dev.off()












table.cll <- subset(table, project_code == "CLLE-ES")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "clle-es_wgs_n95.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "clle-es_wgs_n95.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "BRCA-EU")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "brca-eu_wgs_n78.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "brca-eu_wgs_n78.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "BRCA-US")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "brca-us_wgs_n89.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "brca-us_wgs_n89.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "BRCA-UK")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "brca-uk_wgs_n44.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "brca-uk_wgs_n44.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "ESAD-UK")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "esad-uk_wgs_n88.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "esad-uk_wgs_n88.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "LUAD-US")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "luad-us_wgs_n37.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "luad-us_wgs_n37.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "LUSC-US")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "lusc-us_wgs_n48.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "lusc-us_wgs_n48.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "HNSC-US")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "hnsc-us_wgs_n44.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "hnsc-us_wgs_n44.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "ORCA-IN")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "orca-in_wgs_n13.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "orca-in_wgs_n13.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "OV-AU")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "ov-au_wgs_n70.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "ov-au_wgs_n70.txt"),  colnames=T, rownames=F, sep="\t")

table.cll <- subset(table, project_code == "OV-US")
table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "ov-us_wgs_n42.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "ov-us_wgs_n42.txt"),  colnames=T, rownames=F, sep="\t")

#table.cll <- subset(table, project_code == "UCEC-US")
#table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
#writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "ucec-us_wgs_n49.list"), colnames=F, rownames=F, sep="\t")
#writeTable(table.cll,            file.path(wd.ngs, "ucec-us_wgs_n49.txt"),  colnames=T, rownames=F, sep="\t")

#table.cll <- subset(table, project_code == "THCA-US")
#table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
#writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "thca-us_wgs_n48.list"), colnames=F, rownames=F, sep="\t")
#writeTable(table.cll,            file.path(wd.ngs, "thca-us_wgs_n48.txt"),  colnames=T, rownames=F, sep="\t")

###
##
projects <- unique(table$project_code)
removes <- c("CLLE-ES", "BRCA-EU", "BRCA-US", "BRCA-UK", "ESAD-UK", "LUAD-US", "LUSC-US", "HNSC-US", "ORCA-IN", "OV-AU", "OV-US")
diffs <- setdiff(projects, removes)
table.diffs <- subset(table, project_code %in% diffs)
sort(table(table.diffs$histology_abbreviation), decreasing=T)

table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "icgc_wgs_n2096.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "icgc_wgs_n2096.txt"),  colnames=T, rownames=F, sep="\t")

writeTable(table.cll[1:100,     c(5,1,2)], file.path(wd.ngs, "icgc_wgs_n1-100.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll[101:1000,  c(5,1,2)], file.path(wd.ngs, "icgc_wgs_n101-1000.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll[1001:2096, c(5,1,2)], file.path(wd.ngs, "icgc_wgs_n1001-2096.list"), colnames=F, rownames=F, sep="\t")

##
table.removes <- subset(table, project_code %in% removes)
sort(table(table.removes$histology_abbreviation), decreasing=T)
table.cll <- table.removes

table.cll <- table.cll[sort(table.cll$specimen_id, decreasing=F),]
writeTable(table.cll[,c(5,1,2)], file.path(wd.ngs, "icgc_wgs_n648.list"), colnames=F, rownames=F, sep="\t")
writeTable(table.cll,            file.path(wd.ngs, "icgc_wgs_n648.txt"),  colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------




samples.brca.us <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/BRCA-US/ngs/WGS/brca-us_wgs_n89.txt", header=T, rownames=T, sep="")
samples.esad.uk <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ESAD-UK/ngs/WGS/esad-uk_wgs_n88.txt", header=T, rownames=T, sep="")
samples.clle.es <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/CLLE-ES/ngs/WGS/clle-es_wgs_n95.txt", header=T, rownames=T, sep="")

samples.sclc <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/BRCA-US/ngs/WGS/brca-us_wgs_m2_n89.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ESAD-UK/ngs/WGS/esad-uk_wgs_m2_n88.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/CLLE-ES/ngs/WGS/clle-es_wgs_m2_n95.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_brca+esad+clle_Spearman's.pdf"), height=6, width=6)
ymax <- 0.8   #max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("", "", ""), ylim=c(ymin, ymax), ylab="", main="Tumour read depth vs. RT", yaxt="n", cex.axis=1.8, cex.lab=1.9, cex.main=2.1)
abline(h=0, lty=5, lwd=2)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.9)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col=blue, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col=blue.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col=red.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col=red, pch=19, cex=1, add=T)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Spearman's rho", side=2, line=2.7, cex=1.9)
#mtext("", cex=1.2, line=0.3)
mtext(text=c("BRCA", "ESAD", "CLLE"), side=1, cex=1.9, line=1.3, at=c(1,2,3))
mtext(text=c("n=89", "n=88", "n=95"), side=1, cex=1.8, line=3, at=c(1,2,3))
dev.off()


# -----------------------------------------------------------------------------
# CLL
# Last Modified: 20/02/22
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"   ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

wd.ngs <- file.path(wd, "ICGC", "CLLE-ES", "ngs/WGS")
wd.ngs.data <- file.path(wd.ngs, "data")

##
clle <- readTable(file.path(wd.ngs, "clle-es_wgs_n95.txt"),  header=T, rownames=2, sep="")
cll  <- readTable(file.path(wd.ngs, "cll_wgs_m2q4_n96.txt"), header=T, rownames=1, sep="")
overlaps <- intersect(rownames(cll), rownames(clle))

cll <- cll[overlaps,]
cll$M2 <- cll$M2 + 1
writeTable(cll, file.path(wd.ngs, "cll_wgs_m2_n91.txt"), colnames=T, rownames=F, sep="\t")

##
clle <- readTable(file.path(wd.ngs, "clle-es_wgs_n95.txt"),  header=T, rownames=2, sep="")
cll  <- readTable(file.path(wd.ngs, "cll_wgs_m2q4_n47.txt"), header=T, rownames=1, sep="")
overlaps <- intersect(rownames(cll), rownames(clle))

cll <- cll[overlaps,]
cll$M2 <- cll$M2 + 1
writeTable(cll, file.path(wd.ngs, "cll_wgs_q4_n47.txt"), colnames=T, rownames=F, sep="\t")

# f935c981-4aef-d6a4-e040-11ac0c481062	SP116919	Myeloid-AML	Included
# f940cea2-7e79-e422-e040-11ac0d483224	SP116919	Myeloid-AML	Included

# -----------------------------------------------------------------------------
# Check ICGC partitioning
# Last Modified: 24/02/22
# -----------------------------------------------------------------------------
for (c in 1:22) {
   chr <- chrs[c]  
   bed.gc.chr <- subset(bed.gc, CHR == chr)

   bed.gc.chr$TEST <- mapply(x = 1:nrow(bed.gc.chr), function(x) (bed.gc.chr$START[x+1] - bed.gc.chr$START[x]))
   print(which(bed.gc.chr$TEST < 0))
   bed.gc.chr$TEST <- mapply(x = 1:nrow(bed.gc.chr), function(x) (bed.gc.chr$END[x] - bed.gc.chr$START[x]))
   print(which(bed.gc.chr$TEST < 0))
}

# -----------------------------------------------------------------------------
# Plot RD and RT (see ReplicationTiming.R)
# Last Modified: 28/05/19; 14/02/19; 10/01/19; 31/08/18; 13/06/17
# -----------------------------------------------------------------------------
nrds <- getLog2ScaledRT(wd.rt.data, base, method, PAIR1, PAIR0, n1, n0, chrs, bed.gc)
save(nrds, file=file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "t-n", ".RData")))
#load(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d.rt.log2s_", "t-n", ".RData")))
# > nrow(nrds)
# [1] 2656170 - 22
nrds.sclc <- nrds

ymax <- 0.6
ymin <- 0.15
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   lcl.rt.chr <- subset(lcl.rt, CHR == chr)   ## Koren 2012
 
   ## Plot RT
   main.text <- paste0(BASE, " T/N read depth ratio between tumour (n=", n1, ") and normal (n=", n0, ") samples")  
   file.name <- file.path(wd.rt.plots, paste0("RT_", base, "_", method, ".d.rt.log2s_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ""))   
   plotRT(file.name, main.text, chr, NA, NA, nrds.chr, bed.gc.chr, c("red", "blue"), c("Tumour", "Normal"), c("lightcoral", "lightskyblue3"), c("T", "N"), "png", width=10, peaks=c(), ylim=c(ymin, ymax), lcl.rt.chr)
}

# -----------------------------------------------------------------------------
# RD vs RT (RDS and SPR)
# Last Modified: 11/07/19; 31/05/19
# -----------------------------------------------------------------------------
sprs <- getSPR(nrds, bed.gc)
save(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.RData")))
writeTable(sprs, file=file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.txt")), colnames=T, rownames=F, sep="\t")
#load(file.path(wd.rt.data, paste0("rd-vs-rt_", base, "-t-n_spline_spearman.RData")))

for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   nrds.chr <- nrds[intersect(nrds$BED, rownames(bed.gc.chr)),]
   nrds.chr.T  <- setSpline(nrds.chr, bed.gc.chr, "T")
   nrds.chr.N  <- setSpline(nrds.chr, bed.gc.chr, "N")
   nrds.chr.RT <- setSpline(nrds.chr, bed.gc.chr, "RT")
 
   main.text <- c(paste0("SCLC read depth correlation (", "Chr", c, ")"), paste0("rho = ", round0(sprs$cor[c], digits=2), " (T vs. N)"))
   xlab.text <- "T/N read depth ratio [log2]"
   ylab.text <- "Read depth [RPKM]"
   file.name <- file.path(wd.rt.plots, "chrs", paste0("RD-vs-RT_SCLC-T-N_chr", c, "_spline_spearman"))
   plotRD2vsRT(nrds.chr.T$SPLINE, nrds.chr.N$SPLINE, nrds.chr.RT$SPLINE, file.name, main.text, ylab.text, xlab.text, c("red", "blue"), c("T", "N"), method="spearman")
}

## S-phase progression rate (SPR)
ylab.text <- "SPR"
file.name <- file.path(wd.rt.plots, "SPR_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " S-phase progression rate"), "SPR = (E-L)/(E+L)")
plotSPR(sprs, file.name, main.text, c(13, 17), digits=3, unit=5, ylab.text)

## SPR vs Read depth correlation
file.name <- file.path(wd.rt.plots, "SPR-RDC_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Read depths correlation"), "")
xlab.text <- "T vs. N [rho]"
plotSPRRDC(sprs$spr, sprs$cor, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

## SPR vs Woodfine 2004
file.name <- file.path(wd.rt.plots, "SPR-Woodfine_SCLC-T-N_spline_spearman")
main.text <- c(paste0(BASE, " SPR vs. Woodfine et al. 2004"), "")
xlab.text <- "Mean replication timing ratio"
plotSPRRDC(sprs$spr, lcl.mean$Mean, file.name, main.text, c(4, 13, 17, 19, 22), xlab.text, unit=5, ylab.text)

# -----------------------------------------------------------------------------
# RT vs LCL S/G1
# Last Modified: 27/05/19
# -----------------------------------------------------------------------------
nrds.tmp <- nrds
load(file.path(wd, "LCL/analysis/replication/lcl-wgs-rt/data/lcl_rpkm.gc.cn.d.rt.log2s_s-g1.RData"))
nrds.lcl <- nrds
nrds <- nrds.tmp

cors <- getRTvsRT(nrds, nrds.lcl, bed.gc)
save(cors, file=file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-t-n-vs-lcl-s-g1_spline_spearman.RData")))
#load(file.path(wd.rt.data, paste0("rt-vs-rt_", base, "-t-n-vs-lcl-s-g1_spline_spearman.RData")))

ylab.text <- "Spearman's rho"
xlab.text <- "Chromosome"
file.name <- file.path(wd.rt.plots, "RT-vs-RT_SCLC-T-N-vs-LCL-S-G1_spline_spearman")
main.text <- paste0("SCLC T/N vs. LCL S/G1")
ymin <- 0.25
ymax <- 1.05
plotRTvsRTALL(cors, file.name, main.text, ylab.text, xlab.text, ymin, ymax, col="black", c=2, pos=1)

##
file.name <- file.path(wd.rt.plots, "RTD-vs-RT_SCLC-T-N-vs-LCL-S-G1_spline_spearman")
ymin <- -1.1
ymax <- 1.1
plotRD3vsRTALL(cors, file.name, main.text, ymin, ymax, cols=c("red", "blue", "black"), c("T", "N", "T/N"), c=NA, isRT=T)

# -----------------------------------------------------------------------------
# SCLC T vs LCL S/G1
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# Last Modified: 05/06/19; 20/04/19; 06/03/19
# -----------------------------------------------------------------------------
cors.samples <- getSAMPLEvsRT(wd.rt.data, samples1)
save(cors.samples, file=file.path(wd.rt.data, paste0("samples-vs-rt_", base, "-vs-lcl_spline_spearman.RData")))
# > min(cors.samples[,-c(1:4)])
# [1] -0.8457817
# > max(cors.samples[,-c(1:4)])
# [1] 0.8088429

#load(file.path(wd.rt.data, paste0("samples-vs-rt_sclc-vs-lcl_spline_spearman.RData")))
file.name <- file.path(wd.rt.plots, "SAMPLES-vs-RT_SCLC-vs-LCL_spline_spearman")
main.text <- c("SCLC read depth vs. RT", "")
ymin <- -0.8773492
ymax <- 0.8392611
plotSAMPLEvsRT(cors.samples, samples1, file.name, main.text, ymin, ymax)

# -----------------------------------------------------------------------------
# Last Modified: 30/11/20; 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
nrds <- getOverallRD(wd.rt.data, base, method, PAIR1, n1)
nrow(nrds)
# > nrow(nrds)
# [1] 2684771
nrds.m <- subset(nrds, MEDIAN != 0)   ## ADD 30/11/20
# > nrow(nrds.m)
# [1] 2673828
nrds.d <- getOverallDetectedRD(wd.rt.data, base, method, PAIR1, n1, samples1)
nrow(nrds.d)
# > nrow(nrds.d)
# [1] 2674140

# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 18/11/19; 16/06/19; 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.sclc <- setSamplesQ4(wd.rt.data, samples1)
writeTable(samples.sclc, file.path(wd.ngs, "sclc_wgs_n101.txt"), colnames=T, rownames=F, sep="\t")
#         0%        25%        50%        75%       100% 
# -0.7438588 -0.6840244 -0.6474195 -0.5555826  0.6806205 

writeTable(subset(samples.sclc, Q4 %in% c(4,1)), file.path(wd.ngs, "sclc_wgs_q4_n51.txt"), colnames=T, rownames=F, sep="\t")
writeTable(subset(samples.sclc, Q4 %in% c(3,1)), file.path(wd.ngs, "sclc_wgs_q3_n51.txt"), colnames=T, rownames=F, sep="\t")

## Random 25/26
m2.25 <- sort(rownames(subset(samples.sclc, M2 == 1))[sample(1:50, round(50/2), replace=F)])
m1.26 <- sort(rownames(subset(samples.sclc, M2 == 0))[sample(1:51, round(51/2), replace=F)])

random51 <- sort(c(m2.25, m1.26))
writeTable(samples.sclc[random51,], file.path(wd.ngs, "sclc_wgs_n101-random51.txt"), colnames=T, rownames=F, sep="\t")

random50 <- sort(setdiff(rownames(samples.sclc), random51))
writeTable(samples.sclc[random50,], file.path(wd.ngs, "sclc_wgs_n101-random50.txt"), colnames=T, rownames=F, sep="\t")

## Random 12/13
m2 <- sort(rownames(subset(samples.sclc, M2 == 1)))
m1 <- sort(rownames(subset(samples.sclc, M2 == 0)))

m2.12 <- sort(m2[sample(1:50, round(50/4), replace=F)])
m1.13 <- sort(m1[sample(1:51, round(51/4), replace=F)])
random3 <- sort(c(m2.12, m1.13))
writeTable(samples.sclc[random3,], file.path(wd.ngs, "sclc_wgs_n101-random3.txt"), colnames=T, rownames=F, sep="\t")

m2 <- sort(setdiff(m2, m2.12))
m1 <- sort(setdiff(m1, m1.13))

m2.12 <- sort(m2[sample(1:38, round(50/4), replace=F)])
m1.13 <- sort(m1[sample(1:38, round(51/4), replace=F)])
random4 <- sort(c(m2.12, m1.13))
writeTable(samples.sclc[random4,], file.path(wd.ngs, "sclc_wgs_n101-random4.txt"), colnames=T, rownames=F, sep="\t")

## Random 6/7
m2 <- sort(rownames(subset(samples.sclc, M2 == 1)))
m1 <- sort(rownames(subset(samples.sclc, M2 == 0)))

m2.6 <- sort(m2[sample(1:50, round(50/8), replace=F)])
m1.6 <- sort(m1[sample(1:51, round(51/8), replace=F)])
random5 <- sort(c(m2.6, m1.6))
writeTable(samples.sclc[random5,], file.path(wd.ngs, "sclc_wgs_n101-random5.txt"), colnames=T, rownames=F, sep="\t")

m2 <- sort(setdiff(m2, m2.6))
m1 <- sort(setdiff(m1, m1.6))

m2.6 <- sort(m2[sample(1:44, round(50/8), replace=F)])
m1.6 <- sort(m1[sample(1:45, round(51/8), replace=F)])
random6 <- sort(c(m2.6, m1.6))
writeTable(samples.sclc[random6,], file.path(wd.ngs, "sclc_wgs_n101-random6.txt"), colnames=T, rownames=F, sep="\t")




# -----------------------------------------------------------------------------
# PCA
# Last Modified: 04/06/19; 21/04/19
# -----------------------------------------------------------------------------
## Refer to cmd-rt_2a_nrd.gc.cn.d_sample.R (commandline mode)
nrds.T.chr.d.all <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", "chr1", "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
for (c in 2:22) {
   chr <- chrs[c]
 
   ## Read depth
   #nrds.T.chr.d <- pipeGetDetectedRD(wd.ngs.data, BASE, chr, PAIR1, method)
   nrds.T.chr.d <- readTable(file.path(wd.rt.data, paste0(base, "_", method, ".gc.cn.d_", chr, "_", PAIR1, "_n", n1, ".txt.gz")), header=T, rownames=F, sep="\t")
 
   nrds.T.chr.d.all <- rbind(nrds.T.chr.d.all, nrds.T.chr.d)
}

##
test <- nrds.T.chr.d.all[, -1]   ## BUG 2019/10/14: Remove column BED
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc_chrs.RData")))

#load(file.path(wd.rt.data, paste0("pca_sclc_chrs.RData")))
file.main <- c("SCLC", "")
trait <- as.numeric(samples.sclc$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_SCLC", size=6, file.main, "bottomright", c(red, red.lighter, blue.lighter, blue), NULL, flip.x=1, flip.y=1, legend.title=NA)

## SG1
#trait <- samples.sclc.sg1$SG1
#plotPCA(1, 2, pca.de, trait, wd.rt.plots, "pca_sclc_T_chrs_spline_spearman_SG1", size=6, file.main, "bottomright", c("red", "lightgray", "blue"), NULL, flip.x=1, flip.y=1, legend.title="Consist. CM in all chrs")

# -----------------------------------------------------------------------------
# Beeswarm plots
# Last Modified: 21/04/19
# -----------------------------------------------------------------------------
samples.sclc <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4)

#install.packages('beeswarm')
library(beeswarm)

pdf(file.path(wd.rt.plots, "beeswarm_sclc+nbl+cll_Spearman's_BIGER_tumor.pdf"), height=6, width=6)
ymax <- 0.8   #max(samples$COR)
ymin <- -ymax
boxplot(COR ~ CANCER, data=samples, outline=F, names=c("", "", ""), ylim=c(ymin, ymax), ylab="", main="Tumor read depth vs. RT", yaxt="n", cex.axis=1.8, cex.lab=1.9, cex.main=2.1)
abline(h=0, lty=5, lwd=2)

legend("topright", legend=c("Q4", "Q3", "Q2", "Q1"), pch=19, pt.cex=2.5, col=c(red, red.lighter, blue.lighter, blue), cex=1.9)

beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 1), col=blue, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 2), col=blue.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 3), col=red.lighter, pch=19, cex=1, add=T)
beeswarm(COR ~ CANCER, data=subset(samples, Q4 == 4), col=red, pch=19, cex=1, add=T)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Spearman's rho", side=2, line=2.7, cex=1.9)
#mtext("", cex=1.2, line=0.3)
mtext(text=c("SCLC", "NBL", "CLL"), side=1, cex=1.9, line=1.3, at=c(1,2,3))
mtext(text=c("n=101", "n=56", "n=96"), side=1, cex=1.8, line=3, at=c(1,2,3))
dev.off()

# -----------------------------------------------------------------------------
# PCA for all
# Last Modified: 11/03/20
# -----------------------------------------------------------------------------
# > dim(test.sclc)
# [1] 2657164     102
# > dim(test.nbl)
# [1] 2659570      57
# > dim(test.cll)
# [1] 2653842      97
# > dim(test.sclc.nl)
# [1] 2658679      93

overlaps <- intersect(intersect(intersect(rownames(test.sclc), rownames(test.nbl)), rownames(test.cll)), rownames(test.sclc.nl))
# > length(overlaps)
# [1] 2653842
test.sclc <- test.sclc[overlaps, -1]
test.nbl  <- test.nbl[overlaps, -1]
test.cll  <- test.cll[overlaps, -1]
test.sclc.nl <- test.sclc.nl[overlaps, -1]

test <- cbind(test.sclc, test.nbl, test.cll, test.sclc.nl)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.rt.data, paste0("pca_sclc+nbl_cll+sclc.nl.RData")))

# -----------------------------------------------------------------------------
# PCA ALL+WB
# Last Modified: 12/03/20
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0("pca_sclc+nbl+cll+sclc.nl+nbl.wb+cll.wb_pca.de.RData")))

samples.sclc <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
samples.sclc.nl <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n92_NL.txt", header=T, rownames=T, sep="")
samples.nbl.wb  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1_N.txt", header=T, rownames=T, sep="")
samples.cll.wb  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96_N.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)
n.sclc.nl <- nrow(samples.sclc.nl)
n.nbl.wb  <- nrow(samples.nbl.wb)
n.cll.wb  <- nrow(samples.cll.wb)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)
rownames(samples) <- c(rownames(samples.sclc), rownames(samples.nbl), rownames(samples.cll), paste0(rownames(samples.sclc.nl), "-NL"), paste0(rownames(samples.nbl.wb), "-WB"), paste0(rownames(samples.cll.wb), "-WB"))
samples$SAMPLE_ID <- rownames(samples)

## ALL
q <- quantile(as.numeric(samples$COR))
print(q)

samples.q4 <- list()
samples.q4[[4]] <- rownames(subset(samples, COR > as.numeric(q[4])))
samples.q4[[3]] <- rownames(subset(subset(samples, COR > as.numeric(q[3])), COR <= as.numeric(q[4])))
samples.q4[[2]] <- rownames(subset(subset(samples, COR > as.numeric(q[2])), COR <= as.numeric(q[3])))
samples.q4[[1]] <- rownames(subset(samples, COR <= as.numeric(q[2])))

samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[4]])] <- 4
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[3]])] <- 3
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[2]])] <- 2
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[1]])] <- 1

file.main <- c("Total (n=497) tumour-normal profiles ", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL+WB", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "")

###
## SCLC
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 0),]$Q4 <- NA
file.main <- c("SCLC (n=101) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC")

## NBL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 1),]$Q4 <- NA
file.main <- c("NBL (n=56) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_NBL", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "NBL")

## CLL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 2),]$Q4 <- NA
file.main <- c("CLL (n=96) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_CLL", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "CLL")

## SCLC-NL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 3),]$Q4 <- NA
file.main <- c("SCLC-NL (n=92) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC-NL", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC-NL")

## NBL-WB
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 4),]$Q4 <- NA
file.main <- c("NBL-WB (n=56) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_NBL-WB", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "NBL-WB")

## CLL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb)] <- 4
samples$CANCER[(1+n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb):(n.sclc+n.nbl+n.cll+n.sclc.nl+n.nbl.wb+n.cll.wb)] <- 5
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR, samples.nbl.wb$COR, samples.cll.wb$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4, samples.nbl.wb$Q4, samples.cll.wb$Q4)

samples[which(samples$CANCER != 5),]$Q4 <- NA
file.main <- c("CLL-WB (n=96) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 3, pca.de, trait, wd.rt.plots, "PCA_ALL_CLL-WB", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "CLL-WB")











# -----------------------------------------------------------------------------
# PCA ALL
# Last Modified: 11/03/20
# -----------------------------------------------------------------------------
#load(file.path(wd.rt.data, paste0("pca_sclc+nbl+cll+sclc.nl_pca.de.RData")))

samples.sclc <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n101.txt", header=T, rownames=T, sep="")
samples.nbl  <- readTable("/projects/cangen/tyang2/NBL/ngs/WGS/nbl_wgs_n57-1.txt", header=T, rownames=T, sep="")
samples.cll  <- readTable("/projects/cangen/tyang2/CLL/ngs/WGS/cll_wgs_n96.txt", header=T, rownames=T, sep="")
samples.sclc.nl <- readTable("/projects/cangen/tyang2/SCLC/ngs/WGS/sclc_wgs_n92_NL.txt", header=T, rownames=T, sep="")
n.sclc <- nrow(samples.sclc)
n.nbl  <- nrow(samples.nbl)
n.cll  <- nrow(samples.cll)
n.sclc.nl <- nrow(samples.sclc.nl)

samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)
rownames(samples) <- c(rownames(samples.sclc), rownames(samples.nbl), rownames(samples.cll), rownames(samples.sclc.nl))
samples$SAMPLE_ID <- rownames(samples)
 
## ALL
q <- quantile(as.numeric(samples$COR))
print(q)

samples.q4 <- list()
samples.q4[[4]] <- rownames(subset(samples, COR > as.numeric(q[4])))
samples.q4[[3]] <- rownames(subset(subset(samples, COR > as.numeric(q[3])), COR <= as.numeric(q[4])))
samples.q4[[2]] <- rownames(subset(subset(samples, COR > as.numeric(q[2])), COR <= as.numeric(q[3])))
samples.q4[[1]] <- rownames(subset(samples, COR <= as.numeric(q[2])))

samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[4]])] <- 4
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[3]])] <- 3
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[2]])] <- 2
samples$Q4[which(samples$SAMPLE_ID %in% samples.q4[[1]])] <- 1

file.main <- c("SCLC-NL, SCLC, NBL and CLL (n=345)", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA)

###
## SCLC
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 0),]$Q4 <- NA
file.main <- c("SCLC (n=101) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC")

## NBL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 1),]$Q4 <- NA
file.main <- c("NBL (n=56) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_NBL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "NBL")

## CLL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 2),]$Q4 <- NA
file.main <- c("CLL (n=96) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_CLL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "CLL")

## SCLC-NL
samples <- toTable(0, 3, n.sclc+n.nbl+n.cll+n.sclc.nl, c("CANCER", "COR", "Q4"))
samples$CANCER[1:n.sclc] <- 0
samples$CANCER[(1+n.sclc):(n.sclc+n.nbl)] <- 1
samples$CANCER[(1+n.sclc+n.nbl):(n.sclc+n.nbl+n.cll)] <- 2
samples$CANCER[(1+n.sclc+n.nbl+n.cll):(n.sclc+n.nbl+n.cll+n.sclc.nl)] <- 3
samples$COR <- c(samples.sclc$COR, samples.nbl$COR, samples.cll$COR, samples.sclc.nl$COR)
samples$Q4  <- c(samples.sclc$Q4, samples.nbl$Q4, samples.cll$Q4, samples.sclc.nl$Q4)

samples[which(samples$CANCER != 3),]$Q4 <- NA
file.main <- c("SCLC-NL (n=92) read depth profiles", "")
trait <- as.numeric(samples$Q4)
trait[which(trait == 4)] <- "Q4"
trait[which(trait == 3)] <- "Q3"
trait[which(trait == 2)] <- "Q2"
trait[which(trait == 1)] <- "Q1"
plotPCA(1, 2, pca.de, trait, wd.rt.plots, "PCA_ALL_SCLC-NL_Q4", size=6, file.main, "topleft", c("red", "lightpink1", "lightskyblue2", "blue", "lightgray"), NULL, flip.x=-1, flip.y=1, legend.title=NA, "SCLC-NL")






# -----------------------------------------------------------------------------
# CM2 and CQ4
# Last Modified: 19/05/19
# -----------------------------------------------------------------------------
cors <- t(cors.samples[, -c(1:4)])
#cors <- cors[overlaps,]   ## Overlapping between WGS and RNA (n=70); after sclc-tpm.de.R
colnames(cors) <- paste0("chr", c(1:22))
cq4 <- cors
cm2 <- cors

for (c in 1:22) {
 cors.chr <- cors[,c]
 q <- quantile(as.numeric(cors.chr))
 
 cq4[which(cors[, c] > q[4]), c] <- 4
 cq4[intersect(as.vector(which(cors[, c] > q[3])), as.vector(which(cors[, c] <= q[4]))), c] <- 3
 cq4[intersect(as.vector(which(cors[, c] > q[2])), as.vector(which(cors[, c] <= q[3]))), c] <- 2
 cq4[which(cors[, c] <= q[2]), c] <- 1
 
 cm2[which(cors[, c] > q[3]), c] <- 1
 cm2[which(cors[, c] <= q[3]), c] <- 0
}
cq4 <- as.data.frame(cq4)
cq4$SAMPLE_ID <- ""
cq4$SAMPLE_ID <- rownames(cors)
cm2 <- as.data.frame(cm2)
cm2$SAMPLE_ID <- ""
cm2$SAMPLE_ID <- rownames(cors)
writeTable(cq4[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "sclc_wgs_n101.cq4"), colnames=T, rownames=F, sep="\t")
writeTable(cm2[, c("SAMPLE_ID", paste0("chr", c(1:22)))], file.path(wd.ngs, "sclc_wgs_n101.cm2"), colnames=T, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# Find S-like and G1-like tumour samples
# Last Modified: 04/06/19; 06/03/19
# -----------------------------------------------------------------------------
samples.sclc.sg1 <- setSamplesSG1(wd.rt.data, samples1, cors.samples)
writeTable(samples.sclc.sg1, file.path(wd.ngs, "sclc_wgs_n101.sg1"), colnames=T, rownames=F, sep="\t")
# > length(s_likes)
# [1] 24
# > length(g1_likes)
# [1] 14
# > s_likes
# [1] "S00339" "S00825" "S00827" "S00833" "S00837" "S00838" "S00938" "S00941" "S01020" "S01022" "S01366" "S01524" "S01864" "S02120" "S02219" "S02285" "S02286" "S02287" "S02289" "S02290" "S02291" "S02292"
# [23] "S02295" "S02296"
# > g1_likes
# [1] "S00050" "S00472" "S00831" "S00944" "S01861" "S02139" "S02237" "S02241" "S02248" "S02342" "S02344" "S02352" "S02353" "S02378"
