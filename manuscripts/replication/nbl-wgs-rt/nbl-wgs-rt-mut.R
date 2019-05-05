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
handbooks  <- c("Commons.R", "Mutation.R", "Survival.R")
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
BASE1 <- "NBL"
PAIR1 <- "T"
BASE0 <- "NBL"
PAIR0 <- "N"
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

wd.meta  <- file.path(wd, BASE1, "metadata", "Peifer 2015")
wd.anlys <- file.path(wd, BASE1, "analysis")
wd.rt       <- file.path(wd.anlys, "replication", paste0(base1, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs <- file.path(wd, BASE1, "ngs/WGS")
samples <- readTable(file.path(wd.ngs, "nbl_wgs_n57-1.txt"), header=T, rownames=T, sep="")
maps    <- readTable(file.path(wd.meta, "Samples_NEW_NUMBERS.txt"), header=T, rownames=T, sep="")
samples$SAMPLE_NEW <- maps[samples$SAMPLE_ID,]$NEW
samples$MUTS <- 0

# -----------------------------------------------------------------------------
# Mutation burdens between T29 and T33
# Last Modified: 15/03/19
# -----------------------------------------------------------------------------
#muts <- c()
#for (s in 1:nrow(samples)) {
#   sample <- samples$SAMPLE_ID[s]
#   txt <- read.peiflyne.muts.txt(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_muts.txt")), type="SNM")
# 
#   muts <- unique(c(muts, txt$Gene_Hugo))
#}
txt <- readTable(file.path(wd.meta, "nature14980-s2.txt"), header=T, rownames=F, sep="")
muts <- unique(txt$Gene)

## Find only expressed genes
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
genes <- getGenes(muts)
overlaps <- intersect(genes$ensembl_gene_id, rownames(tpm.gene))
# > length(muts)
# [1] 855
# > nrow(genes)
# [1] 837
# > length(overlaps)
# [1] 715
muts <- genes[overlaps,]$external_gene_name
# > length(muts)
# [1] 715
muts <- unique(genes[overlaps,]$external_gene_name)
# > length(muts)
# [1] 714

samples.rts <- list(rownames(subset(samples, Q4 == 4)), rownames(subset(samples, Q4 == 3)), rownames(subset(samples, Q4 == 2)), rownames(subset(samples, Q4 == 1)))
txts  <- toTable(0, nrow(samples), length(muts), c(rownames(subset(samples, Q4 == 4)), rownames(subset(samples, Q4 == 3)), rownames(subset(samples, Q4 == 2)), rownames(subset(samples, Q4 == 1))))
txts2 <- toTable(0, nrow(samples), length(muts), c(rownames(subset(samples, Q4 == 4)), rownames(subset(samples, Q4 == 3)), rownames(subset(samples, Q4 == 2)), rownames(subset(samples, Q4 == 1))))
rownames(txts)  <- muts
rownames(txts2) <- muts
for (q in 4:1) {
   samples.rt <- samples.rts[[q]]

   for (s in 1:length(samples.rt)) {
      sample <- samples.rt[s]
      sample_new <- samples[sample,]$SAMPLE_NEW
      #txt <- read.peiflyne.muts.txt(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_muts.txt")), type="SNM")
      txt.sample <- subset(txt, Sample == sample_new)
      if (nrow(txt.sample) != 0) {
         mut <- data.frame(table(txt.sample$Gene))
         rownames(mut) <- mut$Var1
         mut.o <- intersect(mut$Var1, muts)
         mut <- mut[mut.o,]

         txts[rownames(mut), sample]  <- mut$Freq
         txts2[rownames(mut), sample] <- 1
      }
      
      samples[sample,]$MUTS <- nrow(txt.sample)
   }
}

txts.out <- toTable(0, 8, length(muts), c("Q4_N", "Q4_S", "Q3_N", "Q3_S","Q2_N", "Q2_S", "Q1_N", "Q1_S"))
rownames(txts.out)  <- muts
for (q in 4:1) {
   txts.q  <- txts[,samples.rts[[q]]]
   txts2.q <- txts2[,samples.rts[[q]]]
   
   txts.out[, paste0("Q", q, "_N")] <- mapply(x = 1:nrow(txts.q), function(x) sum(txts.q[x,]))
   txts.out[, paste0("Q", q, "_S")] <- mapply(x = 1:nrow(txts2.q), function(x) sum(txts2.q[x,]))   
}
writeTable(txts.out, file.path(wd.rt.data, "muts_q4_expressed_ST2.txt"), colnames=T, rownames=T, sep="\t")






sig <- readTable(file.path(wd.meta, "nature14664-s1_ST6.txt"), header=F, rownames=F, sep="")
sig <- unique(sig)
writeTable(txts.out[sig,], file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6.txt"), colnames=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Kaplan-Meier survival analysis (OS~kmeans)
# Last Modified: 02/05/19
# -----------------------------------------------------------------------------
phenos <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
phenos <- phenos[samples$SAMPLE_ID,]
extension <- ""
#phenos <- subset(phenos, sex == "male")
#extension <- "_sex_male_chemotherapy_no"
phenos <- subset(phenos, chemotherapy..yes.no. == "no")
extension <- "_chemotherapy_no"
isRT <- T

## OS_censor
phenos.sub <- phenos[!is.na(phenos$overall_survival..months.),]
phenos.sub$OS_month <- phenos.sub$overall_survival..months.
phenos.sub$Groups <- samples[rownames(phenos.sub),]$Q4
if (isRT)
   phenos.sub$Groups <- samples[rownames(phenos.sub),]$RT

phenos.sub$OS_censor <- phenos.sub$Status..at.time.of.last.follow.up.
phenos.sub$OS_censor <- gsub("dead", 0, phenos.sub$OS_censor)
phenos.sub$OS_censor <- gsub("alive", 1, phenos.sub$OS_censor)
phenos.sub$OS_censor <- as.numeric(phenos.sub$OS_censor)

##
fit <- survfit(Surv(OS_month, OS_censor) ~ Groups, data=phenos.sub) 
file.name <- file.path(wd.rt.plots, "survfit", paste0("survfit_", base, "_OS~Groups", extension))
main.text <- paste0(BASE, " overall survival (n=", nrow(phenos.sub), ")")
plotSurvfitRT(fit, file.name, main.text, isRT)











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
