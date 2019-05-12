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
BASE1 <- "SCLC"
PAIR1 <- "T"
BASE0 <- "SCLC"
PAIR0 <- "N"
base1 <- tolower(BASE1)
base0 <- tolower(BASE0)

wd.wgs  <- file.path(wd, BASE1, "ngs/WGS")
wd.anlys <- file.path(wd, BASE1, "analysis")
wd.meta  <- file.path(wd, BASE1, "metadata", "George 2015")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base1, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

samples <- readTable(file.path(wd.wgs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
phenos  <- readTable(file.path(wd.meta, "nature14664-s1_ST1.txt"), header=T, rownames=T, sep="\t")
phenos  <- survSCLC(phenos, samples, isCensored=F)

purities <- readTable(file.path(wd.meta, "nature14664-s1_ST2.txt"), header=T, rownames=T, sep="\t")
phenos <- cbind(phenos, purities[rownames(phenos),])

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
txt <- readTable(file.path(wd.meta, "nature14664-s1_ST3.txt.gz"), header=T, rownames=F, sep="")
txt <- subset(txt, PAT_ID %in% rownames(samples))   ## ADD 05/05/19
muts <- unique(txt$Gene_Hugo)

## Find only expressed genes
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
genes <- getGenes(muts)
overlaps <- intersect(genes$ensembl_gene_id, rownames(tpm.gene))
# > length(muts)
# [1] 12199
# > nrow(genes)
# [1] 11912
# > length(overlaps)
# [1] 10093
muts <- genes[overlaps,]$external_gene_name
# > length(muts)
# [1] 10093
muts <- unique(genes[overlaps,]$external_gene_name)
# > length(muts)
# [1] 10076

samples.rts <- list(rownames(subset(samples, Q4 == 4)), rownames(subset(samples, Q4 == 3)), rownames(subset(samples, Q4 == 2)), rownames(subset(samples, Q4 == 1)))
txts  <- toTable(0, nrow(samples), length(muts), c(rownames(subset(samples, Q4 == 4)), rownames(subset(samples, Q4 == 3)), rownames(subset(samples, Q4 == 2)), rownames(subset(samples, Q4 == 1))))
txts2 <- toTable(0, nrow(samples), length(muts), c(rownames(subset(samples, Q4 == 4)), rownames(subset(samples, Q4 == 3)), rownames(subset(samples, Q4 == 2)), rownames(subset(samples, Q4 == 1))))
rownames(txts)  <- muts
rownames(txts2) <- muts
for (q in 4:1) {
   samples.rt <- samples.rts[[q]]

   for (s in 1:length(samples.rt)) {
      sample <- samples.rt[s]
      #txt <- read.peiflyne.muts.txt(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_muts.txt")), type="SNM")
      txt.sample <- subset(txt, PAT_ID == sample)
      mut <- data.frame(table(txt.sample$Gene_Hugo))
      rownames(mut) <- mut$Var1
      mut.o <- intersect(mut$Var1, muts)
      mut <- mut[mut.o,]

      txts[rownames(mut), sample]  <- mut$Freq
      txts2[rownames(mut), sample] <- 1
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
writeTable(txts.out, file.path(wd.rt.data, "muts_q4_expressed_ST3.txt"), colnames=T, rownames=T, sep="\t")

sig <- readTable(file.path(wd.meta, "nature14664-s1_ST6.txt"), header=F, rownames=F, sep="")
sig <- unique(sig)
writeTable(txts.out[sig,], file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6_n101.txt"), colnames=T, rownames=T, sep="\t")

mut <- txts.out[sig,]
#mut <- txts.out
mut$M2 <- mut[, "Q4_S"] + mut[, "Q3_S"]
mut$M1 <- mut[, "Q2_S"] + mut[, "Q1_S"]
m2 <- length(samples.rts[[4]]) + length(samples.rts[[3]])
m1 <- length(samples.rts[[2]]) + length(samples.rts[[1]])
q4 <- length(samples.rts[[4]])
q3 <- length(samples.rts[[3]])
q2 <- length(samples.rts[[2]])
q1 <- length(samples.rts[[1]])

###
## Fisher's test for each mutations
results.gene <- toTable(0, 2, nrow(mut), c("P", "BH"))
rownames(results.gene) <- rownames(mut)
for (g in 1:nrow(mut)) {
   test <- toTable(0, 2, 2, c("WT", "MUT"))
   rownames(test) <- c("M2", "M1")
   
   test[1, 2] <- mut[g, "Q4_S"] + mut[g, "Q3_S"]
   test[1, 1] <- m2 - test[1, 2] 
   test[2, 2] <- mut[g, "Q2_S"] + mut[g, "Q1_S"]
   test[2, 1] <- m1 - test[2, 2]
 
   results.gene[g, 1] <- fisher.test(test)[[1]]
}
results.gene <- results.gene[order(results.gene$P),]
results.gene <- cbind(mut[rownames(results.gene), c("M2", "M1")], results.gene)
results.gene$TOTAL <-0
results.gene$TOTAL <- results.gene$M2 + results.gene$M1
results.gene <- subset(results.gene, TOTAL > 10)
results.gene$BH <- p.adjust(results.gene$P, "BH")
writeTable(results.gene, file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6_M2-M1_10genes.txt"), colnames=T, rownames=T, sep="\t")
genes10 <- rownames(results.gene)

###
## Fisher's test for each mutations
results.gene <- toTable(0, 2, nrow(mut), c("P", "BH"))
rownames(results.gene) <- rownames(mut)
for (g in 1:nrow(mut)) {
   test <- toTable(0, 2, 2, c("WT", "MUT"))
   rownames(test) <- c("Q4", "Q1")
 
   test[1, 2] <- mut[g, "Q4_S"]
   test[1, 1] <- q4 - test[1, 2] 
   test[2, 2] <- mut[g, "Q1_S"]
   test[2, 1] <- q1 - test[2, 2]
 
   results.gene[g, 1] <- fisher.test(test)[[1]]
}
results.gene <- results.gene[order(results.gene$P),]
results.gene <- results.gene[genes10,]
results.gene$BH <- p.adjust(results.gene$P, "BH")

results.gene <- cbind(mut[rownames(results.gene), c("Q4_S", "Q1_S")], results.gene)
writeTable(results.gene, file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6_Q4-Q1_10genes.txt"), colnames=T, rownames=T, sep="\t")

###
## Fisher's test for each mutations
results.gene <- toTable(0, 2, nrow(mut), c("P", "BH"))
rownames(results.gene) <- rownames(mut)
for (g in 1:nrow(mut)) {
   test <- toTable(0, 2, 2, c("WT", "MUT"))
   rownames(test) <- c("Q41_S", "Q23_S")
 
   test[1, 2] <- mut[g, "Q4_S"] + mut[g, "Q1_S"]
   test[1, 1] <- q4 + q1 - test[1, 2] 
   test[2, 2] <- mut[g, "Q2_S"] + mut[g, "Q3_S"]
   test[2, 1] <- q2 + q3 - test[2, 2]
 
   results.gene[g, 1] <- fisher.test(test)[[1]]
}
results.gene <- results.gene[order(results.gene$P),]
results.gene <- results.gene[genes10,]
results.gene$BH <- p.adjust(results.gene$P, "BH")

results.gene <- cbind(mut[rownames(results.gene), c("Q4_S", "Q1_S", "Q2_S", "Q3_S")], results.gene)
writeTable(results.gene, file.path(wd.rt.data, "muts_q4_expressed_ST3_ST6_Q41-Q23_10genes.txt"), colnames=T, rownames=T, sep="\t")

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
# Purities
# Last Modified: 02/05/19
# -----------------------------------------------------------------------------
purities <- readTable(file.path(wd.meta, "nature14664-s1_ST2.txt"), header=T, rownames=T, sep="\t")
phenos <- cbind(phenos[samples$SAMPLE_ID,], purities[samples$SAMPLE_ID,])







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
