# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/replication/icgc-wgs-rt/icgc-wgs.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 19/02/22
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"   ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"        ## tyang2@gauss
wd.src <- "/Users/ty2/Work/dev/R"            ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Survival.R", "Transcription.R", "survminer/ggforest.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 18/02/22
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"            ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                 ## tyang2@gauss
#wd <- "/Users/ty2/Work/uni-koeln/tyang2"   ## tpyang@localhost
wd <- "/Users/ty2/Work/uni-koeln/tyang2"    ## tp2@localhost
BASE  <- "ICGC"
base  <- tolower(BASE)

#wd.ngs <- file.path(wd, BASE, "ngs/WGS")
wd.meta  <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
#wd.rt.plots <- file.path(wd.rt, "plots", "PanImmune_2024")
wd.rt.plots <- file.path(wd.rt, "plots")

# -----------------------------------------------------------------------------
# 
# Last Modified: 29/03/22
# -----------------------------------------------------------------------------
raws <- readTable(file.path(wd.meta, "PCAWG.raw.data.copy_number.converted_data.list"), header=F, rownames=F, sep="\t")
length(raws)
# [1] 2951
segs <- readTable(file.path(wd.meta, "PCAWG.raw.data.copy_number.sclust_final_copy_number_analysis_files.list"), header=F, rownames=F, sep="\t")
length(segs)
# [1] 2778
overlaps <- intersect(segs, gsub("_CONVERTED", "", raws))
length(overlaps)
# [1] 2777

totals <- readTable(file.path(wd.meta, "PCAWG.wgs_id_histology_appreviations.txt"), header=T, rownames=T, sep="\t")
nrow(totals)
# [1] 2883

overlaps.2 <- intersect(rownames(totals), overlaps)
totals <- totals[overlaps.2,]
totals$project_code <- ""
nrow(totals)
# [1] 2747

mappings <- readTable(file.path(wd.meta, "PCAWG.wgs_id_histology_table_mached.txt"), header=T, rownames=F, sep="\t")
nrow(mappings)
# [1] 3911
idx <- which(totals$wgs_id %in% mappings$pcawg_wgs_id)
totals <- totals[idx,]
nrow(totals)
# [1] 2744

idx <- which(mappings$pcawg_wgs_id %in% totals$wgs_id)
mappings <- mappings[idx,]
nrow(mappings)
# [1] 3744
mappings <- subset(mappings, specimen_library_strategy != "RNA-Seq")
nrow(mappings)
# [1] 2744
rownames(mappings) <- mappings$pcawg_wgs_id
mappings <- mappings[rownames(totals),]

dim(subset(as.data.frame(table(samples$icgc_donor_id)), Freq > 1))
# [1] 53  2
sum(subset(as.data.frame(table(samples$icgc_donor_id)), Freq > 1)$Freq)
# [1] 163
dim(subset(as.data.frame(table(samples.surv$icgc_donor_id)), Freq > 1))
# [1] 38  2
sum(subset(as.data.frame(table(samples.surv$icgc_donor_id)), Freq > 1)$Freq)
# [1] 131
dim(subset(as.data.frame(table(samples.surv.h$icgc_donor_id)), Freq > 1))
# [1] 33  2
sum(subset(as.data.frame(table(samples.surv.h$icgc_donor_id)), Freq > 1)$Freq)
# [1] 120
                
###
## CAUTION: There are multiple tumor_wgs_aliquot_ids
release <- readTable(file.path(wd.meta, "data_release", "release_may2016.v1.4.tsv"), header=T, rownames=F, sep="")
rownames(release) <- release$tumor_wgs_aliquot_id
nrow(release)
# [1] 2834
#idx <- which(table$wgs_id %in% release$tumor_wgs_aliquot_id)   ## f942b732-7c0b-7ec6-e040-11ac0c483f86,f8f024e5-7096-3047-e040-11ac0d481c0f,f8e61a02-6e5e-c8e2-e040-11ac0d481b70
#table <- table[idx,]
#rownames(table) <- table$specimen_id
#nrow(table)
# [1] 2566

for (s in 1:nrow(totals)) {
   codes <- mappings[which(mappings$pcawg_wgs_id == totals$wgs_id[s]),]$project_code
   if (length(codes) > 1) {
      if (codes[1] != codes[2])
         print(codes)
      totals$project_code[s] <- codes[1]
   } else {
      totals$project_code[s] <- codes
   }
}
#sort(table(totals$project_code), decreasing=T)
sort(table(totals$histology_abbreviation), decreasing=T)

###
##
clinicals <- readTable(file.path("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/metadata/DCC_DATA_RELEASE", "pcawg_donor_clinical_August2016_v9.txt"), header=T, rownames=F, sep="\t")
nrow(clinicals)
# [1] 2834
idx2 <- which(clinicals$icgc_donor_id %in% mapping$icgc_donor_id)
clinicals <- clinicals[idx2,]
nrow(clinicals)
# [1] 2757

clinicals$icgc_specimen_id <- NA
for (s in 1:nrow(clinicals)) {
   ids <- mapping[which(mapping$icgc_donor_id == clinicals$icgc_donor_id[s]),]$icgc_specimen_id
   if (length(ids) > 1) {
      clinicals$icgc_specimen_id[s] <- paste0(paste(unique(ids), collapse=","), ",")
   } else {
      clinicals$icgc_specimen_id[s] <- paste0(ids, ",")
   }
}
rownames(clinicals) <- clinicals$icgc_donor_id
writeTable(clinicals, "/Users/tpyang/Work/uni-koeln/tyang2/ICGC/metadata/DCC_DATA_RELEASE/pcawg_donor_clinical_August2016_v9_tyang2.txt", colnames=T, rownames=F, sep="\t")

save(raws, segs, totals, mappings, clinicals, release, file=file.path(wd.rt.data, paste0("icgc_wgs.RData")), version=2)

# -----------------------------------------------------------------------------
# samples
# Last Modified: 06/02/24 (re-run); 21/04/19
# -----------------------------------------------------------------------------
hists <- as.vector(subset(as.data.frame(sort(table(totals$histology_abbreviation))), Freq >= 20)$Var1)
totals.hist <- subset(totals, histology_abbreviation %in% hists)

icgc <- toTable(0, 2, length(hists), c("MEDIAN", "N"))
rownames(icgc) <- hists
for (h in 1:length(hists)) {
	  samples.hist <- subset(totals.hist, histology_abbreviation == hists[h])
	  icgc$N[h] <- nrow(samples.hist)
	
	  sum <- c()
	  for (s in 1:nrow(samples.hist)) {
		    sample <- samples.hist$specimen_id[s]
		    load(paste0("/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/samples/rd-vs-rt_", sample, "-vs-lcl_spearman.RData"))
		
		    #sum <- c(sum, as.numeric(cor))
		    sum <- c(sum, as.numeric(cors$cor[13]))
	  }
	  icgc$MEDIAN[h] <- median(sum)
}
icgc <- icgc[order(icgc$MEDIAN, decreasing=T),]
writeTable(icgc, "/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/icgc_histology_abbreviation>20_cor13.txt", colnames=T, rownames=T, sep="\t")

###
##
donors <- c("project_code", "icgc_donor_id", "donor_sex", "donor_vital_status", "donor_age_at_diagnosis", "donor_survival_time")
samples <- toTable(0, 12, sum(icgc$N), c("CANCER", "COR", "Q4", "M2", "icgc_specimen_id", "histology_abbreviation", donors))
idx.cor    <- 0
idx.sample <- 0
for (h in 1:nrow(icgc)) {
	  samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
  	q4 <- setSamplesM2(wd.rt.data, samples.hist$specimen_id)
	
	  for (s in 1:nrow(samples.hist)) {
		    sample <- samples.hist$specimen_id[s]
		    load(paste0(wd.rt.data, "/samples/rd-vs-rt_", sample, "-vs-lcl_spearman.RData"))
		
		    samples$COR[idx.sample + s] <- as.numeric(cors$cor[13])
		    samples$icgc_specimen_id[idx.sample + s] <- sample
		    samples$histology_abbreviation[idx.sample + s] <- rownames(icgc)[h]
		
		    samples[idx.sample + s, c(donors)] <- clinicals[grep(paste0(sample, ","), clinicals$icgc_specimen_id), c(donors)]
	  }
	
	  samples$CANCER[(idx.sample+1):(idx.sample+nrow(samples.hist))] <- idx.cor
	  samples$Q4[(idx.sample+1):(idx.sample+nrow(samples.hist))] <- q4$Q4
	  samples$M2[(idx.sample+1):(idx.sample+nrow(samples.hist))] <- q4$M2
  	#samples$Q35[(idx.sample+1):(idx.sample+nrow(samples.hist))] <- q4$Q35
	
	  idx.cor <- idx.cor + 1
	  idx.sample <- idx.sample + nrow(samples.hist)
}
rownames(samples) <- samples$icgc_specimen_id
#samples.h <- samples
samples.v <- samples

samples.tmp <- samples.mut
samples.mut <- samples[intersect(rownames(samples), rownames(samples.mut)),]
samples.mut$tumor_wgs_aliquot_id <- samples.tmp[rownames(samples.mut),]$tumor_wgs_aliquot_id

#save(icgc, samples, samples.v, samples.mut, file=file.path(wd.rt.data, "icgc_wgs_samples_n2612_re-run_cor13.RData"))

# -----------------------------------------------------------------------------
# Stripchart (black; Resting to proliferating)
# Last Modified: 08/02/24; 22/11/22; 22/08/22; 18/08/22; 27/07/22; 21/04/19
# -----------------------------------------------------------------------------
samples.h1 <- toTable(0, 2, 0, c("CANCER", "COR"))
labels1 <- c()
idx.cancer <- 0
for (h in 26:1) {
	  q4 <- subset(samples, histology_abbreviation == rownames(icgc)[h])
	  q4$CANCER <- idx.cancer
	  samples.h1 <- rbind(samples.h1, q4)
	  labels1 <- c(labels1, rownames(icgc)[h])
	  idx.cancer <- idx.cancer + 1
}

###
##
file.name <- "stripchart_ICGC_SCP_n=2612_cor13"
main.text <- c(expression(bolditalic('In silico')~bold("sorting of PCAWG primary tumours (n=2,612)")), "")
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.55, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.9, cex.lab=2, cex.main=2.2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
stripchart(COR ~ CANCER, data=samples.h1, method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:26))

text(labels=labels1[1:26],   x=1:26,   y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.9)
mtext("SCF index", side=2, line=3.5, cex=2)
#legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()

###
##
file.name <- "stripchart_ICGC_SCP_n=2612_cor13_Lymph-BNHL"
main.text <- c(expression(bolditalic('In silico')~bold("sorting of PCAWG primary tumours (n=2,612)")), "")
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.55, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.9, cex.lab=2, cex.main=2.2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
stripchart(COR ~ CANCER, data=subset(samples.h1, histology_abbreviation != "Lymph-BNHL"), method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:9,11:26))
stripchart(COR ~ CANCER, data=subset(subset(samples.h1, histology_abbreviation == "Lymph-BNHL"), M2 == 2), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=10)
stripchart(COR ~ CANCER, data=subset(subset(samples.h1, histology_abbreviation == "Lymph-BNHL"), M2 == 1), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=10)

text(labels=labels1[1:9],   x=1:9,    y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)
text(labels=labels1[10],    x=10,     y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2, font=2)
text(labels=labels1[11:26], x=11:26,  y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.9)
mtext("SCF index", side=2, line=3.5, cex=2)
legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()

###
##
file.name <- "stripchart_ICGC_SCP_n=2612_cor13_Skin-Melanoma_Proliferative"
main.text <- c(expression(bolditalic('In silico')~bold("sorting of PCAWG primary tumours (n=2,612)")), "")
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.55, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.9, cex.lab=2, cex.main=2.2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
stripchart(COR ~ CANCER, data=subset(samples.h1, histology_abbreviation != "Skin-Melanoma"), method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:18,20:26))
stripchart(COR ~ CANCER, data=subset(subset(samples.h1, histology_abbreviation == "Skin-Melanoma"), M2 == 2), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=19)
stripchart(COR ~ CANCER, data=subset(subset(samples.h1, histology_abbreviation == "Skin-Melanoma"), M2 == 1), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=19)

text(labels=labels1[1:18],  x=1:18,   y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)
text(labels=labels1[19],    x=19,     y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2, font=2)
text(labels=labels1[20:26], x=20:26,  y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.9)
mtext("SCF index", side=2, line=3.5, cex=2)
legend("topleft", legend=c("Proliferative (S)", "Resting (G1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()

###
##
file.name <- "stripchart_ICGC_SCP_n=2612_cor13_ColoRect-AdenoCA"
main.text <- c(expression(bolditalic('In silico')~bold("sorting of PCAWG primary tumours (n=2,612)")), "")
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.55, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.9, cex.lab=2, cex.main=2.2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
stripchart(COR ~ CANCER, data=subset(samples.h1, histology_abbreviation != "ColoRect-AdenoCA"), method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:16,18:26))
stripchart(COR ~ CANCER, data=subset(subset(samples.h1, histology_abbreviation == "ColoRect-AdenoCA"), M2 == 2), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=17)
stripchart(COR ~ CANCER, data=subset(subset(samples.h1, histology_abbreviation == "ColoRect-AdenoCA"), M2 == 1), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=17)

text(labels=labels1[1:16],  x=1:16,   y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)
text(labels=labels1[17],    x=17,     y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2, font=2)
text(labels=labels1[18:26], x=18:26,  y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=2)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.9)
mtext("SCF index", side=2, line=3.5, cex=2)
legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()

# -----------------------------------------------------------------------------
# PanImmune
# Last Modified: 06/10/23; 16/03/23
# -----------------------------------------------------------------------------
immunes <- readTable(file.path("/Users/ty2/Work/uni-koeln/tyang2/ICGC/metadata/Thorsson 2018", "PanImmune.txt"), header=T, rownames=T, sep="\t")
nrow(immunes)
# [1] 11080

tcgas <- intersect(rownames(immunes), mappings$submitted_donor_id)
length(tcgas)
# [1] 860

mappings.tcga <- subset(mappings, submitted_donor_id %in% tcgas)
mappings.tcga.rna <- subset(mappings.tcga, specimen_library_strategy == "RNA-Seq")
mappings.tcga.wgs <- subset(mappings.tcga, specimen_library_strategy == "WGS")
nrow(mappings.tcga.rna)
# [1] 792
nrow(mappings.tcga.wgs)
# [1] 860
overlaps <- intersect(mappings.tcga.wgs$icgc_specimen_id, mappings.tcga.rna$icgc_specimen_id)
length(overlaps)
# [1] 792

rownames(mappings.tcga.rna) <- mappings.tcga.rna$submitted_specimen_id
rownames(mappings.tcga.wgs) <- mappings.tcga.wgs$submitted_specimen_id
overlaps.2 <- intersect(rownames(immunes), mappings.tcga.wgs$submitted_donor_id)
length(overlaps.2)
# [1] 860

immunes.tcga <- cbind(immunes[overlaps.2,], mappings.tcga.wgs[overlaps.2,])
rownames(immunes.tcga) <- immunes.tcga$icgc_specimen_id

##
load(file=file.path(wd.driver.data, "icgc_wgs_samples_n2542.RData"))

#overlaps.3 <- intersect(rownames(samples), rownames(immunes.tcga))   ## n=2612 (Total samples with WGSs)
#samples.mut.tcga <- cbind(samples[overlaps.3,], immunes.tcga[overlaps.3,])
#nrow(samples.mut.tcga)
# [1] 736

overlaps.3 <- intersect(rownames(samples.mut), rownames(immunes.tcga))   ## n=2542 (Total samples with SNVs)
samples.mut.tcga <- cbind(samples.mut[overlaps.3,], immunes.tcga[overlaps.3,])
nrow(samples.mut.tcga)
# [1] 710

# -----------------------------------------------------------------------------
# Using PanImmune (Global) to find chrs
# Last Modified: 27/04/24
# -----------------------------------------------------------------------------
insilico <- toTable(0, 4, 22, c("CHR", "Proliferation", "Wound.Healing", "Th17.Cells"))
insilico$CHR <- sprs.order$chr

for (chr in 1:nrow(sprs.order)) {
	  for (c in c(22:23,45)) {
		    col <- colnames(samples.mut.tcga)[c]
		
		    samples.mut.tcga.col <- samples.mut.tcga[!is.na(samples.mut.tcga[, c]),]
		    x=samples.mut.tcga.col[, c]
		    
		    for (s in 1:nrow(samples.mut.tcga.col)) {
		    	  sample <- rownames(samples.mut.tcga.col)[s]
		    	  load(file.path(wd.rt.data, "samples", paste0("rd-vs-rt_", sample, "-vs-lcl_spearman.RData")))
		    	
		    	  samples.mut.tcga.col$COR[s] <- cors$cor[chr]
		    }
		    y=samples.mut.tcga.col$COR
		    
		    cor <- cor.test(y, x, method="spearman", exact=F)
		    if (c == 45)
		    	  insilico[chr, 4] <- cor[[4]]
		    else
		       insilico[chr, c-20] <- cor[[4]]
	  }
}

##
ylim <- c(0.1, 0.35)  #, dns.u$Freq)))
file.name <- file.path(wd.rt.plots, paste0("DNS_ICGC_COR_Proliferation+Th17.Cells_n=710_lwd=4"))
#main.text <- c(expression(bold("NB-CL")~bolditalic('in silico')~bold("SCF estimation")), "")
main.text <- c("PCAWG SCF correlation", "")
ylab.text <- "Correlation to PanImmune"
xlab.text <- "Chromosome"
legends <- c("Proliferation", "Wound Healing", "Th17 Cells (Inverted)")
#cols <- c(red, blue, "dimgray")
cols <- c("gray", "dimgray", "black")

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)

points(insilico$Proliferation ~ rownames(insilico), col=cols[1], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$Proliferation, lty=5, lwd=2.5, col=cols[1])

points(insilico$Wound.Healing ~ rownames(insilico), col=cols[2], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$Wound.Healing, lty=5, lwd=2.5, col=cols[2])

points(insilico$Th17.Cells*-1 ~ rownames(insilico), col=cols[3], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$Th17.Cells*-1, lty=5, lwd=2.5, col=cols[3])

abline(v=13, lty=3, lwd=4, col=red)

axis(side=1, at=seq(1, 22, by=2), labels=insilico$CHR[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=insilico$CHR[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomright", legend=legends[3:1], col=cols[3:1], pch=15, lty=5, lwd=3, pt.cex=3, cex=1.9)
dev.off()

# -----------------------------------------------------------------------------
# samples
# Last Modified: 27/04/24 (re-run)
# -----------------------------------------------------------------------------
samples$COR <- NA
samples$Q4  <- NA
samples$M2  <- NA
for (s in 1:nrow(samples)) {
	  sample <- rownames(samples)[s]
	  load(paste0("/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/samples/rd-vs-rt_", sample, "-vs-lcl_spearman.RData"))
	
	  #samples$COR[s] <- as.numeric(cor)
	  samples$COR[s] <- as.numeric(cors$cor[13])
}

for (h in 1:nrow(icgc)) {
	  samples.hist <- subset(samples, histology_abbreviation == rownames(icgc)[h])
	  median <- median(samples.hist$COR)
	  samples[rownames(subset(samples.hist, COR >= median)),]$M2 <- 2
	  samples[rownames(subset(samples.hist, COR  < median)),]$M2 <- 1
}

samples.mut$COR <- NA
samples.mut$Q4  <- NA
samples.mut$M2  <- NA
for (s in 1:nrow(samples.mut)) {
	  sample <- rownames(samples.mut)[s]
	  load(paste0("/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/samples/rd-vs-rt_", sample, "-vs-lcl_spearman.RData"))
	
	  #samples.mut$COR[s] <- as.numeric(cor)
	  samples.mut$COR[s] <- as.numeric(cors$cor[13])
}

for (h in 1:nrow(icgc)) {
	  samples.hist <- subset(samples.mut, histology_abbreviation == rownames(icgc)[h])
	  median <- median(samples.hist$COR)
	  samples.mut[rownames(subset(samples.hist, COR >= median)),]$M2 <- 2
	  samples.mut[rownames(subset(samples.hist, COR  < median)),]$M2 <- 1
}

save(icgc, samples, samples.mut, file=file.path(wd.driver.data, "icgc_wgs_samples_n2612_re-run_cor13_samples.mut.RData"))

# -----------------------------------------------------------------------------
# PanImmune (Global)
# Last Modified: 06/10/23; 16/03/23
# -----------------------------------------------------------------------------
cols <- colnames(samples.mut.tcga)[18:77]
test <- toTable(0, 3, length(cols), c("PanImmune", "rho", "P"))
for (c in 18:77) {
	  col <- colnames(samples.mut.tcga)[c]
	  col2 <- paste0(unlist(strsplit(col, "\\.")), collapse=" ")
	  if (col == "Homologous.Recombination.Defects")
		    col2 <- "HR Defects"
	
	  samples.mut.tcga.col <- samples.mut.tcga[!is.na(samples.mut.tcga[, c]),]
	  x=samples.mut.tcga.col[, c]
	  y=samples.mut.tcga.col$COR
	  cor <- cor.test(y, x, method="spearman", exact=F)
	  test$PanImmune[c-17] <- col
	  test$rho[c-17] <- cor[[4]]
	  test$P[c-17]   <- cor[[3]]
	
  	file.name <- file.path(wd.rt.plots, "PanImmune_n=2542", paste0("Cor_SCF-vs-", col, ""))
	  main.text <- c("", "")
	  xlab.text <- col2
	  ylab.text <- "SCF index"
	  plotCorrelation(file.name, main.text, xlab.text, ylab.text, x=x, y=y, pos="topright", cols=c("darkgray", "black"), size=5)
}
test <- test[order(test$rho, decreasing=T),]
test$Q <- qvalue(test$P)$qvalue
writeTable(test, file.path(wd.rt.plots, "Cor_SCF-vs-PanImmune_n=2542_cor13.txt"), colnames=T, rownames=F, sep="\t")








# -----------------------------------------------------------------------------
# Overall correlation with LCL S/G1
# Last Modified: 27/04/24 (re-run); 20/11/23
# -----------------------------------------------------------------------------
cols <- colnames(samples.mut.tcga)[22:23]
insilico <- toTable(0, 3, 22, c("CHR", cols))
insilico$CHR <- sprs.order$chr
for (c in 1:nrow(sprs.order)) {
  	cors.all <- toTable(0, 4, nrow(samples.mut.tcga), c("SAMPLE_ID", "COR", "Q4", "M2"))
	  rownames(cors.all) <- rownames(samples.mut.tcga)
	  cors.all$SAMPLE_ID <- rownames(samples.mut.tcga)
	  for (s in 1:nrow(samples.mut.tcga)) {
		    sample <- rownames(samples.mut.tcga)[s]
		    load(paste0("/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/samples/rd-vs-rt_", sample, "-vs-lcl_spline_spearman.RData"))
		    
		    cors.all$COR[s] <- cors$cor[c]
	  }
	
  	cor <- cor.test(samples.mut.tcga[, 22], cors.all$COR, method="spearman", exact=F)
  	insilico$Proliferation[c] <- cor[[4]]
	
  	cor <- cor.test(samples.mut.tcga[, 23], cors.all$COR, method="spearman", exact=F)
  	insilico$Wound.Healing[c] <- cor[[4]]
}

##
ylim <- c(-0.4, 0.9)  #, dns.u$Freq)))
file.name <- file.path(wd.rt.plots, paste0("DNS_NB-CL_COR_-1_7_flowjo_2.5_red_FACS_22_3_3"))
#main.text <- c(expression(bold("NB-CL")~bolditalic('in silico')~bold("SCF estimation")), "")
main.text <- c("PCAWG SCF correlation", "")
ylab.text <- "Correlation to PanImmune"
xlab.text <- "Chromosome"
legends <- c("Wound Healing", "Proliferation")
#cols <- c(red, blue, "dimgray")
cols <- c(flowjo.grey, flowjo.red)

pdf(paste0(file.name, ".pdf"), height=4.5, width=9.8)
par(mar=c(5.1, 4.6, 4.1, 1.5))
plot(NULL, xlim=c(1, 22), ylim=ylim, xlab=xlab.text, ylab=ylab.text, main=main.text, yaxt="n", xaxt="n", pch=19, cex.axis=1.8, cex.lab=1.9, cex.main=2)

points(insilico$Wound.Healing ~ rownames(insilico), col=cols[1], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$Wound.Healing, lty=5, lwd=2.5, col=cols[1])

points(insilico$Proliferation ~ rownames(insilico), col=cols[2], pch=15, cex=2.5)
lines(rownames(insilico), y=insilico$Proliferation, lty=5, lwd=2.5, col=cols[2])

abline(v=11, lty=3, lwd=3, col=red)

axis(side=2, at=seq(-0.4, 0.8, by=0.2), labels=c(-0.4, "", 0, "", 0.4, "", 0.8), cex.axis=1.8)	  
axis(side=2, at=0.4, cex.axis=1.8)	
axis(side=2, at=0.8, cex.axis=1.8)	
axis(side=1, at=seq(1, 22, by=2), labels=insilico$CHR[seq(1, 22, by=2)], cex.axis=1.8)
axis(side=1, at=seq(2, 22, by=2), labels=insilico$CHR[seq(2, 22, by=2)], cex.axis=1.8)
legend("bottomright", legend=legends, col=cols, pch=15, lty=5, lwd=3, pt.cex=3, cex=1.9)
dev.off()





# -----------------------------------------------------------------------------
# samples
# Last Modified: 06/02/24 (re-run); 21/04/19
# -----------------------------------------------------------------------------
hists <- as.vector(subset(as.data.frame(sort(table(totals$histology_abbreviation))), Freq >= 20)$Var1)
totals.hist <- subset(totals, histology_abbreviation %in% hists)

icgc <- toTable(0, 2, length(hists), c("MEDIAN", "N"))
rownames(icgc) <- hists
for (h in 1:length(hists)) {
   samples.hist   <- subset(totals.hist, histology_abbreviation == hists[h])
   icgc$N[h] <- nrow(samples.hist)
   
   sum <- c()
   for (s in 1:nrow(samples.hist)) {
      sample <- samples.hist$specimen_id[s]
      load(paste0("/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/samples/rd-vs-rt_", sample, "-vs-lcl_spearman.RData"))
      
      sum <- c(sum, as.numeric(cors$cor[13]))
   }
   icgc$MEDIAN[h] <- median(sum)
}
icgc <- icgc[order(icgc$MEDIAN, decreasing=T),]
writeTable(icgc, "/Users/ty2/Work/uni-koeln/tyang2/ICGC/analysis/replication/icgc-wgs-rt/data/icgc_histology_abbreviation>20_cor13.txt", colnames=T, rownames=T, sep="\t")


# -----------------------------------------------------------------------------
# Stripchart (black; Resting to proliferating)
# Last Modified: 08/02/24; 22/11/22; 22/08/22; 18/08/22; 27/07/22; 21/04/19
# -----------------------------------------------------------------------------
samples.h1 <- toTable(0, 2, 0, c("CANCER", "COR"))
labels1 <- c()
idx.cancer <- 0
for (h in 26:1) {
   q4 <- subset(samples, histology_abbreviation == rownames(icgc)[h])
	  q4$CANCER <- idx.cancer
	  samples.h1 <- rbind(samples.h1, q4)
	  labels1 <- c(labels1, rownames(icgc)[h])
	  idx.cancer <- idx.cancer + 1
}

###
##
file.name <- "stripchart_ICGC_SCP_n=2612_cor22"
main.text <- c("", "")
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
stripchart(COR ~ CANCER, data=samples.h1, method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:26))
text(labels=labels1[1:26],   x=1:26,   y=par("usr")[3] - 0.05, srt=45, adj=0.965, xpd=NA, cex=1.9)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("SCF index", side=2, line=3.5, cex=1.9)
legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()

save(icgc, samples, samples.v, samples.mut, samples.h1, file=file.path(wd.rt.data, "icgc_wgs_samples_n2612_re-run_cor13.RData"))



















##
file.name <- "stripchart_ICGC_NBL_SCLC_colour_bold_NEW_medcol=red_medlwd=5_SCLC_NB_n57-1_plain_bold"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,769 primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
text(labels=labels1[1],     x=1,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[2:8],   x=2:8,   y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[9],     x=9,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[10:19], x=10:19, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[20],    x=20,    y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[21:28], x=21:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

#stripchart(COR ~ CANCER, data=samples.h1[1:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:28))
#stripchart(COR ~ CANCER, data=samples.h1[1:1889,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:19))
#stripchart(COR ~ CANCER, data=samples.h1[1991:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(21:28))
stripchart(COR ~ CANCER, data=samples.h1[1:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:28))
#stripchart(subset(samples.sclc, M2 == 2)$COR ~ rep(20, time=50), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=20)
#stripchart(subset(samples.sclc, M2 == 1)$COR ~ rep(20, time=51), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=20)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("S-phase cell fraction index", side=2, line=3.5, cex=1.9)
#legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()





















# -----------------------------------------------------------------------------
# PanImmune (Hist)
# Last Modified: 06/10/23; 16/03/23
# -----------------------------------------------------------------------------
for (h in 1:nrow(icgc)) {
	hist <- rownames(icgc)[h]
	samples.mut.tcga.hist <- subset(samples.mut.tcga$histology_abbreviation == )
}














###
## 16/03/23 PanImmune
immunes <- readTable(file.path("/Users/ty2/Work/uni-koeln/tyang2/ICGC/metadata/Thorsson 2018", "PanImmune.txt"), header=T, rownames=T, sep="\t")
nrow(immunes)
# [1] 11080

tcgas <- intersect(rownames(immunes), mappings$submitted_donor_id)
length(tcgas)
# [1] 814

mappings.tcga <- subset(mappings, submitted_donor_id %in% tcgas)
rownames(mappings.tcga) <- mappings.tcga$submitted_donor_id
immunes.tcga <- immunes[mappings.tcga$submitted_donor_id,]
for (h in 1:length(hists)) {
	mappings.tcga.hist <- subset(mappings.tcga, histology_abbreviation == hists[h])
	immunes.tcga.hist  <- immunes[rownames(mappings.tcga.hist),]
	
	##
	immunes.tcga.hist.nona <- immunes.tcga.hist[!is.na(immunes.tcga.hist$Leukocyte.Fraction),]
	samples.tcga <- samples[mappings.tcga.hist[rownames(immunes.tcga.hist.nona),]$icgc_specimen_id,]
	if (nrow(samples.tcga) > 1) {
		file.names <- file.path(wd.rt.plots, paste0("Leukocyte.Fraction_", hists[h], ".pdf"))
		plotCorrelation(file.names, paste0(hists[h], " (n=", nrow(samples.tcga), ")"), "Leukocyte fraction", "S-phase cell fraction index", immunes.tcga.hist.nona$Leukocyte.Fraction, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)
	}
	
	##
	immunes.tcga.hist.nona <- immunes.tcga.hist[!is.na(immunes.tcga.hist$Stromal.Fraction),]
	samples.tcga <- samples[mappings.tcga.hist[rownames(immunes.tcga.hist.nona),]$icgc_specimen_id,]
	if (nrow(samples.tcga) > 1) {
		file.names <- file.path(wd.rt.plots, paste0("Stromal.Fraction_", hists[h], ".pdf"))
		plotCorrelation(file.names, paste0(hists[h], " (n=", nrow(samples.tcga), ")"), "Stromal fraction", "S-phase cell fraction index", immunes.tcga.hist.nona$Stromal.Fraction, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)
	}
	
	##
	immunes.tcga.hist.nona <- immunes.tcga.hist[!is.na(immunes.tcga.hist$Proliferation),]
	samples.tcga <- samples[mappings.tcga.hist[rownames(immunes.tcga.hist.nona),]$icgc_specimen_id,]
	if (nrow(samples.tcga) > 1) {
		file.names <- file.path(wd.rt.plots, paste0("Proliferation_", hists[h], ".pdf"))
		plotCorrelation(file.names, paste0(hists[h], " (n=", nrow(samples.tcga), ")"), "Proliferation", "S-phase cell fraction index", immunes.tcga.hist.nona$Proliferation, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)
	}
}

##
for (c in 5:22) {
	col <- colnames(immunes.tcga)[c]
	col.name <- unlist(strsplit(col, "\\."))
	col.name <- "HR defects"
	if (length(col.name) != 1) {
		for (n in 2:length(col.name))
			col.name[n] <- tolower(col.name[n])
		col.name <- paste(col.name, collapse=" ")
	}
	
	immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga[, col]),]
	samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
	
	file.names <- file.path(wd.rt.plots, paste0("PanImmune_", col, ".pdf"))
	plotCorrelation5(file.names, "", paste0(col.name, " (n=", nrow(samples.tcga), ")"), "SCF index", immunes.tcga.nona[, col], samples.tcga$COR, pos="topright", cols=c("gray", "black"), size=5)
}

##
immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga$Leukocyte.Fraction),]
samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
file.names <- file.path(wd.rt.plots, paste0("Leukocyte.Fraction.pdf"))
plotCorrelation(file.names, paste0("CIBERSORT (n=", nrow(samples.tcga), ")"), "Leukocyte fraction", "S-phase cell fraction index", immunes.tcga.nona$Leukocyte.Fraction, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)

##
immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga$Stromal.Fraction),]
samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
file.names <- file.path(wd.rt.plots, paste0("Stromal.Fraction.pdf"))
plotCorrelation(file.names, paste0("CIBERSORT (n=", nrow(samples.tcga), ")"), "Stromal fraction", "S-phase cell fraction index", immunes.tcga.nona$Stromal.Fraction, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)

##
immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga$Proliferation),]
samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
file.names <- file.path(wd.rt.plots, paste0("Proliferation.pdf"))
plotCorrelation(file.names, paste0("CIBERSORT (n=", nrow(samples.tcga), ")"), "Proliferation", "S-phase cell fraction index", immunes.tcga.nona$Proliferation, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)

##
immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga$Wound.Healing),]
samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
file.names <- file.path(wd.rt.plots, paste0("Wound.Healing.pdf"))
plotCorrelation(file.names, paste0("CIBERSORT (n=", nrow(samples.tcga), ")"), "Wound healing", "S-phase cell fraction index", immunes.tcga.nona$Wound.Healing, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)

##
immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga$IFN.gamma.Response),]
samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
file.names <- file.path(wd.rt.plots, paste0("IFN.gamma.Response.pdf"))
plotCorrelation(file.names, paste0("CIBERSORT (n=", nrow(samples.tcga), ")"), "IFN gamma response", "S-phase cell fraction index", immunes.tcga.nona$IFN.gamma.Response, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)

##
immunes.tcga.nona <- immunes.tcga[!is.na(immunes.tcga$Lymphocyte.Infiltration.Signature.Score),]
samples.tcga <- samples[mappings.tcga[rownames(immunes.tcga.nona),]$icgc_specimen_id,]
file.names <- file.path(wd.rt.plots, paste0("Lymphocyte.Infiltration.Signature.Score.pdf"))
plotCorrelation(file.names, paste0("CIBERSORT (n=", nrow(samples.tcga), ")"), "Lymphocyte infiltration", "S-phase cell fraction index", immunes.tcga.nona$Lymphocyte.Infiltration.Signature.Score, samples.tcga$COR, pos="bottomright", cols=c("dimgray", "black"), size=6)





# -----------------------------------------------------------------------------
# Supplementary Table 1
# Last Modified: 14/04/23; 18/08/22; 27/07/22; 21/04/19
# -----------------------------------------------------------------------------
samples.h2 <- toTable(0, 3, 0, c("histology_abbreviation", "specimen_id", "SCF"))
for (h in 1:8) {
	  samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
	
	  q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
	  colnames(q4) <- c("specimen_id", "SCF")
	  q4$histology_abbreviation <- samples.hist$histology_abbreviation[1]
	  q4 <- q4[, c("histology_abbreviation", "specimen_id", "SCF")]
	  samples.h2 <- rbind(samples.h2, q4)
}

q4 <- samples.sclc[,1:2]
colnames(q4) <- c("specimen_id", "SCF")
q4$histology_abbreviation <- "Lung-SCLC"
q4 <- q4[, c("histology_abbreviation", "specimen_id", "SCF")]
samples.h2 <- rbind(samples.h2, q4)

for (h in 9:18) {
	  samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
	
	  q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
	  colnames(q4) <- c("specimen_id", "SCF")
	  q4$histology_abbreviation <- samples.hist$histology_abbreviation[1]
	  q4 <- q4[, c("histology_abbreviation", "specimen_id", "SCF")]
  	samples.h2 <- rbind(samples.h2, q4)
}

q4 <- samples.nbl[,1:2]
colnames(q4) <- c("specimen_id", "SCF")
q4$histology_abbreviation <- "Neuroblastoma"
q4 <- q4[, c("histology_abbreviation", "specimen_id", "SCF")]
samples.h2 <- rbind(samples.h2, q4)

for (h in 19:26) {
	  samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
	
	  q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
	  colnames(q4) <- c("specimen_id", "SCF")
	  q4$histology_abbreviation <- samples.hist$histology_abbreviation[1]
	  q4 <- q4[, c("histology_abbreviation", "specimen_id", "SCF")]
	  samples.h2 <- rbind(samples.h2, q4)
}

writeTable(samples.h2, file.path("/Users/ty2/Work/uni-koeln/tyang2/ICGC/metadata/", "Supplementary Table 1.txt"), colnames=T, rownames=F, sep="\t")





# -----------------------------------------------------------------------------
# Stripchart (black)
# Last Modified: 18/08/22; 27/07/22; 21/04/19
# -----------------------------------------------------------------------------
samples.h2 <- toTable(0, 3, 0, c("CANCER", "COR"))
labels2 <- c()
idx.cancer <- 0
for (h in 1:8) {
   samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
   
   q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
   colnames(q4) <- c("CANCER", "COR")
   q4$CANCER <- idx.cancer
   samples.h2 <- rbind(samples.h2, q4)
   labels2 <- c(labels2, rownames(icgc)[h])
   idx.cancer <- idx.cancer + 1
}

q4 <- samples.sclc[,1:2]
colnames(q4) <- c("CANCER", "COR")
q4$CANCER <- idx.cancer
samples.h2 <- rbind(samples.h2, q4)
labels2 <- c(labels2, "Lung-SCLC")
idx.cancer <- idx.cancer + 1

for (h in 9:18) {
   samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
   
   q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
   colnames(q4) <- c("CANCER", "COR")
   q4$CANCER <- idx.cancer
   samples.h2 <- rbind(samples.h2, q4)
   labels2 <- c(labels2, rownames(icgc)[h])
   idx.cancer <- idx.cancer + 1
}

q4 <- samples.nbl[,1:2]
colnames(q4) <- c("CANCER", "COR")
q4$CANCER <- idx.cancer
samples.h2 <- rbind(samples.h2, q4)
labels2 <- c(labels2, "Neuroblastoma")
idx.cancer <- idx.cancer + 1

for (h in 19:26) {
   samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
 
   q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
   colnames(q4) <- c("CANCER", "COR")
   q4$CANCER <- idx.cancer
   samples.h2 <- rbind(samples.h2, q4)
   labels2 <- c(labels2, rownames(icgc)[h])
   idx.cancer <- idx.cancer + 1
}

###
##
file.name <- "stripchart_ICGC_SCLC_NBL_colour_bold_n56_legend_violet"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,769 primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h2, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h$CANCER) + c(1.1, 2.8))
text(labels=labels2[1:19], x=1:19, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels2[20], x=20, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels2[21:28], x=21:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

stripchart(COR ~ CANCER, data=samples.h2[1:1832,], method="jitter", cex=1.5, pch=19, col=adjustcolor.gray, vertical=T, add=T, at=c(1:19))
stripchart(COR ~ CANCER, data=samples.h2[1889:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.gray, vertical=T, add=T, at=c(21:28))
stripchart(subset(samples.nbl, GROUP_ID3 == "LR")$COR ~ rep(20.2, time=17), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=20.2)
stripchart(subset(samples.nbl, GROUP_ID3 == "HR")$COR ~ rep(19.8, time=18), method="jitter", cex=1.5, pch=19, col="violet", vertical=T, add=T, at=19.8)
stripchart(subset(samples.nbl, GROUP_ID3 == "TERT")$COR ~ rep(20.066666, time=11), method="jitter", cex=1.5, pch=19, col=yellow, vertical=T, add=T, at=20.066666)
stripchart(subset(samples.nbl, GROUP_ID3 == "MYCN")$COR ~ rep(19.933333, time=10), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=19.933333)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Proliferation rate", side=2, line=3.5, cex=1.9)
legend("topright", legend=c("HR (n=18)", "MYCN (n=10)", "TERT (n=11)", "LR (n=17)"), pch=19, pt.cex=2.5, col=c("violet", red, yellow, blue), cex=1.7)
dev.off()

# -----------------------------------------------------------------------------
# Stripchart (black; Resting to proliferating)
# Last Modified: 22/11/22; 22/08/22; 18/08/22; 27/07/22; 21/04/19
# -----------------------------------------------------------------------------
samples.h1 <- toTable(0, 3, 0, c("CANCER", "COR"))
labels1 <- c()
idx.cancer <- 0
for (h in 26:19) {
   samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
 
   q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
   colnames(q4) <- c("CANCER", "COR")
   q4$CANCER <- idx.cancer
   samples.h1 <- rbind(samples.h1, q4)
   labels1 <- c(labels1, rownames(icgc)[h])
   idx.cancer <- idx.cancer + 1
}

q4 <- samples.nbl[,1:2]
colnames(q4) <- c("CANCER", "COR")
q4$CANCER <- idx.cancer
samples.h1 <- rbind(samples.h1, q4)
labels1 <- c(labels1, "Neuroblastoma")
idx.cancer <- idx.cancer + 1

for (h in 18:9) {
   samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
 
   q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
   colnames(q4) <- c("CANCER", "COR")
   q4$CANCER <- idx.cancer
   samples.h1 <- rbind(samples.h1, q4)
   labels1 <- c(labels1, rownames(icgc)[h])
   idx.cancer <- idx.cancer + 1
}

q4 <- samples.sclc[,1:2]
colnames(q4) <- c("CANCER", "COR")
q4$CANCER <- idx.cancer
samples.h1 <- rbind(samples.h1, q4)
labels1 <- c(labels1, "Lung-SCLC")
idx.cancer <- idx.cancer + 1

for (h in 8:1) {
   samples.hist <- subset(totals.hist, histology_abbreviation == rownames(icgc)[h])
 
   q4 <- setSamplesQ4(wd.rt.data, samples.hist$specimen_id)[,1:2]
   colnames(q4) <- c("CANCER", "COR")
   q4$CANCER <- idx.cancer
   samples.h1 <- rbind(samples.h1, q4)
   labels1 <- c(labels1, rownames(icgc)[h])
   idx.cancer <- idx.cancer + 1
}

###
##
file.name <- "stripchart_ICGC_NBL_SCLC_colour_bold_NEW_medcol=red_medlwd=5_SCLC_NB_n57-1_test"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,769 primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
text(labels=labels1[1],     x=1,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[2:8],   x=2:8,   y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[9],     x=9,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[10:19], x=10:19, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[20],    x=20,    y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[21:28], x=21:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

#stripchart(COR ~ CANCER, data=samples.h1[1:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:28))
stripchart(COR ~ CANCER, data=samples.h1[1:1889,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:19))
stripchart(COR ~ CANCER, data=samples.h1[1991:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(21:28))
stripchart(subset(samples.sclc, M2 == 2)$COR ~ rep(20, time=50), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=20)
stripchart(subset(samples.sclc, M2 == 1)$COR ~ rep(20, time=51), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=20)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("S-phase cell fraction index", side=2, line=3.5, cex=1.9)
legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()

##
file.name <- "stripchart_ICGC_NBL_SCLC_colour_bold_NEW_medcol=red_medlwd=5_SCLC_NB_n57-1_plain_bold"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,769 primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6), medcol=red, medlwd=5)
text(labels=labels1[1],     x=1,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[2:8],   x=2:8,   y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[9],     x=9,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[10:19], x=10:19, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[20],    x=20,    y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[21:28], x=21:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

#stripchart(COR ~ CANCER, data=samples.h1[1:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:28))
#stripchart(COR ~ CANCER, data=samples.h1[1:1889,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:19))
#stripchart(COR ~ CANCER, data=samples.h1[1991:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(21:28))
stripchart(COR ~ CANCER, data=samples.h1[1:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:28))
#stripchart(subset(samples.sclc, M2 == 2)$COR ~ rep(20, time=50), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=20)
#stripchart(subset(samples.sclc, M2 == 1)$COR ~ rep(20, time=51), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=20)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("S-phase cell fraction index", side=2, line=3.5, cex=1.9)
#legend("topleft", legend=c("Second median (M2)", "First median (M1)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.9)
dev.off()





###
##
file.name <- "stripchart_ICGC_NBL_SCLC_colour_bold_n56_legend_NEW_black"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 56 neuroblastoma (NB) primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6))
text(labels=labels1[1:8],   x=1:8,   y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[9],     x=9,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[10:28], x=10:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

stripchart(COR ~ CANCER, data=samples.h1[1:881,], method="jitter", cex=1.5, pch=19, col=adjustcolor.darkgray, vertical=T, add=T, at=c(1:8))
stripchart(COR ~ CANCER, data=samples.h1[938:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.darkgray, vertical=T, add=T, at=c(10:28))
stripchart(subset(samples.nbl, GROUP_ID3 == "LR")$COR ~ rep(8.7, time=17), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=8.7)
stripchart(subset(samples.nbl, GROUP_ID3 == "HR")$COR ~ rep(9.3, time=18), method="jitter", cex=1.5, pch=19, col="black", vertical=T, add=T, at=9.3)
stripchart(subset(samples.nbl, GROUP_ID3 == "TERT")$COR ~ rep(9.1, time=11), method="jitter", cex=1.5, pch=19, col=yellow, vertical=T, add=T, at=9.1)
stripchart(subset(samples.nbl, GROUP_ID3 == "MYCN")$COR ~ rep(8.9, time=10), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=8.9)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("S-phase cell fraction index", side=2, line=3.5, cex=1.9)
legend("topleft", legend=c("HR (n=18)", "TERT (n=11)", "MYCN (n=10)", "LR (n=17)"), pch=19, pt.cex=2.5, col=c("black", yellow, red, blue), cex=1.7)
dev.off()

###
##
file.name <- "stripchart_ICGC_NBL_SCLC_colour_bold_NEW_black_NB"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 56 neuroblastoma (NB) primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6))
text(labels=labels1[1:8],   x=1:8,   y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[9],     x=9,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[10:28], x=10:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

stripchart(COR ~ CANCER, data=samples.h1[1:881,], method="jitter", cex=1.5, pch=19, col=adjustcolor.darkgray, vertical=T, add=T, at=c(1:8))
stripchart(COR ~ CANCER, data=samples.h1[938:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.darkgray, vertical=T, add=T, at=c(10:28))
stripchart(subset(samples.nbl, RISK == "high")$COR ~ rep(9.1, time=39), method="jitter", cex=1.5, pch=19, col=red, vertical=T, add=T, at=9.1)
stripchart(subset(samples.nbl, RISK == "low")$COR ~ rep(8.9, time=17), method="jitter", cex=1.5, pch=19, col=blue, vertical=T, add=T, at=8.9)

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("S-phase cell fraction index", side=2, line=3.5, cex=1.9)
legend("topleft", legend=c("All high-risk (All-HR; n=39)", "Low-risk (LR; n=17)"), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.7)
dev.off()

###
##
file.name <- "stripchart_ICGC_NBL_SCLC_colour_bold_NEW_black"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,769 primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
#par(mar = c(11.2, 5, 4, 2))
par(mar=c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h1, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h1$CANCER) + c(1.4, 0.6))
text(labels=labels1[1],     x=1,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[2:8],   x=2:8,   y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[9],     x=9,     y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[10:19], x=10:19, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)
text(labels=labels1[20],    x=20,    y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9, font=2)
text(labels=labels1[21:28], x=21:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

stripchart(COR ~ CANCER, data=samples.h1[1:2769,], method="jitter", cex=1.5, pch=19, col=adjustcolor.black, vertical=T, add=T, at=c(1:28))

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("S-phase cell fraction index", side=2, line=3.5, cex=1.9)
#legend("topleft", legend=c("HR (n=18)", "TERT (n=11)", "MYCN (n=10)", "LR (n=17)"), pch=19, pt.cex=2.5, col=c("black", yellow, red, blue), cex=1.7)
dev.off()





###
##
file.name <- "stripchart_ICGC_SCLC_NBL_n57-1"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,769 primary tumour samples"))
adjustcolor.gray <- adjustcolor("black", alpha.f=0.25)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
par(mar = c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples.h2, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2, xlim=range(samples.h$CANCER) + c(1.1, 2.8))
text(labels=labels2, x=1:28, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

stripchart(COR ~ CANCER, data=samples.h2, method="jitter", cex=1.5, pch=19, col=adjustcolor.gray, vertical=T, add=T, at=c(1:28))

axis(side=2, at=seq(-0.8, 0.8, by=0.4), labels=c(-0.8, -0.4, 0, 0.4, 0.8), cex.axis=1.8)
mtext("Proliferation rate", side=2, line=3.5, cex=1.9)
dev.off()

save(samples.sclc, samples.nbl, samples.h2, labels2, file=file.path(wd.rt.data, paste0("icgc_wgs_samples_n2612+sclc+nbl.RData")), version=2)
#save(hists, totals.hist, icgc, samples.h, samples.v, samples, file=file.path(wd.rt.data, paste0("icgc_wgs_samples_n2612.RData")), version=2)






# -----------------------------------------------------------------------------
# Stripchart (sorting)
# Last Modified: 31/03/22
# -----------------------------------------------------------------------------
samples <- samples.h
samples$SORTING <- "G1"
idx <- which(samples$COR >= sample$COR)
if (length(idx) != 0)
   samples[idx,]$SORTING <- "S"
samples$SORTING <- as.factor(samples$SORTING)

save(hists, totals.hist, icgc, samples.h, samples.v, samples, file=file.path(wd.rt.data, paste0("icgc_wgs_samples_n2612.RData")), version=2)

##
file.name <- "stripchart_ICGC_samples_0.5_1.7"
main.text <- expression(bold(~bolditalic('In silico')~"sorting of 2,612 ICGC PCAWG samples"))
#labels <- paste0(labels=rownames(icgc), " (n=", icgc$N, ")")
labels <- rownames(icgc)
adjustcolor.red  <- adjustcolor(red, alpha.f=0.5)
adjustcolor.blue <- adjustcolor(blue, alpha.f=0.5)
cols <- c(adjustcolor.blue, adjustcolor.red)

pdf(file.path(wd.rt.plots, paste0(file.name, ".pdf")), height=7, width=20)
par(mar = c(11.2, 5, 4, 2))
boxplot(COR ~ CANCER, data=samples, yaxt="n", xaxt="n", ylab="", xlab="", main=main.text, col="white", outline=F, cex.axis=1.8, cex.lab=1.9, cex.main=2)
text(labels=labels, x=1:26, y=par("usr")[3] - 0.1, srt=45, adj=0.965, xpd=NA, cex=1.9)

abline(h=sample$COR, lty=5, lwd=5, col=green)

#stripchart(COR ~ CANCER, data=samples.h, method="jitter", cex=1.5, pch=19, col=adjustcolor.gray, vertical=T, add=T, at=c(1:26))
stripchart(COR ~ CANCER, data=subset(samples, SORTING == "G1"), method="jitter", cex=1.5, pch=19, col=cols[1], vertical=T, add=T, at=c(1:6, 8:26))
stripchart(COR ~ CANCER, data=subset(samples, SORTING == "S"),  method="jitter", cex=1.5, pch=19, col=cols[2], vertical=T, add=T, at=c(1:26))

axis(side=2, at=seq(-0.4, 0.8, by=0.4), labels=c(-0.4, 0, 0.4, 0.8), cex.axis=1.8)
axis(side=2, at=sample$COR, labels=round(sample$COR, 2), cex.axis=1.8, col.axis=green, font.axis=2)
mtext("Spearman's rho", side=2, line=3.5, cex=1.9)

legend("topright", legend=c(paste0("Proliferative (n=", separator(nrow(subset(samples, SORTING == "S"))), ")"), paste0("Non-proliferative (n=", separator(nrow(subset(samples, SORTING == "G1"))), ")")), pch=19, pt.cex=2.5, col=c(red, blue), cex=1.7)
dev.off()

# -----------------------------------------------------------------------------
# Purities
# Last Modified: 28/03/22
# ----------------------------------------------------------------------------
purities <- readTable("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/consensus/consensus_cnv/consensus.20170218.purity.ploidy.txt", header=T, rownames=F, sep="\t")
nrow(purities)
# [1] 2778

purities$icgc_specimen_id <- NA
for (s in 1:nrow(purities)) {
   ids <- unique(mappings[which(mappings$pcawg_wgs_id == purities$samplename[s]),]$icgc_specimen_id)
   if (length(ids) > 1) {
      purities$icgc_specimen_id[s] <- paste0(paste(unique(ids), collapse=","), ",")
   } else if (length(ids) == 1) {
      purities$icgc_specimen_id[s] <- ids
   }
}
purities <- purities[!is.na(purities$icgc_specimen_id),]
rownames(purities) <- purities$icgc_specimen_id
nrow(purities)
# [1] 2745

###
##
overlaps <- intersect(rownames(samples), rownames(purities))
length(overlaps)
# [1] 2612

x <- purities[overlaps,]$purity
y <- samples[overlaps,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612"))
plotCorrelation(file.name, "PCAWG (n=2,612)", "Purity", txt.In.silico, x, y, "topleft")

x <- purities[overlaps,]$ploidy
y <- samples[overlaps,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "ploidy", "_n2612"))
plotCorrelation(file.name, "PCAWG (n=2,612)", "Ploidy", txt.In.silico, x, y, "topright")

##
quantile(samples$COR)
# 0%        25%        50%        75%       100% 
# -0.7861177 -0.7335852 -0.7081717 -0.4922725  0.7518429 
(-0.7335852 + -0.7081717)/2
# [1] -0.7208785
sample$COR
# [1] -0.7233351

##
s  <- rownames(subset(samples, SORTING == "S"))
g1 <- rownames(subset(samples, SORTING == "G1"))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_S_n1629"))
plotCorrelation(file.name, "Proliferative (n=1,629)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_G1_n983"))
plotCorrelation(file.name, "Non-proliferative (n=983)", "Purity", txt.In.silico, x, y, "topleft", col=blue)

##
s  <- rownames(subset(samples, M2 == 2))
g1 <- rownames(subset(samples, M2 == 1))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_M2_n1297"))
plotCorrelation(file.name, "M2 (n=1,297)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_M1_n1315"))
plotCorrelation(file.name, "M1 (n=1,315)", "Purity", txt.In.silico, x, y, "topleft", col=blue)

##
m  <- quantile(samples$COR)[3]
s  <- rownames(subset(samples, COR >= m))
g1 <- rownames(subset(samples, COR <  m))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_global_M2_n1306"))
plotCorrelation(file.name, "Global M2 (n=1,306)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_global_M1_n1306"))
plotCorrelation(file.name, "Global M1 (n=1,306)", "Purity", txt.In.silico, x, y, "topleft", col=blue)

##
s  <- rownames(subset(samples.surv, SORTING == "S"))
g1 <- rownames(subset(samples.surv, SORTING == "G1"))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n1580_S_n909"))
plotCorrelation(file.name, "Proliferative (n=909)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n1580_G1_n671"))
plotCorrelation(file.name, "Non-proliferative (n=671)", "Purity", txt.In.silico, x, y, "topleft", col=blue)

##
s  <- rownames(subset(samples.surv.hist, SORTING == "S"))
g1 <- rownames(subset(samples.surv.hist, SORTING == "G1"))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n1580_hists_S_n612"))
plotCorrelation(file.name, "Proliferative (n=612)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n1580_hists_G1_n968"))
plotCorrelation(file.name, "Non-proliferative (n=968)", "Purity", txt.In.silico, x, y, "topleft", col=blue)

(-0.7081717-0.4922725)/2
##
q3.5  <- (quantile(samples$COR)[3] + quantile(samples$COR)[4])/2
s  <- rownames(subset(samples, COR >= q3.5))
g1 <- rownames(subset(samples, COR <  q3.5))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_global_q3.5_M2_846"))
plotCorrelation(file.name, "Global Q3.5 ~ Q4 (n=846)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_global_q3.5_M1_n1766"))
plotCorrelation(file.name, "Global Q1 ~ Q3.5 (n=1,766)", "Purity", txt.In.silico, x, y, "topleft", col=blue)

##
s  <- rownames(subset(samples, Q35 == 3.5))
g1 <- rownames(subset(samples, Q35 == 1))

x <- purities[s,]$purity
y <- samples[s,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_q3.5_M2_936"))
plotCorrelation(file.name, "Q3.5 ~ Q4 (n=936)", "Purity", txt.In.silico, x, y, "topleft", col=red)

x <- purities[g1,]$purity
y <- samples[g1,]$COR
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "purity", "_n2612_q3.5_M1_n1676"))
plotCorrelation(file.name, "Q1 ~ Q3.5 (n=1,676)", "Purity", txt.In.silico, x, y, "topleft", col=blue)




###
##
idx.na <- which(is.na(samples$donor_age_at_diagnosis))
x <- samples$donor_age_at_diagnosis[-idx.na]
y <- samples$COR[-idx.na]
file.name <- file.path(wd.rt.plots, paste0("correlation_in-silico_vs_", "age", "_n2548_5*5"))
plotCorrelation(file.name, "PCAWG (n=2,548)", "Age at diagnosis", txt.In.silico, x, y, "topleft")

##
file.name <- file.path(wd.rt.plots, "boxplot_in-silico_F-vs-M")
plotBox(file.name, subset(samples, donor_sex == "male")$COR, subset(samples, donor_sex == "female")$COR, "PCAWG (n=2,612)", names=c("Male", "Female"), cols=c("white", "white"))





































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

# -----------------------------------------------------------------------------
# Pan-cancer survival analysis
# Last Modified: 16/03/22
# -----------------------------------------------------------------------------
samples.surv <- survICGC(samples)

#icgc$HIST_P_OS    <- NA
#icgc$HIST_P_HAZARD <- NA
for (h in 1:nrow(icgc)) {
   hist <- rownames(icgc)[h]
   samples.hist <- subset(samples.surv, histology_abbreviation == hist)
   samples.hist <- samples.hist[order(samples.hist$COR),]
   
   if (nrow(samples.hist) > 20) {
      pvals <- c()
      rhos  <- c()
      
      for (s in 1:nrow(samples.hist)) {
         samples.hist$SG1 <- "G1"
         idx <- which(samples.hist$COR >= samples.hist$COR[s])
         if (length(idx) != 0)
            samples.hist[idx,]$SG1 <- "S"
         samples.hist$SG1 <- as.factor(samples.hist$SG1)
       
         ##
         fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.hist)
         if (!is.na(surv_pvalue(fit)$pval)) {
            pvals <- c(pvals, surv_pvalue(fit)$pval)
            rhos <- c(rhos, samples.hist$COR[s])
            #file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_in-silico_", h, "_", hist))
            #main.text <- c("PCAWG", hist)   ##, expression(italic('in silico')~"sorting"))
            #plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))
            #icgc$COR_P_OS[h] <- surv_pvalue(fit)$pval
          
            #res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR + donor_sex + donor_age_at_diagnosis, data=samples.hist.surv)
            #res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR, data=samples.hist.surv)
            #pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/hazard_", hist, ".pdf"), height=3, width=5)
            #ggforest(res.cox, data=samples.hist.surv, main=paste0("Hazard ratio in ", hist), cpositions = c(0.02, 0.22, 0.4), fontsize=0.7)
            #dev.off()
         }
      }
      
      if (-log10(min(pvals)) > 3) {
         s <- which(pvals == min(pvals))
         samples.hist$SG1 <- "G1"
         idx <- which(samples.hist$COR > samples.hist$COR[s])
         if (length(idx) != 0)
            samples.hist[idx,]$SG1 <- "S"
         samples.hist$SG1 <- as.factor(samples.hist$SG1)
         
         if (as.numeric(summary(samples.hist$SG1)[1]) >= 10 && as.numeric(summary(samples.hist$SG1)[2] >= 10)) {
            x <- rhos
            y <- -log10(pvals)
            file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/correlation_in-silico_P_KM_OS_", hist)
            plotJust(file.name, hist, expression(italic('in silico')~"sorting"), "-log10(p-value)", x, y, line=2.7)
          
            fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.hist)
            file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_in-silico_", hist))
            main.text <- c(hist, paste0("Divided by rho = ", round0(rhos[s], 2)))   ##, expression(italic('in silico')~"sorting"))
            plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))
            
            res.cox <- coxph(Surv(OS_month, OS_censor) ~ SG1 + COR + SEX + AGE, data=samples.hist)
            pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/hazard_in-silico_", hist, ".pdf"), height=2.5, width=5)
            ggforest(res.cox, data=samples.hist, main=hist, cpositions = c(0.02, 0.22, 0.4), fontsize=2)
            dev.off()
         }
      }
   }
}

samples.surv <- survICGC(samples)
samples.surv$COR_P_OS     <- NA
samples.surv$COR_P_HAZARD <- NA
for (s in 1:nrow(samples.surv)) {
   sample <- samples.surv[s,]
   
   samples.surv$SG1 <- "G1"
   idx <- which(samples.surv$COR >= samples.surv$COR[s])
   if (length(idx) != 0)
      samples.surv[idx,]$SG1 <- "S"
   samples.surv$SG1 <- as.factor(samples.surv$SG1)
 
   ##
   fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv)
   #file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_in-silico_", s, "_", sample$icgc_specimen_id, "_80_0.07_2_3_1.6"))
   #main.text <- c("Overall survival (OS)", paste0("Divided by rho = ", round0(sample$COR, 2)))   ##, expression(italic('in silico')~"sorting"))
   #plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))
   samples.surv$COR_P_OS[s] <- surv_pvalue(fit)$pval
 
   #res.cox <- coxph(Surv(OS_month, OS_censor) ~ SG1 + COR + SEX + AGE, data=samples.surv)
   #pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/hazard_", sample$icgc_specimen_id, "_2.5.pdf"), height=2.5, width=5)
   #ggforest(res.cox, data=samples.surv, main=paste0("Hazard ratio divided by rho = ", round0(sample$COR, 2)), cpositions = c(0.02, 0.22, 0.4), fontsize=2)
   #dev.off()
}
samples.surv <- samples.surv[!is.na(samples.surv$COR_P_OS),]
spot <- which(samples.surv$COR_P_OS == min(samples.surv$COR_P_OS))

x <- samples.surv$COR
y <- -log10(samples.surv$COR_P_OS)
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/correlation_in-silico_P_KM_OS_1679_P_bold")
#plotJust(file.name, "Kaplan-Meier OS", expression(italic('in silico')~"sorting"), "-log10(p-value)", x, y, line=2.7, ymax=19.5)
plotJust(file.name, "OS of 1,679 patients", "Spearman's rho", expression('-log10('*italic('P')*')'), x, y, line=2.53)

##
test <- toTable(0, 2, 2, c("G1", "S"))
rownames(test) <- c("M", "F")
test[1, 1] <- nrow(subset(subset(samples.surv, donor_sex == "male"), SG1 == "G1"))
test[1, 2] <- nrow(subset(subset(samples.surv, donor_sex == "male"), SG1 == "S"))
test[2, 1] <- nrow(subset(subset(samples.surv, donor_sex == "female"), SG1 == "G1"))
test[2, 2] <- nrow(subset(subset(samples.surv, donor_sex == "female"), SG1 == "S"))
fisher.test(test)[[1]]
# [1] 0.5063105

file.name <- paste0("boxplot_", sample$icgc_specimen_id, "_M1-vs-M2_Age")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/", file.name, subset(samples.surv, SG1 == "G1")$donor_age_at_diagnosis, subset(samples.surv, SG1 == "S")$donor_age_at_diagnosis, hist, names=c("G1-like", "S-like"), cols=c(blue, red))








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
for (h in 1:nrow(icgc)) {
   hist <- rownames(icgc)[h]
   samples.hist <- subset(samples, histology_abbreviation == hist)
   samples.hist <- survICGC(samples.hist)
   
   if (nrow(samples.hist) > 0) {
      x <- samples.hist[, colname]
      y <- samples.hist$COR
      file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/correlation_in-silico_vs_", colname, "_", hist)
      #plotCorrelation(file.name, hist, "Survival time", expression(italic('in silico')~"sorting"), x, y, "topright", line=2.4)
  
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
#library(survival)
#library(survminer)
#library(broom)
library(ggalt)
library(grid)
library(dplyr)

for (h in 1:nrow(icgc)) {
   hist <- rownames(icgc)[h]
   samples.hist <- subset(samples, histology_abbreviation == hist)
   samples.hist.surv <- survICGC(samples.hist)
   
   if (nrow(samples.hist.surv) > 10) {
       ## Cox regression model
     
       #print(hist)
       #res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR + donor_sex + donor_age_at_diagnosis, data=samples.hist.surv)
       #res.cox <- coxph(Surv(OS_month, OS_censor) ~ COR, data=samples.hist.surv)
       #print(res.cox)
       #print("----------------------------------------------------------------")
       #pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/hazard_", hist, ".pdf"), height=3, width=5)
       #ggforest(res.cox, data=samples.hist.surv, main=paste0("Hazard ratio in ", hist), cpositions = c(0.02, 0.22, 0.4), fontsize=0.7)
       #dev.off()
       #ggforest(res.cox)
     
       ##
       #fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.hist.surv)
       #file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/M2/survfit_", hist, "_in-silico"))
       #main.text <- c(hist, expression(italic('in silico')~"sorting"))
       #plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))
     
       ##
       fit <- survfit(Surv(OS_month, OS_censor) ~ M2, data=samples.hist.surv)
       file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_", hist, "_M2"))
       main.text <- c(hist, expression(italic('in silico')~"sample sorting"))
       plotSurvfit(fit, file.name, main.text, c("M1", "M2"), c(blue, red))
     
       ##
       #fit <- survfit(Surv(OS_month, OS_censor) ~ Q4, data=samples.hist.surv)
       #file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/Q4/survfit_", hist, "_Q4"))
       #main.text <- c(hist, expression(bolditalic('in silico')~"sorting"))
       #plotSurvfit(fit, file.name, main.text, c("Q1", "Q2", "Q3", "Q4"), c(blue, blue.lighter, red.lighter, red))
     
       ##
       #fit <- survfit(Surv(OS_month, OS_censor) ~ Sex, data=samples.hist.surv)
       #file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_", hist, "_SEX"))
       #main.text <- c(hist, expression(bolditalic('in silico')~"sorting"))
       #plotSurvfit(fit, file.name, main.text, c("M", "F"), c(blue.lighter, red.lighter))
   }
}

###
## Wow
samples.surv <- survICGC(samples)
res.cox <- coxph(Surv(OS_month, OS_censor) ~ SG1 + COR + SEX + AGE, data=samples.surv)

#samples.surv <- samples.surv.neg
colname <- "donor_survival_time"
x <- samples.surv[, colname]
y <- samples.surv$COR
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/correlation_in-silico_vs_", colname, "_", "ALL")
plotCorrelation(file.name, "1,683 PCAWG samples", "Survival time", expression(italic('in silico')~"sorting"), x, y, "topright", line=2.4)

colname <- "donor_age_at_diagnosis"
x <- samples.surv[, colname]
y <- samples.surv$COR
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/correlation_in-silico_vs_", colname, "_", "ALL_flip")
plotCorrelation(file.name, "1,683 PCAWG samples", expression(italic('in silico')~"sorting"), "Age at diagnosis", y, x, "topright", line=2.4)

x <- samples.surv.pos[, colname]
y <- samples.surv.pos$COR
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/correlation_in-silico_vs_", colname, "_", "ALL_RHO>0_flip")
plotCorrelation(file.name, "879 PCAWG samples (RHO > 0)", expression(italic('in silico')~"sorting"), "Age at diagnosis", y, x, "topright", line=2.4)

x <- samples.surv.neg[, colname]
y <- samples.surv.neg$COR
file.name <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/correlation_in-silico_vs_", colname, "_", "ALL_RHO<0_flip")
plotCorrelation(file.name, "764 PCAWG samples (RHO < 0)",  expression(italic('in silico')~"sorting"), "Age at diagnosis", y, x, "topright", line=2.4)

fit <- survfit(Surv(OS_month, OS_censor) ~ donor_sex, data=samples.surv)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_ALL_Sex"))
main.text <- c("Sex", "")
plotSurvfit(fit, file.name, main.text, c("F", "M"), c(red.lighter, blue.lighter))

fit <- survfit(Surv(OS_month, OS_censor) ~ donor_sex, data=samples.surv.pos)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_ALL_Sex_pos"))
main.text <- c("Sex", "Slope > 0")
plotSurvfit(fit, file.name, main.text, c("F", "M"), c(red.lighter, blue.lighter))

fit <- survfit(Surv(OS_month, OS_censor) ~ donor_sex, data=samples.surv.neg)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_ALL_Sex_neg"))
main.text <- c("Sex", "Slope < 0")
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
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_", "ALL", "_S-vs-G1"))
main.text <- c("Cell cycle statue", "")
plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

icgc.pos <- subset(icgc, survival_rho > 0)
icgc.neg <- subset(icgc, survival_rho < 0)
samples.surv.pos <- subset(samples.surv, histology_abbreviation %in% rownames(icgc.pos))
samples.surv.neg <- subset(samples.surv, histology_abbreviation %in% rownames(icgc.neg))

##
res.cox <- coxph(Surv(OS_month, OS_censor) ~ In_silico + G1_vs_S + Sex + Age, data=samples.surv)
fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv.pos)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_", "ALL", "_S-vs-G1_Slope>0"))
main.text <- c("Cell cycle statue", "Slope > 0")
plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/hazard_", "ALL", ".pdf"), height=3, width=5)
ggforest(res.cox, data=samples.surv, main=paste0("Hazard ratio in ", "ALL"), cpositions = c(0.02, 0.22, 0.4), fontsize=2)
dev.off()
#ggforest(res.cox)

##
#res.cox <- coxph(Surv(OS_month, OS_censor) ~ In_silico + G1_vs_S + Sex + Age, data=samples.surv.pos)
#fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv.pos)
#file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_", "ALL", "_S-vs-G1_Slope>0"))
#main.text <- c("Cell cycle statue", "Slope > 0")
#plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/hazard_", "ALL_RHO>0", "_label.pdf"), height=3, width=5)
ggforest(res.cox, data=samples.surv.pos, main=paste0("Hazard ratio in ", "RHO > 0"), cpositions = c(0.02, 0.22, 0.4), fontsize=2)
dev.off()
#ggforest(res.cox)

##
#res.cox <- coxph(Surv(OS_month, OS_censor) ~ In_silico + G1_vs_S + Sex + Age, data=samples.surv.neg)
fit <- survfit(Surv(OS_month, OS_censor) ~ SG1, data=samples.surv.neg)
file.name <- file.path(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/plots/survfit_", "ALL", "_S-vs-G1_Slope<0"))
main.text <- c("Cell cycle statue", "Slope < 0")
plotSurvfit(fit, file.name, main.text, c("G1-like", "S-like"), c(blue, red))

pdf(paste0("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/hazard_", "ALL_RHO<0", "_label.pdf"), height=3, width=5)
ggforest(res.cox, data=samples.surv.neg, main=paste0("Hazard ratio in ", "RHO < 0"), cpositions = c(0.02, 0.22, 0.4), fontsize=2)
dev.off()

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
test[1, 1] <- nrow(subset(subset(samples.hist.surv, donor_sex == "male"), SG1 == "G1"))
test[1, 2] <- nrow(subset(subset(samples.hist.surv, donor_sex == "male"), SG1 == "S"))
test[2, 1] <- nrow(subset(subset(samples.hist.surv, donor_sex == "female"), SG1 == "G1"))
test[2, 2] <- nrow(subset(subset(samples.hist.surv, donor_sex == "female"), SG1 == "S"))
# > fisher.test(test)[[1]]
# [1] 0.4319473

## Wow
file.name <- paste0("boxplot_", hist, "_F-vs-M_COR")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/", file.name, subset(samples.hist.surv, donor_sex == "male")$COR, subset(samples.hist.surv, donor_sex == "female")$COR, hist, names=c("Male", "Female"), cols=c(blue.lighter, red.lighter))

file.name <- paste0("boxplot_", hist, "_F-vs-M_Age")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/", file.name, subset(samples.hist.surv, donor_sex == "male")$donor_age_at_diagnosis, subset(samples.hist.surv, donor_sex == "female")$donor_age_at_diagnosis, hist, names=c("Male", "Female"), cols=c(blue.lighter, red.lighter))

file.name <- paste0("boxplot_", hist, "_M1-vs-M2_Age")
plotBox02("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/ICGC/", file.name, subset(samples.hist.surv, M2 == 1)$donor_age_at_diagnosis, subset(samples.hist.surv, M2 == 2)$donor_age_at_diagnosis, hist, names=c("M1", "M2"), cols=c(blue, red))





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
