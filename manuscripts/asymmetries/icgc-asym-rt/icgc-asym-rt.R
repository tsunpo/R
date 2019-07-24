# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/asymmetries/cll-asym-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 20/06/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Asymmetry.R", "Mutation.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.1kb.gc.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 10/07/18
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.icgc  <- file.path(wd, BASE, "consensus")
wd.anlys <- file.path(wd, BASE, "analysis")

wd.icgc.vcf <- file.path(wd.icgc, "snv_mnv", "Converted_Data")

release <- readTable(file.path(wd.icgc, "release_may2016.v1.4.tsv"), header=T, rownames=F, sep="")
list <- readTable(file.path(wd.icgc, "Sclust_Mar_2017_S2703.txt"), header=T, rownames=T, sep="")
samples <- readTable(file.path(wd.icgc, "Sclust_Mar_2017_final_summary.txt.gz"), header=T, rownames=F, sep="")
rownames(samples) <- samples$samplename
# > dim(samples)
# [1] 2703   8

icgcs <- as.data.frame(table(samples$cancer_type))
icgcs <- icgcs[order(icgcs$Freq, decreasing=T),]
icgcs$icgc <- tolower(icgcs$Var1)
writeTable(icgcs, file.path(wd.icgc, "icgc_t40_n2703.txt"), colnames=F, rownames=F, sep="\t")
# > dim(icgcs)
# [1] 40  2

samples2 <- samples[0,]
icgcs2   <- icgcs[0,]
for (t in 1:nrow(icgcs)) {
   BASE <- as.vector(icgcs[t,]$Var1)
   samples.code <- subset(samples, cancer_type == BASE)
   
   bases <- release[grep(BASE, release$dcc_project_code),]
   bases <- subset(bases, tumor_wgs_aliquot_id %in% samples.code$samplename)
   
   codes <- unique(bases$dcc_project_code)
   for (b in 1:length(codes)) {
      code <- codes[b]
      bases.code <- subset(bases, dcc_project_code == code)
      
   }
}
writeTable(icgcs2, file.path(wd.icgc, "icgc_t48_n2703.txt"), colnames=F, rownames=F, sep="\t")
readTable(samples2, file.path(wd.icgc, "Sclust_Mar_2017_final_summary_t48.txt.gz"), colnames=T, rownames=F, sep="")

# for (s in 1:nrow(list)) {
#    system(paste0("gzip ", file.path(wd.icgc.vcf, paste0(rownames(list[s,]), ".consensus.20160830.somatic.snv_mnv.vcf"))))
# }

# -----------------------------------------------------------------------------
# Step 1: Load all mutations
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
for (t in 1:nrow(icgcs)) {
   BASE <- as.vector(icgcs[t,]$Var1)
   base <- tolower(BASE)
   
   system(paste0("mkdir -p ",  file.path(wd.anlys, "asymmetries", paste0(base, "-asym-rt"), "data", "logs")))
   system(paste0("mkdir -p ",  file.path(wd.anlys, "asymmetries", paste0(base, "-asym-rt"), "plots")))
   wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-rt"))
   wd.asym.data  <- file.path(wd.asym, "data")
   wd.asym.plots <- file.path(wd.asym, "plots")
   
   samples.code <- subset(samples, cancer_type == BASE)
   colnames <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO", "SAMPLE")
   snvs <- toTable("", length(colnames), 0, colnames)
   for (s in 1:nrow(samples.code)) {
      sample <- rownames(samples.code[s,])
      vcf <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.icgc.vcf, paste0(sample, ".consensus.20160830.somatic.snv_mnv.vcf.gz")), pass=T, rs=F)
      vcf$SAMPLE <- sample
   
      ## Seperate SNVs
      vcf.snv <- subset(vcf,     REF %in% c("A", "T", "C", "G"))
      vcf.snv <- subset(vcf.snv, ALT %in% c("A", "T", "C", "G"))   
   
      #writeTable(vcf.snv, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.snv.gz"))), colnames=T, rownames=F, sep="\t")
      snvs <- rbind(snvs, vcf.snv)
   }
   colnames <- c("SAMPLE", "CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")
   snvs <- snvs[, colnames]
   save(snvs, file=file.path(wd.asym.data, paste0(base, "_mut_snvs.RData")))
}

# -----------------------------------------------------------------------------
# Step 2: Build S6 table
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
# > nrow(snvs)
# [1] 5268240

isWithinWindow <- function(bed.gc.chr.rt, POS) {
   bed.gc.chr.rt.start <- subset(bed.gc.chr.rt, START <= POS)
   bed.gc.chr.rt.start.end <- subset(bed.gc.chr.rt.start, END >= POS)
   
   if (nrow(bed.gc.chr.rt.start.end) == 1)
      return(rownames(bed.gc.chr.rt.start.end))
   else
      return("")
}

muts <- toTable(0, 9, 22, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
muts$CHR <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")
   rpkms.chr.rt <- cbind(rpkms.chr.rt[rownames(rpkms.chr.rt.RT),], rpkms.chr.rt.RT[, "SPLINE"])
   colnames(rpkms.chr.rt) <- c("BED", "T", "N", "RT", "SPLINE")
   
   bed.gc.chr.rt <- bed.gc.chr[rownames(rpkms.chr.rt),]
   snvs.chr <- subset(snvs, CHROM == chr)
   muts$LENGTH[c] <- nrow(rpkms.chr.rt) * 1000/1000000
   
   ## Keep only SNVs on overlapping 1kb windows
   snvs.chr$BED <- mapply(x = 1:nrow(snvs.chr), function(x) isWithinWindow(bed.gc.chr.rt, snvs.chr[x,]$POS))
   
   muts$TOTAL[c]  <- nrow(subset(snvs.chr, BED != ""))
   snvs.chr <- subset(snvs.chr, BED != "")
   snvs.chr <- cbind(snvs.chr, rpkms.chr.rt[snvs.chr$BED, 2:5])
   save(snvs.chr, file=file.path(wd.asym.data, paste0(base, "_mut_snvs_chr", c,".RData")))
   
   snvs.chr.s6 <- getTableS6SNV(snvs.chr[, c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL")])
   for (s in 1:6) {
      muts[c, 3+s] <- nrow(snvs.chr.s6[[s]][[1]]) + nrow(snvs.chr.s6[[s]][[2]])
   }
}
save(muts, file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6.RData")))

###
## After copying files back from cheops
muts <- toTable(0, 9, 0, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
for (c in 1:22) {
   load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_chr", c,".RData")))
   muts <- rbind(muts, muts.chr)
}
save(muts, file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6.RData")))

# -----------------------------------------------------------------------------
# SCLC RD vs RD
# Last Modified: 07/07/19
# -----------------------------------------------------------------------------
load(file=file.path(wd.rt.data, paste0("rds-vs-rt_", base, "-m2-m1_spline_spearman.RData")))   ## Load skews
load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6.RData")))
#load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_early.RData")))
#load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_late.RData")))
time <- "EARLY"

s6 <- c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C")
for (s in 1:6) {
   cors <- toTable(0, 2, 22, c("chr", "cor"))
   cors$chr <- 1:22
   idx <- idxs[s]
   
   for (c in 1:22)
      cors$cor[c] <- muts[c, 3+s] / muts[c, "LENGTH"]

   file.name <- file.path(wd.asym.plots, paste0(time, "_RDS-vs-MUT-", gsub("/", "_", s6[s]), "_", BASE))
   main.text <- c(paste0(s6[s], ""), BASE)
   xlab.text <- "SNVs/Mb"
   plotSPRRDC(cors, skews, file.name, main.text, c(4, 13, 17, 19, 21, 22), xlab.text, unit=4.5)
   
   #file.name <- file.path(wd.asym.plots, paste0(time, "_LENGTH-vs-MUT-", gsub("/", "_", s6[s]), "_", BASE))
   #main.text <- c(paste0(s6[s], ""), "")
   #plotMUTLength(cors, muts, file.name, main.text, c(1, 4, 5, 13, 17, 19, 22), "SNVs/Mb")
   
   #file.name <- file.path(wd.asym.plots, paste0(time, "_LENGTH-vs-SNV-", gsub("/", "_", s6[s]), "_", BASE))
   #main.text <- c(paste0(s6[s], ""), "")
   #plotSNVLength(muts[, 3+s], muts, file.name, main.text, c(1, 4, 5, 13, 17, 19, 22), "SNVs")
}

# -----------------------------------------------------------------------------
# Step 2.1: Build S6 table for early/late replicated region
# Last Modified: 08/07/19
# -----------------------------------------------------------------------------
muts <- toTable(0, 9, 22, c("CHR", "LENGTH", "TOTAL", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
muts$CHR <- 1:22
for (c in 1:22) {
   chr <- chrs[c]
   bed.gc.chr <- subset(bed.gc, CHR == chr)
 
   rpkms.chr.rt <- readTable(file.path(wd.rt.data, paste0(base, "_rpkm.corr.gc.d.rt_", chr, "_", PAIR1, "-", PAIR0, "_n", n1, "-", n0, ".txt.gz")), header=T, rownames=T, sep="\t")
   rpkms.chr.rt <- setScaledRT(rpkms.chr.rt, pseudocount=0.01, recaliRT=T, scaledRT=T) 
   rpkms.chr.rt.RT <- setSpline(rpkms.chr.rt, bed.gc.chr, "RT")

   rpkms.chr.rt.RT <- subset(rpkms.chr.rt.RT, SPLINE > 0)
   muts$LENGTH[c] <- nrow(rpkms.chr.rt.RT) * 1000/1000000
   
   ##
   load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_chr", c,".RData")))   ## Load snvs.chr
   snvs.chr <- subset(snvs.chr, SPLINE > 0)
   muts$TOTAL[c] <- nrow(snvs.chr)

   snvs.chr.s6 <- getTableS6SNV(snvs.chr[, c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL")])
   for (s in 1:6) {
      muts[c, 3+s] <- nrow(snvs.chr.s6[[s]][[1]]) + nrow(snvs.chr.s6[[s]][[2]])
   }
}
save(muts, file=file.path(wd.asym.data, paste0(base, "_mut_snvs_s6_early.RData")))

# -----------------------------------------------------------------------------
# Step ?: Vocanoplot
# Last Modified: 22/07/19
# -----------------------------------------------------------------------------
plotVolcano <- function(de, pvalue, genes, file.de, file.main) {
   #pvalue <- fdrToP(fdr, de)
   #fdr <- pvalueToFDR(pvalue, de)
   de.sig <- subset(de, P <= pvalue)
   de.sig$log10P <- -log10(de.sig$P)
   genes <- de$external_gene_name
   
   de$log10P <- -log10(de$P)
   xmax <- max(0.98)
   ymax <- max(-log10(7.614703e-14))
 
   pdf(file.de, height=8, width=3.5)
   plot(de$LOG2_FC, de$log10P, pch=16, xlim=c(-xmax, xmax), ylim=c(0, ymax), xlab="SPR correlation [rho]", ylab="-log10(p-value)", col="purple", main=file.main[1])
 
   abline(v=0, lty=5)
   #text(xmax*-1 + 2*xmax/40, -log10(pvalue) + ymax/42, paste0("rho=±", rho), cex=0.85)
 
   de.up   <- subset(de.sig, LOG2_FC > 0)
   points(de.up$LOG2_FC, de.up$log10P, pch=16, col="purple")
   de.down <- subset(de.sig, LOG2_FC < 0)
   points(de.down$LOG2_FC, de.down$log10P, pch=16, col="purple")
 
   for (g in 1:length(genes)) {
      gene <- subset(de, external_gene_name == genes[g])

      #if (gene$LOG2_FC > 0)
               text(gene$LOG2_FC, gene$log10P, genes[g], col="black", adj=c(0, -0.5), cex=0.75)
               #else
               #text(gene$LOG2_FC, gene$log10P, genes[g], col="black", adj=c(1, -0.5), cex=0.75)

   }
 
   mtext(file.main[2], cex=1.2, line=0.3)
   #legend("topright", legend=c("Inverse", "Positive"), col=c("dodgerblue", "red"), pch=19)
   dev.off()
}

##
s6 <- c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C")
for (s in 1:6) {
   rhos.s1 <- toTable(0, 2, nrow(icgcs), c("rho", "P"))
   rownames(rhos.s1) <- icgcs$Var1
 
   for (t in 1:nrow(icgcs)) {
      BASE <- as.vector(icgcs[t,]$Var1)
      base <- tolower(BASE)
  
      wd.asym       <- file.path(wd.anlys, "asymmetries", paste0(base, "-asym-rt"))
      wd.asym.data  <- file.path(wd.asym, "data")
      wd.asym.plots <- file.path(wd.asym, "plots")
  
      load(file=file.path(wd.asym.data, paste0(base, "_mut_snvs_rho_s6.RData")))
      rhos.s1$rho[t] <- rhos$rho[s]
      rhos.s1$P[t]   <- rhos$P[s]
   }
   rhos.s1$external_gene_name <- rownames(rhos.s1)
   #rhos.s1$BH <- qvalue(rhos.s1$P)$qvalue
   #rhos.s1 <- rhos.s1[order(rhos.s1$P),]
 
   ##
   plot.main <- s6[s]
   plot.de <- file.path(wd.anlys, "asymmetries", "volcanoplot_icgc_")
   file.main <- c(plot.main, "SPR vs. Mutation rates")
   file.de <- paste0(plot.de, unlist(strsplit(s6[s], "/"))[1], ".pdf")
 
   rhos.s1$LOG2_FC <- rhos.s1$rho
   plotVolcano(rhos.s1, pvalue=1e-4, NULL, file.de, file.main)
}














# -----------------------------------------------------------------------------
# Step 1: Finding mutations locate within Ensembl genes
# Last Modified: 23/01/18
# -----------------------------------------------------------------------------
for (s in 1:length(samples)) {
   sample <- samples[s]
   
   vcf <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.ngs, "peiflyne", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_mutcall_filtered.vcf")), pass=T, rs=F)
   vcf.gene <- getSNVinEnsGene(sample, vcf, ensGene)
   
   writeTable(vcf.gene, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.vcf.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 2: Keep only genes that are transcribed in this cancer type 
# Last Modified: 24/01/18
# -----------------------------------------------------------------------------   
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)   ## CHANGE 02/04/18

for (s in 1:length(samples)) {
   #sample <- SAMPLE 
   sample <- samples[s]
   
   vcf <- readTable(file.path(wd.asym.files, paste0(sample, "_mut.ens.vcf.gz")), header=T, rownames=F, sep="")
   ens.tx <- intersect(unique(vcf$ensembl_gene_id), rownames(tpm.gene.input))   ## All SNVs on "expressed" genes (regardless protein coding or not)
   vcf.tx <- subset(vcf, ensembl_gene_id %in% ens.tx)
   writeTable(vcf.tx, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.vcf.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Seperate SNVs
   vcf.tx.snv <- subset(vcf.tx,     REF %in% c("A", "T", "C", "G"))
   vcf.tx.snv <- subset(vcf.tx.snv, ALT %in% c("A", "T", "C", "G"))   
   writeTable(vcf.tx.snv, gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.snv.vcf.gz"))), colnames=T, rownames=F, sep="\t")
   
   ## Seperate Indels
   vcf.tx$TMP <- paste0(vcf.tx$REF, vcf.tx$ALT)
   vcf.tx.indel <- vcf.tx[which(nchar(vcf.tx$TMP) != 2),]
   writeTable(vcf.tx.indel[,-16], gzfile(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.indel.vcf.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Step 3.1: Read in all the SNVs
# Last Modified: 02/02/18
# -----------------------------------------------------------------------------
tx.snv <- initMergedReport(F)
for (s in 1:length(samples)) {
   sample <- samples[s]
 
   vcf <- readTable(file.path(wd.asym.files, paste0(sample, "_mut.ens.tx.snv.vcf.gz")), header=T, rownames=F, sep="")
   vcf <- subset(vcf, CHROM %in% paste0("chr", c(1:22)))   ## Keep only autosomes
   snv <- getMergedReport(sample, vcf)
   tx.snv <- rbind(tx.snv, snv)
}
save(tx.snv, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv.RData")))   ## All SNVs on "expressed" genes (regardless protein coding or not)
# > nrow(tx.snv)
# [1] 1441303

###
## Build up S6 table
tx.snv.s6 <- getTableS6(tx.snv, isExon=F)
save(tx.snv.s6, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6.RData")))   ## All SNVs on "expressed" genes (regardless protein coding or not)

# -----------------------------------------------------------------------------
# Step 4.1: Keep only SNVs on expressed, "non-redundant" protein coding genes
# Last Modified: 06/10/18
# -----------------------------------------------------------------------------
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
#tpm.gene <- tpm.gene[setdiff(rownames(tpm.gene), outliers1.5),]   ## ADD 26/06/18
#tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=F, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)   ## n=19131
tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=F, proteinCodingNonRedundantOnly=F)   ## n=18440
#tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=F)   ## n=16410
#tpm.gene.input      <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)   ## n=10604

tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

ens.tx.snv.input <- intersect(rownames(tpm.gene.input), unique(tx.snv$ensembl_gene_id))   ## Needed in Step 5; ADD 02/02/18
tx.snv.input <- subset(tx.snv, ensembl_gene_id %in% ens.tx.snv.input)
save(ens.tx.snv.input, tx.snv.input, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_input.RData")))
# > nrow(tx.snv.input)   ## All SNVs on expressed
# [1] 1441303
# > nrow(tx.snv.input)   ## All SNVs on expressed, "all" protein coding genes
# [1] 1414777
# > nrow(tx.snv.input)   ## All SNVs on expressed, "non-redundant" protein coding genes
# [1] 1021827
# > nrow(tx.snv.input)   ## After removing CNTNAP2, PTPRD, LSAMP              ## ADD 26/06/18
# [1] 992324
# > nrow(tx.snv.input)   ## After removing CNTNAP2, PTPRD, LSAMP and RBFOX1   ## ADD 26/06/18
# [1] 988794

###
## Keep only genes with at least one SNV                   ## REMOVED 20/06/18
#tx.q4.fix <- getTxQ4(ens.tx.snv.input, tpm.gene.input.log2)   ## ADD 02/04/18; Divide Q4 based on genes with at least one SNV (not all expressed genes)
#save(tx.q4, file=file.path(wd.asym.data, paste0(base, "_asym_tx_q4.RData")))
# for (q in 1:4)
#    print(length(intersect(ens.tx.snv.input, tx.q4[[q]])))
# [1] 2587
# [1] 2587
# [1] 2586
# [1] 2587
tx.q4.fix <- getTxQ4(tpm.gene.input.log2, NA)
# for (q in 1:4)
#    print(length(tx.q4.fix[[q]]))
# [1] 4610
# [1] 4610
# [1] 4610
# [1] 4610

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
for (q in 1:4)
   print(length(tx.q4[[q]]))
# [1] 4939
# [1] 4939
# [1] 4939
# [1] 4939

# -----------------------------------------------------------------------------
# Distribution of G2-M in 4 quantiles amongst SCLC
# Last Modified: 24/08/18
# -----------------------------------------------------------------------------
tx.q4.length <- initLength(tx.q4[[1]], 1)[0,]
tx.q4.cycle <- list(list(), list(), list(), list())

for (q in 1:4) {
  #genes <- intersect(ens.tx.snv.input, tx.q4[[q]])
   genes <- tx.q4[[q]]
   
   tx.q4.length <- rbind(tx.q4.length, initLength(genes, q))
   
   genes.g1s <- intersect(genes, unique(c(core.G1S, genes.G1S)))
   genes.g2m <- intersect(genes, unique(c(core.G2M, genes.G2M)))
   tx.q4.cycle[[q]][[1]] <- length(genes.g1s)
   tx.q4.cycle[[q]][[2]] <- length(genes.g2m)
   tx.q4.cycle[[q]][[3]] <- length(setdiff(genes, c(genes.g1s, genes.g2m)))
}

##
tx.q4.length$Group <- as.factor(tx.q4.length$Group)
tx.q4.length$Length <- log10(tx.q4.length$Length)
colnames <- c("Q1", "Q2", "Q3", "Q4")
rownames <- c("G1-S", "G2-M", "Others")

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_length_all.pdf"))
pdf(file.name, height=6, width=4)
ymin <- min(tx.q4.length$Length)
ymax <- max(tx.q4.length$Length)
boxplot(Length ~ Group, data=tx.q4.length, outline=T, names=colnames, ylim=c(ymin, ymax), ylab="Gene length (log10)", main="SCLC expression")
dev.off()

##
data <- toTable(0, length(colnames), 3, colnames)
rownames(data) <- rownames
for (q in 1:4)
   for (r in 1:3)
      data[r, q] <- tx.q4.cycle[[q]][[r]]
writeTable(data, file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_length_all.txt")), colnames=T, rownames=T, sep="\t")
data <- as.matrix(data)

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_tx_q4_cycle_all.pdf"))
cols <- c("lightgray", "gray", "white")
pdf(file.name, height=6, width=4)
barplot(data, ylab="Proportion of cell cycle genes", col=cols, main="SCLC expression")
legend("topleft", legend=rownames, fill=cols)
dev.off()

# -----------------------------------------------------------------------------
# Step 5.1: SNV asymetrey (Figure S1)
# Last Modified: 20/06/18
# -----------------------------------------------------------------------------
ens.asyms <- list()   ## ADD 20/06/18
asyms     <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   
   txs.onboth.fw <- intersect(unique(tx.snv.s6[[i]][[1]][[1]]$ensembl_gene_id), unique(tx.snv.s6[[i]][[2]][[1]]$ensembl_gene_id))   ## ADD 27/06/18
   txs.onboth.re <- intersect(unique(tx.snv.s6[[i]][[1]][[2]]$ensembl_gene_id), unique(tx.snv.s6[[i]][[2]][[2]]$ensembl_gene_id))   
   ens.asyms[[i]] <- intersect(rownames(tpm.gene.input), c(txs.onboth.fw, txs.onboth.re))   ## CHANGE 27/06/18; ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   
   asym <- as.matrix(toTable(0, 2, 2, c("Tx(+)", "Tx(-)")))
   rownames(asym) <- c(paste0(REFS[idx], ">", ALTS[idx]), paste0(REFS[idx+1], ">", ALTS[idx+1]))
   for (j in 1:2)
      for (k in 1:2) {
         #ens.asym   <- c(ens.asym, getMutGenes(tx.snv.s6[[i]][[j]][[k]], tpm.gene.input))   ## ADD 20/06/18
         asym[j, k] <- getMutPerMb(tx.snv.s6[[i]][[j]][[k]], tpm.gene.input[ens.asyms[[i]],])   ## Only look at SNVs on expressed genes, which is with at least one SNV on each strand
               ## E.g. getMutPerMb(tx.snv.s6[[1]][[1]][[1]])   ## E.g. C>A Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[1]][[2]])   ##      C>A Tx(-)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[1]])   ##      G>T Tx(+)
               ##      getMutPerMb(tx.snv.s6[[1]][[2]][[2]])   ##      G>T Tx(-)
      }

   #ens.asyms[[i]] <- unique(ens.asym)   ## BUG 27/06/18
   asyms[[i]] <- asym
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()
save(tx.snv.s6, ens.asyms, asyms, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6.RData")))

###
## Refining the plot
max(asyms[[1]])
# [1] 306.3362   ## After removing genes with only SNVs on one strand (ALL genes) 06/10/18
# [1] 313.1942
# [1] 326.1996   ## After removing genes with only SNVs on one strand
# [1] 318.4453   ## After removing genes with only SNVs on one strand, and genes longer than 1.5MB (CNTNAP2, PTPRD, LSAMP and RBFOX1)
# [1] 305.9242   ## After removing genes longer than 2MB (CNTNAP2, PTPRD and LSAMP)
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_ylim306.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   asym <- asyms[[i]]
   idx <- idxs[i]
   
   barplot(asym, ylab="SNVs/Mb", main=getMain(rownames(asym)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown", "sandybrown", "lightskyblue"), ylim=c(0, 306.3362))
   mtext(paste(c(paste0(REFS[idx], ":", REFS[idx+1]), paste0(REFS[idx+1], ":", REFS[idx])), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5.2: Transcription-associated SNV asymmetry (Figure S2)
# Last Modified: 22/06/18
# -----------------------------------------------------------------------------
q4s <- list()
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   idx <- idxs[i]
   ## Cntx:Gtx  =   tx.snv.s6[[1]][[1]][[1]]) C>A Tx(+)  +  tx.snv.s6[[1]][[2]][[2]] G>T Tx(-)
   CntxGtx <- rbind(tx.snv.s6[[i]][[1]][[1]],               tx.snv.s6[[i]][[2]][[2]])
   ## Gntx:Ctx  =   tx.snv.s6[[1]][[2]][[1]]) G>T Tx(+)  +  tx.snv.s6[[1]][[1]][[2]] C>A Tx(-)       
   GntxCtx <- rbind(tx.snv.s6[[i]][[2]][[1]],               tx.snv.s6[[i]][[1]][[2]])
   
   ens.asym <- ens.asyms[[i]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
   ens.asym.q4 <- getTxQ4(tpm.gene.input.log2, ens.asym) 
   
   q4 <- as.matrix(toTable(0, 4, 2, c("n=", "n=", "n=", "n=")))
   rownames(q4) <- c(paste0(REFS[idx], "ntx:", REFS[idx+1], "tx"), paste0(REFS[idx+1], "ntx:", REFS[idx], "tx"))
   for (q in 1:4) {
      txs <- ens.asym.q4[[q]]   ## ADD 20/06/18; Divide Q4 more precisely based on genes with specific base changes
      
      CntxGtx.cg <- subset(CntxGtx, ensembl_gene_id %in% txs)
      txs.cg <- intersect(txs, unique(CntxGtx.cg$ensembl_gene_id))   ## SWAP txs order BUG BUG BUG 29/06/18
      q4[1, q] <- getMutPerMbTxs(CntxGtx.cg, txs.cg)   ## ADD 22/06/18: Also implemented txs.cg in getMutPerMbTxs() now

      GntxCtx.gc <- subset(GntxCtx, ensembl_gene_id %in% txs)
      txs.gc <- intersect(txs, unique(GntxCtx.gc$ensembl_gene_id))   ## SWAP txs order BUG BUG BUG 29/06/18
      q4[2, q] <- getMutPerMbTxs(GntxCtx.gc, txs.gc)   ## ADD 22/06/18: Also implemented txs.gc in getMutPerMbTxs() now
      
      colnames(q4)[q] <- paste0(colnames(q4)[q], length(unique(c(txs.cg, txs.gc))), " (", length(txs.cg), "+", length(txs.gc), ")")
   }
   
   q4s[[i]] <- q4
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(q4)), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))
}
dev.off()
save(tx.snv.s6, ens.asyms, asyms, q4s, file=file.path(wd.asym.data, paste0(base, "_asym_tx_snv_s6_q4s.RData")))

## ADD 16/03/18
q4s.rt <- list()
q4s.rt[[1]] <- q4s

###
## Refining the plot
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "            Tx", "", "")
   #for (c in 1:ncol(q4))
   #   colnames(q4)[c] <- gsub("n=", "", unlist(strsplit(colnames(q4)[c], " "))[1])
   
   barplot(q4, ylab="SNVs/Mb", main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=c("lightskyblue", "sandybrown"))   #, xlab="Number of genes")
   mtext(paste(rownames(q4), collapse=" vs "), cex=0.55, font=3, line=0.5)
}
dev.off()

# -----------------------------------------------------------------------------
# Step 5.3: Log2 transcription-associated SNV asymmetry (Figure S3)
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
pdf(file.path(wd.asym.plots, paste0(base, "_asym_tx_snv_s6_q4s_log2_ylim1.75.pdf")), height=4.5, width=7)
par(mfrow=c(2, 3))
for (i in 1:length(idxs)) {
   q4 <- q4s[[i]]
   colnames(q4) <- c("", "           Tx", "", "")
   
   barplot(-log2(q4[1,]/q4[2,]), ylab="TCR efficiency", ylim=c(-0.5, 1.75), main=getMain(rownames(asyms[[i]])), beside=TRUE, width=.3, col=getLog2Colours(q4))   ## ADD 08/03/18
   mtext(paste0("-log2(", paste(rownames(q4), collapse="/"), ")"), cex=0.55, font=3, line=0.5)
}
dev.off()











# -----------------------------------------------------------------------------
# Test: EZH2 ~ CNTNAP2 ?
# Last Modified: 26/01/18
# -----------------------------------------------------------------------------
file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2-EZH2.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2~EZH2, ylab="CNTNAP2", xlab="EZH2", main="Gene expression")
dev.off()

CNTNAP2.muts <- table[overlaps,]$Freq

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2-CNTNAP2.muts.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2~CNTNAP2.muts, ylab="CNTNAP2", xlab="CNTNAP2 mutation", main="Gene expression")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_EZH2-CNTNAP2.muts.pdf"))
pdf(file.name, height=6, width=6)
plot(EZH2~CNTNAP2.muts, ylab="EZH2", xlab="CNTNAP2 mutation", main="Gene expression")
dev.off()

##
CNTNAP2 <- as.numeric(tpm.gene.input.log2["ENSG00000174469", overlaps])

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2-CNTNAP2_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2.ratios~CNTNAP2, ylab="CNTNAP2 TCR ratio (log2)", xlab="CNTNAP2 log2(TPM + 0.01)", main="SCLC (n=70)")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2_mut-CNTNAP2_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2.ratios~CNTNAP2.muts, ylab="CNTNAP2 TCR ratio (log2)", xlab="CNTNAP2 mutation", main="SCLC (n=70)")
dev.off()

##
RB1 <- as.numeric(tpm.gene.input.log2["ENSG00000139687", overlaps])

file.name <- file.path(wd.asym.plots, paste0(base, "_RB1-RB1_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(RB1.ratios~RB1, ylab="RB1 TCR ratio (log2)", xlab="RB1 log2(TPM + 0.01)", main="SCLC (n=70)")
dev.off()

file.name <- file.path(wd.asym.plots, paste0(base, "_CNTNAP2_mut-CNTNAP2_ratio_y.pdf"))
pdf(file.name, height=6, width=6)
plot(CNTNAP2.ratios~CNTNAP2.muts, ylab="CNTNAP2 TCR ratio (log2)", xlab="CNTNAP2 mutation", main="SCLC (n=70)")
dev.off()

###
## 34 R(L)-Tx(+)/Q3/TCD genes
i <- 1
ratios <- log2(as.numeric(q4s.s6.rt.st.mut.cg[[1]][[1]][[1]][[3]])/as.numeric(q4s.s6.rt.st.mut.gc[[1]][[1]][[1]][[3]]))
idx <- which(ratios > 0)
genes <- q4s.s6.rt.st.gen[[1]][[1]][[1]][[3]][idx]

CntxGtx      <- tx.snv.s6[[i]][[1]][[1]]      ## ADD 23/06/18
GntxCtx      <- tx.snv.s6[[i]][[2]][[1]]
n <- c()
for (g in 1:length(genes)) {
   gene <- genes[g]

   CntxGtx.gene <- subset(CntxGtx, ensembl_gene_id == gene)
   GntxCtx.gene <- subset(GntxCtx, ensembl_gene_id == gene)
   s1 <- rbind(CntxGtx.gene, GntxCtx.gene)
   samples <- unique(s1$SAMPLE)
   n <- c(n, length(samples))
}

de$SRC_RHO <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[4]])
de$SRC_P   <- mapply(x = 1:nrow(expr.pheno.log2), function(x) cor.test(as.numeric(expr.pheno.log2[x,]), pheno.all$RB1_RATE, method="spearman", exact=F)[[3]])

