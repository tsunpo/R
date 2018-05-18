# =============================================================================
# Manuscript   : The dangeous case of DNA replication
# Chapter      : Gene-to-gene distance is indicative of replication origin
# Name         : manuscripts/expression/lcnec-tpm-g2g.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/05/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"            ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"                 ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"                  ## tpyang@localhost

wd.src.handbook <- file.path(wd.src, "handbook-of")   ## Required handbooks/libraries for the manuscript
handbooks <- c("Common.R", "Asymmetry.R", "ReplicationOrigin.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.handbook, x))))

wd.src.guide <- file.path(wd.src, "guide-to-the")     ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.guide, "hg19.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
BASE <- "LCNEC"
#wd     <- paste0("/ngs/cangen/tyang2/", BASE, "/analysis/")                   ## tyang2@gauss
#wd.ngs <- paste0("/ngs/cangen/tyang2/", BASE, "/ngs/RNA/")
wd     <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/")   ## tpyang@localhost
wd.ngs <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/ngs/RNA/")

wd.asym       <- paste0(wd, "asymmetries/", tolower(BASE), "-asym-tx/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "lcnec_rna_n69.list"), header=F, rownames=F, sep="")[,1]

## Phenotypic data
load("/Users/tpyang/Work/uni-koeln/tyang2/LCNEC/analysis/expression/FPKM/lcnec_expr.RData")
pheno.expr <- pheno[samples,]

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
BASE <- "LUSQ"
#wd     <- paste0("/ngs/cangen/tyang2/", BASE, "/analysis/")                   ## tyang2@gauss
#wd.ngs <- paste0("/ngs/cangen/tyang2/", BASE, "/ngs/RNA/")
wd     <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/analysis/")   ## tpyang@localhost
wd.ngs <- paste0("/Users/tpyang/Work/uni-koeln/tyang2/", BASE, "/ngs/RNA/")

wd.asym       <- paste0(wd, "asymmetries/", tolower(BASE), "-asym-tx/")
wd.asym.data  <- paste0(wd.asym, "data/")
wd.asym.plots <- paste0(wd.asym, "plots/")
setwd(wd.asym)

samples <- readTable(paste0(wd.ngs, "lusq_rna_n21.list"), header=F, rownames=F, sep="")[,1]

# =============================================================================
# Step 1: Gene-to-gene minmum distance 
# Last Modified: 18/05/18
# =============================================================================
load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/data/all_kallisto_0.43.1_tpm.gene_r5_p47.RData")
tpm.gene <- tpm.gene[,samples]
#tpm.gene <- tpm.gene[,rownames(subset(pheno.expr, RB1NEW == 1))]
#tpm.gene <- tpm.gene[,rownames(subset(pheno.expr, RB1NEW == 0))]
tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2633
# [1] 2633
# [1] 2633
# [1] 2633

g2g.q4 <- getG2GQ4(tx.q4)
p3 <- testT(g2g.q4[[3]], g2g.q4[[4]])
p2 <- testT(g2g.q4[[2]], g2g.q4[[4]])
p1 <- testT(g2g.q4[[1]], g2g.q4[[4]])
c(p1, p2, p3)
# [1] 4.417072e-06 2.653005e-04 2.932690e-02
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 1.325122e-05 7.959016e-04 8.798070e-02

distances <- c()
quantiles <- c()
for (q in 1:4) {
   distances <- c(distances, as.numeric(g2g.q4[[q]]))
   quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
}
pdf(paste0(wd.asym.plots, "lcnec_asym_tx_g2g.pdf"), height=6, width=4)
boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=c("25%", "50%", "75%", "100%"), ylim=c(0.3, ymax), main="LCNEC (n=69)")
dev.off()

###
## RB1
c(p1, p2, p3)
# [1] 2.788324e-06 2.237995e-05 1.335054e-02
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 8.364973e-06 6.713984e-05 4.005161e-02

distances <- c()
quantiles <- c()
for (q in 1:4) {
   distances <- c(distances, as.numeric(g2g.q4[[q]]))
   quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
}
pdf(paste0(wd.asym.plots, "lcnec_asym_tx_g2g_RB1.pdf"), height=6, width=4)
boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=c("25%", "50%", "75%", "100%"), ylim=c(0.3, ymax), main="LCNEC (RB1; n=20)")
dev.off()

###
## RB1-WT
c(p1, p2, p3)
# [1] 2.702025e-05 7.955730e-04 5.846751e-02
p.adjust(c(p1, p2, p3), method="bonferroni")
# [1] 8.106074e-05 2.386719e-03 1.754025e-01

ymax <- log10(max(as.numeric(g2g.q4[[2]])))
distances <- c()
quantiles <- c()
for (q in 1:4) {
   distances <- c(distances, as.numeric(g2g.q4[[q]]))
   quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
}
pdf(paste0(wd.asym.plots, "lcnec_asym_tx_g2g_RB1-WT.pdf"), height=6, width=4)
boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=c("25%", "50%", "75%", "100%"), ylim=c(0.3, ymax), main="LCNEC (RB1-WT; n=34)")
dev.off()







# =============================================================================
# Step X: Min gene to gene distance 
# Last Modified: 14/05/18
# =============================================================================
getTSS <- function(ensGene.genes) {
   ensGene.genes$TSS <- 0
 
   for (g in 1:nrow(ensGene.genes)) {
      if (ensGene.genes$strand[g] > 0)
         ensGene.genes$TSS[g] <- ensGene.genes$start_position[g]
      else
         ensGene.genes$TSS[g] <- ensGene.genes$end_position[g]
   }
 
   return(ensGene.genes)
}

###
##
load("/Users/tpyang/Work/uni-koeln/tyang2/ALL/analysis/expression/kallisto/luad-lcnec-sclc-rnaseq-de/data/all_kallisto_0.43.1_tpm.gene_r5_p47.RData")
tpm.gene <- tpm.gene[,samples]
#tpm.gene <- tpm.gene[,rownames(subset(pheno.expr, RB1NEW == 1))]
#tpm.gene <- tpm.gene[,rownames(subset(pheno.expr, RB1NEW == 0))]

tpm.gene.input <- getEnsGeneFiltered(tpm.gene, ensGene, autosomeOnly=T, proteinCodingOnly=T, proteinCodingNonRedundantOnly=T)
tpm.gene.input.log2 <- getLog2andMedian(tpm.gene.input)

tx.q4 <- getTxQ4(NA, tpm.gene.input.log2)
# > for (q in 1:4)
#  + print(length(tx.q4[[q]]))
# [1] 2633
# [1] 2633
# [1] 2633
# [1] 2633

g2g.q4 <- list()
for (q in 1:4) {
   genes <- tx.q4[[q]]
   ensGene.genes <- ensGene[genes,]
   ensGene.genes <- getTSS(ensGene.genes)
 
   g2g <- list()
   for (g in 1:length(genes)) {
      ensGene.gene <- ensGene.genes[genes[g],]
      ensGene.genes.chr <- subset(ensGene.genes, chromosome_name == ensGene.gene$chromosome_name)
      ensGene.genes.chr <- subset(ensGene.genes.chr, ensembl_gene_id != ensGene.gene$ensembl_gene_id)
    
      g2g[[g]] <- min(abs(ensGene.gene$TSS - ensGene.genes.chr$TSS))
   }
   g2g.q4[[q]] <- g2g
}
## RB1-MUT
# > t.test(as.numeric(g2g.q4[[3]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 0.01335054
# > t.test(as.numeric(g2g.q4[[2]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 2.237995e-05
# > t.test(as.numeric(g2g.q4[[1]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 2.788324e-06
# > p.adjust(c(0.01335054, 2.237995e-05, 2.788324e-06), method="BH")
# [1] 1.335054e-02 3.356993e-05 8.364972e-06

## RB1-WT
# > t.test(as.numeric(g2g.q4[[3]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 0.05846751
# > t.test(as.numeric(g2g.q4[[2]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 0.000795573
# > t.test(as.numeric(g2g.q4[[1]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 2.702025e-05
# > p.adjust(c(0.05846751, 0.000795573, 2.702025e-05), method="BH")
# [1] 5.846751e-02 1.193360e-03 8.106075e-05

## LCNEC
# > t.test(as.numeric(g2g.q4[[3]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 0.0293269
# > t.test(as.numeric(g2g.q4[[2]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 0.0002653005
# > t.test(as.numeric(g2g.q4[[1]]), as.numeric(g2g.q4[[4]]))$p.value
# [1] 4.417072e-06
# > p.adjust(c(0.0293269, 0.0002653005, 4.417072e-06), method="BH")
# [1] 2.932690e-02 3.979508e-04 1.325122e-05

###
##
distances <- c()
quantiles <- c()
for (q in 1:4) {
   distances <- c(distances, as.numeric(g2g.q4[[q]]))
   quantiles <- c(quantiles, mapply(x = 1:length(g2g.q4[[q]]), function(x) paste0("q", q)))
}

#ymax <- log10(max(as.numeric(g2g.q4[[2]])))
pdf(paste0(wd.asym.plots, "lcnec_asym_tx_g2g.pdf"), height=6, width=4)
boxplot(log10(distances)~quantiles, ylab="Gene-to-gene min dist. (log10)", xlab="Expression", names=c("25%", "50%", "75%", "100%"), ylim=c(0.3, ymax), main="LCNEC")
dev.off()
