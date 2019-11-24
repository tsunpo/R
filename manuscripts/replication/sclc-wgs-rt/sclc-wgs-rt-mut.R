# =============================================================================
# Manuscript   : 
# Chapter      : Chromosome replication timing of the human genome
# Name         : manuscripts/replication/sclc-wgs-rt.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/10/19; 05/03/19; 25/02/19; 30/01/18
# =============================================================================
wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
#wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "DifferentialExpression.R", "ReplicationTiming.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Step 0: Set working directory
# Last Modified: 30/01/18
# -----------------------------------------------------------------------------
wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
#wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE  <- "SCLC"
PAIR1 <- "T"
PAIR0 <- "N"
base  <- tolower(BASE)
method <- "rpkm"

wd.ngs   <- file.path(wd, BASE, "ngs/WGS")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "George 2015")

wd.rt       <- file.path(wd.anlys, "replication", paste0(base, "-wgs-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

wd.ngs.data <- file.path(wd.ngs, "data")
samples <- readTable(file.path(wd.ngs, "sclc_wgs_n101.txt"), header=T, rownames=T, sep="")
genes <- readTable(file.path(wd.meta, "SCLC_Mutationtable_ordered.txt"), header=T, rownames=T, sep="")[, -1]
genes <- genes[rownames(samples),]

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
plotBox <- function(gene, wd.de.plots, tpm.gene.log2, de.tpm.gene, samples, type, ylim=NULL) {
   ids <- subset(ensGene, external_gene_name == gene)   ## E.g. RBL1 (ENSG00000080839 and ENSG00000269846)
   chr <- ids$chromosome_name[1]
 
   for (i in 1:nrow(ids)) {
      id <- ids$ensembl_gene_id[i]
 
      file.name <- paste0("boxplot_", tolower(type), ".gene.log2_", gene)
      if (nrow(ids) != 1)
         file.name <- paste0(file.name, "_", id)
  
      gene.tpm <- cbind(t(tpm.gene.log2[id, rownames(samples)]), samples)
      colnames(gene.tpm)[1] <- "MEDIAN"
  
      pdf(file.path(wd.de.plots, paste0(file.name, ".pdf")), height=6, width=3.5)
      boxplot(gene.tpm$MEDIAN ~ gene.tpm[, chr], outline=T, names=c("CM1", "CM2"), ylab=paste0("log2(", type, "+0.01)"), main=paste0(gene, " (", id, ")"))
  
      mtext(paste0("p-value = ", scientific(de.tpm.gene[id,]$P)), cex=1.2, line=0.3)
      dev.off()
   }
}

##
pvalues <- toTable(0, 2, ncol(genes), c("GENE", "P"))
for (g in 1:ncol(genes)) {
   samples$MUT <- as.factor((genes[, g]))
   gene <- colnames(genes)[g]
   pvalue <- testU(subset(samples, MUT == 0)$COR, subset(samples, MUT == 1)$COR)
   pvalues$GENE[g] <- gene
   pvalues$P[g] <- pvalue
   
   pdf(file.path(wd.rt.plots, paste0("COR~", gene, ".pdf")), height=6, width=3.5)
   boxplot(samples$COR ~ samples$MUT, outline=T, names=c("WT", "MUT"), ylab="Spearman's rho", main=paste0(gene, " vs. cell-cycle status"))
 
   mtext(paste0("p-value = ", scientific(pvalue)), cex=1.2, line=0.3)
   dev.off()
}
pvalues <- pvalues[order(pvalues$P),]
pvalues$Q <- qvalue(pvalues$P)$qvalue
pvalues$BH <- p.adjust(pvalues$P, "BH")
writeTable(pvalues, file=file.path(wd.rt.plots, paste0("COR~MUTs.txt")), colnames=T, rownames=F, sep="\t")
