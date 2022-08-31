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
handbooks  <- c("Commons.R", "Asymmetry.R", "Mutation.R", "Survival.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
#wd <- "/projects/cangen/tyang2"              ## tyang2@cheops
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost
BASE <- "NBL"
base <- tolower(BASE)

wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata", "data_release")
wd.vcf   <- file.path(wd, BASE, "ngs/WGS/peiflyne")

wd.driver <- file.path(wd.anlys, "driver", paste0(base, "-driver"))
wd.driver.data  <- file.path(wd.driver, "data")
wd.driver.plots <- file.path(wd.driver, "plots")

# -----------------------------------------------------------------------------
# Load data
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
for (s in 1:nrow(samples.nbl)) {
   sample <- samples.nbl[s,]$SAMPLE_ID
 
   colnames <- c("ensembl_gene_id", "CN")
   ensGene.cna <- toTable("", length(colnames), 0, colnames)
   segs <- read.peiflyne.cn.seg(file.path("/Users/tpyang/Work/uni-koeln/tyang2/NBL/ngs/WGS/peiflyne/2015", sample, paste0(sample, "_ANALYSIS"), paste0(sample, "_iCN.seg")))   ## See ReplicationTiming.R / Inner Class / PeifLyne File Reader
   for (g in 1:nrow(segs)) {
      seg <- segs[g,]
      ensGene.chr.end.start <- getEnsGenesFromSegment(seg)
  
      if (nrow(ensGene.chr.end.start) > 0) {
         ens <- toTable("", length(colnames), nrow(ensGene.chr.end.start), colnames)
         ens$ensembl_gene_id <- ensGene.chr.end.start$ensembl_gene_id
         ens$CN              <- seg$CN
   
         ensGene.cna <- rbind(ensGene.cna, ens)
      }
   }
 
   writeTable(ensGene.cna, gzfile(file.path(wd.driver.data, "copy_number_alterations", "2015", paste0(sample, "_ens.iCN.seg.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Differential copy numbers (CNA)
# Last Modified: 25/06/22
# -----------------------------------------------------------------------------
#samples.nbl$SORTING <- "G1"
#idx <- which(samples.nbl$COR >= rho.nbl)
#samples.nbl[idx,]$SORTING <- "S"
#samples.nbl$SORTING <- as.factor(samples.nbl$SORTING)

cna.gene <- getGeneCNV(file.path(wd.driver.data, "copy_number_alterations", "2015"), samples.nbl)
cna.gene.nona <- cna.gene[removeMissing(cna.gene),]
save(samples.nbl, cna.gene, cna.gene.nona, file=file.path(wd.driver.data, "copy_number_alterations", "2015", "nbl-driver-cna_iCN.seg_NONA_n57.RData"), version=2)
#cna.gene.no0  <- cna.gene[getExpressed(cna.gene),]
#cna.gene.nona.no0  <- cna.gene[getExpressed(cna.gene.nona),]
dim(cna.gene)
# [1] 57773    56
dim(cna.gene.nona)
# [1] 57736    56

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.RData")))
cna.gene.nona.tpm <- cna.gene.nona[intersect(rownames(cna.gene.nona), rownames(tpm.gene)),]
dim(cna.gene.nona.tpm)
# [1] 34895    56
load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene.median0.RData")))
cna.gene.nona.tpm.median0 <- cna.gene.nona[intersect(rownames(cna.gene.nona), rownames(tpm.gene)),]
dim(cna.gene.nona.tpm.median0)
# [1] 22886    56

# -----------------------------------------------------------------------------
# WGD
# Last Modified: 28/06/22; 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
wgds <- readTable(file.path(wd.ngs, "nbl_wgs_n57_WGD.txt"), header=T, rownames=T, sep="\t")
wgds <- wgds[rownames(samples.nbl),]
wgds[which(wgds$wgd_predicted == 0), ]$wgd_predicted <- "non-WGD"
wgds[which(wgds$wgd_predicted == 1), ]$wgd_predicted <- "WGD"

##
test <- cna.gene.nona.tpm.median0[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.iCN <- getPCA(t(test))
save(pca.iCN, file=file.path(wd.driver.data, "2015", paste0("PCA-NBL-iCN-TPM-MEDIAN0_n57_g22886.RData")))

trait <- wgds$wgd_predicted
file.main <- c("NB gene CNA", "")
plotPCA(1, 2, pca.iCN, trait, wd.driver.plots, "PCA-NBL-iCN-TPM-MEDIAN0_WGD_n57_g22886", size=5, file.main, "topleft", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=1, flip.y=-1)

scores.iCN <- pcaScores(pca.iCN)
x <- scores.iCN$PC1
y <- wgds$decision_value
file.name <- file.path(wd.driver.plots, "2015", paste0("correlation_PCA-NBL-iCN-TPM-MEDIAN0_PC1-vs-WGD"))
plotCorrelation(file.name, "NBL gene CNA", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "WGD", x, y, size=5)

x <- scores.iCN$PC1
y <- samples.nbl$COR
file.name <- file.path(wd.driver.plots, "2015", paste0("correlation_PCA-NBL-iCN-TPM-MEDIAN0_PC1-vs-SORTING"))
plotCorrelation(file.name, "NBL gene CNA", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Proliferation rate", x, y, size=5)

x <- scores.iCN$PC1
y <- samples.wgs$purity2
file.name <- file.path(wd.driver.plots, "2015", paste0("correlation_PCA-NBL-iCN-TPM-MEDIAN0_PC1-vs-purity_n57"))
plotCorrelation(file.name, "NBL gene CNA", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Purity", x, y, size=5)

x <- scores.iCN$PC1
y <- samples.wgs$ploidy2
file.name <- file.path(wd.driver.plots, "2015", paste0("correlation_PCA-NBL-iCN-TPM-MEDIAN0_PC1-vs-ploidy_n57"))
plotCorrelation(file.name, "NBL gene CNA", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Ploidy", x, y, size=5)










##
test <- cna.gene.nona[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.iCN <- getPCA(t(test))
save(pca.iCN, file=file.path(wd.driver.data, paste0("PCA-NBL-iCN_n56_g57736.RData")))

trait <- wgds$wgd_predicted
file.main <- c("NB gene CNA (57,736 genes)", "")
plotPCA(1, 2, pca.iCN, trait, wd.driver.plots, "PCA-NBL-iCN_WGD_n56_g57736", size=6, file.main, "topleft", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=-1, flip.y=-1)

scores.iCN <- pcaScores(pca.iCN)
x <- scores.iCN$PC1 * -1
y <- wgds$decision_value
file.name <- file.path(wd.driver.plots, paste0("correlation_PCA-NBL-iCN_PC1-vs-WGD"))
plotCorrelation(file.name, "NB gene CNA (57,736 genes)", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "WGD", x, y, size=6)

x <- scores.iCN$PC1 * -1
y <- samples.nbl$COR
file.name <- file.path(wd.driver.plots, paste0("correlation_PCA-NBL-iCN_PC1-vs-SORTING"))
plotCorrelation(file.name, "NB gene CNA (57,736 genes)", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Proliferation", x, y, size=6)

##
#test <- cna.gene.nona.tpm[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#pca.iCN <- getPCA(t(test))
#save(pca.iCN, file=file.path(wd.driver.data, paste0("PCA-SCLC-iCN-TPM_n56_g34895.RData")))
#
#trait <- wgds$wgd_predicted
#file.main <- c("SCLC gene CNA (34,895 genes)", "")
#plotPCA(1, 2, pca.iCN, trait, wd.driver.plots, "PCA-SCLC-iCN-TPM_WGD_n101_g34895", size=6, file.main, "topleft", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=1, flip.y=1)
#
#scores.iCN <- pcaScores(pca.iCN)
#x <- scores.iCN$PC1
#y <- wgds$decision_value
#file.name <- file.path(wd.driver.plots, paste0("correlation_PCA-SCLC-iCN-TPM_PC1-vs-WGD"))
#plotCorrelation(file.name, "SCLC gene CNA (34,895 genes)", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "WGD", x, y, size=6)
#
#x <- scores.iCN$PC1
#y <- samples.sclc$COR
#file.name <- file.path(wd.driver.plots, paste0("correlation_PCA-SCLC-iCN-TPM_PC1-vs-SORTING"))
#plotCorrelation(file.name, "SCLC gene CNA (34,895 genes)", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), expression(italic("In silico") ~ "sorting"), x, y, size=6)

##
test <- cna.gene.nona.tpm.median0[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.iCN <- getPCA(t(test))
#save(pca.iCN, file=file.path(wd.driver.data, paste0("PCA-NBL-iCN-TPM-MEDIAN0_n56_g22886.RData")))

trait <- wgds$wgd_predicted
file.main <- c("NB gene CNA", "")
plotPCA(1, 2, pca.iCN, trait, wd.driver.plots, "PCA-NBL-iCN-TPM-MEDIAN0_WGD_n56_g22886", size=6, file.main, "topleft", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=-1, flip.y=1)

scores.iCN <- pcaScores(pca.iCN)
x <- scores.iCN$PC1 * -1
y <- wgds$decision_value
file.name <- file.path(wd.driver.plots, paste0("correlation_PCA-NBL-iCN-TPM-MEDIAN0_PC1-vs-WGD"))
plotCorrelation(file.name, "NB gene CNA", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "WGD", x, y, size=6)

x <- scores.iCN$PC1 * -1
y <- samples.nbl$COR
file.name <- file.path(wd.driver.plots, paste0("correlation_PCA-NBL-iCN-TPM-MEDIAN0_PC1-vs-SORTING"))
plotCorrelation(file.name, "NB gene CNA", paste0("PC", 1, " (", pcaProportionofVariance(pca.iCN, 1), "%)"), "Proliferation", x, y, size=6)






# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/06/22; 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
#test <- cna.gene.nona[, rownames(samples.nbl)]
#pca.de <- getPCA(t(test))
#save(pca.iCN, file=file.path(wd.driver.data, paste0("PCA_NBL_iCN_n57_g45906.RData")))

test <- cna.gene.nona.tpm[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.iCN <- getPCA(t(test))
save(pca.iCN, file=file.path(wd.driver.data, paste0("PCA_NBL_iCN_median0_n55_g22886.RData")))

samples.nbl$GROUP_ID <- "LR"
samples.nbl[rownames(subset(samples.nbl, RISK == "high")),]$GROUP_ID <- "HR"

samples.nbl$GROUP_ID2 <- "LR"
samples.nbl[rownames(subset(subset(samples.nbl, RISK == "high"), SORTING == "G1")),]$GROUP_ID2 <- "G1-HR"
samples.nbl[rownames(subset(subset(samples.nbl, RISK == "high"), SORTING == "S")),]$GROUP_ID2  <- "S-HR"

trait <- samples.nbl$GROUP_ID
file.main <- c("NB gene CNA (57,736 genes)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA-NBL-iCN_LR-HR_n56_g57736", size=6, file.main, "topright", c("HR", "LR"), c(yellow, "skyblue"), flip.x=1, flip.y=1)

trait <- samples.nbl$GROUP_ID2
file.main <- c("NB gene CNA (57,736 genes)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA-NBL-iCN_LR-HR-S-G1_n56_g57736", size=6, file.main, "topright", c("S-HR", "G1-HR", "LR"), c(red, yellow, blue), flip.x=1, flip.y=1)

trait <- samples.nbl$GROUP_ID2
file.main <- c("NB gene CNA (57,736 genes)", "")
plotPCA(2, 3, pca.de, trait, wd.driver.plots, "PCA-NBL-iCN_LR-HR-S-G1_n56_g57736", size=6, file.main, "topright", c("S-HR", "G1-HR", "LR"), c(red, yellow, blue), flip.x=1, flip.y=1)

trait <- as.vector(samples.nbl$SORTING)
file.main <- c("NB gene CNA (57,736 genes)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA-NBL-iCN_G1-S_n56_g57736", size=6, file.main, "topright", c("S", "G1"), c(red, blue), flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# WGD
# Last Modified: 28/06/22; 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
load(file=file.path(wd.driver.data, paste0("PCA_NBL_iCN_n57_g57736.RData")))

scores <- pcaScores(pca.iCN)
overlaps <- intersect(rownames(scores), rownames(samples.nbl))
wgds <- readTable(file.path(wd.ngs, "nbl_wgs_n57_WGD.txt"), header=T, rownames=T, sep="\t")
wgds <- wgds[overlaps,]
wgds[which(wgds$wgd_predicted == 0), ]$wgd_predicted <- "non-WGD"
wgds[which(wgds$wgd_predicted == 1), ]$wgd_predicted <- "WGD"

y <- wgds$decision_value
x <- scores$PC1
file.name <- file.path(wd.driver.plots, paste0("correlation_WGD_vs_PC1"))
plotCorrelation(file.name, "NBL samples", "PC1 (55%)", "WGD", x, y, size=6)

y <- wgds$decision_value
x <- scores$PC2
file.name <- file.path(wd.driver.plots, paste0("correlation_WGD_vs_PC2"))
plotCorrelation(file.name, "NBL samples", "PC2 (24.9%)", "WGD", x, y, size=6)

y <- wgds$decision_value
x <- samples.nbl$COR
file.name <- file.path(wd.driver.plots, paste0("correlation_WGD_vs_SORTING"))
plotCorrelation(file.name, "NBL samples", expression(italic("In silico") ~ "sorting"), "WGD", x, y, size=6)

##
trait <- wgds$wgd_predicted
file.main <- c("NBL gene CNA (57,736 genes)", "")
plotPCA(1, 2, pca.iCN, trait, wd.driver.plots, "PCA-NBL-iCN_WGD_n56_g57736", size=6, file.main, "topright", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=1, flip.y=1)
plotPCA(2, 3, pca.iCN, trait, wd.driver.plots, "PCA-NBL-iCN_WGD_n56_g57736", size=6, file.main, "topright", c("WGD", "non-WGD"), c("black", "lightgray"), flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# 
# Last Modified: 26/06/22
# -----------------------------------------------------------------------------
scores <- pcaScores(pca.de)

y <- samples.nbl$COR
x <- scores$PC1
file.name <- file.path(wd.driver.plots, paste0("correlation_SORTING_vs_PC1_g57736"))
plotCorrelation(file.name, "NB samples (57,736 genes)", "PC1 (55%)", expression(italic("In silico") ~ "sorting"), x, y, size=6)

y <- samples.nbl$COR
x <- scores$PC2
file.name <- file.path(wd.driver.plots, paste0("correlation_SORTING_vs_PC2_g57736"))
plotCorrelation(file.name, "NB samples (57,736 genes)", "PC2 (24.9%)", expression(italic("In silico") ~ "sorting"), x, y, size=6)

y <- samples.nbl$COR
x <- scores$PC3
file.name <- file.path(wd.driver.plots, paste0("correlation_SORTING_vs_PC3_g57736"))
plotCorrelation(file.name, "NB samples (57,736 genes)", "PC3 (6.6%)", expression(italic("In silico") ~ "sorting"), x, y, size=6)

# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/06/22; 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
test <- cna.gene.nona.tpm[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.driver.data, paste0("PCA_NBL_n56_g18751.RData")))

samples.nbl$GROUP_ID <- "LR"
samples.nbl[rownames(subset(samples.nbl, RISK == "high")),]$GROUP_ID <- "HR"

samples.nbl$GROUP_ID2 <- "LR"
samples.nbl[rownames(subset(subset(samples.nbl, RISK == "high"), SORTING == "G1")),]$GROUP_ID2 <- "G1-HR"
samples.nbl[rownames(subset(subset(samples.nbl, RISK == "high"), SORTING == "S")),]$GROUP_ID2  <- "S-HR"

trait <- samples.nbl$GROUP_ID
file.main <- c("NB gene CNA (18,751 expressed)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR_n56_g18751", size=6, file.main, "topright", c("HR", "LR"), c(yellow, "skyblue"), flip.x=1, flip.y=1)

trait <- samples.nbl$GROUP_ID2
file.main <- c("NB gene CNA (18,751 expressed)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR-S-G1_n56_g18751", size=6, file.main, "topright", c("S-HR", "G1-HR", "LR"), c(red, yellow, blue), flip.x=1, flip.y=1)

trait <- samples.nbl$GROUP_ID2
file.main <- c("NB gene CNA (18,751 expressed)", "")
plotPCA(1, 3, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR-S-G1_n56_g18751", size=6, file.main, "topright", c("S-HR", "G1-HR", "LR"), c(red, yellow, blue), flip.x=1, flip.y=1)

trait <- as.vector(samples.nbl$SORTING)
file.main <- c("NB gene CNA (18,751 expressed)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA_NBL_G1-S_n56_g18751", size=6, file.main, "topright", c("S", "G1"), c(red, blue), flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# 
# Last Modified: 26/06/22
# -----------------------------------------------------------------------------
scores <- pcaScores(pca.de)

y <- samples.nbl$COR
x <- scores$PC1
file.name <- file.path(wd.driver.plots, paste0("correlation_SORTING_vs_PC1_g18751"))
plotCorrelation(file.name, "NB samples (18,751 expressed)", "PC1 (55.7%)", expression(italic("In silico") ~ "sorting"), x, y, size=6)

y <- samples.nbl$COR
x <- scores$PC2
file.name <- file.path(wd.driver.plots, paste0("correlation_SORTING_vs_PC2_g18751"))
plotCorrelation(file.name, "NB samples (18,751 expressed)", "PC2 (24.6%)", expression(italic("In silico") ~ "sorting"), x, y, size=6)

y <- samples.nbl$COR
x <- scores$PC3
file.name <- file.path(wd.driver.plots, paste0("correlation_SORTING_vs_PC3_g18751"))
plotCorrelation(file.name, "NB samples (18,751 expressed)", "PC3 (6.4%)", expression(italic("In silico") ~ "sorting"), x, y, size=6)






trait <- samples.nbl$GROUP_ID
file.main <- c("NB CNA on genes (n=39,998)", "")
plotPCA(1, 3, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR_n56_g39998", size=6, file.main, "topright", c("HR", "LR"), c(yellow, "skyblue"), flip.x=1, flip.y=1)

trait <- samples.nbl$GROUP_ID2
file.main <- c("NB CNA on genes (n=39,998)", "")
plotPCA(1, 3, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR-S-G1_n56_g39998", size=6, file.main, "topright", c("S-HR", "G1-HR", "LR"), c(red, yellow, blue), flip.x=1, flip.y=1)

trait <- as.vector(samples.nbl$SORTING)
file.main <- c("NB CNA on genes (n=39,998)", "")
plotPCA(1, 3, pca.de, trait, wd.driver.plots, "PCA_NBL_G1-S_n56_g39998", size=6, file.main, "topright", c("S", "G1"), c(red, blue), flip.x=1, flip.y=1)








# -----------------------------------------------------------------------------
# PCA
# Last Modified: 26/03/21; 23/08/20
# -----------------------------------------------------------------------------
test <- cna.gene.tpm.no0[, rownames(samples.nbl)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
pca.de <- getPCA(t(test))
save(pca.de, file=file.path(wd.driver.data, paste0("PCA_NBL_iCN_n56_g13529.RData")))

samples.nbl$GROUP_ID <- "LR"
samples.nbl[rownames(subset(samples.nbl, RISK == "high")),]$GROUP_ID <- "HR"

samples.nbl$GROUP_ID2 <- "LR"
samples.nbl[rownames(subset(subset(samples.nbl, RISK == "high"), SORTING == "G1")),]$GROUP_ID2 <- "G1-HR"
samples.nbl[rownames(subset(subset(samples.nbl, RISK == "high"), SORTING == "S")),]$GROUP_ID2  <- "S-HR"

trait <- samples.nbl$GROUP_ID
file.main <- c("NB CNA on expressed genes (n=13,529)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR_n56_g13529", size=6, file.main, "topright", c("HR", "LR"), c(yellow, "skyblue"), flip.x=1, flip.y=1)

trait <- samples.nbl$GROUP_ID2
file.main <- c("NB CNA on expressed genes (n=13,529)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA_NBL_LR-HR-S-G1_n56_g13529", size=6, file.main, "topright", c("S-HR", "G1-HR", "LR"), c(red, yellow, blue), flip.x=1, flip.y=1)

trait <- as.vector(samples.nbl$SORTING)
file.main <- c("NB CNA on expressed genes (n=13,529)", "")
plotPCA(1, 2, pca.de, trait, wd.driver.plots, "PCA_NBL_G1-S_n56_g13529", size=6, file.main, "topright", c("S", "G1"), c(red, blue), flip.x=1, flip.y=1)




# -----------------------------------------------------------------------------
# RMRs (RISK: Low vs. High)
# Last Modified: 05/06/22
# -----------------------------------------------------------------------------


s=0
mut=0

pmr.cna.gene <- toTable(0, 12, nrow(cna.gene), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_L", "MUT_H", "WT_L", "WT_H", "P", "FDR", "LMR", "HMR", "RMR"))
pmr.cna.gene$MUT <- rownames(cna.gene)
rownames(pmr.cna.gene) <- pmr.cna.gene$MUT
 
pmr.mut.gene$MUT_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene[x,])))
pmr.mut.gene$SAMPLE_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) length(which(mut.gene[x,] != 0)))
pmr.mut.gene <- subset(pmr.mut.gene, SAMPLE_FREQ > s)
 
for (g in 1:nrow(pmr.mut.gene)) {
   gene <- pmr.mut.gene$MUT[g]
   ids  <- colnames(mut.gene)[which(mut.gene[gene,] != 0)]

   samples.mut.mut <- subset(samples.nbl, SAMPLE_ID %in% ids)
   samples.mut.wt  <- samples.nbl[setdiff(rownames(samples.nbl), rownames(samples.mut.mut)),]
  
   mut.g1 <- rownames(subset(samples.mut.mut, RISK == "low"))
   mut.s  <- rownames(subset(samples.mut.mut, RISK == "high"))
   wt.g1  <- rownames(subset(samples.mut.wt,  RISK == "low"))
   wt.s   <- rownames(subset(samples.mut.wt,  RISK == "high"))
  
   if ((length(mut.g1) > mut) && (length(mut.s) > mut) && (length(wt.g1) > mut) && (length(wt.s) > mut)) {
      test <- toTable(0, 2, 2, c("L", "H"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- length(mut.g1)
      test[1, 2] <- length(mut.s)
      test[2, 1] <- length(wt.g1)
      test[2, 2] <- length(wt.s)
   
      pmr.mut.gene$MUT_L[g] <- length(mut.g1)
      pmr.mut.gene$MUT_H[g] <- length(mut.s)
      pmr.mut.gene$WT_L[g]  <- length(wt.g1)
      pmr.mut.gene$WT_H[g]  <- length(wt.s)
      #pmr.mut.gene$P[g]    <- fisher.test(test)[[1]]
      pmr.mut.gene$P[g]     <- chisq.test(test)[[3]]
   
      pmr.mut.gene$LMR[g] <- pmr.mut.gene$MUT_L[g] / pmr.mut.gene$WT_L[g]
      pmr.mut.gene$HMR[g] <- pmr.mut.gene$MUT_H[g] / pmr.mut.gene$WT_H[g]
      pmr.mut.gene$RMR[g] <- pmr.mut.gene$HMR[g]   / pmr.mut.gene$LMR[g]
   }
}
pmr.mut.gene <- subset(pmr.mut.gene, P != 0)
pmr.mut.gene <- pmr.mut.gene[order(pmr.mut.gene$P, decreasing=F),]

annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
pmr.mut.gene <- cbind(annot[rownames(pmr.mut.gene),], pmr.mut.gene[, -1])
pmr.mut.gene$FDR <- testFDR(pmr.mut.gene$P, "Q")

file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>0_mut>0_g=", nrow(pmr.mut.gene), "_chisq_RISK_Low-vs-High")
save(pmr.mut.gene, file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(pmr.mut.gene, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

rmr.mut.gene.s10 <- subset(pmr.mut.gene, SAMPLE_FREQ >= 10)
rmr.mut.gene.s10$FDR <- testFDR(rmr.mut.gene.s10$P, "Q")
file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>10_mut>0_g=", nrow(pmr.mut.gene.s10), "_chisq_RISK_Low-vs-High")
writeTable(pmr.mut.gene.s10, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

#colnames(rmr.mut.gene.s10)[18] <- "PMR"
#file.name <- file.path(wd.driver.plots, paste0("RMR-RISK_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>10_mut>0_g=", nrow(rmr.mut.gene.s10), "_chisq_GRID2"))
#plotPMR(file.name, "Risk mutational rate", "RMR", text.Log10.P, rmr.mut.gene.s10, c("GRID2", "ESR1", "CEP112", "ATRX"), size=5, cols=c(adjustcolor.skyblue, adjustcolor.yellow), xmax=28, ymax=5)

#colnames(rmr.mut.gene.s10)[18] <- "PMR"
#file.name <- file.path(wd.driver.plots, paste0("RMR-RISK_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>10_mut>0_g=", nrow(rmr.mut.gene.s10), "_chisq_CDH12"))
#plotPMR(file.name, "Risk mutational rate", "RMR", text.Log10.P, rmr.mut.gene.s10, c("GRID2", "CDH12", "CEP112", "ATRX", "PDE4DIP"), size=5, cols=c(adjustcolor.skyblue, adjustcolor.yellow), xmax=28, ymax=5)

# -----------------------------------------------------------------------------
# PMRs (SORTING: G1 vs. S)
# Last Modified: 05/06/22
# -----------------------------------------------------------------------------
s=0
mut=0

pmr.mut.gene <- toTable(0, 12, nrow(mut.gene), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "P", "FDR", "G1MR", "SMR", "PMR"))
pmr.mut.gene$MUT <- rownames(mut.gene)
rownames(pmr.mut.gene) <- pmr.mut.gene$MUT

pmr.mut.gene$MUT_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene[x,])))
pmr.mut.gene$SAMPLE_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) length(which(mut.gene[x,] != 0)))
pmr.mut.gene <- subset(pmr.mut.gene, SAMPLE_FREQ > s)

for (g in 1:nrow(pmr.mut.gene)) {
   gene <- pmr.mut.gene$MUT[g]
   ids  <- colnames(mut.gene)[which(mut.gene[gene,] != 0)]
 
   samples.mut.mut <- subset(samples.nbl, SAMPLE_ID %in% ids)
   samples.mut.wt  <- samples.nbl[setdiff(rownames(samples.nbl), rownames(samples.mut.mut)),]
 
   mut.g1 <- rownames(subset(samples.mut.mut, SORTING == "G1"))
   mut.s  <- rownames(subset(samples.mut.mut, SORTING == "S"))
   wt.g1  <- rownames(subset(samples.mut.wt,  SORTING == "G1"))
   wt.s   <- rownames(subset(samples.mut.wt,  SORTING == "S"))
 
   if ((length(mut.g1) > mut) && (length(mut.s) > mut) && (length(wt.g1) > mut) && (length(wt.s) > mut)) {
      test <- toTable(0, 2, 2, c("G1", "S"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- length(mut.g1)
      test[1, 2] <- length(mut.s)
      test[2, 1] <- length(wt.g1)
      test[2, 2] <- length(wt.s)
  
      pmr.mut.gene$MUT_G1[g] <- length(mut.g1)
      pmr.mut.gene$MUT_S[g] <- length(mut.s)
      pmr.mut.gene$WT_G1[g]  <- length(wt.g1)
      pmr.mut.gene$WT_S[g]  <- length(wt.s)
      #pmr.mut.gene$P[g]    <- fisher.test(test)[[1]]
      pmr.mut.gene$P[g]     <- chisq.test(test)[[3]]
  
      pmr.mut.gene$G1MR[g] <- pmr.mut.gene$MUT_G1[g] / pmr.mut.gene$WT_G1[g]
      pmr.mut.gene$SMR[g] <- pmr.mut.gene$MUT_S[g] / pmr.mut.gene$WT_S[g]
      pmr.mut.gene$PMR[g] <- pmr.mut.gene$SMR[g]   / pmr.mut.gene$G1MR[g]
   }
}
pmr.mut.gene <- subset(pmr.mut.gene, P != 0)
pmr.mut.gene <- pmr.mut.gene[order(pmr.mut.gene$P, decreasing=F),]

annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
pmr.mut.gene <- cbind(annot[rownames(pmr.mut.gene),], pmr.mut.gene[, -1])
pmr.mut.gene$FDR <- testFDR(pmr.mut.gene$P, "Q")

file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>0_mut>0_g=", nrow(pmr.mut.gene), "_chisq_SORTING_G1-vs-S")
save(pmr.mut.gene, file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(pmr.mut.gene, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

##
pmr.mut.gene.s10 <- subset(pmr.mut.gene, SAMPLE_FREQ >= 10)
pmr.mut.gene.s10$FDR <- testFDR(pmr.mut.gene.s10$P, "Q")
file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>10_mut>0_g=", nrow(pmr.mut.gene.s10), "_chisq_SORTING_G1-vs-S")
writeTable(pmr.mut.gene.s10, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

file.name <- file.path(wd.driver.plots, paste0("PMR-SORTING_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>10_mut>0_g=", nrow(pmr.mut.gene.s10), "_chisq_ESR1"))
plotPMR(file.name, "Proliferative mutation rate", "PMR", text.Log10.P, pmr.mut.gene.s10, c("ESR1", "ATRX"), size=5, xmax=49, ymax=6.7)

##
#pmr.mut.gene.s20 <- subset(pmr.mut.gene, SAMPLE_FREQ >= 20)
#pmr.mut.gene.s20$FDR <- testFDR(pmr.mut.gene.s20$P, "Q")
#file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>20_mut>0_g=", nrow(pmr.mut.gene.s20), "_chisq_SORTING_G1-vs-S")
#writeTable(pmr.mut.gene.s20, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
#
#file.name <- file.path(wd.driver.plots, paste0("PMR-SORTING_", base, "_mut.gene_n=", nrow(samples.nbl), "_s>20_mut>0_g=", nrow(pmr.mut.gene.s20), "_chisq_ESR1+NCAM2"))
#plotPMR(file.name, "Proliferational mutation rate", "PMR", text.Log10.P, pmr.mut.gene.s20, c("ESR1", "NCAM2", "ATRX"), size=5, xmax=47, ymax=6.5)

# -----------------------------------------------------------------------------
# PMRs (SORTING:  vs. )
# Last Modified: 05/06/22
# -----------------------------------------------------------------------------
samples.nbl.high <- subset(samples.nbl, RISK == "high")
mut.gene.high <- mut.gene[, rownames(samples.nbl.high)]

s=0
mut=0

pmr.mut.gene <- toTable(0, 12, nrow(mut.gene), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "P", "FDR", "G1MR", "SMR", "PMR"))
pmr.mut.gene$MUT <- rownames(mut.gene)
rownames(pmr.mut.gene) <- pmr.mut.gene$MUT

pmr.mut.gene$MUT_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene.high[x,])))
pmr.mut.gene$SAMPLE_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) length(which(mut.gene.high[x,] != 0)))
pmr.mut.gene <- subset(pmr.mut.gene, SAMPLE_FREQ > s)

for (g in 1:nrow(pmr.mut.gene)) {
   gene <- pmr.mut.gene$MUT[g]
   ids  <- colnames(mut.gene.high)[which(mut.gene.high[gene,] != 0)]
 
   samples.mut.mut <- subset(samples.nbl.high, SAMPLE_ID %in% ids)
   samples.mut.wt  <- samples.nbl.high[setdiff(rownames(samples.nbl.high), rownames(samples.mut.mut)),]
 
   mut.g1 <- rownames(subset(samples.mut.mut, SORTING == "G1"))
   mut.s  <- rownames(subset(samples.mut.mut, SORTING == "S"))
   wt.g1  <- rownames(subset(samples.mut.wt,  SORTING == "G1"))
   wt.s   <- rownames(subset(samples.mut.wt,  SORTING == "S"))
 
   if ((length(mut.g1) > mut) && (length(mut.s) > mut) && (length(wt.g1) > mut) && (length(wt.s) > mut)) {
      test <- toTable(0, 2, 2, c("G1", "S"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- length(mut.g1)
      test[1, 2] <- length(mut.s)
      test[2, 1] <- length(wt.g1)
      test[2, 2] <- length(wt.s)
  
      pmr.mut.gene$MUT_G1[g] <- length(mut.g1)
      pmr.mut.gene$MUT_S[g] <- length(mut.s)
      pmr.mut.gene$WT_G1[g]  <- length(wt.g1)
      pmr.mut.gene$WT_S[g]  <- length(wt.s)
      #pmr.mut.gene$P[g]    <- fisher.test(test)[[1]]
      pmr.mut.gene$P[g]     <- chisq.test(test)[[3]]
  
      pmr.mut.gene$G1MR[g] <- pmr.mut.gene$MUT_G1[g] / pmr.mut.gene$WT_G1[g]
      pmr.mut.gene$SMR[g] <- pmr.mut.gene$MUT_S[g] / pmr.mut.gene$WT_S[g]
      pmr.mut.gene$PMR[g] <- pmr.mut.gene$SMR[g]   / pmr.mut.gene$G1MR[g]
   }
}
pmr.mut.gene <- subset(pmr.mut.gene, P != 0)
pmr.mut.gene <- pmr.mut.gene[order(pmr.mut.gene$P, decreasing=F),]

annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
pmr.mut.gene <- cbind(annot[rownames(pmr.mut.gene),], pmr.mut.gene[, -1])
pmr.mut.gene$FDR <- testFDR(pmr.mut.gene$P, "Q")

file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl.high), "_s>0_mut>0_g=", nrow(pmr.mut.gene), "_chisq_SORTING_G1-vs-S")
save(pmr.mut.gene, file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(pmr.mut.gene, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

##
pmr.mut.gene.s10 <- subset(pmr.mut.gene, SAMPLE_FREQ >= 10)
pmr.mut.gene.s10$FDR <- testFDR(pmr.mut.gene.s10$P, "Q")
file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl.high), "_s>10_mut>0_g=", nrow(pmr.mut.gene.s10), "_chisq_SORTING_G1-vs-S")
writeTable(pmr.mut.gene.s10, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

file.name <- file.path(wd.driver.plots, paste0("PMR-SORTING_", base, "_mut.gene_n=", nrow(samples.nbl.high), "_s>10_mut>0_g=", nrow(pmr.mut.gene.s10), "_chisq_ESR1"))
plotPMR(file.name, "Proliferative mutation rate", "PMR", text.Log10.P, pmr.mut.gene.s10, c("ESR1", "ATRX"), size=5, xmax=31.2, ymax=3.73)

##
#pmr.mut.gene.s20 <- subset(pmr.mut.gene, SAMPLE_FREQ >= 20)
#pmr.mut.gene.s20$FDR <- testFDR(pmr.mut.gene.s20$P, "Q")
#file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl.high), "_s>20_mut>0_g=", nrow(pmr.mut.gene.s20), "_chisq_SORTING_G1-vs-S")
#writeTable(pmr.mut.gene.s20, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
#
#file.name <- file.path(wd.driver.plots, paste0("PMR-SORTING_", base, "_mut.gene_n=", nrow(samples.nbl.high), "_s>20_mut>0_g=", nrow(pmr.mut.gene.s20), "_chisq_ESR1+NCAM2"))
#plotPMR(file.name, "Proliferational mutation rate", "PMR", text.Log10.P, pmr.mut.gene.s20, c("ESR1", "NCAM2", "ATRX"), size=5, xmax=47, ymax=6.5)






# -----------------------------------------------------------------------------
# PMRs (SORTING:  vs. )
# Last Modified: 07/06/22
# -----------------------------------------------------------------------------
samples.nbl.G1 <- subset(samples.nbl, SORTING == "G1")
mut.gene.G1 <- mut.gene[, rownames(samples.nbl.G1)]

s=0
mut=0

pmr.mut.gene <- toTable(0, 12, nrow(mut.gene), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_L", "MUT_H", "WT_L", "WT_H", "P", "FDR", "LMR", "HMR", "RMR"))
pmr.mut.gene$MUT <- rownames(mut.gene)
rownames(pmr.mut.gene) <- pmr.mut.gene$MUT

pmr.mut.gene$MUT_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) sum(as.numeric(mut.gene.G1[x,])))
pmr.mut.gene$SAMPLE_FREQ <- mapply(x = 1:nrow(pmr.mut.gene), function(x) length(which(mut.gene.G1[x,] != 0)))
pmr.mut.gene <- subset(pmr.mut.gene, SAMPLE_FREQ > s)

for (g in 1:nrow(pmr.mut.gene)) {
   gene <- pmr.mut.gene$MUT[g]
   ids  <- colnames(mut.gene.G1)[which(mut.gene.G1[gene,] != 0)]
 
   samples.mut.mut <- subset(samples.nbl.G1, SAMPLE_ID %in% ids)
   samples.mut.wt  <- samples.nbl.G1[setdiff(rownames(samples.nbl.G1), rownames(samples.mut.mut)),]
 
   mut.g1 <- rownames(subset(samples.mut.mut, RISK == "low"))
   mut.s  <- rownames(subset(samples.mut.mut, RISK == "high"))
   wt.g1  <- rownames(subset(samples.mut.wt,  RISK == "low"))
   wt.s   <- rownames(subset(samples.mut.wt,  RISK == "high"))
 
   if ((length(mut.g1) > mut) && (length(mut.s) > mut) && (length(wt.g1) > mut) && (length(wt.s) > mut)) {
      test <- toTable(0, 2, 2, c("L", "H"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- length(mut.g1)
      test[1, 2] <- length(mut.s)
      test[2, 1] <- length(wt.g1)
      test[2, 2] <- length(wt.s)
  
      pmr.mut.gene$MUT_L[g] <- length(mut.g1)
      pmr.mut.gene$MUT_H[g] <- length(mut.s)
      pmr.mut.gene$WT_L[g]  <- length(wt.g1)
      pmr.mut.gene$WT_H[g]  <- length(wt.s)
      #pmr.mut.gene$P[g]    <- fisher.test(test)[[1]]
      pmr.mut.gene$P[g]     <- chisq.test(test)[[3]]
  
      pmr.mut.gene$LMR[g] <- pmr.mut.gene$MUT_L[g] / pmr.mut.gene$WT_L[g]
      pmr.mut.gene$HMR[g] <- pmr.mut.gene$MUT_H[g] / pmr.mut.gene$WT_H[g]
      pmr.mut.gene$RMR[g] <- pmr.mut.gene$HMR[g]   / pmr.mut.gene$LMR[g]
   }
}
pmr.mut.gene <- subset(pmr.mut.gene, P != 0)
pmr.mut.gene <- pmr.mut.gene[order(pmr.mut.gene$P, decreasing=F),]

annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
pmr.mut.gene <- cbind(annot[rownames(pmr.mut.gene),], pmr.mut.gene[, -1])
pmr.mut.gene$FDR <- testFDR(pmr.mut.gene$P, "Q")

file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl.G1), "_s>0_mut>0_g=", nrow(pmr.mut.gene), "_chisq_RISK_L-vs-H_in_G1")
save(pmr.mut.gene, file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(pmr.mut.gene, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

pmr.mut.gene.s10 <- subset(pmr.mut.gene, SAMPLE_FREQ >= 10)
pmr.mut.gene.s10$FDR <- testFDR(pmr.mut.gene.s10$P, "Q")
file.name <- paste0("pmr_", base, "_mut.gene_n=", nrow(samples.nbl.G1), "_s>10_mut>0_g=", nrow(pmr.mut.gene.s10), "_chisq_RISK_L-vs-H_in_G1")
writeTable(pmr.mut.gene.s10, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")










# =============================================================================
# Reference    : HGNC (hg19)
# Last Modified: 04/05/22
# =============================================================================
vcfs <- c()
for (s in 1:nrow(samples.wes)) {
   sample <- paste0(samples.wes[s,]$Patient_ID, "-B")
 
   muts <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.vcf, sample, paste0(sample, "_mutcall_filtered.vcf")), pass=T, rs=F)
   vcfs <- c(vcfs, mapply(x = 1:nrow(muts), function(x) unlist(strsplit(muts$INFO[x], ";TG="))[2]))
}

hgnc.ref.ens

ensGene.hgnc.overlaps <- subset(ensGene.hgnc, external_gene_name %in% hgncs)
hgncs.diff <- setdiff(hgncs, ensGene.overlaps$external_gene_name)

hgnc <- subset(hgnc, ensembl_gene_id %in% ensGene$ensembl_gene_id)







#overlaps <- intersect(rownames(release), list)
#length(overlaps)
# [1] 2521
overlaps <- intersect(totals$wgs_id, list)
length(overlaps)
# [1] 2674

#nrow(samples.mut)
# [1] 2373
overlaps <- intersect(rownames(samples), totals[overlaps,]$specimen_id)
length(overlaps)
# [1] 2542
samples.mut <- samples[overlaps,]
rownames(totals) <- totals$specimen_id
samples.mut$tumor_wgs_aliquot_id <- totals[overlaps,]$wgs_id

# -----------------------------------------------------------------------------
# Assigning Ensembl genes to each SNVs
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
for (s in 1:nrow(samples.mut)) {
   sample.mut <- samples.mut[s,]
 
   muts <- read.peiflyne.mutcall.filtered.vcf(file.path(wd.icgc.vcf, paste0(sample.mut$tumor_wgs_aliquot_id, "_mutcall_filtered.vcf.gz")), pass=T, rs=F)
   muts.gene <- getSNVinEnsGene(muts, ensGene)
 
   writeTable(muts.gene, gzfile(file.path(wd.driver.data, "point_mutations", paste0(sample.mut$icgc_specimen_id, "_mutcall_filtered_ens.vcf.gz"))), colnames=T, rownames=F, sep="\t")
}

# -----------------------------------------------------------------------------
# Differential point mutations (SNV)
# Last Modified: 23/03/22
# -----------------------------------------------------------------------------
colnames <- rownames(samples.mut)
rownames <- rownames(ensGene)
mut.gene <- toTable(0, length(colnames), length(rownames), colnames)
rownames(mut.gene) <- rownames

for (s in 1:nrow(samples.mut)) {
   sample.mut <- samples.mut[s,]
 
   mut.dupl <- readTable(file.path(wd.driver.data, "point_mutations", paste0(sample.mut$icgc_specimen_id, "_mutcall_filtered_ens.vcf.gz")), header=T, rownames=F, sep="")
   mut <- as.data.frame(sort(table(mut.dupl$ensembl_gene_id), decreasing=T))
   rownames(mut) <- mut$Var1
   
   overlaps <- intersect(rownames, rownames(mut))
   mut.gene[overlaps, s] <- mut[overlaps,]$Freq
}

save(samples.mut, mut.gene, file=file.path(wd.driver.data, "icgc-driver-mut.RData"), version=2)

# -----------------------------------------------------------------------------
# Differential point mutations (SNV)
# Last Modified: 11/04/22; 23/03/22
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.RData")
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.log2.m.RData")
tpm.gene.m1 <- tpm.gene[rownames(subset(tpm.gene.m, MEDIAN >= 1)),]
nrow(tpm.gene.m1)
# [1] 20720   ## TPM > 0
# [1] 13475   ## TPM > 1
expressed <- intersect(rownames(mut.gene), rownames(tpm.gene.m1))

overlaps <- intersect(rownames(samples.surv.hist), rownames(samples.mut))
length(overlaps)
# [1] 1545
samples.surv.hist.mut <- samples.surv.hist[overlaps,]

#for (h in c(3, 6, 7, 13, 11)) {
#for (h in c(8, 12, 14)) { 
for (h in c(3, 6, 4, 7, 11, 13, 8, 12, 14)) {
   hist <- hists.surv[h]
   samples.test <- subset(samples.surv.hist.mut, histology_abbreviation == hist)
   mut.gene.test <- mut.gene[expressed, rownames(samples.test)]
 
   ##
   de.mut.gene.test <- getPMR(mut.gene.test, samples.test)
   file.name <- paste0("pmr_", base, "_mut.gene_chisq_s>20_tpm>1_g", nrow(de.mut.gene.test), "_", hist)
   save(de.mut.gene.test,  file=file.path(wd.driver.data, "hists", paste0(file.name, ".RData")))
   writeTable(de.mut.gene.test, file.path(wd.driver.data, "hists", paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
 
   ##
   file.name <- file.path(wd.driver.plots, "hists", paste0("PMR_", base, "_mut.gene_chisq_s>20_tpm>1_g", nrow(de.mut.gene.test), "_", hist))
   x <- de.mut.gene.test$PMR
   y <- -log10(de.mut.gene.test$P)
   plotPMR(file.name, hist, "PMR", text.Log10.P, x, y)  
}

# -----------------------------------------------------------------------------
# Differential point mutations (SNV)
# Last Modified: 11/04/22; 23/03/22
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.RData")
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.log2.m.RData")

overlaps <- intersect(rownames(samples.surv.hist), rownames(samples.mut))
length(overlaps)
# [1] 1545
samples.surv.hist.mut <- samples.surv.hist[overlaps,]

for (h in c(3, 6, 7, 13, 11)) {
   hist <- hists.surv[h]
   samples.test <- subset(samples.surv.hist.mut, histology_abbreviation == hist)
   
   ##
   overlaps <- intersect(rownames(samples.test), colnames(tpm.gene))
   tpm.gene.h   <- tpm.gene[, overlaps]   ## VERY VERY VERY IMPORTANT!!!
   tpm.gene.h.m <- tpm.gene.h
   tpm.gene.h.m$MEDIAN <- mapply(x = 1:nrow(tpm.gene.h), function(x) median(as.numeric(tpm.gene.h[x,])))
   
   tpm.gene.h.m <- subset(tpm.gene.h.m, MEDIAN > 0)
   overlaps2 <- intersect(rownames(mut.gene), rownames(tpm.gene.h.m))
   mut.gene.test <- mut.gene[, rownames(samples.test)]
  
   
   ##
   de.mut.gene.test <- getPMR(mut.gene.test, samples.test)
   file.name <- paste0("pmr_", base, "_mut.gene_chisq_s>20_tpm>1_g", nrow(de.mut.gene.test), "_", hist)
   save(de.mut.gene.test,  file=file.path(wd.driver.data, "hists", paste0(file.name, ".RData")))
   writeTable(de.mut.gene.test, file.path(wd.driver.data, "hists", paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
   
   ##
   file.name <- file.path(wd.driver.plots, "hists", paste0("PMR_", base, "_mut.gene_chisq_s>20_tpm>1_g", nrow(de.mut.gene.test), "_", hist))
   x <- de.mut.gene.test$PMR
   y <- -log10(de.mut.gene.test$P)
   plotPMR(file.name, hist, "PMR", text.Log10.P, x, y)  
}

# -----------------------------------------------------------------------------
# Wilcoxon rank sum test (non-parametric)
# Last Modified: 12/04/22; 12/03/22
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Mann–Whitney/U/wilcox.test
## FDR : Q
argv      <- data.frame(predictor="SORTING", predictor.wt="G1", test="Wilcoxon", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", hist, "_tpm-gene-median0_SG1_wilcox_q_n904")
file.main <- paste0("", BASE)

de <- differentialAnalysis(tpm.gene.log2, samples.tpm, argv$predictor, argv$predictor.wt, argv$test, argv$test.fdr)
# Samples with MUT SG1: 594   ## N=904
# Samples with  WT SG1: 310
# Samples with MUT SG1: 247   ## N=395
# Samples with  WT SG1: 148

## Ensembl gene annotations
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.tpm.gene <- cbind(annot[rownames(de),], de)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")), version=2)
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
nrow(de.tpm.gene)


# -----------------------------------------------------------------------------
# Differential point mutations (SNV) in each cancer type
# Last Modified: 05/04/22
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.RData")
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.log2.m.RData")

tpm.gene <- tpm.gene[rownames(subset(tpm.gene.m, MEDIAN >= 5)),]
nrow(tpm.gene)
# [1] 20720   ## TPM > 0
# [1] 13475   ## TPM > 1
# [1] 10260   ## TPM > 5

colnames <- intersect(colnames(mut.gene), colnames(tpm.gene))   ## 902 patients
mut.gene.tpm <- mut.gene[, colnames]
tpm.gene.mut <- tpm.gene[, colnames]
samples.mut.tpm <- samples.mut[colnames,]

rownames <- intersect(rownames(mut.gene), rownames(tpm.gene))   ## 18523 genes
mut.gene.tpm <- mut.gene.tpm[rownames,]
tpm.gene.mut <- tpm.gene.mut[rownames,]
dim(mut.gene.tpm)
# [1] 20720   904
dim(tpm.gene.mut)
# [1] 20720   904
dim(samples.mut.tpm)
# [1] 904  13

hists.tpm <- as.data.frame(sort(table(samples.mut.tpm$histology_abbreviation), decreasing=T))
hists.tpm <- subset(hists.tpm, Freq > 20)
for (h in 1:nrow(hists.tpm)) {
   hist <- hists.tpm$Var1[h]
   samples.mut.tpm.hist <- subset(samples.mut.tpm, histology_abbreviation == hist)
   
   colnames <- rownames(samples.mut.tpm.hist)
   mut.gene.tpm.hist <- mut.gene.tpm[, colnames]
   tpm.gene.mut.hist <- tpm.gene.mut[, colnames]

   de.mut.tpm.hist <- getPMR(mut.gene.tpm.hist, samples.mut.tpm.hist)
   file.name <- paste0("pmr_", base, "_mut-tpm_chisq_q_tpm>0_s>20_g", nrow(de.mut.tpm.hist), "_mean_", hist)
   save(de.mut.tpm.hist,  file=file.path(wd.driver.data, "hists", paste0(file.name, ".RData")))
   writeTable(de.mut.tpm.hist, file.path(wd.driver.data, "hists", paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
   
   ###
   ##
   file.name <- file.path(wd.driver.plots, "hists", paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>0_s>20_g", nrow(de.mut.tpm.hist), "_median_", hist, "_EPMR"))
   x <- de.mut.tpm.hist$PMR_2
   y <- -log10(de.mut.tpm.hist$P_2)
   plotPMR(file.name, "", "EPMR", text.Log10.P, x, y)
   
   file.name <- file.path(wd.driver.plots, "hists", paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>0_s>20_g", nrow(de.mut.tpm.hist), "_median_", hist))
   x <- de.mut.tpm.hist$PMR
   y <- -log10(de.mut.tpm.hist$P)
   plotPMR(file.name, "", "PMR", text.Log10.P, x, y)
}

# -----------------------------------------------------------------------------
# Differential point mutations (SNV)
# Last Modified: 23/03/22
# -----------------------------------------------------------------------------
de.mut.gene <- toTable(0, 12, nrow(table(rownames(mut.gene))), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "G1MR", "SMR", "PMR", "P", "FDR"))
de.mut.gene$MUT <- rownames(mut.gene)
rownames(de.mut.gene) <- de.mut.gene$MUT
de.mut.gene$MUT_FREQ <- mapply(x = 1:nrow(de.mut.gene), function(x) sum(as.numeric(mut.gene[x,])))
de.mut.gene$SAMPLE_FREQ <- mapply(x = 1:nrow(de.mut.gene), function(x) length(which(mut.gene[x,] != 0)))
de.mut.gene.tmp <- de.mut.gene
de.mut.gene <- subset(de.mut.gene, SAMPLE_FREQ >= 20)
nrow(de.mut.gene)
# [1] 31420
file.name <- paste0("de_", base, "_mut-gene_chisq_q_n2542_s>20_n31420")
save(de.mut.gene, file=file.path(wd.driver.data, paste0(file.name, ".RData")))

mut.gene.test <- mut.gene[de.mut.tpm$MUT,]
## Fishers' test
for (g in 1:nrow(mut.gene.test)) {
   ids <- colnames(mut.gene.test[,which(mut.gene.test[g,] != 0)])
 
   samples.mut.mut <- subset(samples.mut, icgc_specimen_id %in% ids)
   samples.mut.wt  <- samples.mut[setdiff(rownames(samples.mut), rownames(samples.mut.mut)),]
  
   mut.g1 <- rownames(subset(samples.mut.mut, SORTING == "G1"))
   mut.s  <- rownames(subset(samples.mut.mut, SORTING == "S"))
   wt.g1  <- rownames(subset(samples.mut.wt,  SORTING == "G1"))
   wt.s   <- rownames(subset(samples.mut.wt,  SORTING == "S"))
      
   if ((length(mut.g1) > 10) && (length(mut.s) > 10) && (length(wt.g1) > 10) && (length(wt.s) > 10)) {
      test <- toTable(0, 2, 2, c("G1", "S"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- length(mut.g1)
      test[1, 2] <- length(mut.s)
      test[2, 1] <- length(wt.g1)
      test[2, 2] <- length(wt.s)
       
      de.mut.gene$MUT_G1[g] <- length(mut.g1)
      de.mut.gene$MUT_S[g]  <- length(mut.s)
      de.mut.gene$WT_G1[g]  <- length(wt.g1)
      de.mut.gene$WT_S[g]   <- length(wt.s)
      #de.mut.tpm$P[g]     <- fisher.test(test)[[1]]
      de.mut.gene$P[g]      <- chisq.test(test)[[3]]
       
      de.mut.gene$G1MR[g] <- de.mut.gene$MUT_G1[g] / de.mut.gene$WT_G1[g]
      de.mut.gene$SMR[g]  <- de.mut.gene$MUT_S[g]  / de.mut.gene$WT_S[g]
      de.mut.gene$PMR[g]  <- de.mut.gene$SMR[g]    / de.mut.gene$G1MR[g]
   }
}
de.mut.gene <- de.mut.gene[order(de.mut.gene$P),]
de.mut.gene <- subset(de.mut.gene, P != 0)
nrow(de.mut.gene)
# [1] 

## Ensembl gene annotations
de <- de.mut.gene
rownames(de) <- de$MUT
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.mut.gene <- cbind(annot[rownames(de),], de[, -1])

#de.mut.gene <- de.mut.gene[intersect(rownames(de.mut.gene), rownames(tpm.gene.log2)),]
#de.mut.gene$FDR <- testFDR(de.mut.gene$P, "Q")

###
##
file.name <- paste0("de_", base, "_mut-gene_chisq_q_n2542")
save(de.mut.gene,  file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(de.mut.gene, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

de.mut.gene <- subset(subset(de.mut.gene, MUT_G1 > 10), MUT_S > 10)

###
##
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_n2343_g", nrow(de.mut.tpm)))
x <- de.mut.tpm$PMR
y <- -log10(de.mut.tpm$P)
plotPMR(file.name, "", "PMR", text.Log10.P, x, y)

expressed <- intersect(rownames(de.mut.gene), rownames(tpm.gene.log2))
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_n2343_g", length(expressed)))
x <- de.mut.tpm[expressed,]$PMR
y <- -log10(de.mut.tpm[expressed,]$P)
plotPMR(file.name, "", "PMR", text.Log10.P, x, y)

unexpressed <- setdiff(rownames(tpm.gene.log2), expressed)
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_n2343_g", length(unexpressed)))
x <- de.mut.tpm[unexpressed,]$PMR
y <- -log10(de.mut.tpm[unexpressed,]$P)
plotPMR(file.name, "", "PMR", text.Log10.P, x, y)




##
overlaps <- intersect(rownames(de.mut.gene), rownames(de.tpm.gene))
x <- -log10(de.tpm.gene[overlaps,]$P)
y <- -log10(de.mut.gene[overlaps,]$P)
file.name <- file.path(wd.driver.data, paste0("correlation_", base, "_de-mut-gene_vs_de-tpm-gene_SG1_fishers_q_n2373-n902_g18392"))
plotCorrelation(file.name, "Differental mutational expression", "Differental expression", "Differental mutation", x, y, pos="topright", line=2.75)

sig.mut <- rownames(subset(de.mut.gene, P < 1E-16))
sig.mut.tpm <- rownames(subset(de.tpm.gene[sig.mut,], P < 1E-16))
length(sig.mut.tpm)
de.mut.tpm[sig.mut.tpm,]


sig.mut <- rownames(subset(de.mut.gene, P <= 4.37075393850223E-07))
sig.mut.tpm <- rownames(subset(de.tpm.gene[sig.mut,], P <= 3.93766312113366E-17))

sig.mut.tpm.o <- intersect(rownames(de.mut.tpm), sig.mut.tpm)
writeTable(de.mut.tpm[sig.mut.tpm.o,], file.path(wd.driver.data, paste0(paste0("de_", base, "_mut-tpm-gene_SG1_fishers_q_n267"), ".txt")), colnames=T, rownames=F, sep="\t")


de.mut.tpm["ENSG00000125730",]



# -----------------------------------------------------------------------------
# Differential point mutations (SNV)
# Last Modified: 23/03/22
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.RData")
load("/Users/tpyang/Work/uni-koeln/tyang2/ICGC/analysis/expression/icgc-tpm-de/data/icgc_kallisto_0.42.1_tpm.gene.log2.m.RData")

tpm.gene <- tpm.gene[rownames(subset(tpm.gene.m, MEDIAN >= 1)),]
nrow(tpm.gene)
# [1] 13475   ## TPM > 1
# [1] 10260   ## TPM > 5

colnames <- intersect(colnames(mut.gene), colnames(tpm.gene))   ## 902 patients
mut.gene.tpm <- mut.gene[, colnames]
tpm.gene.mut <- tpm.gene[, colnames]
samples.mut.tpm <- samples.mut[colnames,]

rownames <- intersect(rownames(mut.gene), rownames(tpm.gene))   ## 18523 genes
mut.gene.tpm <- mut.gene.tpm[rownames,]
tpm.gene.mut <- tpm.gene.mut[rownames,]
dim(mut.gene.tpm)
# [1] 13475   904
dim(tpm.gene.mut)
# [1] 13475   904
dim(samples.mut.tpm)
# [1] 904  13

de.mut.tpm <- toTable(0, 21, length(rownames), c("MUT", "MUT_FREQ", "SAMPLE_FREQ", "MUT_G1", "MUT_S", "WT_G1", "WT_S", "G1MR", "SMR", "PMR", "P", "FDR", "MUT_G1_2", "MUT_S_2", "WT_G1_2", "WT_S_2", "G1MR_2", "SMR_2", "PMR_2", "P_2", "FDR_2"))
de.mut.tpm$MUT <- rownames
rownames(de.mut.tpm) <- de.mut.tpm$MUT
de.mut.tpm$MUT_FREQ <- mapply(x = 1:nrow(de.mut.tpm), function(x) sum(as.numeric(mut.gene.tpm[x,])))
de.mut.tpm.tmp <- de.mut.tpm
de.mut.tpm$SAMPLE_FREQ <- mapply(x = 1:nrow(de.mut.tpm), function(x) length(which(mut.gene.tpm[x,] != 0)))
de.mut.tpm.tmp <- de.mut.tpm
de.mut.tpm <- de.mut.tpm.tmp
de.mut.tpm <- subset(de.mut.tpm, SAMPLE_FREQ >= 20)
nrow(de.mut.tpm)
# [1] 12089   ## TPM > 1
# [1] 10260
# [1] 9257    ## TPM > 5

rownames <- de.mut.tpm$MUT   ## 13,475 genes
mut.gene.tpm <- mut.gene.tpm[rownames,]
tpm.gene.mut <- tpm.gene.mut[rownames,]

## Fishers' test
for (g in 1:nrow(de.mut.tpm)) {
   ids <- colnames(mut.gene.tpm[,which(mut.gene.tpm[g,] != 0)])
   
   samples.mut.mut <- subset(samples.mut.tpm, icgc_specimen_id %in% ids)
   samples.mut.wt  <- samples.mut.tpm[setdiff(rownames(samples.mut.tpm), rownames(samples.mut.mut)),]
  
   mut.g1 <- rownames(subset(samples.mut.mut, SORTING == "G1"))
   mut.s  <- rownames(subset(samples.mut.mut, SORTING == "S"))
   wt.g1  <- rownames(subset(samples.mut.wt,  SORTING == "G1"))
   wt.s   <- rownames(subset(samples.mut.wt,  SORTING == "S"))
      
   if ((length(mut.g1) > 10) && (length(mut.s) > 10) && (length(wt.g1) > 0) && (length(wt.s) > 0)) {
      test <- toTable(0, 2, 2, c("G1", "S"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- length(mut.g1)
      test[1, 2] <- length(mut.s)
      test[2, 1] <- length(wt.g1)
      test[2, 2] <- length(wt.s)
    
      de.mut.tpm$MUT_G1[g] <- length(mut.g1)
      de.mut.tpm$MUT_S[g]  <- length(mut.s)
      de.mut.tpm$WT_G1[g]  <- length(wt.g1)
      de.mut.tpm$WT_S[g]   <- length(wt.s)
      #de.mut.tpm$P[g]     <- fisher.test(test)[[1]]
      de.mut.tpm$P[g]      <- chisq.test(test)[[3]]
       
      de.mut.tpm$G1MR[g] <- de.mut.tpm$MUT_G1[g] / de.mut.tpm$WT_G1[g]
      de.mut.tpm$SMR[g]  <- de.mut.tpm$MUT_S[g]  / de.mut.tpm$WT_S[g]
      de.mut.tpm$PMR[g]  <- de.mut.tpm$SMR[g]    / de.mut.tpm$G1MR[g]
         
      ##
      test <- toTable(0, 2, 2, c("G1", "S"))
      rownames(test) <- c("MUT", "WT")
      test[1, 1] <- mean(as.numeric(tpm.gene.mut[g, mut.g1]))
      test[1, 2] <- mean(as.numeric(tpm.gene.mut[g, mut.s]))
      test[2, 1] <- mean(as.numeric(tpm.gene.mut[g, wt.g1]))
      test[2, 2] <- mean(as.numeric(tpm.gene.mut[g, wt.s]))
  
      de.mut.tpm$MUT_G1_2[g] <- mean(as.numeric(tpm.gene.mut[g, mut.g1]))
      de.mut.tpm$MUT_S_2[g]  <- mean(as.numeric(tpm.gene.mut[g, mut.s]))
      de.mut.tpm$WT_G1_2[g]  <- mean(as.numeric(tpm.gene.mut[g, wt.g1]))
      de.mut.tpm$WT_S_2[g]   <- mean(as.numeric(tpm.gene.mut[g, wt.s]))
      #de.mut.tpm$P[g]      <- fisher.test(test)[[1]]
      de.mut.tpm$P_2[g]      <- chisq.test(test)[[3]]
         
      de.mut.tpm$G1MR_2[g] <- de.mut.tpm$MUT_G1_2[g] / de.mut.tpm$WT_G1_2[g]
      de.mut.tpm$SMR_2[g]  <- de.mut.tpm$MUT_S_2[g]  / de.mut.tpm$WT_S_2[g]
      de.mut.tpm$PMR_2[g]  <- de.mut.tpm$SMR_2[g]    / de.mut.tpm$G1MR_2[g]
   }
}
de.mut.tpm <- de.mut.tpm[order(de.mut.tpm$PMR_2, decreasing=F),]
de.mut.tpm <- subset(de.mut.tpm, P   != 0)
de.mut.tpm <- subset(de.mut.tpm, P_2 != 0)
nrow(de.mut.tpm)
# [1] 12071   ## tpm > 0, s > 20
# [1] 8323    ## tpm > 1, s > 20, mut > 10
# [1] 6605    ## tpm > 5, s > 20, mut > 10

# [1] 17553   ## Unfiltered
# [1] 16320   ## s=20_1
# [1] 14606   ## s=5
# [1] 8279   ## s=20
# [1]    ## s=35

## Ensembl gene annotations
de <- de.mut.tpm
rownames(de) <- de$MUT
annot <- ensGene[,c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position", "gene_biotype")]
de.mut.tpm <- cbind(annot[rownames(de),], de[, -1])
#de.mut.tpm$FDR <- testFDR(de.mut.tpm$P, "Q")

file.name <- paste0("de_", base, "_mut-tpm_chisq_q_tpm>5_s>20_mut>10_g6605_mean")
save(de.mut.tpm,  file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(de.mut.tpm, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

###
## 27/04/22
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_MUT_TP53_EP300_MYC"))
plotPMR(file.name, "Proliferative mutational rate", "PMR", text.Log10.P, de.mut.tpm, c("TP53", "PIK3CA", "RB1", "NOTCH1", "EP300", "CREBBP", "TP73", "MYC", "MYCN", "MYCL"))

overlaps <- intersect(rownames(subset(de.mut.tpm, PMR >= 1)), rownames(de.tpm.gene))
file.name <- file.path(wd.driver.plots, paste0("Correlation_PMR_", base, "_mut.gene_chisq_tpm>5_s>20_mut>10_g6605_DE_S-G1_log2FC"))
x <- de.mut.tpm[overlaps,]$PMR
y <- de.tpm.gene[overlaps,]$LOG2_FC
plotCorrelation(file.name, "", "PMR", expression("Log" * ""[2] * " fold change"), x, y, pos="bottomright", col="dimgray", size=6)

file.name <- file.path(wd.driver.plots, paste0("Correlation_PMR_", base, "_mut.gene_chisq_tpm>5_s>20_mut>10_g6605_DE_S-G1_P"))
x <- de.mut.tpm[overlaps,]$PMR
y <- -log10(de.tpm.gene[overlaps,]$P)
plotCorrelation(file.name, "", "PMR", text.Log10.P, x, y, pos="bottomright", col="dimgray", size=6)

## Up
overlaps <- intersect(rownames(subset(de.mut.tpm, PMR >= 1)), rownames(subset(de.tpm.gene, LOG2_FC < 0)))
file.name <- file.path(wd.driver.plots, paste0("Correlation_PMR_", base, "_mut.gene_chisq_tpm>5_s>20_mut>10_g6605_DE_S-G1_log2FC_DOWN"))
x <- de.mut.tpm[overlaps,]$PMR
y <- de.tpm.gene[overlaps,]$LOG2_FC
plotCorrelation(file.name, "", "PMR", expression("Log" * ""[2] * " fold change"), x, y, pos="bottomright", col="dimgray", size=6)

file.name <- file.path(wd.driver.plots, paste0("Correlation_PMR_", base, "_mut.gene_chisq_tpm>5_s>20_mut>10_g6605_DE_S-G1_P_DOWN"))
x <- de.mut.tpm[overlaps,]$PMR
y <- -log10(de.tpm.gene[overlaps,]$P)
plotCorrelation(file.name, "", "PMR", text.Log10.P, x, y, pos="bottomright", col="dimgray", size=6)






###
##
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_TPM"))
x <- de.mut.tpm$PMR_2
y <- -log10(de.mut.tpm$P_2)
plotPMR(file.name, "", "EPMR", text.Log10.P, x, y)

file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_MUT_TP53_EP300_MYC"))
plotPMR(file.name, "Proliferative mutational rate", "PMR", text.Log10.P, de.mut.tpm, c("TP53", "PIK3CA", "RB1", "NOTCH1", "EP300", "CREBBP", "TP73", "MYC", "MYCN", "MYCL"))

file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_MUT_NPMR"))
x <- de.mut.tpm$G1MR/de.mut.tpm$SMR
y <- -log10(de.mut.tpm$P)
plotPMR(file.name, "", "NPMR", text.Log10.P, x, y)







###
##
overlaps <- intersect(rownames(de.mut.tpm), rownames(de.tpm.gene))
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene_chisq_tpm>5_s>20_mut>10_g6605_DE_log2FC"))
x <- de.tpm.gene[overlaps,]$LOG2_FC
y <- -log10(de.mut.tpm[overlaps,]$P)
plotCorrelation(file.name, "", expression("Log" * ""[2] * " fold change"), text.Log10.P, x, y)

file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene_chisq_tpm>5_s>20__mut>10_g8699_DE_P"))
x <- -log10(de.tpm.gene[overlaps,]$P)
y <- -log10(de.mut.tpm[overlaps,]$P)
plotCorrelation(file.name, "", text.Log10.P, text.Log10.P, x, y, pos="topright")





file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_TPM-N-of-sample"))
x <- de.mut.tpm$PMR_2
y <- de.mut.tpm$SAMPLE_FREQ
plotPMR(file.name, "", "Expressional PMR", "N of samples", x, y)

file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_N-of-sample_P"))
x <- -log10(de.mut.tpm$P)
y <- de.mut.tpm$SAMPLE_FREQ
plotPMR(file.name, "", text.Log10.P, "N of samples", x, y)






###
##
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene_chisq_s>20_tpm>1_mut>10_g8699"))
x <- de.mut.tpm$PMR
y <- -log10(de.mut.tpm$P)
plotCorrelation(file.name, "", "PMR", text.Log10.P, x, y, pos="topright", line=2.75)

##
overlaps <- intersect(rownames(de.mut.tpm), rownames(de.tpm.gene))
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene_chisq_s>20_tpm>1_mut>10_g8699_DE_log2FC"))
x <- de.tpm.gene[overlaps,]$LOG2_FC
y <- de.mut.tpm[overlaps,]$PMR
plotCorrelation(file.name, "", expression("Log" * ""[2] * " fold change"), "PMR [ratio]", x, y, pos="topright", line=2.75)

sig.mut <- rownames(subset(de.mut.tpm, PMR >= 2))
sig.mut.tpm <- rownames(subset(de.tpm.gene[sig.mut,], LOG2_FC >= 2))

file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_s>20_tpm>1_mut>10_g8699_sum"))
x <- de.mut.tpm$PMR_2
y <- -log10(de.mut.tpm$P_2)
plotCorrelation(file.name, "", "PMR", text.Log10.P, x, y, pos="topright", line=2.75)


##
de.mut.tpm.sig.tpm <- subset(de.mut.tpm.sig, MUT_G1_2 >= 3)
de.mut.tpm.sig.tpm <- subset(de.mut.tpm.sig.tpm, MUT_S_2 >= 3)
nrow(de.mut.tpm.sig.tpm)
# [1] 8602 TPM > 1
# [1] 7830 TPM > 2
# [1] 7259 TPM > 3





file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_mean_PMR-vs-PMR2"))
x <- de.mut.tpm$PMR_2
y <- de.mut.tpm$PMR
plotCorrelation(file.name, "", "Expressional PMR", "Mutational PMR", x, y, pos="topright", line=2.4)


overlaps <- intersect(rownames(de.mut.tpm), rownames(de.tpm.gene))
file.name <- file.path(wd.driver.plots, paste0("PMR_", base, "_mut.gene.tpm_chisq_tpm>5_s>20_mut>10_g6605_DE_log2FC"))
x <- de.tpm.gene[overlaps,]$LOG2_FC
y <- de.mut.tpm[overlaps,]$PMR
plotCorrelation(file.name, "", expression("Log" * ""[2] * " fold change"), "PMR [ratio]", x, y, pos="topright", line=2.75)







##
file.name <- paste0("de_", base, "_mut-tpm_SG1_chisq_q_n902_unfiltered")
load(file=file.path(wd.driver.data, paste0(file.name, ".RData")))
de.mut.tpm <- subset(de.mut.tpm, SAMPLE_FREQ >= 65)                                            ## Number of patients >= 65
nrow(de.mut.tpm)
# [1] 10927

de.mut.tpm <- filtered(de.mut.tpm, c("MUT_G1", "MUT_S", "WT_G1", "WT_S"), cutoff=10)           ## At least 10 patients in each test
nrow(de.mut.tpm)
# [1] 10735

#de.mut.tpm <- filtered(de.mut.tpm, c("MUT_G1_2", "MUT_S_2", "WT_G1_2", "WT_S_2"), cutoff=1)   ## Expressed genes (TPM >= 5)
overlaps <- intersect(rownames(de.mut.tpm), genes.tpm1)
de.mut.tpm <- de.mut.tpm[overlaps,]
nrow(de.mut.tpm)
# [1] 8142
de.mut.tpm$FDR   <- testFDR(de.mut.tpm$P,   "Q")
de.mut.tpm$FDR_2 <- testFDR(de.mut.tpm$P_2, "Q")
nrow(de.mut.tpm)
# [1] 10927   ## N > 65
# [1] 5497    ## N > 65; TPM > 5

file.name <- paste0("de_", base, "_mut-tpm_SG1_chisq_q_n902_p>65_p>10_tpm>1_g8142")
save(de.mut.tpm,  file=file.path(wd.driver.data, paste0(file.name, ".RData")))
writeTable(de.mut.tpm, file.path(wd.driver.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")

x <- -log10(de.mut.tpm$P_2)
y <- -log10(de.mut.tpm$P)
file.name <- file.path(wd.driver.data, paste0("correlation_", base, "_mut-tpm_SG1_chisq_q_n902_p>65_p>10_tpm>1_g8142"))
plotCorrelation(file.name, "Differental mutational expression", "Differental expression", "Differental mutation", x, y, pos="topright", line=2.75)

##
overlaps <- intersect(rownames(de.mut.tpm), rownames(de.mut.gene))
x <- -log10(de.mut.tpm[overlaps,]$P_2)
y <- -log10(de.mut.gene[overlaps,]$P)
file.name <- file.path(wd.driver.data, paste0("correlation_", base, "_mut-tpm_SG1_chisq_q_n902_p>65_g10927_n2373mut"))
plotCorrelation(file.name, "Differental mutational expression", "Differental expression", "Differental mutation", x, y, pos="topright", line=2.75)







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

