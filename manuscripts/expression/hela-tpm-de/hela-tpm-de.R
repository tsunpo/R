# =============================================================================
# Manuscript   : 
# Chapter I    : RB1-loss differential gene expression in neuroendocrine tumours
# Name         : manuscripts/expression/hela-tpm-de.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/05/19; 09/08/18
# =============================================================================
#wd.src <- "/projects/cangen/tyang2/dev/R"        ## tyang2@cheops
#wd.src <- "/ngs/cangen/tyang2/dev/R"             ## tyang2@gauss
wd.src <- "/Users/tpyang/Work/dev/R"              ## tpyang@localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for the manuscript
handbooks  <- c("Common.R", "DifferentialExpression.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Human Genome
load(file.path(wd.src.ref, "hg19.RData"))

# -----------------------------------------------------------------------------
# Set working directory
# -----------------------------------------------------------------------------
#wd <- "/ngs/cangen/tyang2"                   ## tyang2@gauss
wd <- "/Users/tpyang/Work/uni-koeln/tyang2"   ## tpyang@localhost

BASE <- "HeLa"
base <- tolower(BASE)
wd.rna   <- file.path(wd, BASE, "ngs/RNA")
wd.anlys <- file.path(wd, BASE, "analysis")
wd.meta  <- file.path(wd, BASE, "metadata/Dominguez 2016")

wd.de       <- file.path(wd.anlys, "expression/kallisto", paste0(base, "-tpm-de"))
wd.de.data  <- file.path(wd.de, "data")
wd.de.plots <- file.path(wd.de, "plots")
wd.de.path  <- file.path(wd.de, "pathway")

samples <- readTable(file.path(wd.rna, "hela_rna_n14.txt"), header=T, rownames=T, sep="\t")
rownames(samples) <- gsub("D", "T", rownames(samples))
samples$SAMPLE_ID <- gsub("D", "T", samples$SAMPLE_ID)

load(file.path(wd, base, "analysis/expression/kallisto", paste0(base, "-tpm-de/data/", base, "_kallisto_0.43.1_tpm.gene_r5p47.RData")))
colnames(tpm.gene) <-  gsub("D", "T", colnames(tpm.gene))
tpm.gene <- tpm.gene[,rownames(samples)]

# -----------------------------------------------------------------------------
# PCA on G1-S/G2-M genes
# Last Modified: 18/05/19
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

plotCellCyclePCA <- function(wd.de.plots, BASE, tpm.gene, samples, trait, genes.G1S, genes.G2M, isID, cols, flip.x, flip.y) {
   genes.G1S <- intersect(genes.G1S, rownames(tpm.gene))
   genes.G2M <- intersect(genes.G2M, rownames(tpm.gene))
   traits <- samples[,trait]
   ids <- NA
   if (isID) ids <- rownames(samples)
   
   ##
   test <- tpm.gene[genes.G1S, rownames(samples)]
   pca.de <- getPCA(t(test))
  
   file.main <- paste0(BASE, " on ", length(genes.G1S), "/304 G1-S genes")
   plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G1-S"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G1-S"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G1-S"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   
   ##
   test <- tpm.gene[genes.G2M, rownames(samples)]
   pca.de <- getPCA(t(test))
   
   file.main <- paste0(BASE, " on ", length(genes.G2M), "/876 G2-M genes")
   plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G2-M"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G2-M"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
   plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", tolower(BASE), "_G2-M"), size=6.5, file.main, "topleft", cols, ids, flip.x, flip.y)
}


# -----------------------------------------------------------------------------
# PCA of HeLa on G1-S/G2-M genes
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
plotCellCyclePCA(wd.de.plots, "HeLa", tpm.gene.hela, samples.hela, "CYCLE_PCA_2", genes.G1S, genes.G2M, isID=T, c("purple3", "forestgreen", "gold"), 1, 1)

## SCLC
overlaps <- intersect(colnames(tpm.gene.sclc), rownames(samples.sclc))
samples.sclc <- samples.sclc[overlaps,]
plotCellCyclePCA(wd.de.plots, "SCLC", tpm.gene.sclc, samples.sclc, "G1S", genes.G1S, genes.G2M, isID=F, c("blue", "gray", "red"), 1, 1)

## LCNEC
samples.lcnec$PCA <- NA
samples.lcnec$PCA[which(samples.lcnec$RB1_MUT == 0)] <- "WT"
samples.lcnec$PCA[which(samples.lcnec$RB1_MUT == 1)] <- "RB1"
plotCellCyclePCA(wd.de.plots, "LCNEC", tpm.gene.lcnec, samples.lcnec, "PCA", genes.G1S, genes.G2M, isID=F, c("gray", "red", "blue"), 1, 1)

## NBL
plotCellCyclePCA(wd.de.plots, "NBL", tpm.gene.nbl, samples.nbl, "G1S", genes.G1S, genes.G2M, isID=F, c("blue", "gray", "red"), 1, 1)

###
## Combine HeLa and SCLC
samples.hela$PCA <- samples.hela$CYCLE_PCA_2
samples.sclc$PCA <- paste0("SCLC (", samples.sclc$G1S, ")")
samples <- rbind(samples.hela[,c("SAMPLE_ID", "PCA")], samples.sclc[,c("SAMPLE_ID", "PCA")])

overlaps <- intersect(rownames(tpm.gene.hela), rownames(tpm.gene.sclc))
tpm.gene <- cbind(tpm.gene.hela[overlaps,], tpm.gene.sclc[overlaps,])
tpm.gene <- tpm.gene[,samples$SAMPLE_ID]
#plotCellCyclePCA(wd.de.plots, "SCLC and HeLa", tpm.gene.pca, samples.pca, "PCA", genes.G1S, genes.G2M, isID=F, c("purple3", "forestgreen", "gold", "blue", "gray", "red"), -1, 1)

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene))
traits <- samples$PCA
cols <- c("purple3", "forestgreen", "gold", "blue", "gray", "red")
BASE2 <- "SCLC and HeLa"
base2 <- "sclc+hela"

##
test <- tpm.gene[genes.G1S, rownames(samples)]
pca.de <- getPCA(t(test))

file.main <- paste0(BASE2, " on ", length(genes.G1S), "/304 G1-S genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "topleft", cols, NA, -1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "topleft", cols, NA, -1, -1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "topleft", cols, NA, 1, -1)
plotPCA(4, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "topleft", cols, NA, 1, -1)

##
test <- tpm.gene[genes.G2M, rownames(samples)]
pca.de <- getPCA(t(test))

file.main <- paste0(BASE2, " on ", length(genes.G2M), "/876 G2-M genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, 1, 1)

###
## Combine HeLa and LCNEC
samples.lcnec$PCA <- paste0("LCNEC (", samples.lcnec$PCA, ")")
samples <- rbind(samples.hela[,c("SAMPLE_ID", "PCA")], samples.lcnec[,c("SAMPLE_ID", "PCA")])

overlaps <- intersect(rownames(tpm.gene.hela), rownames(tpm.gene.lcnec))
tpm.gene <- cbind(tpm.gene.hela[overlaps,], tpm.gene.lcnec[overlaps,])
tpm.gene <- tpm.gene[,samples$SAMPLE_ID]
#plotCellCyclePCA(wd.de.plots, "SCLC and HeLa", tpm.gene.pca, samples.pca, "PCA", genes.G1S, genes.G2M, isID=F, c("purple3", "forestgreen", "gold", "blue", "gray", "red"), -1, 1)

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene))
traits <- samples$PCA
cols <- c("purple3", "gray", "red", "blue", "forestgreen", "gold")
BASE2 <- "LCNEC and HeLa"
base2 <- "lcnec+hela"

##
test <- tpm.gene[genes.G1S, rownames(samples)]
pca.de <- getPCA(t(test))

file.main <- paste0(BASE2, " on ", length(genes.G1S), "/304 G1-S genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "bottomleft", cols, NA, 1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "bottomleft", cols, NA, 1, -1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "bottomleft", cols, NA, 1, -1)

##
test <- tpm.gene[genes.G2M, rownames(samples)]
pca.de <- getPCA(t(test))

file.main <- paste0(BASE2, " on ", length(genes.G2M), "/876 G2-M genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, 1, 1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, 1, 1)

###
## Combine HeLa and NBL
samples.nbl$PCA <- paste0("NBL (", samples.nbl$G1S, ")")
samples.nbl$SAMPLE_ID <- samples.nbl$V1
samples <- rbind(samples.hela[,c("SAMPLE_ID", "PCA")], samples.nbl[,c("SAMPLE_ID", "PCA")])

overlaps <- intersect(rownames(tpm.gene.hela), rownames(tpm.gene.nbl))
tpm.gene <- cbind(tpm.gene.hela[overlaps,], tpm.gene.nbl[overlaps,])
tpm.gene <- tpm.gene[,samples$SAMPLE_ID]
#plotCellCyclePCA(wd.de.plots, "SCLC and HeLa", tpm.gene.pca, samples.pca, "PCA", genes.G1S, genes.G2M, isID=F, c("purple3", "forestgreen", "gold", "blue", "gray", "red"), -1, 1)

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene))
traits <- samples$PCA
cols <- c("purple3", "blue", "gray", "red", "forestgreen", "gold")
BASE2 <- "NBL and HeLa"
base2 <- "nbl+hela"

##
test <- tpm.gene[genes.G1S, rownames(samples)]
pca.de <- getPCA(t(test))

file.main <- paste0(BASE2, " on ", length(genes.G1S), "/304 G1-S genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "bottomleft", cols, NA, 1, 1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "bottomleft", cols, NA, 1, 1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G1-S"), size=6.5, file.main, "bottomleft", cols, NA, 1, 1)

##
test <- tpm.gene[genes.G2M, rownames(samples)]
pca.de <- getPCA(t(test))

file.main <- paste0(BASE2, " on ", length(genes.G2M), "/876 G2-M genes")
plotPCA(1, 2, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, -1, -1)
plotPCA(1, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, -1, 1)
plotPCA(2, 3, pca.de, traits, wd.de.plots, paste0("pca_", base2, "_G2-M"), size=6.5, file.main, "topleft", cols, NA, -1, 1)








# -----------------------------------------------------------------------------
# PCA (on G1-S/G2-M genes)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.hela))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.hela))
# > length(genes.G1S)
# [1] 279   ## Original 304 from Dominguez et al.
# > length(genes.G2M)
# [1] 838   ## Original 876 from Dominguez et al.

test <- tpm.gene[genes.G1S,]
#test <- tpm.gene[genes.G2M,]
pca.de <- getPCA(t(test))
trait <- samples.hela$CYCLE_PCA_2

##
file.main <- "HeLa (n=14) on 279/304 G1-S genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_279G1-S", size=6.5, file.main, "topleft", , rownames(samples), flip.x=1, flip.y=1)

file.main <- "HeLa (n=14) on 838/876 G2-M genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_838G2-M", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# PCA of SCLC on RT (Q4 to Q1)
# Last Modified: 15/05/19
# -----------------------------------------------------------------------------
## Go to sclc-tpm-de.R
tpm.gene.sclc <- tpm.gene
samples.sclc <- samples

## Come back to hela-tpm-de.R
overlaps <- intersect(rownames(tpm.gene.hela), rownames(tpm.gene.sclc))
tpm.gene <- cbind(tpm.gene.hela[overlaps,], tpm.gene.sclc[overlaps,])

samples.hela$PCA <- samples.hela[,"CYCLE_PCA_2"]
#samples.sclc$PCA <- paste0("Q", samples.sclc$Q4)
#samples.sclc$PCA <- as.numeric(samples.sclc[,"RT"])
#samples.sclc$PCA[which(samples.sclc$PCA == 1)] <- "M1"
#samples.sclc$PCA[which(samples.sclc$PCA == 2)] <- "M2"
samples.sclc$PCA <- samples.sclc$Q4
samples.sclc$PCA[which(samples.sclc$PCA == 1)] <- "SCLC (non-G1)"
samples.sclc$PCA[which(samples.sclc$PCA == 2)] <- "SCLC (G1)"
samples.sclc$PCA[which(samples.sclc$PCA == 3)] <- "SCLC (G1)"
samples.sclc$PCA[which(samples.sclc$PCA == 4)] <- "SCLC (non-G1)"

samples.pca <- rbind(samples.hela[,c("SAMPLE_ID", "PCA")], samples.sclc[,c("SAMPLE_ID", "PCA")])
tpm.gene <- tpm.gene[,samples.pca$SAMPLE_ID]

##
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene))
# > length(genes.G1S)
# [1] 277   ## Original 304 from Dominguez et al.
# > length(genes.G2M)
# [1] 830   ## Original 876 from Dominguez et al.

genes.G1S.0.6 <- rownames(subset(ensGene[genes.G1S,], chromosome_name %in% cors.sclc.0.6$chr))
genes.G2M.0.6 <- rownames(subset(ensGene[genes.G2M,], chromosome_name %in% cors.sclc.0.6$chr))
# > length(genes.G1S.0.6)
# [1] 110
# > length(genes.G2M.0.6)
# [1] 354

test <- tpm.gene[genes.G1S.0.6,]
#test <- tpm.gene[genes.G2M.0.6,]
pca.de <- getPCA(t(test))
trait <- samples.pca$PCA

##
file.main <- "Hela and SCLC on 110/277/304 G1-S genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_hela+sclc_0.6_G1_110G1-S", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold", "blue", "red"), NA, flip.x=1, flip.y=1)

file.main <- "Hela and SCLC on 354/830/876 G2-M genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_hela+sclc_0.6_G1_354G2-M", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold", "blue", "red"), NA, flip.x=1, flip.y=1)

##
file.main <- "Hela and SCLC on 277/304 G1-S genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela+sclc_M2_277G1-S", size=6.5, file.main, "topleft", c("purple3", "blue", "red", "forestgreen", "gold"), NA, flip.x=1, flip.y=1)

file.main <- "Hela and SCLC on 830/876 G2-M genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela+sclc_M2_830G2-M", size=6.5, file.main, "topleft", c("purple3", "blue", "red", "forestgreen", "gold"), NA, flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# PCA of SCLC on RT (G1 vs. non-G1)
# Last Modified: 17/05/19
# -----------------------------------------------------------------------------
load("/Users/tpyang/Work/uni-koeln/tyang2/SCLC/analysis/replication/sclc-wgs-rt/data/rt-vs-rt_SCLC-M2-M1-vs-LCL-S-G1_cors-pearson.RData")
cors$chr <- paste0("chr", cors$chr)
cors.sclc <- cors
cors.sclc.0.6 <- subset(cors, cor >= 0.6)

table.g1s.sclc <- as.data.frame(table(ensGene[genes.G1S,]$chromosome_name))
table.g1s.sclc <- table.g1s.sclc[order(table.g1s.sclc$Freq, decreasing=T),]

genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.sclc))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.sclc))

genes.G1S.0.6 <- rownames(subset(ensGene[genes.G1S,], chromosome_name %in% cors.sclc.0.6$chr))
genes.G2M.0.6 <- rownames(subset(ensGene[genes.G2M,], chromosome_name %in% cors.sclc.0.6$chr))
# > length(genes.G1S.0.6)
# [1] 111
# > length(genes.G2M.0.6)
# [1] 357

test <- tpm.gene.sclc[genes.G1S.0.6, samples.sclc$SAMPLE_ID]
#test <- tpm.gene.sclc[genes.G2M.0.6, samples.sclc$SAMPLE_ID]
pca.de <- getPCA(t(test))
#trait <- samples.nbl$PCA
#trait <- samples.sclc$Q4
#trait[which(trait == 1)] <- "non-G1"
#trait[which(trait == 2)] <- "G1"
#trait[which(trait == 3)] <- "G1"
#trait[which(trait == 4)] <- "non-G1"
trait <- as.numeric(samples.nbl[,"RT"])
trait[which(trait == 1)] <- "M1"
trait[which(trait == 2)] <- "M2"

##
file.main <- "SCLC on 111/279/304 G1-S genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_sclc_cor0.6_RT_111G1-S", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)

file.main <- "SCLC on 357/834/876 G2-M genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_sclc_cor0.6_RT_357G2-M", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)



# -----------------------------------------------------------------------------
# PCA of SCLC on RT (G1 vs. non-G1)
# Last Modified: 17/05/19
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.sclc))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.sclc))
# > length(genes.G1S)
# [1] 279   ## Original 304 from Dominguez et al.
# > length(genes.G2M)
# [1] 834   ## Original 876 from Dominguez et al.

test <- tpm.gene.sclc[genes.G1S, samples.sclc$SAMPLE_ID]
#test <- tpm.gene.sclc[genes.G2M, samples.sclc$SAMPLE_ID]
pca.de <- getPCA(t(test))
#trait <- samples.nbl$PCA
trait <- as.numeric(samples.nbl[,"RT"])
trait[which(trait == 1)] <- "M1"
trait[which(trait == 2)] <- "M2"
#trait <- samples.sclc$Q4
#trait[which(trait == 1)] <- "non-G1"
#trait[which(trait == 2)] <- "G1"
#trait[which(trait == 3)] <- "G1"
#trait[which(trait == 4)] <- "non-G1"

##
file.main <- "SCLC on 279/304 G1-S genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_sclc_g1-like_279G1-S", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)

file.main <- "SCLC on 834/876 G2-M genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_sclc_g1-like_834G2-M", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# PCA SCLC on LCNEC RB1 D.E. genes (Q < 0.05/0.1)
# Last Modified: 11/05/19
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.q0.1.sclc  <- intersect(genes.rb1.q0.1, rownames(tpm.gene.sclc))   ## Do NOT use filtered tpm.gene
genes.rb1.q0.05.sclc <- intersect(genes.rb1.q0.05, rownames(tpm.gene.sclc))
## > length(genes.rb1.q0.1.sclc)
## [1] 621   ## Original 639 from lcnec-tpm-de-pca.R
## > length(genes.rb1.q0.05.sclc)
## [1] 143   ## Original 145 from lcnec-tpm-de-pca.R

## NBL on RB1 D.E genes
test <- tpm.gene.sclc[genes.rb1.q0.05.sclc, samples.sclc$SAMPLE_ID]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene.sclc[genes.rb1.q0.1.sclc, samples.sclc$SAMPLE_ID]
pca.de <- getPCA(t(test))
trait <- samples.sclc$Q4
trait[which(trait == 1)] <- "non-G1"
trait[which(trait == 2)] <- "G1"
trait[which(trait == 3)] <- "G1"
trait[which(trait == 4)] <- "non-G1"

##
file.main <- "SCLC on 139/145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_sclc_RB1_Q0.05_143DE", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)

file.main <- "SCLC on 621/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(1, 3, pca.de, trait, wd.de.plots, "pca_sclc_RB1_Q0.1_621DE", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)








# -----------------------------------------------------------------------------
# PCA of NBL on RT (Q4 to Q1)
# Last Modified: 17/05/19
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.nbl))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.nbl))
# > length(genes.G1S)
# [1] 276   ## Original 304 from Dominguez et al.
# > length(genes.G2M)
# [1] 826   ## Original 876 from Dominguez et al.

#samples.nbl$PCA <- paste0("Q", samples.nbl$Q4)
test <- tpm.gene.nbl[genes.G1S, samples.nbl$V1]
#test <- tpm.gene.nbl[genes.G2M, samples.nbl$V1]
pca.de <- getPCA(t(test))
#trait <- samples.nbl$PCA
trait <- as.numeric(samples.nbl[,"RT"])
trait[which(trait == 1)] <- "M1"
trait[which(trait == 2)] <- "M2"

##
file.main <- "NBL on 276/304 G1-S genes"
plotPCA(1, 3, pca.de, trait, wd.de.plots, "pca_nbl_m2_276G1-S", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)

file.main <- "NBL on 826/876 G2-M genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_nbl_m2_826G2-M", size=6.5, file.main, "topleft", c("blue", "red"), NA, flip.x=1, flip.y=1)

# -----------------------------------------------------------------------------
# PCA of LCNEC on RT (RB1 and WT)
# Last Modified: 17/05/19
# -----------------------------------------------------------------------------
genes.G1S <- readTable(file.path(wd.meta, "Dominguez_G1-S.list"), header=F, rownames=F, sep="")
genes.G2M <- readTable(file.path(wd.meta, "Dominguez_G2-M.list"), header=F, rownames=F, sep="")

genes.G1S <- intersect(genes.G1S, rownames(tpm.gene.lcnec))
genes.G2M <- intersect(genes.G2M, rownames(tpm.gene.lcnec))
# > length(genes.G1S)
# [1] 279   ## Original 304 from Dominguez et al.
# > length(genes.G2M)
# [1] 833   ## Original 876 from Dominguez et al.

test <- tpm.gene.lcnec[genes.G1S, samples.lcnec$SAMPLE_ID]
#test <- tpm.gene.lcnec[genes.G2M,]
pca.de <- getPCA(t(test))
trait <- as.numeric(samples.lcnec[,"RB1_MUT"])
trait[which(trait == 0)] <- "LCNEC (WT)"
trait[which(trait == 1)] <- "LCNEC (RB1)"
trait[which(is.na(trait))] <- "LCNEC (NA)"

##
file.main <- "LCNEC on 279/304 G1-S genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_lcnec_279G1-S", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NA, flip.x=1, flip.y=1)

file.main <- "LCNEC on 833/876 G2-M genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_lcnec_833G2-M", size=6.5, file.main, "topleft", c("gray", "red", "dodgerblue"), NA, flip.x=1, flip.y=1)







# -----------------------------------------------------------------------------
# PCA (on NBL RB1 D.E. genes; Q < 0.05/0.1)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.q0.1  <- intersect(genes.rb1.q0.1, rownames(tpm.gene.nbl))   ## Do NOT use filtered tpm.gene
genes.rb1.q0.05 <- intersect(genes.rb1.q0.05, rownames(tpm.gene.nbl))
## > length(genes.rb1.q0.1)
## [1] 601   ## Original 639 from lcnec-tpm-de-pca.R
## > length(genes.rb1.q0.05)
## [1] 139   ## Original 145 from lcnec-tpm-de-pca.R

## NBL on RB1 D.E genes
test <- tpm.gene.nbl[genes.rb1.q0.05,]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene.nbl[genes.rb1.q0.1,]
pca.de <- getPCA(t(test))
trait <- samples.nbl$PCA

##
file.main <- "NBL on 139/145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(1, 3, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.05_135DE", size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NA, flip.x=1, flip.y=1)

file.main <- "NBL on 601/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.1_578DE", size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NA, flip.x=1, flip.y=-1)








## HeLa on 
test <- tpm.gene
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

##
file.main <- "HeLa (n=14) on 135/145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.05_135DE", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=1)

file.main <- "HeLa (n=14) on 578/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.1_578DE", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=-1)

##
## > length(genes.rb1.q0.1)
## [1] 573   ## Original 639 from lcnec-tpm-de-pca.R
## > length(genes.rb1.q0.05)
## [1] 133   ## Original 145 from lcnec-tpm-de-pca.R

## HeLa on RB1 D.E genes
test <- tpm.gene.sclc[genes.rb1.q0.05, rownames(samples.sclc)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene.sclc[genes.rb1.q0.1, rownames(samples.sclc)]
pca.de <- getPCA(t(test))
trait <- samples.sclc$PCA

##
file.main <- "SCLC on 133/145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_sclc_RB1_Q0.05_133DE", size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NA, flip.x=1, flip.y=1)

file.main <- "SCLC on 573/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(2, 3, pca.de, trait, wd.de.plots, "pca_sclc_RB1_Q0.1_573DE", size=6.5, file.main, "topleft", c("blue", "skyblue3", "lightcoral", "red"), NA, flip.x=1, flip.y=-1)








# -----------------------------------------------------------------------------
# PCA (on LCNEC RB1 D.E. genes; Q < 0.05/0.1)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
base.lcnec <- "LCNEC"
load(file.path(wd, base.lcnec, "analysis/expression/kallisto", paste0(base.lcnec, "-tpm-de/data/", "de_", base.lcnec, "_tpm-gene-r5p47_rb1_wilcox_q_n54.RData")))

genes.rb1.q0.1  <- rownames(subset(de.tpm.gene, FDR <= 0.1))
genes.rb1.q0.05 <- rownames(subset(de.tpm.gene, FDR <= 0.05))
genes.rb1.q0.1  <- intersect(genes.rb1.q0.1, rownames(tpm.gene))   ## Do NOT use filtered tpm.gene
genes.rb1.q0.05 <- intersect(genes.rb1.q0.05, rownames(tpm.gene))
## > length(genes.rb1.q0.1)
## [1] 578   ## Original 639 from lcnec-tpm-de-pca.R
## > length(genes.rb1.q0.05)
## [1] 135   ## Original 145 from lcnec-tpm-de-pca.R

## HeLa on RB1 D.E genes
test <- tpm.gene[genes.rb1.q0.05, rownames(samples)]   ## BUG FIX 13/02/17: Perform PCA using normalised data (NOT log2-transformed)
#test <- tpm.gene[genes.rb1.q0.1, rownames(samples)]
pca.de <- getPCA(t(test))
trait <- samples$CYCLE_PCA

##
file.main <- "HeLa (n=14) on 135/145 D.E. (LCNEC RB1; Q < 0.05) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.05_135DE", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=1)

file.main <- "HeLa (n=14) on 578/639 D.E. (LCNEC RB1; Q < 0.1) genes"
plotPCA(1, 2, pca.de, trait, wd.de.plots, "pca_hela_RB1_Q0.1_578DE", size=6.5, file.main, "topleft", c("purple3", "forestgreen", "gold"), rownames(samples), flip.x=1, flip.y=-1)

# -----------------------------------------------------------------------------
# Gene length
# Last Modified: 10/08/18
# -----------------------------------------------------------------------------
initLength <- function(genes, group) {
   ens.genes <- ensGene[genes,]
   ens.genes$Length <- mapply(x = 1:nrow(ens.genes), function(x) getLengthTx(genes[x]))
   ens.genes$Group  <- group
   
   return(ens.genes)
}

plotBox <- function(gene, wd.de, expr.pheno.log2, pheno.all) {
   ensembl_gene_id <- subset(ensGene, external_gene_name == gene)$ensembl_gene_id[1]   ## To avoid ENSG00000269846 (one of RBL1)
   gene.tpms <- cbind(t(expr.pheno.log2)[rownames(pheno.all), ensembl_gene_id], pheno.all)
   colnames(gene.tpms)[1] <- "LOG2_TPM"
 
   pdf(paste0(wd.de, "plots/boxplot/boxplot_tpm.gene.log2_", gene, ".pdf"), height=6, width=4)
   ymin <- min(gene.tpms$LOG2_TPM)
   ymax <- max(gene.tpms$LOG2_TPM)
   boxplot(LOG2_TPM ~ Cancer_Type, data=gene.tpms, outline=T, names=c("LUAD", "LCNEC", "SCLC"), ylim=c(ymin, ymax), ylab="log2(TPM + 0.01)", main=paste0(gene, " (", ensembl_gene_id, ")"))
 
   dev.off()
}

ens.genes.rb1.q0.1  <- initLength(genes.rb1.q0.1, 0)
ens.genes.rb1.q0.05 <- initLength(genes.rb1.q0.05, 1)
ens.genes.ALL <- initLength(rownames(tpm.gene), 2)   ## ADD 17/08/18
ens.genes.G1S <- initLength(genes.G1S, 3)
ens.genes.G2M <- initLength(genes.G2M, 4)
ens.genes <- rbind(ens.genes.rb1.q0.1, ens.genes.rb1.q0.05, ens.genes.ALL, ens.genes.G1S, ens.genes.G2M)
ens.genes$Group <- as.factor(ens.genes$Group)
ens.genes$Length <- log10(ens.genes$Length)

file.name <- file.path(wd.de.plots, paste0("boxplot_", base, "_genes_RB1-G1S-G2M-ALL_length.pdf"))
pdf(file.name, height=6, width=4.5)
ymin <- min(ens.genes$Length)
ymax <- max(ens.genes$Length)
boxplot(Length ~ Group, data=ens.genes, outline=T, names=c("RB1", "RB1", "All", "G1-S", "G2-M"), col=c("white", "white", "gray", "gray", "gray"), ylim=c(ymin, ymax), ylab="Gene length (log10)", main=c("LCNEC RB1 D.E. and HeLa cell-cycle", "gene lists"))
dev.off()

## G1-S vs G2-M
testW(ens.genes.G1S$Length, ens.genes.ALL$Length)
# [1] 0.831461
testW(ens.genes.G2M$Length, ens.genes.ALL$Length)
# [1] 9.266777e-32
testW(ens.genes.G1S$Length, ens.genes.G2M$Length)
# [1] 1.012973e-11

## RB1 (q<0.1) vs G1-S and G2-M
testW(ens.genes.rb1.q0.1$Length, ens.genes.ALL$Length)
# [1] 9.892606e-05
testW(ens.genes.rb1.q0.1$Length, ens.genes.G1S$Length)
# [1] 3.853642e-03
testW(ens.genes.rb1.q0.1$Length, ens.genes.G2M$Length)
# [1] 5.253057e-08

## RB1 (q<0.05) vs G1-S and G2-M
testW(ens.genes.rb1.q0.05$Length, ens.genes.ALL$Length)
# [1] 0.53314
testW(ens.genes.rb1.q0.05$Length, ens.genes.G1S$Length)
# [1] 0.3546051
testW(ens.genes.rb1.q0.05$Length, ens.genes.G2M$Length)
# [1] 4.708969e-06

# -----------------------------------------------------------------------------
# D.E. using non-parametric test: 5 First cycle vs 4 Second cycle (on D.E. LCNCE RB1 genes)
# Last Modified: 09/08/18
# -----------------------------------------------------------------------------
## Test: Wilcoxon/Wilcox/U/Students/ttest
## FDR : Q/BH
## D.E.: RB1_MUT (1) vs RB1_WT(0) as factor
samples$CYCLE_DE[13] <- NA   ## Preclude T13 which is similar to G1-S

argv      <- data.frame(predictor="CYCLE_DE", predictor.wt=0, test="Wilcox", test.fdr="Q", stringsAsFactors=F)
file.name <- paste0("de_", base, "_tpm-gene-135rb1_cycle2_wilcox_q_n9")
#file.name <- paste0("de_", base, "_tpm-gene-578rb1_cycle2_wilcox_q_n9")

de.tpm.gene <- pipeDE(tpm.gene.log2[genes.rb1.q0.05, rownames(samples)], samples, argv, ensGene)
#de.tpm.gene <- pipeDE(tpm.gene.log2[genes.rb1.q0.1, rownames(samples)], samples, argv, ensGene)

save(de.tpm.gene, file=file.path(wd.de.data, paste0(file.name, ".RData")))
writeTable(de.tpm.gene, file.path(wd.de.data, paste0(file.name, ".txt")), colnames=T, rownames=F, sep="\t")
