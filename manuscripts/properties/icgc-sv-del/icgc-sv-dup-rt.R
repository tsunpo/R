# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/driver/icgc-driver-cna.R
# Author       : Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 22/03/22
# =============================================================================
#wd.src <- "/nfs/users/nfs_t/ty2/dev/R"           ## @nfs
wd.src <- "/Users/ty2/Work/dev/R"                 ## @localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Asymmetry.R", "Mutation.R", "Survival.R", "Transcription.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
#wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## @lustre
wd <- "/Users/ty2/Work/sanger/ty2"               ## @localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.icgc     <- file.path(wd, BASE, "consensus")
wd.icgc.sv  <- file.path(wd.icgc, "sv")
wd.icgc.sv.plots <- file.path(wd.icgc.sv, "dup", "plots")

wd.meta     <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "properties", paste0(base, "-sv-dup"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

# -----------------------------------------------------------------------------
# Quantile skew plot
# Last Modified: 17/06/23
# -----------------------------------------------------------------------------
dups.all <- dels[0, 1:14]
for (h in 1:nrow(table.sv)) {
	  hist <- as.vector(table.sv$Var1[h])
	  load(file=file.path(wd.rt.data, "beds", paste0(hist, ".sv.dup.rt", 0, ".beds.RData")))

	  dups.all <- rbind(dups.all, dups[, 1:14])
}
#dups.all.auto <- subset(subset(dups.all, chrom1 != "X"), chrom1 != "Y")
#dups.all.auto <- subset(subset(dups.all, chrom2 != "X"), chrom2 != "Y")
#dups.all.nona <- dups.all[!is.na(dups.all$BED),]
dups.all$RT <- nrds.RT[dups.all$BED,]$RT
dups.all.nona <- dups.all[!is.na(dups.all$RT),]
dups.all.nona.2 <- subset(subset(dups.all.nona, RT >= -2), RT <= 2)
dups.all.nona.2.cut <- dups.all.nona[setdiff(rownames(dups.all.nona), rownames(dups.all.nona.2)),]
(67616 - 67330) / 67616
# [1] 0.004229768
(67330 - 66111) / 67330
# [1] 0.01810486
nrow(subset(dups.all.nona.2, RT > 0))
# [1] 40800
nrow(subset(dups.all.nona.2, RT < 0))
# [1] 25311
nrow(subset(dups.all.nona.2.cut, RT > 0))
# [1] 40
nrow(subset(dups.all.nona.2.cut, RT < 0))
# [1] 1179

##
overlaps <- intersect(dups.all$BED, nrds.RT.NRFD.lcl$BED)
dups.all.lcl$RT <- nrds.RT.NRFD.lcl[dups.all.lcl$BED,]$RT
dups.all.lcl.nona <- dups.all.lcl[!is.na(dups.all.lcl$RT),]
dups.all.lcl.nona.2 <- subset(subset(dups.all.lcl.nona, RT >= -2), RT <= 2)
dups.all.lcl.nona.2.cut <- dups.all.lcl.nona[setdiff(rownames(dups.all.lcl.nona), rownames(dups.all.lcl.nona.2)),]

## PCAWG deletions (ALL)
tests <- dups.all.nona.2$RT
randoms <- replicate(1000, getMonteCarloSimulations(nrds.RT.2, length(tests)))

#file.name <- file.path(wd.rt.plots, paste0("QS_RT_PCAWG.pdf"))
main.text <- "PCAWG DUP"
#xlab.text <- ""
#plotMonteCarloSimulation(tests, randoms, file.name, c(red, blue), main.text, xlab.text)

file.name <- file.path(wd.rt.plots, paste0("Density_RT_PCAWG-DUP.pdf"))
plotPropertyDensity(tests, randoms, file.name, c(red, blue), main.text, xlab.text)

## PCAWG deletions (< 10 kb)
tests <- subset(dups.all.nona.2, size < 50000)$RT
randoms <- replicate(1000, getMonteCarloSimulations(nrds.RT.2, length(tests)))

#file.name <- file.path(wd.rt.plots, paste0("QS_RT_PCAWG_<10kb.pdf"))
main.text <- "PCAWG DUP < 50 kb"
#xlab.text <- ""
#plotMonteCarloSimulation(tests, randoms, file.name, c(red, blue), main.text, xlab.text)

file.name <- file.path(wd.rt.plots, paste0("Density_RT_PCAWG-DUP<50kn.pdf"))
plotPropertyDensity(tests, randoms, file.name, c(red, blue), main.text, xlab.text)

## PCAWG deletions (> 10 kb)
tests <- subset(dups.all.nona.2, size > 50000)$RT
randoms <- replicate(1000, getMonteCarloSimulations(nrds.RT.2, length(tests)))

#file.name <- file.path(wd.rt.plots, paste0("QS_RT_PCAWG_>10kb.pdf"))
main.text <- "PCAWG DUP > 50 kb"
#xlab.text <- ""
#plotMonteCarloSimulation(tests, randoms, file.name, c(red, blue), main.text, xlab.text)

file.name <- file.path(wd.rt.plots, paste0("Density_RT_PCAWG-DUP>50kn.pdf"))
plotPropertyDensity(tests, randoms, file.name, c(red, blue), main.text, xlab.text)
