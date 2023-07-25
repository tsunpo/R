# =============================================================================
# Manuscript   : 
# Chapter      : 
# Name         : manuscripts/properties/icgc-hapmap-rt_chr.R
# Author       : Tsun-Po Yang (t2@sanger.ac.uk)
# Last Modified: 24/07/23
# =============================================================================
wd.src <- "/nfs/users/nfs_t/ty2/dev/R"           ## @nfs
#wd.src <- "/Users/ty2/Work/dev/R"                 ## @localhost

wd.src.lib <- file.path(wd.src, "handbook-of")    ## Required handbooks/libraries for this manuscript
handbooks  <- c("Commons.R", "Graphics.R", "GenomicProperty.R", "ReplicationForkDirectionality.R")
invisible(sapply(handbooks, function(x) source(file.path(wd.src.lib, x))))

wd.src.ref <- file.path(wd.src, "guide-to-the")   ## The Bioinformatician's Guide to the Genome
load(file.path(wd.src.ref, "hg19.RData"))
load(file.path(wd.src.ref, "hg19.bed.gc.1kb.RData"))

# -----------------------------------------------------------------------------
# Set up working directory
# Last Modified: 22/03/22
# -----------------------------------------------------------------------------
wd <- "/lustre/scratch127/casm/team294rr/ty2"   ## @lustre
#wd <- "/Users/ty2/Work/sanger/ty2"             ## @localhost
BASE <- "ICGC"
base <- tolower(BASE)

wd.meta <- file.path(wd, BASE, "metadata")

wd.anlys  <- file.path(wd, BASE, "analysis")
wd.rt <- file.path(wd.anlys, "properties", paste0(base, "-hapmap-rt"))
wd.rt.data  <- file.path(wd.rt, "data")
wd.rt.plots <- file.path(wd.rt, "plots")

# -----------------------------------------------------------------------------
# After running cmd_icgc-hapmap-rt_chr.R
# Last Modified: 24/07/23
# -----------------------------------------------------------------------------
map <- get(load(file.path(wd.rt.data, paste0("genetic_map_GRCh37_", chrs[1], ".RData"))))
for (c in 2:22) {
	  load(file.path(wd.rt.data, paste0("genetic_map_GRCh37_", chrs[c], ".RData")))
	  map <- rbind(map, map.chr)
}
nrow(map)
# [1] 3303900
map.nona <- map[!is.na(map$BED),]
nrow(map.nona)
# [1] 3303898

# -----------------------------------------------------------------------------
# RT
# Last Modified: 24/07/23
# -----------------------------------------------------------------------------
kb <- 20
rfd <- 0.9
load(file.path(wd.meta, "bstrps", paste0("nrds.RT.2_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))

map.nona$RT <- nrds.RT[map.nona$BED,]$RT
map.nona.e <- subset(map.nona, RT >= 0)
map.nona.l <- subset(map.nona, RT < 0)
nrow(map.nona.e)
nrow(map.nona.l)

##
FILE <- c("RATE", "MAP")
xlab.text <- c("Rate (cM/Mb)", "Map (cM)")
min <- c(0, NA)
max <- c(3, NA)

legends <- c("Early", "Late")
cols <- c(red, blue)
for (p in 1:2) {
	  main.text <- "Recombination (HapMap II)"
	  file.name <- file.path(wd.rt.plots, paste0("Density_Recombination_SCLC-NL_RT_", FILE[p], ".pdf"))
	  plotDensity2(map.nona.e[, FILE[p]], map.nona.l[, FILE[p]], file.name, cols, legends, main.text, xlab.text[p], min=min[p], max=max[p])
	  
	  #file.name <- file.path(wd.rt.plots, paste0("Cor_Recombination_", FILE[p], "_E.pdf"))
	  #plotCorrelation5(file.name, main.text, FILE[p], "RT (Early)", map.nona.e[, FILE[p]], map.nona.e$RT, pos="bottomright", cols=c(red.lighter, red))

	  #file.name <- file.path(wd.rt.plots, paste0("Cor_Recombination_", FILE[p], "_L.pdf"))
	  #plotCorrelation5(file.name, main.text, FILE[p], "RT (Late)", map.nona.l[, FILE[p]], map.nona.l$RT, pos="bottomright", cols=c(blue.lighter, blue))
}

# -----------------------------------------------------------------------------
# RFD
# Last Modified: 24/07/23
# -----------------------------------------------------------------------------
kb <- 20
load(file=file.path(wd.meta, "bstrps", paste0("sclc-nl_rpkm.gc.cn.d.rt.log2s.nrfd.", kb, "kb_", "m2-m1", ".RData")))
overlaps <- intersect(map.nona$BED, rownames(nrds.RT.NRFD))

map.nona$RT     <- nrds.RT.NRFD[map.nona$BED,]$RT
map.nona$SPLINE <- nrds.RT.NRFD[map.nona$BED,]$SPLINE
map.nona$SLOPE  <- nrds.RT.NRFD[map.nona$BED,]$SLOPE
map.nona$RFD    <- nrds.RT.NRFD[map.nona$BED,]$RFD
map.nona$NRFD   <- nrds.RT.NRFD[map.nona$BED,]$NRFD

save(map.nona, file=file.path(wd.rt.data, "Recombination_nrds.RT.NRFD_SCLC-NL.RData"))

##
load(file=file.path(wd.rt.data, "Recombination_nrds.RT.NRFD_SCLC-NL.RData"))
map.nona <- map.nona[!is.na(map.nona$RT), ]
nrow(map.nona)
# [1] 3299054
map.nona <- map.nona[!is.na(map.nona$RFD), ]
nrow(map.nona)
# [1] 3299054

nrds.ttr <- getBootstrapTTR(map.nona, 0.9)
nrds.ctr <- getBootstrapCTR(map.nona, 0.9)
nrds.ctr.iz <- subset(nrds.ctr, NRFD >= 0)
nrds.ctr.tz <- subset(nrds.ctr, NRFD < 0)

##
FILE <- c("RATE", "MAP")
xlab.text <- c("Rate (cM/Mb)", "Map (cM)")
min <- c(0, NA)
max <- c(3, NA)

legends <- c("TTR", "TZ", "IZ")
cols <- c("black", blue, red)
for (p in 1:2) {
	  main.text <- "Recombination (HapMap II)"
	  file.name <- file.path(wd.rt.plots, paste0("Density_Recombination_SCLC-NL_RFD_", FILE[p], ".pdf"))
	  plotDensity3(nrds.ttr[, FILE[p]], nrds.ctr.tz[, FILE[p]], nrds.ctr.iz[, FILE[p]], file.name, cols, legends, main.text, xlab.text[p], min=min[p], max=max[p])
	
	  #file.name <- file.path(wd.rt.plots, paste0("Cor_Recombination_RFD_RT-", FILE[p], "_TTR.pdf"))
	  #plotCorrelation5(file.name, main.text, FILE[p], "RT (TTR)", nrds.ttr[, FILE[p]], nrds.ttr$RT, pos="bottomright", cols=c("dimgray", "black"))
	
	  #file.name <- file.path(wd.rt.plots, paste0("Cor_Recombination_RFD_RT-", FILE[p], "_TZ.pdf"))
	  #plotCorrelation5(file.name, main.text, FILE[p], "RT (TZ)", nrds.ctr.tz[, FILE[p]], nrds.ctr.tz$RT, pos="bottomright", cols=c(blue.lighter, blue))
	  
	  #file.name <- file.path(wd.rt.plots, paste0("Cor_Recombination_RFD_RT-", FILE[p], "_IZ.pdf"))
	  #plotCorrelation5(file.name, main.text, FILE[p], "RT (IZ)", nrds.ctr.iz[, FILE[p]], nrds.ctr.iz$RT, pos="bottomright", cols=c(red.lighter, red))
}
