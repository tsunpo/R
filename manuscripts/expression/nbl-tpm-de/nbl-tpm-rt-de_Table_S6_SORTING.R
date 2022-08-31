convert <- function(x) {
   if (x == "SYMPATHOBLAST")
      x <- "Sympathoblasts"
   if (x == "CHROMAFFIN CELL")
      x <- "Chromaffin cells"
   if (x == "PRE-EMT NCC II")
      x <- "Pre-EMT NCC II"
   if (x == "PRE-EMT NCC I")
      x <- "Pre-EMT NCC I"
   #if (x == "EMT NCC")
   #   x <- "EMT neural crest cells"
   if (x == "NCC")
      x <- "Neural crest cells"
   if (x == "NT")
      x <- "Neural tube cells"
   if (x == "SCP")   
      x <- "Schwann cell precursors"
   if (x == "SN") 
      x <- "Sensory neurons"
 
   return(x)
}
#x <- "Schwann cell precursors"
#x <- "Sensory neurons"

table.S6 <- readTable(file.path(wd.de.gsea, "Table_S6.txt"), header=T, rownames=F, sep="\t")
writeTable(subset(table.S6, Cell.type == "Sympathoblasts")$Gene, file.path(wd.de.gsea, "Table_S6_Sympathoblasts.txt"), colnames=F, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# ALL
# Last Modified: 18/07/22
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "ALL_SORTING", "Table_S6", "gsea_report_for_na_pos_1661861532377.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[7:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "ALL_SORTING", "Table_S6_POS.pdf")
pdf(file.de, height=2.5, width=4)
par(mar=c(2,10.7,2,2.8))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "ALL_SORTING", "Table_S6", "gsea_report_for_na_neg_1661861532377.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[2:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "ALL_SORTING", "Table_S6_NEG.pdf")
pdf(file.de, height=0.85, width=4)
par(mar=c(2,4.2,0,9.3))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-6, 0, by=2), labels=F)
axis(side=1, at=seq(-6, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# LR
# Last Modified: 17/07/22
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "LR_SORTING", "Table_S6", "gsea_report_for_na_pos_1661863013854.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[2:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "LR_SORTING", "Table_S6_POS.pdf")
pdf(file.de, height=1.25, width=4)
par(mar=c(2,10.5,2,3))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "LR_SORTING", "Table_S6", "gsea_report_for_na_neg_1661863013854.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[7:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "LR_SORTING", "Table_S6_NEG.pdf")
pdf(file.de, height=2.15, width=4)
par(mar=c(2,2.8,0,9.7))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-5, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# HR
# Last Modified: 18/07/22
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "HR_SORTING", "Table_S8", "gsea_report_for_na_pos_1661352241204.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[10:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "HR_SORTING", "Table_S8_POS.pdf")
pdf(file.de, height=3.2, width=4)
par(mar=c(2,9,2,4.5))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 6, by=2), labels=F)
axis(side=1, at=seq(0, 6, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()










## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "HR", "Table_S8", "gsea_report_for_na_neg_1661291304495.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[1:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "HR", "Table_S8_NEG.pdf")
pdf(file.de, height=0.6, width=4)
par(mar=c(2,2.5,0,11))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()
