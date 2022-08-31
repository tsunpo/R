convert <- function(x) {
   if (x == "SYMPATHOBLAST")
      x <- "Sympathoblast"
   if (x == "CHROMAFFIN CELL")
      x <- "Chromaffin cells"
   if (x == "PRE-EMT NCC II")
      x <- "Pre-EMT NCC II"
   if (x == "PRE-EMT NCC I")
      x <- "Pre-EMT NCC I"
   if (x == "EMT NCC")
      x <- "EMT NCC"
   #if (x == "NCC")
   #   x <- "NCC"
   #if (x == "NT")
   #   x <- "Neural tube cells"
   #if (x == "SCP")   
   #   x <- "SCP"
   #if (x == "SN") 
   #   x <- "SN"
   return(x)
}
#x <- "Schwann cell precursors"
#x <- "Sensory neurons"

table.S6 <- readTable(file.path(wd.de.gsea, "Table_S6.txt"), header=T, rownames=F, sep="\t")
writeTable(subset(table.S6, Cell.type == "Sympathoblasts")$Gene, file.path(wd.de.gsea, "Table_S6_Sympathoblasts.txt"), colnames=F, rownames=F, sep="\t")

# -----------------------------------------------------------------------------
# RISK
# Last Modified: 18/07/22
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "RISK", "Table_S8", "gsea_report_for_na_pos_1661329853361.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[6:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "RISK", "Table_S8_POS.pdf")
pdf(file.de, height=2.2, width=4)
par(mar=c(2,11.5,2,2))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "RISK", "Table_S8", "gsea_report_for_na_neg_1661329853361.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[4:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "RISK", "Table_S8_NEG.pdf")
pdf(file.de, height=1.5, width=4)
par(mar=c(2.8,1.7,0,8.6))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-6, 0), xaxt="n", names.arg="", col=blue, xlab="Normalised enrichment score (NES)")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-6, 0, by=2), labels=F)
axis(side=1, at=seq(-6, 0, by=2), pos=0.5, tick=F)
mtext(main.text[2], line=0.3)
mtext("                                           Normalised enrichment score (NES)", side=1, line=1.9)
dev.off()

# -----------------------------------------------------------------------------
# ALL
# Last Modified: 18/07/22
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "ALL_PC2", "Table_S8", "gsea_report_for_na_pos_1661286115548.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[3:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "ALL_PC2", "Table_S8_POS.pdf")
pdf(file.de, height=1.5, width=4)
par(mar=c(2,11,2,2.5))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "ALL_PC2", "Table_S8", "gsea_report_for_na_neg_1661286115548.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[7:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "ALL_PC2", "Table_S8_NEG.pdf")
pdf(file.de, height=2.1, width=4)
par(mar=c(2,1.1,0,9.1))   # increase y-axis margin.

posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-6, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
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
gsea.pos <- readTable(file.path(wd.de.gsea, "LR_PC2", "Table_S8", "gsea_report_for_na_pos_1661286629387.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[5:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "LR_PC2", "Table_S8_POS.pdf")
pdf(file.de, height=2, width=4)
par(mar=c(2,11.3,2,2.2))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "LR_PC2", "Table_S8", "gsea_report_for_na_neg_1661286629387.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[5:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "LR_PC2", "Table_S8_NEG.pdf")
pdf(file.de, height=1.6, width=4)
par(mar=c(2,1.3,0,8.9))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-6, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-6, 0, by=2), labels=F)
axis(side=1, at=seq(-6, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

# -----------------------------------------------------------------------------
# HR
# Last Modified: 18/07/22
# -----------------------------------------------------------------------------
## POS
gsea.pos <- readTable(file.path(wd.de.gsea, "HR_PC2", "Table_S8", "gsea_report_for_na_pos_1661287584936.tsv"), header=T, rownames=F, sep="\t")
#gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) unlist(strsplit(gsea.pos$NAME[x], "_"))[2])
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) gsub("_", " ", gsea.pos$NAME[x]))
gsea.pos <- gsea.pos[3:1, c("NAME", "NES")]
gsea.pos$NAME <- mapply(x = 1:nrow(gsea.pos), function(x) convert(gsea.pos$NAME[x]))

main.text <- c("Dong et al. gene sets NES                        ", "")
file.de <- file.path(wd.de.gsea, "HR_PC2", "Table_S8_POS.pdf")
pdf(file.de, height=1.5, width=4)
par(mar=c(2,11,2,2.5))   # increase y-axis margin.
barplot(gsea.pos$NES, main=main.text[1], las=1, horiz=T, xlim=c(0, 4), xaxt="n", names.arg=gsea.pos$NAME, col=red.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
#abline(v=2, lty=5)
axis(side=1, at=seq(0, 4, by=2), labels=F)
axis(side=1, at=seq(0, 4, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()

## NEG
gsea.neg <- readTable(file.path(wd.de.gsea, "HR_PC2", "Table_S8", "gsea_report_for_na_neg_1661287584936.tsv"), header=T, rownames=F, sep="\t")
#gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) unlist(strsplit(gsea.neg$NAME[x], "HALLMARK_"))[2])
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) gsub("_", " ", gsea.neg$NAME[x]))
gsea.neg <- gsea.neg[7:1, c("NAME", "NES")]
gsea.neg$NAME <- mapply(x = 1:nrow(gsea.neg), function(x) convert(gsea.neg$NAME[x]))

main.text <- c("", "")
file.de <- file.path(wd.de.gsea, "HR_PC2", "Table_S8_NEG.pdf")
pdf(file.de, height=2.1, width=4)
par(mar=c(2,4.3,0,9.2))   # increase y-axis margin.
posbar <- barplot(gsea.neg$NES, main=main.text[1], las=1, horiz=T, xlim=c(-4, 0), xaxt="n", names.arg="", col=blue.lighter, xlab="")   #cex.names=0.8) cex.axis=1.1, cex.lab=1.15, cex.main=1.3
text(y=posbar, x=0, pos=4, labels=gsea.neg$NAME, xpd=NA)
#abline(v=-2, lty=5)
axis(side=1, at=seq(-4, 0, by=2), labels=F)
axis(side=1, at=seq(-4, 0, by=2), pos=0.5, tick=F)
#mtext(main.text[2], line=0.3)
dev.off()
