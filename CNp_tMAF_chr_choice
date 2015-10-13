MAF <- read.table("QS.Normal.Tumor.MAF.ref.alt.txt", header=TRUE, sep="\t")
CP <- read.table("QS2-tumor.bwa.realigned.rmDups.recal.bam_ratio.txt", header=TRUE, sep="\t")
PV <- read.table("QS2-tumor.bwa.realigned.rmDups.recal.bam_CNVs.p.value.txt", header=TRUE, sep="\t")

loh <- data.frame(MAF)
ratio <- data.frame(CP)
pvalue <- data.frame(PV)

ploidy <- 2

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}

for (i in c(17)) {
	ll <- which(loh$chromosome==(paste("chr", i, sep = "")))
	tt <- which(ratio$Chromosome==i)
	pG <- which(pvalue$chr==i & pvalue$status=="gain")
	pL <- which(pvalue$chr==i & pvalue$status=="loss")
	maxPosition <- max(loh$position[ll], ratio$Start[tt])
	
	png(filename = paste("PatientQS_CNp_tLOH_chr", i, ".png", sep = ""), h=2200, w=4000, pointsize=20)
	plot(1:10)
	layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
	par(mar = c(0, 0, 0, 0), oma = c(5, 4.5, 5, 2) + 0.1, xpd = NA)
	opt <- options(scipen = 20)

	if (length(tt)>0) {
		par(xpd=TRUE)
		plot(ratio$Start[tt], ratio$Ratio[tt]*ploidy, xlim=c(0,maxPosition), ylim = c(0,maxLevelToPlot*ploidy), xaxs="i", yaxs="i", yaxt="n", xaxt="n", xlab = "", ylab = "", pch = 20, col = colors()[88])
		axis(2, yaxp = c(1, maxLevelToPlot*ploidy, 5), cex.axis=1.5)
		tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
		points(ratio$Start[tt],ratio$Ratio[tt]*ploidy, xlim=c(0,maxPosition), pch = 20, col = colors()[136])
	
		tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
		points(ratio$Start[tt],ratio$Ratio[tt]*ploidy, xlim=c(0,maxPosition), pch = 20, col = colors()[136], cex=2)
	 
		tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
	 		points(ratio$Start[tt],ratio$Ratio[tt]*ploidy, xlim=c(0,maxPosition), pch = 20, col = colors()[461])
		segments(x0 = pvalue$start[pG], x1 = pvalue$end[pG], y0 = (5.5 + pvalue$copy.number[pG] - pvalue$copy.number[pG]), y1 = (5.5 + pvalue$copy.number[pG] - pvalue$copy.number[pG]), col="red", lwd = 2, lend = 3)
		segments(x0 = pvalue$start[pL], x1 = pvalue$end[pL], y0 = (0.5 + pvalue$copy.number[pL] - pvalue$copy.number[pL]), y1 = (0.5 + pvalue$copy.number[pL] - pvalue$copy.number[pL]), col="blue", lwd = 2, lend = 3)
        par(xpd=NA)
		title(ylab="normalized copy number profile", cex.lab=2)
		par(xpd=TRUE)
		}
	par(xpd=TRUE)
	if (length(ll)>0) {
		###plot(loh$position[ll], loh$nMAFref[ll], col='#008837', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xlab="", ylab="", cex.axis=1.25, pch=".")
		###par(new=T)
		###plot(loh$position[ll], loh$nMAFalt[ll], col='#008837', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
		###par(new=T)
		plot(loh$position[ll], loh$tMAFref[ll], col='#7b3294', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xlab="", ylab="", cex.axis=1.25, pch=".")
		par(new=T)
		plot(loh$position[ll], loh$tMAFalt[ll], col='#7b3294', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", pch=19)
		par(xpd=NA)
		title(main = paste("patient 1 - chr", i, sep =""), outer = TRUE, line = 2, cex.main=2)
		title(xlab="position", cex.lab=2)
		title(ylab="allele frequency", cex.lab=2)
		}
}

dev.off()
