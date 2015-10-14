MAF <- read.table("QS.Normal.Tumor.MAF.ref.alt.txt", header=TRUE, sep="\t")
CP <- read.table("QS2-tumor.bwa.realigned.rmDups.recal.bam_ratio.txt", header=TRUE, sep="\t")
PV <- read.table("QS2-tumor.bwa.realigned.rmDups.recal.bam_CNVs.p.value.txt", header=TRUE, sep="\t")
m <- read.table("5x5matrix.csv", sep=",", header=FALSE)
colnames(m) <- NULL
m <- as.matrix(m)

loh <- data.frame(MAF)
ratio <- data.frame(CP)
pvalue <- data.frame(PV)

pvalue.split  <- split(x = pvalue[,c("chr", "start", "end", "copy.number", "status", "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue")], f = pvalue$status)
pvalue.gain <- pvalue.split[[1]]
pvalue.loss <- pvalue.split[[2]]

ploidy <- 2

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}

png(filename = "Patient1_CNp_LOH_genome.png", width = 1760, height = 1180, units = "px", bg = "white", res = NA)
plot(1:10)
layout(m)
par(mar = c(0, 0.5, 0, 0), oma = c(0, 0, 0, 0), xpd = NA)
	
for (i in c(1:22,'X','Y')) {
	ll <- which(loh$chromosome==(paste("chr", i, sep = "")))
	tt <- which(ratio$Chromosome==i)
	pg <- which(pvalue.gain$chr==i)
	pl <- which(pvalue.loss$chr==i)
	maxPosition <- max(loh$position[ll], ratio$Start[tt])
	
	if (length(tt)>0) {
		plot(ratio$Start[tt], ratio$Ratio[tt]*ploidy, xlim=c(0,maxPosition), ylim = c(0,maxLevelToPlot*ploidy), xaxs="i", yaxs="i", yaxt="n", xaxt="n", xlab = "", ylab = "normalized CN", cex.lab=1.5, pch = ".", col = colors()[88])
		axis(2, yaxp = c(1, maxLevelToPlot*ploidy, 5), cex.axis=1.5)
		tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
		points(ratio$Start[tt],ratio$Ratio[tt]*ploidy, pch = ".", col = colors()[136])
	
		tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
		points(ratio$Start[tt],ratio$Ratio[tt]*ploidy, pch = ".", col = colors()[136], cex=2)
	 
		tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
	 		points(ratio$Start[tt],ratio$Ratio[tt]*ploidy, pch = ".", col = colors()[461])
	 		tt <- which(ratio$Chromosome==i)
		par(xpd=NA)
    	segments(x0 = pvalue.gain$start[pg], x1 = pvalue.gain$end[pg],
            y0 = (5.5 + pvalue.gain$copy.number[pg] - pvalue.gain$copy.number[pg]), y1 = (5.5 + pvalue.gain$copy.number[pg] - pvalue.gain$copy.number[pg]),
            col="red", lwd = 2, lend = 3)
    	segments(x0 = pvalue.loss$start[pl], x1 = pvalue.loss$end[pl],
            y0 = (0.5 + pvalue.loss$copy.number[pl] - pvalue.loss$copy.number[pl]), y1 = (0.5 + pvalue.loss$copy.number[pl] - pvalue.loss$copy.number[pl]),
            col="blue", lwd = 2, lend = 3)
        }
	
	if (length(ll)>0) {
		###plot(loh$position[ll], loh$nMAFref[ll], col='#008837', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xlab="", ylab="", cex.axis=1.25, pch=".")
		###par(new=T)
		###plot(loh$position[ll], loh$nMAFalt[ll], col='#008837', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
		###par(new=T)
		plot(loh$position[ll], loh$tMAFref[ll], col='#7b3294', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xlab="", ylab="", cex.axis=1.25, pch=".")
		par(new=T)
		plot(loh$position[ll], loh$tMAFalt[ll], col='#7b3294', ylim=c(0,100), xlim=c(0,maxPosition), yaxs="i", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
		par(xpd=NA)
		title(xlab= paste("position - chr", i), cex.lab=2)
		title(ylab="allele frequency", cex.lab=1.5)
		}
}
par(xpd=TRUE)
title(main = paste("patient 1"), outer = TRUE, line = -3, cex.main=2)

dev.off()
