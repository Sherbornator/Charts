MAF <- read.table("QS.Normal.Tumor.MAF.ref.alt.txt", header=TRUE, sep="\t")
CP <- read.table("QS2-tumor.bwa.realigned.rmDups.recal.bam_ratio.txt", header=TRUE, sep="\t")
PV <- read.table("QS2-tumor.bwa.realigned.rmDups.recal.bam_CNVs.p.value.txt", header=TRUE, sep="\t")

loh <- data.frame(MAF)
ratio <- data.frame(CP)
pvalue <- data.frame(PV)

chr.order <- c(as.numeric(c(1:22)),'X','Y')
chr.order <- as.matrix(chr.order)

seg.list.cn  <- split(x = ratio[,c("Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber")], f = ratio$Chromosome)
seg.list.cn  <- seg.list.cn[order(order(chr.order))]
seg.max.cn   <- lapply(X = seg.list.cn, FUN = function(x) x[nrow(x), "Start" ])
seg.pos.cn   <- lapply(seg.list.cn, "[", TRUE, c("Start", "Ratio"))
seg.max.cn   <- cumsum(as.numeric(do.call(rbind, seg.max.cn)))
chr.offset <- 0
for (i in 1:length(seg.pos.cn)){
	seg.pos.cn[[i]] <- seg.pos.cn[[i]] + chr.offset
    colnames(seg.pos.cn[[i]]) <- c("abs.start.cn", "abs.ratio.cn")
    chr.offset   <- seg.max.cn[i]
}
seg.max.cn      <- sapply(X = seg.pos.cn, FUN = function(x) x[nrow(x), "abs.start.cn" ])
abs.list.cn     <- mapply(cbind, seg.list.cn, seg.pos.cn, SIMPLIFY = FALSE)
abs.segments.cn <- do.call(rbind, abs.list.cn)

seg.list.pv  <- split(x = pvalue[,c("chr", "start", "end", "copy.number", "status", "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue")], f = pvalue$chr)
seg.list.pv  <- seg.list.pv[order(order(chr.order))]
seg.max.pv   <- lapply(X = seg.list.cn, FUN = function(x) x[nrow(x), "Start" ])
seg.pos.pv   <- lapply(seg.list.pv, "[", TRUE, c("start", "end"))
seg.max.pv   <- cumsum(as.numeric(do.call(rbind, seg.max.pv)))

chr.offset <- 0
for (i in 1:length(seg.pos.pv)){
    seg.pos.pv[[i]] <- seg.pos.pv[[i]] + chr.offset
    colnames(seg.pos.pv[[i]]) <- c("abs.start.pv","abs.end.pv")
    chr.offset   <- seg.max.pv[i]
}

seg.max.pv      <- sapply(X = seg.pos.pv, FUN = function(x) x[nrow(x), "abs.end.pv" ])
abs.list.pv     <- mapply(cbind, seg.list.pv, seg.pos.pv, SIMPLIFY = FALSE)
abs.segments.pv <- do.call(rbind, abs.list.pv)
abs.segments.pv <- na.exclude(abs.segments.pv)

pvalue.split  <- split(x = abs.segments.pv[,c("chr", "start", "end", "copy.number", "status", "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue", "abs.start.pv", "abs.end.pv")], f = abs.segments.pv$status)
pvalue.gain <- pvalue.split[[1]]
pvalue.loss <- pvalue.split[[2]]
pvalue.gain.trim <- subset(pvalue.gain, WilcoxonRankSumTestPvalue<7.03E-65, select=c("chr", "start", "end", "copy.number", "status", "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue", "abs.start.pv", "abs.end.pv"))
pvalue.loss.trim <- subset(pvalue.loss, WilcoxonRankSumTestPvalue<7.03E-65, select=c("chr", "start", "end", "copy.number", "status", "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue", "abs.start.pv", "abs.end.pv"))

seg.list.loh  <- split(x = loh[,c("chromosome", "position", "ref", "alt", "nMAFref", "nMAFalt", "tMAFref", "tMAFalt")], f = loh$chromosome)
seg.list.loh  <- seg.list.loh[order(order(chr.order))]
seg.max.loh   <- lapply(X = seg.list.loh, FUN = function(x) x[nrow(x), "position" ])
seg.pos.loh   <- lapply(seg.list.loh, "[", TRUE, c("position", "nMAFref"))
seg.max.loh   <- cumsum(as.numeric(do.call(rbind, seg.max.loh)))
chr.offset <- 0
for (i in 1:length(seg.pos.loh)){
	seg.pos.loh[[i]] <- seg.pos.loh[[i]] + chr.offset
    colnames(seg.pos.loh[[i]]) <- c("abs.start.loh", "abs.nMAFref.loh")
    chr.offset   <- seg.max.loh[i]
}
seg.max.loh      <- sapply(X = seg.pos.loh, FUN = function(x) x[nrow(x), "abs.start.loh" ])
abs.list.loh     <- mapply(cbind, seg.list.loh, seg.pos.loh, SIMPLIFY = FALSE)
abs.segments.loh <- do.call(rbind, abs.list.loh)

ploidy <- 2

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}



png(filename = "Patient1_CN_Wilp_tLOH_genome_wide.png", width = 2360, height = 1180, units = "px", bg = "white", res = NA)
plot(1:10)
layout(matrix(c(1, 2), 2, 1, byrow = TRUE), heights=c(3,2))
par(mar = c(0, 0, 0, 0), oma = c(5, 4.5, 5, 2) + 0.1)
			
plot(abs.segments.cn$abs.start.cn, abs.segments.cn$Ratio*ploidy, xlim=c(0, max(abs.segments.cn$abs.start.cn)), ylim = c(0,maxLevelToPlot*ploidy), xaxs="i", yaxs="i", yaxt="n", xaxt="n", xlab = "", ylab = "", cex.lab=1.5, pch = ".", col = colors()[88])
axis(2, yaxp = c(1, maxLevelToPlot*ploidy, 5), cex.axis=1.5)
par(xpd=NA)
title(ylab= "normalized copy number ratio", cex.lab=2)
par(xpd=FALSE)
tt <- which(abs.segments.cn$CopyNumber>ploidy)
points(abs.segments.cn$abs.start.cn[tt],abs.segments.cn$Ratio[tt]*ploidy, pch = ".", col = colors()[136])
tt <- which(abs.segments.cn$Ratio==maxLevelToPlot & abs.segments.cn$CopyNumber>ploidy)	
points(abs.segments.cn$abs.start.cn[tt],abs.segments.cn$Ratio[tt]*ploidy, pch = ".", col = colors()[136], cex=2)
tt <- which(abs.segments.cn$CopyNumber<ploidy & abs.segments.cn$CopyNumber!= -1)
points(abs.segments.cn$abs.start.cn[tt],abs.segments.cn$Ratio[tt]*ploidy, pch = ".", col = colors()[461])
segments(x0 = pvalue.gain.trim$abs.start.pv, x1 = pvalue.gain.trim$abs.end.pv,
        y0 = (5.5 + pvalue.gain.trim$copy.number - pvalue.gain.trim$copy.number), y1 = (5.5 + pvalue.gain.trim$copy.number - pvalue.gain.trim$copy.number),
        col="red", lwd = 2, lend = 3)
segments(x0 = pvalue.loss.trim$abs.start.pv, x1 = pvalue.loss.trim$abs.end.pv,
        y0 = (0.5 + pvalue.loss.trim$copy.number - pvalue.loss.trim$copy.number), y1 = (0.5 + pvalue.loss.trim$copy.number - pvalue.loss.trim$copy.number),
        col="blue", lwd = 2, lend = 3)
abline(v = c(0, seg.max.cn), lty = 3)
   for (i in 1:length(abs.list.cn)){
      max.pos <- nrow(abs.list.cn[[i]])
      mtext(chr.order[i], side = 3, line = 0, cex=1.5,
            at = sum(abs.list.cn[[i]]$abs.start.cn[1], abs.list.cn[[i]]$abs.start.cn[max.pos])/2)
    }

par(xpd=TRUE)
###plot(abs.segments.loh$abs.start.loh, abs.segments.loh$nMAFref, xlim=c(0, max(abs.segments.cn$abs.start.cn)), ylim = c(0,100), col='#008837', yaxs="i", xaxs="i", xaxt="n", xlab = "", ylab = "", cex.axis=1.5, pch = ".")
###par(new=T)
###plot(abs.segments.loh$abs.start.loh, abs.segments.loh$nMAFalt, col='#008837', xlim=c(0, max(abs.segments.loh$abs.start.loh)), ylim=c(0,100), yaxs="i", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
###par(new=T)
plot(abs.segments.loh$abs.start.loh, abs.segments.loh$tMAFref, col='#7b3294', ylim=c(0,100), xlim=c(0, max(abs.segments.loh$abs.start.loh)), yaxs="i", xaxs="i", xaxt="n", cex.axis=1.5, xlab="", ylab="", pch=".")
par(new=T)
plot(abs.segments.loh$abs.start.loh, abs.segments.loh$tMAFalt, col='#7b3294', ylim=c(0,100), xlim=c(0, max(abs.segments.loh$abs.start.loh)), yaxs="i", xaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", pch=".")
abline(v = c(0, seg.max.loh), lty = 3)
   for (i in 1:length(abs.list.loh)){
      max.pos <- nrow(abs.list.loh[[i]])
    }
axis(labels = as.character(round(seq(abs.list.cn[[1]]$Start[1]/1e6, abs.list.cn[[1]]$Start[nrow(abs.list.cn[[1]])]/1e6, by = 50), 0)),
        at = seq(abs.list.cn[[1]]$Start[1], abs.list.cn[[1]]$Start[nrow(abs.list.cn[[1]])], by = 5e7), outer = FALSE, cex.axis = 1.5, cex = par("cex.axis")*par("cex"),
        side = 1 , line = 1)
        
par(xpd=NA)
title(main = paste("Patient 1"), outer = TRUE, line = 2, cex.main=2)
title(ylab="allele frequency", cex.lab=2)
mtext(side=1, "position (Mb)", cex=2, line=3.85, adj=0)               

    
dev.off()
