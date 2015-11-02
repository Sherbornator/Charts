SS.allele.frac <- read.table("SS.variant.allele.freq.trim.csv", header=TRUE, sep=",", row.names=1)
attach(SS.allele.frac)

code <- as.factor(SS.allele.frac$code)

col.list <- c("#C4C4C466","#FF0000FF", "#2E8B57FF")
palette(col.list)

png("Allele_freq_scatterplot_labels.png", h=1000, w=1000, pointsize=20)

plot(TumorA_frac, TumorB_frac, xlim=c(-0.02,1.0), ylim=c(-0.02,1.0), yaxs="i", xaxs="i", main="patient 2 somatic variants", xlab="tumor a", ylab="tumor b", cex.lab=1.25, pch=19, col=code)

loc <- par("usr")
text(loc[1], loc[4], "allele\nfrequency", adj = c(0.5,-0.8),  xpd = T)
text(loc[2], loc[3], "allele\nfrequency", adj = c(0.5,2.5), xpd = T)
legend("topright", legend=levels(code), col=col.list, pch=19)

SS.allele.frac.ex.noncoding <- subset(SS.allele.frac, code != "noncoding", select = TumorA_frac:SIFT_code, drop = TRUE)
attach(SS.allele.frac.ex.noncoding)
gene <-as.factor(SS.allele.frac.ex.noncoding$gene)
text(TumorA_frac, TumorB_frac, labels=gene, cex= 0.5, pos=3)


dev.off()
