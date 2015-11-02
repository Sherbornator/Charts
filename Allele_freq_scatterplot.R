P3.allele.frac <- read.table("Patient3.variant.allele.freq.csv", header=TRUE, sep=",", row.names=1)
attach(P3.allele.frac)

cluster <- as.factor(P3.allele.frac$cluster)
col.list <- c("#C4C4C4FF", "#4daf4a", "#e41a1c", "#377eb8", "#984ea3", "#ff7f00", "#ffff33")
palette(col.list)

png("Allele_freq_scatterplot_P3_cluster.png", h=1000, w=1000, pointsize=20)

plot(Pre_frac,Post_frac, xlim=c(-0.02,1.0), ylim=c(-0.02,1.0), yaxs="i", xaxs="i", main="patient 3 somatic variants", xlab="pre", ylab="post", cex.lab=1.25, pch=19, col=cluster)

loc <- par("usr")
text(loc[1], loc[4], "allele\nfrequency", adj = c(0.5,-0.8),  xpd = T)
text(loc[2], loc[3], "allele\nfrequency", adj = c(0.5,2.5), xpd = T)
legend("topright", legend=c("no cluster", "cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5", "cluster 6"), col=col.list, pch=19)

dev.off()
