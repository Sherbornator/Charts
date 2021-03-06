
print(files <- list.files(pattern="all.txt$"))

print(labs <- paste(gsub("\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))
 
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
cov[[i]] <- read.table(files[i])
cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}

library(RColorBrewer)
cols <- c("lightgreen", "seagreen", "lightblue", "cornflowerblue", "blue")

png("Human-exome-coverage-plots.png", h=1000, w=1000, pointsize=20)
 
plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Exome Coverage - Human tumors")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))
 
for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])

legend("topright", legend=labs, col=cols, lty=1, lwd=4)
 
dev.off() 
