path = "/Volumes/HAEMATONC/SHARED/Myeloma/digitalMLPA/CN_array_data/Tumour_split/"

#Takes the names of the files from the segmentation files (which are simply called the 5 digit identifier) "^" identifies the start of a string, and the "$" the end of a line.
sample.names <- dir(path, pattern ="^[0123456789]{5}$")

out.file<-""

#loop through
for(sample in 1:length(sample.names)){
	#This script takes copy number data from a PennCNV file and plots it as one chromosome output
	#Could be modified for other data types
	###N.B. Data needs to be sorted by position

	#read in the file and attach data frame
	nsp <- read.table(paste("nsp.gw6.",sample.names[sample],"_DNANspI.sort", sep=""), header=TRUE, sep="\t")
	sty <- read.table(paste("sty.gw6.",sample.names[sample],"_DNAStyI.sort", sep=""), header=TRUE, sep="\t")
	arraycn <- data.frame(rbind(nsp, sty))

	CN <- read.table(sample.names[sample], header=FALSE, sep=",")
	CNchr_full <- read.table("CN_chr_blocks.csv", header=FALSE, sep=",")
	names(CN) = c("Sample", "CNVlocation", "Chr", "Start", "End", "SNPnumber", "Length", "CN", "State", "PATH", "StartSNP", "EndSNP")
	names(CNchr_full) = c("Sample", "CNVlocation", "Chr", "Start", "End", "SNPnumber", "Length", "CN", "State", "PATH", "StartSNP", "EndSNP")
	CNseg <- data.frame(CN)
	CNchrs <- data.frame(CNchr_full)
	CNseg <- rbind(CNseg, CNchrs)

	#define chromosomes (for males) and order
	#chr <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y')
	#ends <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
	#chr.order <- c(as.numeric(c(1:22)),'X','Y')

	#OR if the sample is female, use one without a Y
	chr <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X')
	ends <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560)
	chr.order <- c(as.numeric(c(1:22)),'X')

	chr.ends <- data.frame(chr, ends)
	chr.order <- as.matrix(chr.order)

	#read in data and calculate max positions for each chromosome
	#split data by chromosome
	seg.list.cn  <- split(x = arraycn[,c("Name", "Chr", "Position", "LogR_ratio", "Allele_freq")], f = arraycn$Chr)
	#order by predefined chr.order
	seg.list.cn  <- seg.list.cn[order(order(chr.order))]
	#work out the max position for each chromosome
	seg.max.cn   <- split(x = chr.ends[,c("ends")], f = chr.ends$chr)
	seg.max.cn  <- seg.max.cn[order(order(chr.order))]
	#define a new matrix with position and logR ratio
	seg.pos.cn   <- lapply(seg.list.cn, "[", TRUE, c("Position", "LogR_ratio"))
	#bind max position per chromosome to this matrix
	seg.max.cn   <- cumsum(as.numeric(do.call(rbind, seg.max.cn)))
	#define chr.offset as zero to begin with
	chr.offset <- 0
	#loop through chromosomes to assign new position based on sum of all previous chromosomes
	for (i in 1:length(seg.pos.cn)){
		seg.pos.cn[[i]] <- seg.pos.cn[[i]] + chr.offset
    	colnames(seg.pos.cn[[i]]) <- c("abs.start.cn", "abs.ratio.cn")
    	chr.offset   <- seg.max.cn[i]
	}
	#define new positions for end of chromosomes based on the sums calculated above
	seg.max.cn      <- sapply(X = seg.pos.cn, FUN = function(x) x[nrow(x), "abs.start.cn" ])
	#bind everything together
	abs.list.cn     <- mapply(cbind, seg.list.cn, seg.pos.cn, SIMPLIFY = FALSE)
	abs.segments.cn <- do.call(rbind, abs.list.cn)



	###For the CN lines
	#read in data and calculate max positions for each chromosome
	#split data by chromosome
	seg.list.line  <- split(x = CNseg[,c("Sample", "CNVlocation", "Chr", "Start", "End", "SNPnumber", "Length", "CN", "State", "PATH", "StartSNP", "EndSNP")], f = CNseg$Chr)
	#order by predefined chr.order
	seg.list.line  <- seg.list.line[order(order(chr.order))]
	#work out the max position for each chromosome
	seg.max.line   <- split(x = chr.ends[,c("ends")], f = chr.ends$chr)
	seg.max.line  <- seg.max.line[order(order(chr.order))]
	#define a new matrix with position and logR ratio
	seg.pos.line   <- lapply(seg.list.line, "[", TRUE, c("Start", "End", "CN"))
	#bind max position per chromosome to this matrix
	seg.max.line   <- cumsum(as.numeric(do.call(rbind, seg.max.line)))
	#define chr.offset as zero to begin with
	chr.offset <- 0
	#loop through chromosomes to assign new position based on sum of all previous chromosomes	
	for (i in 1:length(seg.pos.line)){
		seg.pos.line[[i]] <- seg.pos.line[[i]] + chr.offset
	    colnames(seg.pos.line[[i]]) <- c("abs.start.line", "abs.end.line", "abs.ratio.line")
    	chr.offset   <- seg.max.line[i]
	}
	#define new positions for end of chromosomes based on the sums calculated above
	seg.max.line      <- sapply(X = seg.pos.line, FUN = function(x) x[nrow(x), "abs.start.line" ])
	#bind everything together
	abs.list.line     <- mapply(cbind, seg.list.line, seg.pos.line, SIMPLIFY = FALSE)
	abs.segments.line <- do.call(rbind, abs.list.line)

	#split the CN segments into gains and losses for color purposes
	CN.loss <- subset(abs.segments.line, CN < 2, select=c("Sample", "CNVlocation", "Chr", "Start", "End", "SNPnumber", "Length", "CN", "State", "PATH", "StartSNP", "EndSNP", "abs.start.line", "abs.end.line"))
	CN.loss <- subset(CN.loss, Start != 0, select=c("Sample", "CNVlocation", "Chr", "Start", "End", "SNPnumber", "Length", "CN", "State", "PATH", "StartSNP", "EndSNP", "abs.start.line", "abs.end.line"))
	CN.gain <- subset(abs.segments.line, CN > 2, select=c("Sample", "CNVlocation", "Chr", "Start", "End", "SNPnumber", "Length", "CN", "State", "PATH", "StartSNP", "EndSNP", "abs.start.line", "abs.end.line"))

	#open png
	png(filename = paste(sample.names[sample],".png", sep = ""), width = 2360, height = 960, units = "px", bg = "white", res = NA)
	par(mar = c(0, 0, 0, 0), oma = c(5, 4.5, 5, 5) + 0.1)

	#draw plot with Logr ratio in green points
	plot(abs.segments.cn$abs.start.cn, abs.segments.cn$LogR_ratio, xlim=c(0, max(abs.segments.cn$abs.start.cn)), ylim=c(-1.75, 1.75), xaxs="i", yaxs="i", yaxt="n", xaxt="n", xlab = "", ylab = "", cex.lab=1.5, pch = ".", col = colors()[88])
	axis(2, cex.axis=1.5)
	#allow next title to be outside plot area, add title, then shift back to plot area
	par(xpd=NA)
	title(ylab= "LogR ratio", cex.lab=2)
	par(xpd=FALSE)
	#color points equivalent to copy number > 1.3 as red
	tt <- which(abs.segments.cn$LogR_ratio>(0.37))
	points(abs.segments.cn$abs.start.cn[tt],abs.segments.cn$LogR_ratio[tt], pch = ".", col = colors()[136])
	#color points equivalent to copy number < 0.7 as blue
	tt <- which(abs.segments.cn$LogR_ratio<(-0.51))
	points(abs.segments.cn$abs.start.cn[tt],abs.segments.cn$LogR_ratio[tt], pch = ".", col = colors()[461])

	#add in the regions called as CN gain/loss
	segments(x0 = CN.gain$abs.start.line, x1 = CN.gain$abs.end.line, y0 = (((CN.gain$CN) / 4 * 3) - 1.5), y1 = (((CN.gain$CN) / 4 * 3) - 1.5), col="red", lwd = 2, lend = 3)
	segments(x0 = CN.loss$abs.start.line, x1 = CN.loss$abs.end.line, y0 = (((CN.loss$CN) / 4 * 3) - 1.5), y1 = (((CN.loss$CN) / 4 * 3) - 1.5), col="blue", lwd = 2, lend = 3)
        
	#add in vertical lines at chromosome boundries
	abline(v = c(0, seg.max.cn), lty = 3)
	for (i in 1:length(abs.list.cn)){
		max.pos <- nrow(abs.list.cn[[i]])
		mtext(chr.order[i], side = 3, line = 0, cex=1.5,
            at = sum(abs.list.cn[[i]]$abs.start.cn[1], abs.list.cn[[i]]$abs.start.cn[max.pos])/2)
    }
	par(xpd=TRUE)
	#add in legend
	#NB. have subtracted 700,000 here to make numbers round, as the first snp is not near zero - may need to alter if other chips have earlier snps.
	axis(labels = as.character(round(seq((abs.list.cn[[1]]$Position[1]-700000)/1e6, abs.list.cn[[1]]$Position[nrow(abs.list.cn[[1]])]/1e6, by = 50), 0)),
    	at = seq(abs.list.cn[[1]]$Position[1]-700000, abs.list.cn[[1]]$Position[nrow(abs.list.cn[[1]])]-700000, by = 5e7), cex.axis = 1.15, cex = par("cex.axis")*par("cex"),
		side = 1 , line = 1)
        
	par(xpd=NA)
	axis(side = 4, cex.axis=1.5, labels = c(0, 1, 2, 3, 4), at = c(-1.5, -0.75, 0, 0.75, 1.5))
	mtext(side = 4, 'Copy Number Segments', cex=2,  line = 3.85)
	title(main = paste(sample.names[sample], sep=""), outer = TRUE, line = 2, cex.main=2)
	mtext(side=1, "position (Mb)", cex=2, line=3.85, adj=0)               

    dev.off()
}
