#!/usr/bin/R

path = getwd()

#Takes the names of the files from the allele freq files (which are in the following format "10068_4;14_7_concatenate.txt") "^" identifies the start of a string, and the "$" the end of a line.
sample.names <- dir(path, pattern ="^[0123456789]+_[0123456789]{1,2};[0123456789]{2}_[0123456789]{1,2}_concatenate.txt$")


#loop through
for(sample in 1:length(sample.names)){
	#This script plots allele frequency data for matched presentation and relapse samples
	#Could be modified for other data types
	print(paste("Processing", sample.names[sample], sep=" "))
	names = strsplit(sample.names[sample], "_")
	patientID = names[[1]][1]
	transID = paste(names[[1]][2],"_",names[[1]][3], sep="")

    allele.frac <- read.table(sample.names[sample], header=TRUE, sep="\t")
    attach(allele.frac)
    
    total <- nrow(allele.frac)
    pres_zero <- sum(allele.frac$pres_MAF == 0)
    relapse_zero <- sum(allele.frac$relapse_MAF == 0)
    
    png(paste("Allele_freq_", patientID, "_", transID, "_cosmic.png", sep=""), h=1000, w=1000, pointsize=20)

    plot(pres_MAF,relapse_MAF, xlim=c(-0.02,0.52), ylim=c(-0.02,0.52), yaxs="i", xaxs="i", main=paste("Patient", patientID, transID, "\nsomatic variants (cosmic only)\ntotal variants =", total, sep=" "), xlab="presentation", ylab="relapse", cex.lab=1.25, pch=19)

    loc <- par("usr")
    #label axis
    text(loc[1], loc[4], "allele\nfrequency", adj = c(0.5,-0.8),  xpd = T)
    text(loc[2], loc[3], "allele\nfrequency", adj = c(0.5,2.5), xpd = T)
    #give numbers of variants with zero MAF in pres and relapse
    text(loc[2], loc[3], paste("variants with\nzero MAF in\nrelapse = ", relapse_zero, sep=""), adj = c(1.2,-1),  xpd = T)
    text(loc[1], loc[4], paste("variants with\nzero MAF in\npres = ", pres_zero, sep=""), adj = c(-0.5,1.2),  xpd = T)
    
    detach(allele.frac)

    dev.off()
}
