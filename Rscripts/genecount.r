#!/usr/bin/env R

#
base = dirname(dirname(getwd()))

data <- read.table(paste0(base,"/code/supplemental/geneCount.txt"))
postscript(paste0(base,"/R_graphs/geneMutationCount.ps"), width=6, height=4, paper="special", onefile=F)

par(mar=c(6,6,4,1))
hist(data[,2], breaks = 200, xlim = c(0, 300),  axes=F, xlab="", ylab="" , main="", col=2)
axis(1, at=c(0, 100, 200, 300), cex.axis=1.2)
axis(2)
mtext("Number of occurences", side=1, line=3, cex=1.2)
mtext("Frequency", side=2, line=3, cex=1.3)
mtext("second label", side=2, line=1.8, cex=0.8)
mtext("Silent Mutations Iin COSMIC database", side=3, line=2, cex=1.3)

dev.off()
