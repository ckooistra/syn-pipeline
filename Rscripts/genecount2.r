
data <- read.table("/home/chris/Dropbox/BIN_3005/code/supplemental/geneCount.txt")
postscript("/home/chris/Dropbox/BIN_3005/R_graphs/geneMutationCount.ps", width=5, height=4, paper="special", onefile=F)

par(mar=c(6,6,4,1))
hist(data[,2], breaks = 500, xlim = c(0, 150),  axes=F, xlab="", ylab="" , main="", col="#29ABE2", border="white")
axis(1, at=c(0, 50, 100, 150), cex.axis=1.1)
axis(2, at=c(0, 1000, 2000, 3000), cex.axis=1.1)
mtext("Number of occurences", side=1, line=3, cex=1.1)
mtext("Frequency", side=2, line=2.5, cex=1.1)
mtext("Silent Mutations Iin COSMIC database", side=3, line=2, cex=1.3)

dev.off()
