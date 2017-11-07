#!/usr/bin/env R

base = dirname(dirname(getwd()))
par(mar=c(10,6,2,1))

inf <- read.csv(paste0(base,"/code/supplemental/cosmic_cleanup_bar_graph.txt"))

cols <- colnames(inf)

vals <- as.numeric(inf[1,])

t <- "Data removed at various points within the data clean up"
xa <- "Reason for extraction"
ya <- "Data values remaining"
#barplot(vals, names.arg = cols, xlab = xa, ylab = ya, 
#main = t, col = c("gold"), las = 2, border="white")

postscript(paste0(base,"/R_graphs/figures/cleanup.ps")
b <- barplot(vals, names.arg = cols, axes=F, xlab="", ylab="", col = c("gold"), las = 2, border="white")
axis(2, cex.axis=1.1, las=2)
mtext(ya, side=2, line=4.5, cex=1.2)
dev.off()
