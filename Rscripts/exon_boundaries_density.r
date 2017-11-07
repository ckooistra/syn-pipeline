#!/usr/bin/env R

base = dirname(dirname(getwd()))

data <- read.csv(paste0(base,"/output.csv"), sep=",", header=TRUE)
subs <- scan(paste0(base,"/code/supplemental/subtypes.txt"), what="character")
fp <- paste0(base,"/R_graphs/figures/exon_boundary_density.ps")
s <- c()
c <- c()
nc <- c()
perc <- c()

t <- "Distribution of values close to exon splice sites"
xax <- "Percentage of values close exon/intron splice sites"
yax <- "Number of values"
for (i in 1:length(subs)){
	#Get subtype dataframe to be check
	subt <- data[data$cancerType == subs[i],]
	l = nrow(subt)

	#Only take array if over 1500 entries for subtype
	if (l < 1500) next
	close <-sum(subt$closeToExon == "yes")
	notClose <- sum(subt$closeToExon == "no")
	
	s <- append(s, subs[i])
	c <- append(c, close)
	nc <- append(nc, notClose)
	perc <- append(perc, close/notClose)

}

df = data.frame(s,c,nc,perc)

postscript(fp)
plot(density(df$perc), xlim = c(0,0.6), xlab = xax, ylab = yax)
dev.off()
