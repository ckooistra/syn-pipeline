#!/usr/bin/env R

#Stacked Barchart of relative amount of mutations close to exon boundaries

base = dirname(dirname(getwd()))

exonThreshold = 30

data <- read.csv(paste0(base,"/output_new_30nuc.csv"), sep=",", header=TRUE)
subtypes <- scan(paste0(base,"/code/supplemental/cancer_spots.txt"), what="character")

l <- vector("list", length(subtypes))
for (i in 1:length(subtypes)){

	print(subtypes[i])
	sub <- subset(data, cancerType == subtypes[i])
	num_sub <- nrow(sub)

	no <- subset(sub, closeToExon >= 30)
	notclose = nrow(no)
	print(notclose)

	print("yes")
	yes <- subset(sub, closeToExon <= exonThreshold )
	close = nrow(yes)
	print(close)
	num_total <-close+notclose
	print(num_total)

	l[[i]] <-list(sub = subtypes[i], t = num_total, nc = notclose, c = close )
	

}

nc <- vector("list", length(subtypes))
c <- vector("list", length(subtypes))
for (i in 1:length(subtypes)){
	nc[[i]] <- l[[i]]$nc
	c[[i]] <- l[[i]]$c

}

nc<-unlist(nc)
c<-unlist(c)

m <- matrix(c(nc,c),ncol=length(subtypes), nrow=2, byrow=TRUE)
dimnames(m) = list(c('notClose', 'close'),subtypes)

title <- "Splice Site Proximity For Subtypes in COSMIC Synonymouns\nMutation Data Set"

fn <- ("/home/chris/Dropbox/BIN_3005/R_graphs/figures/exon_prox_stacked_bar_subtypes.ps")
postscript(fn)
#png(fn)
barplot(m, main = title, xlab = "Subtypes", ylab = "Proximity to Splice Sites", col=c("lightblue","lightgreen"), legend = rownames(m))
dev.off()


