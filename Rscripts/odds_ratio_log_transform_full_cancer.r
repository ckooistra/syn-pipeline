#!/usr/bin/env R

base = dirname(dirname(getwd()))

d <- read.csv(paste0(base, "/code/python/test/full_data_set_nCodon.csv"))
c <- as.vector(d[,1])
or <- as.numeric(as.vector(d[,2]))
pv <- as.numeric(as.vector(d[,3]))
lor <- log2(or)

t = "Barplot of under/over expression of codons between random and generated mutations"
dev.new(width = 15, height=4)
barplot(lor, main = t, names.arg=c, cex.axis = 0.05, 
	las = 2, col="#EA61D1", border="white")
