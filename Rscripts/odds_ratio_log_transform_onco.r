#!/usr/bin/env R

onco <- read.csv("/home/chris/Dropbox/BIN_3005/code/python/onco_full_data_set_nCodon.csv")
full <- read.csv("/home/chris/Dropbox/BIN_3005/code/python/full_data_set_nCodon.csv")

onco_m <- read.csv("/home/chris/Dropbox/BIN_3005/code/python/23_oct2017_full_data_set_mCodon.csv")
full_m <- read.csv("/home/chris/Dropbox/BIN_3005/code/python/23_oct2017_onco_full_data_set_mCodon.csv")

c <- as.vector(onco[,1])
or <- as.numeric(as.vector(onco[,2]))
pv <- as.numeric(as.vector(onco[,3]))
lor <- log2(or)

c1 <- as.vector(full[,1])
or1 <- as.numeric(as.vector(full[,2]))
pv1 <- as.numeric(as.vector(full[,3]))
lor1 <- log2(or1)

c2 <- as.vector(full_m[,1])
or2 <- as.numeric(as.vector(full_m[,2]))
pv2 <- as.numeric(as.vector(full_m[,3]))
lori2 <- log2(or2)

c3 <- as.vector(onco_m[,1])
or3 <- as.numeric(as.vector(onco_m[,2]))
pv3 <- as.numeric(as.vector(onco_m[,3]))
lor3 <- log2(or3)

f<-"/home/chris/Dropbox/BIN_3005/R_graphs/figures/overUnder.ps"

#dev.new(width = 15, height= 8)
postscript(f)
par(mfrow=c(4,1))
t = "Barplot of under/over expression of codons between random and generated mutations - Full Dataset \nMutated from"
barplot(lor, main = t, names.arg=c, cex.axis = 0.05, 
	las = 2, col="#EA61D1", border="white")


t1 = "Barplot of under/over expression of codons between random and generated mutations - 'Oncogenes' \nMutated from"
#dev.new(width = 15, height=4)
barplot(lor1, main = t1, names.arg=c, cex.axis = 0.05, 
	las = 2, col="#E66433", border="white")

t2 = "Barplot of under/over expression of codons between random and generated mutations - 'Full Dataset\n Mutated to'"
#dev.new(width = 15, height=4)
barplot(lori2, main = t2, names.arg=c, cex.axis = 0.05, 
	las = 2, col="#EA61D1", border="white")

t3 = "Barplot of under/over expression of codons between random and generated mutations - 'Oncogenes'\n Mutated to"
#dev.new(width = 15, height=4)
barplot(lor3, main = t3, names.arg=c, cex.axis = 0.05, 
	las = 2, col="#E66433", border="white")
dev.off()
