#!/usr/bin/env R


#Get Pie data
args <- commandArgs()
subtype <- args[6] #subtype

data <- read.csv("/home/chris/Dropbox/BIN_3005/output_new_30nuc.csv", sep=",", header=TRUE)

print(nrow(data))

if (!is.na(subtype)){
	data <- subset(data, cancerType == subtype)
}else{
	subtype <- "Full"
}
#Print number of rows in dataFrame
print(nrow(data))


print("No")
no <- subset(data, closeToExon == "no")
notclose = nrow(no)
print(notclose)

print("yes")
yes <- subset(data, closeToExon == "yes")
close = nrow(yes)
print(close)

title <- paste0("Proximity to Splice sites for ", subtype, " Data Set")


pie_vector <- c(notclose,close)
ncl <- paste0("Not close\n",notclose)
cl <- paste0("Close to Exon\n",close)
pie_vector_lables <- c(ncl,cl)

fn <- paste0("/home/chris/Dropbox/BIN_3005/R_graphs/figures/correctedFigures/subtypes/", subtype,".ps")
postscript(fn)
pie(pie_vector, labels=pie_vector_lables, main=title)
dev.off()
