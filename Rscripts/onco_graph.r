#!/usr/bin/env R
library(ggplot2)
library(utils)

#Create graphs for oncogenes i.e. genes foung in cancer database over 350 times
cancer_data <- read.csv("/home/chris/Dropbox/BIN_3005/output.csv", sep=",", header=TRUE)
oncos <- scan("/home/chris/Dropbox/BIN_3005/code/supplemental/oncogenes.txt", what="character")

subtypes <- scan("/home/chris/Dropbox/BIN_3005/code/supplemental/cancer_spots.txt", what="character")

print("Number of rows in full data frame")
print(nrow(cancer_data))

#Mutations that in genes with a count over 300
oncoGenes_data <- subset(cancer_data, gene %in% oncos)
oncoGenes_data <- subset(oncoGenes_data, cancerType %in% subtypes)
print("Number of rows in oncoGene data frame")
print(nrow(oncoGenes_data))

#png("/home/chris/Dropbox/BIN_3005/R_graphs/figures/boxplot_of_distributions.png", width = 1000, height = 1000)
#boxplot(oncoGenes_data$optimality.change, cancer_data$optimality.change, notch=TRUE)
#dev.off()

png("/home/chris/Dropbox/BIN_3005/R_graphs/figures/boxplot_of_onco_subtype_distributions.png", width = 2000, height = 1000)
print(bp)
bp <- ggplot(oncoGenes_data, aes(x=oncoGenes_data$cancerType, y=oncoGenes_data$optimality.change, colour = oncoGenes_data$cancerType )) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()



for (i in 1:length(subtypes)){
	

}
