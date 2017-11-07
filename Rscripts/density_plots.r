#!/usr/bin/env R

base = dirname(dirname(getwd()))
data <- read.csv(paste0(base,"/output_new_30nuc.csv"), sep=",", header=TRUE)


syn_pro_data <- read.csv(paste0(base,"/generatedMutations/generatedSynMutationsWithProb.csv"), sep=",", header=TRUE)
oncos <- scan(paste0(base,"/code/supplemental/oncogenes.txt"), what="character")

onco_df <- subset(data, gene %in% oncos)

fp <- paste0(base,"/R_graphs/figures/optimality_density_plots.ps")
fp1 <- paste0(base,"/R_graphs/figures/optimality_density_plots_oncos.ps")

t <- "Distribution of Cosmic Data Set vs Random \nfor Optimality Change Values"

xax <- "Optimality"
yax <- "Number of Values"

postscript(fp)
plot(density(data$optimality.change, bw = 0.1), col = c("red"), main = t, xlab = xax, ylab = yax, lwd = 3)
lines(density(syn_pro_data$OC, bw = 0.1), col=c("blue"),lwd = 3) 
dev.off()

t <- "Distribution of Oncogenes Cancer Data Set vs Random \nfor optimality change values"

xax <- "Optimality"
yax <- "Number of Values"
postscript(fp1)
plot(density(onco_df$optimality.change, bw = 0.1), col = c("red"), main = t, xlab = xax, ylab = yax,lwd = 3)
lines(density(syn_pro_data$OC, bw = 0.1), col=c("blue"), lwd=3) 
dev.off()
