#!/usr/bin/env R


base = dirname(dirname(getwd()))

data <- read.csv(paste0(base,"/output_new_30nuc.csv"), sep=",", header=TRUE)
cancer_data <- read.csv(paste0(base,"/output.csv"), sep=",", header=TRUE)

syn_pro_data <- read.csv(paste0(base,"/generatedMutations/generatedSynMutationsWithProb.csv"), sep=",", header=TRUE)

syn_ran_data <- read.csv(paste0(base,"/generatedMutations/generatedMutations.csv"), sep=",", header=TRUE)


#png(paste0(base,"/R_graphs/figures/boxplot_of_distributions.png"), width = 1000, height = 1000)
boxplot(cancer_data$optimality.change, syn_ran_data$OC, syn_pro_data$OC, notch=TRUE)
#dev.off()
