#!/usr/bin/env R

base = dirname(dirname(getwd()))

data <- read.csv(paste0(base,"/output_new_30nuc.csv"), sep=",", header=TRUE)

#parser <- ArgumentParser(description='Process some integers')

