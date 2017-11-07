#!/usr/bin/env R

data <- read.csv("/home/chris/Dropbox/BIN_3005/output.csv", sep=",", header=TRUE)

fp = "/home/chris/Dropbox/BIN_3005/R_graphs/figures/mut_pos.png"

png(fp, width = 1000, height = 1000)
hist(data$position_AA, main="This is a test")
dev.off()
