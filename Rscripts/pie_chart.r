#Get Pie data
pie_data <- read.table("/home/chris/Dropbox/BIN_3005/code/supplemental/pie_data.txt", sep=" ", header=TRUE)


pie_vector <- c(pie_data$no,pie_data$yes)
pie_vector_lables <- c("Not Close", "Close to Exon")

print(pie_data)
print(pie_vector)

postscript("/home/chris/Dropbox/BIN_3005/R_graphs/figures/exon_pie.ps")
pie(pie_vector, labels=pie_vector_lables, main="Mutations close to exon junctures")
dev.off()
