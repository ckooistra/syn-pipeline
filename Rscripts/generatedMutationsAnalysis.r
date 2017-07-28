library(ggplot2)

output2 <- utils::read.csv("/home/chris/hd1/COSMIC_V80/generatedMutations.csv")


#Boxplot of transitions versus transversion for Generated Mutations
png("/home/chris/Dropbox/BIN 3005/R graphs/generatedMutations/transitionVtransversionforGeneratedMutations.png", width = 1000, height = 1000)
ggplot(output2, aes(x=output2$`TransitionOrTransversion`, y=output2$OC, fill=output2$TransitionOrTransversion)) + geom_boxplot() + theme_bw()
dev.off()

#Histogram of codon optimality for Generated Mutations
png("~/Dropbox/BIN\ 3005/R\ graphs/generatedMutations/CodonOptimalityDelataforGeneratedMutations.png", width = 1000, height = 1000)
hist(output2$OC, main="Codon Optimality Delta Distribution for Generated Mutations")
dev.off()

# Not Possible
#Mutations in first 10 percent of amino acid chain
#firstTen <- subset(output2, output2$ <= 0.1)


#Boxplot for first ten with transition and transversion
png("/home/chris/Dropbox/BIN 3005/R graphs/firstTenBoxplot.png", width = 1000, height = 1000)
ggplot(firstTen, aes(x=firstTen$TransitionOrTransversion, y=firstTen$OC, fill=firstTen$TransitionOrTransversion)) + geom_boxplot() + theme_bw()
dev.off()