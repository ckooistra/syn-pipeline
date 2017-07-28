library(ggplot2)
library(utils)

#Load in Cosmic Mutation Data into variable output and prepare various analyses

output <- utils::read.csv("~/hd1/COSMIC_V80/output.csv")

#Boxplot of transitions versus transversion 
png("/home/chris/Dropbox/BIN 3005/R graphs/transitionVtransversion.png", width = 1000, height = 1000)
ggplot(output, aes(x=output$TransitionOrTransversion, y=output$optimality.change ),fill=output$TransitionOrTransversion) + geom_boxplot() + theme_bw()
dev.off()

#Histogram of codon optimality
png("/home/chris/Dropbox/BIN 3005/R graphs/CodonOptimalityDelata.png", width = 1000, height = 1000)
hist(output$optimality.change, main="Codon Optimality Delta Distribution")
dev.off()

#Subset of transitions
#transi<- subset(output, TransitionOrTransversion == 'transition' & closeToExon == 'no')

#Subset of transversions
#transver<- subset(output, TransitionOrTransversion == 'transversion' & closeToExon == 'no')

#######################################################################################
#Mutations in first 10 percent of amino acid chain
firstTen <- subset(output, output$position.in.AA.chain.. <= 0.1)

#Transitions in the first 10% of amino acid chain
firstTenTransition <-subset(firstTen, TransitionOrTransversion == 'transition' & closeToExon == 'no')

#Transversions in the first 10% of amino acid chain
firstTenTransversion <-subset(firstTen, TransitionOrTransversion == 'transversion' & closeToExon == 'no')

#Boxplot for first ten with transition and transversion
png("/home/chris/Dropbox/BIN 3005/R graphs/firstTenBoxplot.png", width = 1000, height = 1000)
ggplot(firstTen, aes(x=firstTen$TransitionOrTransversion, y=firstTen$optimality.change, fill=firstTen$TransitionOrTransversion)) + geom_boxplot() + theme_bw()
dev.off()

########################################################################################
#Mutations in last 10% of amino acid chain
lastTen <- subset(output, output$position.in.AA.chain.. >= 0.9)

#Transition in that last 10% of amino acid chain
lastTenTransition <- subset(lastTen, TransitionOrTransversion == 'transition' & closeToExon == 'no')

#Transversions in last 10% of amino acid chain
lastTenTransversion <- subset(lastTen, TransitionOrTransversion == 'transversion' & closeToExon == 'no')

png("/home/chris/Dropbox/BIN 3005/R graphs/lastTenBoxplot.png", width = 1000, height = 1000)
ggplot(lastTen, aes(x=lastTen$TransitionOrTransversion, y=lastTen$position.in.AA.chain.., fill=lastTen$TransitionOrTransversion)) + geom_boxplot() + theme_bw()
dev.off()

###########################################################################################
#List of genes that show up over 300 times
over300 <- scan(file="/home/chris/hd1/COSMIC_V80/genesOver300.txt", what="character")

#Mutations that in genes with a count over 300
mutationsOver300 <- subset(output, gene %in% over300)

#Transitions - Mutations that in genes with a count over 300
mutationsOver300Transition <- subset(mutationsOver300, TransitionOrTransversion == 'transition' & closeToExon == 'no')

#Transversions - Mutations that in genes with a count over 300
mutationsOver300Transversion <- subset(mutationsOver300, TransitionOrTransversion == 'transversion' & closeToExon == 'no')

png("/home/chris/Dropbox/BIN 3005/R graphs/over300Boxplot.png", width = 1000, height = 1000)
ggplot(mutationsOver300, aes(x=mutationsOver300$TransitionOrTransversion, y=mutationsOver300$optimality.change, fill=mutationsOver300$TransitionOrTransversion)) + geom_boxplot() + theme_bw()
dev.off()



