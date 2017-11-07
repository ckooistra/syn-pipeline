
#Get and load data
data <- read.csv("/home/chris/Dropbox/BIN_3005/output.csv")

#Get precompiled genelist
genes <- read.table("/home/chris/Dropbox/BIN_3005/code/supplemental/genesOver350Mutations.txt")[,1]

#Get list of subtypes
subtypes <- read.table("/home/chris/Dropbox/BIN_3005/code/supplemental/subtypes.txt")[,1]

#Create matrix to hold values for heatmap [geneName,subtype]
m <- matrix(0*(length(genes)*length(subtypes)), nrow=length(genes), ncol=length(subtypes))

#Give column and rownames to matrix
colnames(m) <- subtypes
rownames(m) <- genes

print(subtypes)
print(genes)

#axis(1, at=c(0, 100, 200, 300), cex.axis=1.2)
#Iterate over df data and caclulate the instances of a cancer subtype per gene
for (i in 1:length(data[,1])){

	if (data[i,1] %in% genes){
		ge<-toString(data[i,1])
		sub<-toString(data[i,15])

		m[ge,sub]=m[ge,sub]+1
	}

}

print(m)
#postscript("/home/chris/Dropbox/BIN_3005/R_graphs/figures/heatmap.ps", width=6, height=4, paper="special", onefile=F)

png("/home/chris/Dropbox/BIN_3005/R_graphs/figures/heatmap.png", width = 1000, height = 1000)
heatmap(m)
#par(mar=c(6,6,4,1))
#axis(2)
#mtext("Number of occurences", side=1, line=3, cex=1.2) 
#mtext("Genes with >350 mutations and subtype heatmap", side=3, line=2, cex=1.3) 
dev.off()

