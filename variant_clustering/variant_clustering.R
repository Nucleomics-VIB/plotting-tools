#!/usr/bin/Rscript

# create a dendrogram-plot from multigenome vcf data converted with vcf2binary.pl
# usage variant_clustering.R <input.txt>

library(ape)
library(ggplot2)

args=commandArgs(TRUE)
dataFile <- args[1]
mytitle<-sub(".txt$", "", dataFile)

# load data, convert & remove unused objects to free memory
data<-read.table(file=dataFile, header=TRUE, sep="\t")
rownames(data) <- apply(data[,c(1,2,4,5)], 1, function(x) paste0(x, sep="_") )

# remove left columns
data1<-as.matrix(data[,-1:9])
rm(data)

# permute data
perm.data<-aperm(data1)
rm(data1)

# save data for re-use!
write.csv(perm.data, file = "perm_data.csv", row.names = FALSE)

euclidian.distance<-dist(perm.data, method="euclidian")
euclidian.tree<-hclust(euclidian.distance, method="mcquitty")
manhattan.distance<-dist(perm.data, method="manhattan")
manhattan.tree<-hclust(manhattan.distance, method="mcquitty")

# save distances to text file
distFile<-paste("distances-",dataFile,sep="")
sink(distFile)
cat("# pairwise distances between the input genomes\n")
cat("euclidian distance\n")
euclidian.distance
cat("\n\n# manhattan distance\n")
manhattan.distance
sink()

# print multigraph graphic to pdf
# avoid clipping labels at page edges
par(xpd=TRUE)
outFile<-paste("dendogram-",dataFile,".pdf",sep="")
pdf(outFile, bg = "white")
# graph
op<-par(mfrow=c(2,1), cex=0.5)
plot(euclidian.tree,col="blue", main=paste(mytitle," dendrogram (euclidian-distance)",sep=""))
plot(manhattan.tree,col="red",  main=paste(mytitle," dendrogram (manhattan-distance)",sep=""))
dev.off()

# print multigraph graphic to pdf
outFile<-paste("phylogram-",dataFile,".pdf",sep="")
pdf(outFile, bg = "white")
# graph
# op<-par(mfrow=c(2,1), cex=0.95)

# euclidian
plot(as.phylo(euclidian.tree), 
	font=2,
	no.margin = FALSE,
	cex=0.5,
	type = "fan", 
	tip.color = hsv(runif(15, 0.65, 0.95), 1, 1, 0.7), 
	edge.color = hsv(runif(10, 0.65, 0.75), 1, 1, 0.7), 
	edge.width = runif(20, 0.5, 3), 
	use.edge.length = FALSE, 
	col = "gray80")
title(paste(mytitle," (euclidian.tree)",sep=""))

# manhattan
plot(as.phylo(manhattan.tree), 
	font=2,
	no.margin = FALSE,
	cex=0.5,
	type = "fan",
	tip.color = hsv(runif(15, 0.65, 0.95), 1, 1, 0.7), 
	edge.color = hsv(runif(10, 0.65, 0.75), 1, 1, 0.7), 
	edge.width = runif(20, 0.5, 3), 
	use.edge.length = FALSE, 
	col = "gray80")
title(paste(mytitle," (manhattan.tree)",sep=""))
dev.off()

# save workspace for faster edit
save.image(file = "clustered_genomes.RData")

q()
