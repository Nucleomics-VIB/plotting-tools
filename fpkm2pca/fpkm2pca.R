# load data from expXXXX-RNAseqCounts.xlsx
# generate PCA plots for the first 3 components (2 by 2)
# visit our Git: https://github.com/Nucleomics-VIB

# dependencies. Please install them if not yet available on your machine
library("openxlsx"); # read from excel files
library("pheatmap"); # plot pretty heatmaps (different module than used at the NC)
library("gplots"); # for heatmaps
#library(pcaMethods); # for PCA
library("RColorBrewer"); # build nice color palettes

# edit following lines to point to the folder and file
basedir <- ("~/git_repos/Nucleomics-VIB/plotting-tools/fpkm2PCA/data")
data.file <- "expXXXX-RNAseqCounts.xlsx"

log.transform <- TRUE
#log.transform <- FALSE

# load data from file (the first worksheet contains raw counts while the second worksheet contains FPKM data)
setwd(basedir)
excel.data <- read.xlsx(data.file, sheet=2)

# keep only data columns (remove last columns including "Chromosome")
chromosome.col <- which(colnames(excel.data)==as.vector("Chromosome"))
fpkm.data <- excel.data[,c(2,1,3:(chromosome.col-1))]

# remove some columns
row.names(fpkm.data) <- paste(fpkm.data[,1], fpkm.data[,2], sep=":")
fpkm.data <- fpkm.data[,-c(1,2)]

# convert 0 to NA to avoid log-normalization error
fpkm.data[fpkm.data == 0] <- NA

# remove full NA rows
fpkm.data <- fpkm.data[rowSums(is.na(fpkm.data[,2:length(fpkm.data)]))<length(fpkm.data[,2:length(fpkm.data)]),]

# kick useless part of names for samples
colnames(fpkm.data) <- sub("@.*", "", colnames(fpkm.data))

###################################
## to transform or not to transform
###################################

if (log.transform == TRUE) {
  data <- log(fpkm.data,2)
  log.txt=" (log2)"
  log.filename="_log"
} else {
  data <- fpkm.data
  log.txt=""
  log.filename=""
}

###################################
## plot heatmap of sample distances
###################################

# prepare output
filename <- paste("sample_correlation-plot", log.filename, ".pdf", sep="")
outfile <- paste(basedir, filename, sep="/")

# define metrics for clustering
drows <- "euclidean"
dcols <- "euclidean"
clustmet <- "average"

# color palette
hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(100)

# create heatmap
main.title <- paste("sample pairwise correlation", log.txt, sep="")
dist.data <- dist(t(data))
mat.samples <- as.matrix(dist.data)

hm.parameters <- list(mat.samples,
                      color = hmcol,
                      fontsize = 8,
                      cellwidth = 12, cellheight = 12, scale = "none",
                      treeheight_row = 200,
                      kmeans_k = NA,
                      show_rownames = T,
                      show_colnames = T,
                      main = main.title,
                      clustering_method = clustmet,
                      cluster_rows = TRUE,
                      cluster_cols = TRUE,
                      clustering_distance_rows = drows,
                      clustering_distance_cols = dcols)

# To draw the heatmap on screen (comment-out if you run the script from terminal)
do.call("pheatmap", hm.parameters)

# To draw to file (you may want to adapt the sizes)
do.call("pheatmap", c(hm.parameters, filename=outfile, width=10, height=7))


########################
## PCA analysis on FPKM
########################

# FPKM data (or log2-transformed) is transposed and centered
# The resulting matrix is used to compute PCA
# biplots for the first 3 principal components are shown

###############
# on fpkm data
###############
tdata <- t(data)

# compute PCA
p <- prcomp(tdata, scale=TRUE)

# compute PCA summary
summary <- summary(p)

# prepare output
filename <- paste("pca_plot", log.filename, ".pdf", sep="")
outfile <- paste(basedir, filename, sep="/")

pdf(file=outfile, bg="transparent", width=7, height=6)

# group plots
par(mfrow=c(2,2), mar=c(2,3,3,1), oma=c(1,1,1,1), mai=c(0.5,1,0.25,0))
cextitle <- 1
cexlab <- 0.75
cextxt <- 0.5
cexaxis <- 0.75
ptcol <- "blue"

title <- paste("variance across PC#1:5", log.txt, sep="")
plot(summary$importance[2,1:5],
     type="p",
     pch=20,
     cex.title=2,
     cex.lab=cexlab,
     cex=cextxt+0.5,
     cex.axis=cexaxis,
     col=ptcol,
     xlab=NA,
     ylab=NA,
     las=1
     )
mtext(side=1, text="Principal Components", line = 2, cex=cexlab)
mtext(side=2, text="% Variance", line = 3, cex=cexlab)
mtext(side=3, text=title, line = 0.5, cex=cexlab)

# plot first 1st and 2nd PC
plot(p$x[,1],p$x[,2],
     cex=0,
     cex.lab=cexlab,
     cex.axis=cexaxis,
     xlab=NA,
     ylab=NA,
     las=1)
mtext(side=1, text="PCA:1", line = 2, cex=cexlab)
mtext(side=2, text="PCA:2", line = 2.5, cex=cexlab)
text(p$x[,1],p$x[,2],
     labels=rownames(p$x),
     cex=cextxt,
     col=ptcol)

# plot next 2nd and 3rd PC
plot(p$x[,2],p$x[,3],
     cex=0,
     cex.lab=cexlab,
     cex.axis=cexaxis,
     xlab=NA,
     ylab=NA,
     las=1)
mtext(side=1, text="PCA:2", line = 2, cex=cexlab)
mtext(side=2, text="PCA:3", line = 2.5, cex=cexlab)
text(p$x[,2],p$x[,3],
     labels=rownames(p$x),
     cex=cextxt,
     col=ptcol)

# plot 1st and 3rd PC
plot(p$x[,1],p$x[,3],
     cex=0,
     cex.lab=cexlab,
     cex.axis=cexaxis,
     xlab=NA,
     ylab=NA,
     las=1)
mtext(side=1, text="PCA:1", line = 2, cex=cexlab)
mtext(side=2, text="PCA:3", line = 2.5, cex=cexlab)
text(p$x[,1],p$x[,3],
     labels=rownames(p$x),
     cex=cextxt,
     col=ptcol)

dev.off()




