#!/usr/bin/RScript

# plot coverage per chromosome from bedtools grouped-by data and a list of fasta headers
# usage: btcvg2plots.R -b <bedtools stats> -t <fasta header list> -s <stat>
#
# Stephane Plaisance VIB-BITS September-11-2017 v1.0
# version="1.1, 2017_09_19"
# visit our Git: https://github.com/Nucleomics-VIB

# requirements
# asm2cvg.sh (ngs-tools) for generating the data for plotting and the plot titles

# dependencies. Please install them if not yet available on your machine
suppressPackageStartupMessages(library("optparse")); # take care of command-line arguments
suppressPackageStartupMessages(library("ggplot2")); # nice plots
suppressPackageStartupMessages(library("gplots")); #1 table on plot

#####################################
### Handle COMMAND LINE arguments ###
#####################################

# parameters
option_list <- list(
  make_option(c("-b", "--bedstats"), type="character",
			  help="bedtools grouped-by results name [REQUIRED]"),
  make_option(c("-t", "--titles"), type="character",
			  help="fasta headers as titles for the plots [REQUIRED]"),
  make_option(c("-s", "--stat"), type="character",
			  help="stat to be used for the plots [min, mean, median, max, all | all])")
  )

## parse options
opt <- parse_args(OptionParser(usage = "%prog [options] file",
							   option_list=option_list, add_help_option = TRUE))

# check if arguments provided
# check that bam file exists
if ( is.null(opt$bedstats) ) {
  stop("no bedtools stats file provided, run with --help")
}
if( file.access(opt$bedstats) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", opt$bedstats))
}
# check that fasta file exists
if ( is.null(opt$titles) ) {
  stop("no title file provided, run with --help")
}
if( file.access(opt$titles) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", opt$titles))
}

if ( is.null(opt$stat) ) {
	plotstat <- "all"
} else {
	plotstat <- opt$stat
}

# prevent Rplot.pdf creation
pdf(NULL)

# load data in
bedtools.data <- read.delim(opt$bedstats, sep="\t", header=FALSE)
chromosome.titles <- read.delim(opt$titles, sep="\t", header=FALSE, stringsAsFactors=FALSE)$V1
colnames(bedtools.data) <- c("seq", "start", "end", "min", "median", "mean", "max")

# list of contigs/chromosomes
sequences <- as.vector(unique(bedtools.data$seq))

# create output folder
curfolder=getwd()
subfolder <- basename(gsub("\\..+?\\.titles", "", opt$titles))
outfolder=paste("coverage_plots", subfolder ,sep="-")
output_dir <- file.path(curfolder, outfolder)

if (!dir.exists(output_dir)){
dir.create(output_dir, recursive=TRUE)
}

# loop here
for (seqname in sequences) {

oneseq <- subset(bedtools.data, bedtools.data$seq == seqname)
title <- chromosome.titles[grep(seqname, chromosome.titles, fixed=TRUE)]
  
filename <- paste(outfolder, "/", seqname,"_plots.pdf", sep="")
pdf(file=filename)

if (plotstat=="all") {

par(mfrow=c(4,1),
	oma = c(2,2,3,1) + 0.1,
	mar = c(4,4,2,1) + 0.1)
cex=0.5

for (onestat in c("min", "mean", "median", "max")) {
plot(oneseq[,onestat],
	 log="y",
	 main = "",
	 xlab = "window index",
	 ylab = paste(onestat,"-coverage (log10)",sep=""),
	 type="p",
	 pch=20,
	 cex=cex,
	 col = 'blue')
}

# add title for the page
page.title <- paste("Coverage plots for ", title, sep="")
page.title <- paste(strwrap(page.title,80), collapse="\n")
mtext(page.title, outer=TRUE,  cex=0.75, line=-1.5)

} else {

# one plot only with loess line
par(mfrow=c(1,1),
	oma = c(2,2,3,1) + 0.1,
	mar = c(4,4,2,1) + 0.1)
cex=1; #0.3

plot.data <- data.frame(x=oneseq$start+(oneseq$end-oneseq$start)/2, y=oneseq[,plotstat])
suppressWarnings(lw1 <- loess(y ~ x,data=plot.data))

plot(plot.data$x, plot.data$y,
	log="y",
	main = "",
	xlab = "window index",
	ylab = paste(plotstat,"-coverage (log10)",sep=""),
	type="p",
	pch=20,
	cex=cex,
	col = 'blue')
j <- order(plot.data$x)
lines(plot.data$x[j],lw1$fitted[j], col="red", lwd=3)

# add title for the page
page.title <- paste("Coverage plots for ", title, sep="")
page.title <- paste(strwrap(page.title,80), collapse="\n")
mtext(page.title, outer=TRUE,  cex=0.75, line=-1.5)
}
done <- dev.off()
}

# plot coverage stats across all contigs
filename <- paste(outfolder, "/", subfolder, "-coverage_depth.pdf", sep="")
pdf(file=filename)

par(oma = c(2,2,3,1) + 0.1,
	mar = c(4,4,2,1) + 0.1)
cex=1; #0.3

max=4*quantile(bedtools.data[,plotstat], 0.75)
hist(bedtools.data[,plotstat], 
	xlim=c(0,max), 
	breaks=max,
	main = "",
	xlab = "coverage depth",
	ylab = "number of windows at a given depth",
	cex=cex,
	col = 'gray')

# customize quantiles
p <- c(0, 5, 25, 50, 75, 95, 100)/100
qt <- round( quantile(bedtools.data[,plotstat],
	probs = p,
	na.rm = FALSE,
	names = TRUE,
	type = 7),
	3 )

sum.table <- data.frame(q=paste(100*p,"%",sep=''), coverage=qt)

# overlay summary
par(new=TRUE)

textplot(sum.table, halign="right", valign="top", cex=1, show.rownames=FALSE, show.colnames=TRUE)

# add title for the page
page.title <- paste("Coverage depth distribution for ", subfolder, sep="")
page.title <- paste(strwrap(page.title,80), collapse="\n")
mtext(page.title, outer=TRUE,  cex=0.75, line=-1.5)

done <- dev.off()
