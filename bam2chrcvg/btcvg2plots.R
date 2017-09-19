#!/usr/bin/RScript

# plot coverage per chromosome from bedtools grouped-by data and a list of fasta headers
# usage: btcvg2plots.R -b <bedtools stats> -t <fasta header list>
#
# Stephane Plaisance VIB-BITS September-11-2017 v1.0
# visit our Git: https://github.com/Nucleomics-VIB

# requirements
# asm2cvg.sh (ngs-tools) for generating the data for plotting and the plot titles

# dependencies. Please install them if not yet available on your machine
suppressPackageStartupMessages(library("optparse")); # take care of command-line arguments
suppressPackageStartupMessages(library("ggplot2")); # nice plots

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
outfolder="coverage_plots"
output_dir <- file.path(curfolder, outfolder)

if (!dir.exists(output_dir)){
dir.create(output_dir)
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
cex=0.3

for (onestat in c("min", "mean", "median", "max")) {
plot(oneseq[,onestat],
	 log="y",
	 main = "",
	 xlab = "",
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

# one plot only
par(mfrow=c(1,1),
    oma = c(2,2,3,1) + 0.1,
    mar = c(4,4,2,1) + 0.1)
cex=0.3

plot(oneseq[,plotstat],
     log="y",
     main = "",
     xlab = "",
     ylab = paste(plotstat,"-coverage (log10)",sep=""),
     type="p",
     pch=20,
     cex=cex,
     col = 'blue')
     
# add title for the page
page.title <- paste("Coverage plots for ", title, sep="")
page.title <- paste(strwrap(page.title,80), collapse="\n")
mtext(page.title, outer=TRUE,  cex=0.75, line=-1.5)
}
dev.off()
}