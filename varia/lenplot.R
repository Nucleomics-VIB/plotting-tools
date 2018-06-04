#!/usr/bin/env Rscript

# script: lenplot.R
# Aim: plot length distributions for a multi-fasta (avoid top outliers 1%)
# input can be archived fasta files
# requires bioawk (https://github.com/lh3/bioawk) for rapid feature extraction
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2017-06-04 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# depends on
library("readr")
library("ggplot2")

# adapt the path to your bioawk executable
bioawk <- "/opt/biotools/bioawk/bioawk"

# provided arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("A fasta must be supplied\n", call.=FALSE)
}

infile <- args[1]
base <- basename(infile)

# compute XNXX function
Nvalue <- function(lim, x, na.rm = TRUE) {
  # handle NA values
  if(isTRUE(na.rm)){
    x <- x[!is.na(x)]
    }
  cutval <- 100/lim
  # compute LXX and NXX
  sorted <- sort(x, decreasing = TRUE)
  SXX <- sum(x)/cutval
  csum <- cumsum(sorted)
  GTLXX <- as.vector(csum >= SXX)
  LXX=min(which(GTLXX == TRUE))
  NXX <- round(sorted[LXX], 1)
  # eg: get NXX with lst['NXX']
  NXX
  }

# build command
cmd <- sprintf("%s -c fastx \'{print $name, length($seq)}\' \'%s\' > /tmp/sizes.txt", bioawk, infile)

cat("# collecting lengths from file(s)\n")

t1 <- try(system(cmd, intern = TRUE, wait = TRUE))

cat("# loading results and plotting\n")

data <- read_delim("/tmp/sizes.txt",
                   "\t",
                   escape_double = FALSE,
                   col_names = FALSE,
                   trim_ws = TRUE,
                   col_types = cols())

colnames(data)<-c("name","len")

summary(data)
plotfile=paste0(basename(infile), "_size-plot.pdf")
plots <- pdf(file = plotfile, width=5, height=4, onefile = TRUE)

# avoid top outliers (1%)
maxplot <- quantile(data$len,0.99)
selection <- data[data$len<=maxplot,]

# set bwidth
#nc <- nclass.Sturges(selection$len); #low density
nc <- nclass.scott(selection$len); # mid density
#nc <- nclass.FD(selection$len); # high density
breaks <- pretty(range(selection$len), n = nc, min.n = 1)
bwidth <- breaks[2]-breaks[1]

ggplot(data, aes(x=len)) +
  geom_histogram(binwidth=bwidth, colour="white", fill="grey60") +
  scale_x_continuous(limits=c(0,maxplot)) +
  geom_vline(aes(xintercept=mean(len, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dotted", size=0.8) +
  geom_vline(aes(xintercept=Nvalue(50, selection$len)),   # Ignore NA values for median
             color="blue", linetype="dotted", size=0.8) +
  labs(title="", x="Sequence size distribution (bps) [mean=red, N50=blue]", y="count")

close.plots <- dev.off()

cat(paste0("# results were plotted to ", plotfile, "\n"))
