#!/usr/bin/RScript
# load data from provided file-path
# usage: fpkm2heatmap-cmd.R -i expXXXX-RNAseqCounts.xlsx
#
# Stephane Plaisance VIB-BITS October-11-2016 v1.0
# find Chromosome column at runtime (structure has changed); v1.1 2017-05-12
# visit our Git: https://github.com/Nucleomics-VIB

# dependencies. Please install them if not yet available on your machine
suppressPackageStartupMessages(library("optparse")); # take care of command-line arguments
suppressPackageStartupMessages(library("openxlsx")); # read from excel files
suppressPackageStartupMessages(library("pheatmap")); # plot pretty heatmaps (different module than used at the NC)
suppressPackageStartupMessages(library("RColorBrewer")); # build nice color palettes

#####################################
### Handle COMMAND LINE arguments ###
#####################################

# parameters
option_list <- list(
  make_option(c("-i", "--infile"), type="character",
              help="input file name [REQUIRED]"),
  make_option(c("-s", "--signature"), type="character", default=NULL,
              help="signature file name (if absent, all data will be plotted)"),
  make_option(c("-o", "--outfile"), type="character", default="full_heatmap",
              help="base name for output [default: %default]"),
  make_option(c("-f", "--outformat"), type="integer", default=1,
              help="file format for output 1:PNG, 2:PDF [default: %default]")
  )

## parse options
opt <- parse_args(OptionParser(usage = "%prog [options] file",
                               option_list=option_list, add_help_option = TRUE))

# check if arguments provided
if ( is.null(opt$infile) ) {
  stop("no input file provided, run with --help")
}

  # check that infile exists
  if( file.access(opt$infile) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", opt$infile))
  }

  # check if signature argument provided and file exists
  if( !is.null(opt$signature) ) {
    if( file.access(opt$signature) == -1) {
      stop(sprintf("Specified file ( %s ) does not exist", opt$signature))
      }

    # load signature data in
    sig.name <- readLines(opt$signature, n=1)

    # check format or die
    if( startsWith(sig.name, "# ") ) {
      sig.name <- gsub("# ", "", sig.name)
    } else {
      stop("Signature first row should be '# signature_name' (only space between # and name)")
    }

    sig.vect <- read.table(opt$signature, skip=1, sep=",", header=FALSE)
    sig.vect <- as.vector(t(sig.vect))

    if(length(sig.vect)==0) {
      stop("Signature second row should be a comma-separated list of ENSEMBL-IDs")
    }

    # change opt$outfile
    opt$outfile <- paste(sig.name, "-signature-heatmap", sep="")
  }

##### load data in
excel.data <- read.xlsx(opt$infile, sheet=2)
attach(excel.data)

# create output plot file
if (opt$outformat==1){
  # png
  filename <- paste(opt$outfile, ".png", sep="")
} else {
  # pdf
  filename <- paste(opt$outfile, ".pdf", sep="")
}

##################
# compute and plot
# keep only data columns (remove last columns including "Chromosome")
chromosome.col <- which(colnames(excel.data)==as.vector("Chromosome"))
fpkm.data <- excel.data[,c(2,1,3:(chromosome.col-1))]

# remove some columns
row.names(fpkm.data) <- paste(fpkm.data[,1], fpkm.data[,2], sep=":")
fpkm.data <- fpkm.data[,-1]

# kick useless part of names for samples
colnames(fpkm.data) <- sub("@.*", "", colnames(fpkm.data))

# custom signatures or full_data

# custom signatures or full_data
if( !is.null(opt$signature) ) {
  # only signature genes are kept
  # select only signature rows and discard Gene.ID column to keep only FPKM in data.frame
  selection <- fpkm.data[fpkm.data$Gene.ID %in% sig.vect, 2:length(fpkm.data)]
  sig.name <- ifelse(exists("sig.name"), sig.name, "signature")
  main.title <- paste(sig.name,"signature", sep=" ")
} else {
  # full data gets plotted
  selection <- fpkm.data[,2:length(fpkm.data)]
  main.title <- "Full data heatmap"
}

# convert 0 to NA to avoid log-normalization error
selection[selection == 0] <- NA

# remove full NA rows
selection <- selection[rowSums(is.na(selection[,2:length(selection)]))
                       <length(selection[,2:length(selection)]),]

# color pallet
#col.pal <- brewer.pal(9,"Blues")

# define metrics for clustering
drows <- "euclidean"
dcols <- "euclidean"
clustmet <- "average"

# create heatmap
# type "?pheatmap()" for more help
hm.data <- log(selection, 2)

# prevent Rplot.pdf creation
pdf(NULL)

# plot
pheatmap(hm.data,
         filename=filename,
         color = (brewer.pal(9,"Blues")),
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
