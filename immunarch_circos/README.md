[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
==========

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/)**

## Installing and loading the required R libraries

**[BiocManager](https://www.bioconductor.org/install/)** is used here to install packages as it is often moore efficient than R at doing it.

```{r setup}
# Install the required package
# Define the list of required packages
packages <- c("immunarch", "circlize", "stringr")

# Check which packages are not installed
to_install <- packages[!packages %in% installed.packages()[, "Package"]]

# Install only the missing packages
if (length(to_install) > 0) {
  BiocManager::install(to_install)
} else {
  message("All packages are already installed.")
}

# load the corresponding libraries
library("circlize")
library("stringr")
library("immunarch")

# load demo data from immunarch
data(immdata)

# select one dataset
mydat <- immdata$data$`A2-i129`
```

## Custom function

The example below takes the demo data from the **[immunarch](https://github.com/immunomind/immunarch)** package and creates a plot.

The function below should create an object compatible with **circlize** plotting (ref:**[https://github.com/immunomind/immunarch/issues/103](https://github.com/immunomind/immunarch/issues/103)**).

```{r custom function}
# custom function to plot a circos from immdata
# © Stephane Plaisance, VIB Nucleomics Core, 2020/08/24

vis_genus <- function(x, title="", ...) {
# require("stringr", "circlize")
V <- sort(unique(unlist(strsplit(paste(x$V.name, collapse=","), ","))))
J <- sort(unique(unlist(strsplit(paste(x$J.name, collapse=","), ","))))

# create and fill matrix of pairwise with sum(Clones)
mat <- matrix(nrow = length(J), ncol = length(V))
colnames(mat) <- V
rownames(mat) <- J
for (j in seq(1:length(J))) {
  for (v in seq(1:length(V))) {
  sub <- x[(str_detect(x$V.name, V[[v]], negate = FALSE) & 
               str_detect(x$J.name, J[[j]], negate = FALSE)),]
  mat[j,v] <- sum(sub$Clones)
  }
}

# convert to proportions / fraction of total
mat <- mat/sum(mat)

# reorder by V and J decreasing order (left to right)
vmax <- colSums(mat)
jmax <- rowSums(mat) 
mat2 <- mat[order(jmax, decreasing = FALSE), order(vmax, decreasing = TRUE)]

# use circlize
chordDiagram(mat2, annotationTrack = "grid", ...)

# add legends
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 
              CELL_META$ylim[1], 
              CELL_META$sector.index, 
              facing = "clockwise", 
              niceFacing = TRUE,
              adj = c(-0.25, 0.5))
}, bg.border = NA)

# add title
title(title)
}
```

## Plotting

```{r plot, fig.height=8, fig.width=8}

circos.clear()
circos.par("canvas.xlim" = c(-1.25, 1.25), "canvas.ylim" = c(-1.25, 1.25))
vis_genus(mydat, title="A2-i129_sample")
```

<img src="pictures/circus_plot.png?raw=true" alt="circus_plot.png" style="width: 600px;"/>

*[[back-to-top](#top)]*  

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>
