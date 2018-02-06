#!/usr/bin/env Rscript

# Installs/updates the packages necessary to run the dada2 pipeline
source("https://bioconductor.org/biocLite.R")

biocLite("dada2")
biocLite("ShortRead")

install.packages("ggplot2")

