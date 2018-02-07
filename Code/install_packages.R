#!/usr/bin/env Rscript

# Installs/updates the packages necessary to run the dada2 pipeline
source("https://bioconductor.org/biocLite.R")

biocLite("dada2", lib = "/pollard/home/mpittman/apps/R_pkg")
biocLite("ShortRead", lib = "/pollard/home/mpittman/apps/R_pkg")

install.packages("ggplot2", lib = "/pollard/home/mpittman/apps/R_pkg")

