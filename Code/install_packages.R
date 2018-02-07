#!/usr/bin/env Rscript

# Installs/updates the packages necessary to run the dada2 pipeline
source("https://bioconductor.org/biocLite.R")

alt_path = "/pollard/home/mpittman/apps/R_pkg"

install.packages("ggplot2", lib.loc=alt_path, lib=alt_path)
biocLite("dada2", lib.loc=alt_path, lib=alt_path)
biocLite("ShortRead", lib.loc=alt_path, lib=alt_path)



