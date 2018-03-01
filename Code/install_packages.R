#!/usr/bin/env Rscript

# Installs/updates the packages necessary to run the dada2 pipeline
source("http://bioconductor.org/biocLite.R")

install.packages("ggplot2",repos = "http://cran.us.r-project.org")
biocLite("dada2")
biocLite("ShortRead")
biocLite("phyloseq")
biocLite("decontam")

library(devtools)
install_github("benjjneb/decontam")
