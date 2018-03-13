#!/usr/bin/env Rscript

# Installs/updates the packages necessary to run the dada2 pipeline
source("http://bioconductor.org/biocLite.R")

install.packages("ggplot2",repos = "http://cran.us.r-project.org")
install.packages("RAM",repos = "http://cran.us.r-project.org")
biocLite("dada2")
biocLite("ShortRead")
biocLite("phyloseq")
biocLite("sva")
biocLite("DECIPHER")
biocLite("phangorn")

install.packages("devtools", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library(devtools)
install_github("benjjneb/decontam")
