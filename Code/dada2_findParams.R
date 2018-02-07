#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)

# # # # # # # # # # # # MAIN FUNCTION # # # # # # # # # # # # # # # # 
find_params <- function(R1,R2,name){
  # Where R1 and R2 are the file locations for the sequence files for a paired-end read
  setwd("/Users/student/Documents/PollardRotation/dada2/Data/test_img")
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # Look at quality profiles - everything look normal?
  pdf(paste0("QualityProfiles/",name,".pdf"))
  print(plotQualityProfile(R1))
  print(plotQualityProfile(R2))
  dev.off()
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # Filtering - choose a range of truncation lengths and plot
  filt_path <- "/Users/student/Documents/PollardRotation/InputData/filt_16S"
  
  trunc_list = c(75, 100, 125, 150, 200, 250)
  
  for (trunc in trunc_list) {
    
    # Just chose tutorial defaults...maybe play with this?
    filtF = paste0(filt_path,"/",name,"-truncLen",trunc,"F-filt.fastq.gz")
    filtR = paste0(filt_path,"/",name,"-truncLen",trunc,"-R-filt.fastq.gz")
    
    out <- filterAndTrim(R1, filtF, 
                         R2, filtR, 
                         truncLen=c(trunc,trunc),
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
    head(out)
    
    # Save out to a pdf
    pdf(paste0("ErrorRates/",name,"-truncLen",trunc,".pdf"))
    errF <- learnErrors(filtF, multithread=TRUE)
    print(plotErrors(errF, nominalQ=TRUE))
    
    errR <- learnErrors(filtR, multithread=TRUE)
    print(plotErrors(errR, nominalQ=TRUE))
    dev.off()
    
  }
  
}

# # # # # # # # # # # # RUN FUNCTION # # # # # # # # # # # # # # # # 

path = "/Users/student/Documents/PollardRotation/dada2/Data/raw_16S"

#TMM_AOMDSS_2014
R1 = "SRR4004921_pass_1.fastq.gz"
R2 = "SRR4004921_pass_2.fastq.gz"
name = "TMM_AOMDSS_2014"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)

#TMM_AOMDSS_2016
R1 = "SRR4417483_pass_1.fastq.gz"
R2 = "SRR4417483_pass_2.fastq.gz"
name = "TMM_AOMDSS_2016"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)

#Baxter_AOMDSS
R1 = "c1a_1399_d00.R1.fastq.gz"
R2 = "c1a_1399_d00.R2.fastq.gz"
name = "Baxter_AOMDSS"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)

#TMM_DSS
R1 = "SRR4423081_pass_1.fastq.gz"
R2 = "SRR4423081_pass_2.fastq.gz"
name = "TMM_DSS"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)

#UTS_DSS
R1 = "ERR1806597_pass_1.fastq.gz"
R2 = "ERR1806597_pass_2.fastq.gz"
name = "UTS_DSS"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)

# Helm_DSS
R1 = "SRR6127305_pass_1.fastq.gz"
R2 = "SRR6127305_pass_2.fastq.gz"
name = "Helm_DSS"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)

# UMAA_DSS
R1 = "SRR1914841_pass_1.fastq.gz"
R2 = "SRR1914841_pass_2.fastq.gz"
name = "UMAA_DSS"
find_params(paste0(path, "/",R1),paste0(path, "/", R2),name)
