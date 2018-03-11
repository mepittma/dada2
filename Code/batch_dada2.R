#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)

base_path = "/pollard/home/mpittman/dada2/"

# # # # # # # # 1. Get sample files and name # # # # # # # #

paired_read <- function(name){
  path <- paste0(base_path,"Data/raw_16S/", name)
  
  # Forward and reverse filenames have format RUNID_pass_1 for forward, RUNID_pass_2 for reverse
  fnFs <- sort(list.files(path, pattern="_pass_1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="_pass_2.fastq", full.names = TRUE))
  
  # If there are samples with only a forward read, exclude them
  if (length(fnFs) != length(fnRs)){
    
    # Find and select names in F but not in R
    names_f <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    names_r <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)
    only_one <- setdiff(names_f, names_r)
    
    # Add back the filenames so we can remove them from our list
    to_del_f <- sapply(only_one, function(x) paste0(path,"/",x,"_pass_1.fastq.gz"))
    to_del_r <- sapply(only_one, function(x) paste0(path,"/",x,"_pass_2.fastq.gz"))
    
    # Remove any singletons from both vectors
    fnFs <- fnFs[!fnFs %in% to_del_f]
    fnRs <- fnRs[!fnRs %in% to_del_r]
  }
  
  # Save sample names of files that will be analyzed
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  # Return sample.names and fnFs and fnRs
  return(list(fnFs,fnRs,sample.names))
  
}

single_read <- function(name){
  path <- paste0(base_path,"Data/raw_16S/", name)
  fnFs <- sort(list.files(path, pattern="_pass_1.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  return(list(fnFs, sample.names))
}

# # # # # # # # # # # # # # # # # # # # # # # # 
# 2. QC/Filtering

paired_filt <- function(name, sample.names, f_trunc, r_trunc){
  
  # Save an image to summarize quality profiles of the samples
  img_path = paste0(base_path,"Data/test_img")
  
  pdf(paste0(img_path,"/QualityProfiles/",name,"_QP.pdf"))
  for (sample in sample.names){
    print(plotQualityProfile(c(paste0(path,"/", sample ,"_pass_1.fastq.gz"),
                               paste0(path,"/", sample ,"_pass_2.fastq.gz"))
    ))
    }
  dev.off()
  
  # Filter the data
  filt_path <- paste0(base_path,"Data/filt_16S/",name)
  
  filtFs <- file.path(filt_path,"/", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path,"/", paste0(sample.names, "_R_filt.fastq.gz"))
  
  write(paste0("fnFs: ", fnFs),stderr())
  write(paste0("filtFs: ", filtFs), stderr())
  write(paste0("fnRs: ", fnRs), stderr())
  write(paste0("filtRs: ", filtRs), stderr())
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(f_trunc,r_trunc),
                       maxN=0, maxEE=c(f_EE,r_EE), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  
  # Save an image to summarize the error rates of the samples after filtering
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  
  pdf(paste0(img_path,"/ErrorRates/",name,"_ER.pdf"))
  print(plotErrors(errF, nominalQ=TRUE))
  print(plotErrors(errR, nominalQ=TRUE))
  dev.off()
  
  # Return error variables
  return(list(errF, errR))
}

single_filt <- function(name, sample.names, trunc){

  img_path = paste0(base_path,"Data/test_img")
  pdf(paste0(img_path,"/QualityProfiles/",name,"_QP.pdf"))
  for (sample in sample.names){
    print(plotQualityProfile(
      paste0(path,"/", sample ,"_pass_1.fastq.gz")
    ))
  }
  dev.off()
  
  # Filter the data
  filt_path <- paste0(base_path,"Data/filt_16S/",name)
  filtFs <- file.path(filt_path,"/", paste0(sample.names, "_F_filt.fastq.gz"))
  write(paste0("fnFs: ", fnFs),stderr())
  write(paste0("filtFs: ", filtFs), stderr())
  
  out <- filterAndTrim(fnFs, filtFs, truncLen=trunc,
                       maxN=0, maxEE=EE, truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  
  # Save an image to summarize the error rates of the samples after filtering
  errF <- learnErrors(filtFs, multithread=TRUE)
  
  pdf(paste0(img_path,"/ErrorRates/",name,"_ER.pdf"))
  print(plotErrors(errF, nominalQ=TRUE))
  dev.off()
  
  # Return error rates
  return(errF)
}

# # # # # # # # # # # # # # # # # # # # # # # # 
# 3. Dereplication, Sample Inference, Merging (if applicable)

paired_inference <- function(filtFs, filtRs, errF, errR, sample.names){

  # Dereplication
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  #  Sample inference
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  # Merge paired reads
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  
  # return table input
  return(mergers)
}

single_inference <- function(filtFs, errF, sample.names){
  
  # Dereplication
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  names(derepFs) <- sample.names
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  
  return(dadaFs)

}

# # # # # # # # # # # # # # # # # # # # # # # # 
# 4. Create the tables

get_tabs <- function(seqs, track_reads){
  
  out_path = paste0(base_path,"Output")
  silva_path = paste0(base_path,"Data/taxon")
  
  # Create and save out sequence tables
  seqtab <- makeSequenceTable(seqs)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  saveRDS(seqtab.nochim,file = paste0(out_path,"/SeqTables/",name,"_seqtab_nochim.rds"))
  
  # Track reads through pipeline
  getN <- function(x) sum(getUniques(x))
  if (track_reads == "full") {
    track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
    
  } else if (track_reads == "full_single") {
    track <- cbind(out, getN(dadaFs), rowSums(seqtab), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
    
  } else if (track_reads == "part") {
    track <- cbind(sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
    colnames(track) <- c("denoised", "merged", "tabled", "nonchim")
  }
  
  rownames(track) <- sample.names
  saveRDS(track, file=paste0(out_path,"/QC/",name,"_trackedReads.rds"))

  # Taxa tables
  taxa <- assignTaxonomy(seqtab.nochim, paste0(silva_path,"/silva_nr_v128_train_set.fa.gz"), 
                         multithread=TRUE, tryRC = TRUE)
  taxa <- addSpecies(taxa,paste0(silva_path,"/silva_species_assignment_v128.fa.gz"))
  saveRDS(taxa, file = paste0(out_path, "/Taxa/",name,"_taxa_silva_plus.rds"))
  
}

# # # # # # # # # # # # # # # # # # # # # # # # 
# Multi-functions
dada_paired <- function(name, f_trunc, r_trunc, f_EE, r_EE){
  
  # Read in names
  name_list <- paired_read(name)
  fnFs <- name_list[1]
  fnRs <- name_list[2]
  sample.names <- name_list[3]
  
  # Filter data
  error_list <- paired_filt(name, sample.names, f_trunc, r_trunc, f_EE, r_EE)
  errF <- error_list[1]
  errR <- error_list[2]
  
  # Inference and table creation
  seqs <- paired_inference(filtFs, filtRs, errF, errR, sample.names)
  get_tabs(seqs, track_reads = "full")
}

dada_single <- function(name, f_trunc, f_EE){
  
  # Read in names
  name_list <- paired_read(name)
  fnFs <- name_list[1]
  sample.names <- name_list[2]
  
  # Filter data
  errF <- paired_filt(name, sample.names, f_trunc, f_EE)
  
  # Inference and table creation
  seqs <- paired_inference(filtFs, errF, sample.names)
  get_tabs(seqs, track_reads = "full_single")
}

dada_prefilt <- function(name){
  filt_path <- paste0(base_path,"Data/filt_16S/sl_",name)
  img_path = paste0(base_path,"Data/test_img")
  
  # Forward and reverse filenames have format RUNID_pass_1 for forward, RUNID_pass_2 for reverse
  filtFs <- sort(list.files(filt_path, pattern="_F_filt.fastq.gz", full.names = TRUE))
  filtRs <- sort(list.files(filt_path, pattern="_R_filt.fastq.gz", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  
  # Save an image to summarize the error rates of the samples after filtering
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  
  pdf(paste0(img_path,"/ErrorRates/",name,"_ER.pdf"))
  print(plotErrors(errF, nominalQ=TRUE))
  print(plotErrors(errR, nominalQ=TRUE))
  dev.off()
  
  seqs <- paired_inference(filtFs, filtRs, errF, errR, sample.names)
  get_tabs(seqs, track_reads = "part")
  
}

# # # # # # # # COMMANDS # # # # # # # # 

#dada_paired("Baxter_AOMDSS", 190, 170, 2, 2)
dada_paired("Helm_DSS", 200, 100, 2, 2)
#dada_single("UCSD_TNBS",240,2)
#dada_paired("UMAA_DSS", 240, 170, 2, 2)
#dada_single("UTA_TNBS",300,2)
#dada_paired("UTS_DSS", 160, 230, 2, 2)
#dada_paired("TMM_DSS", 290, 200, 2, 5)
#dada_prefilt("UCSF_DNR")
dada_paired("UCSD_IL10", 190,170,2,2)
