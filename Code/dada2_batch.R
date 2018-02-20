#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)

# # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Get sample files and name

base_path = "/pollard/home/mpittman/dada2/"
#base_path = "/Users/student/Documents/PollardRotation/dada2/"

dada2 <- function(name, f_trunc, r_trunc, f_EE, r_EE){
  
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
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 2. QC/Filtering
  
  # Save an image to summarize quality profiles of the samples
  img_path = paste0(base_path,"Data/test_img")
  
  pdf(paste0(img_path,"/QualityProfiles/",name,"_QP.pdf"))
  for (sample in sample.names){
    print(plotQualityProfile(
      c(paste0(path,"/", sample ,"_pass_1.fastq.gz"), 
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
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 3. Dereplication
  
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 4. Sample inference
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 5. Merge paired reads
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 6. Construct tables
  
  out_path = paste0(base_path,"Output")
  silva_path = paste0(base_path,"Data/taxon")
  
  # Create and save out sequence tables
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  saveRDS(seqtab.nochim, 
          file = paste0(out_path,"/SeqTables/",name,"_seqtab_nochim.rds"))
  
  # Track reads through pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
  rownames(track) <- sample.names
  saveRDS(track, file=paste0(out_path,"/QC/",name,"_trackedReads.rds"))
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 7. Assign taxonomy
#  taxa <- assignTaxonomy(seqtab.nochim, 
#                         paste0(silva_path,"/silva_nr_v128_train_set.fa.gz"), multithread=TRUE)
#  taxa <- addSpecies(taxa, 
#                     paste0(silva_path,"/silva_species_assignment_v128.fa.gz"))
#  saveRDS(taxa, file = paste0(out_path, "/Taxa/",name,"_taxa_silva_plus.rds"))
}

# # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # #
dada2_single <- function(name, trunc, EE){
  
  path <- paste0(base_path,"Data/raw_16S/", name)
  
  # Forward and reverse filenames have format RUNID_pass_1 for forward, RUNID_pass_2 for reverse
  fnFs <- sort(list.files(path, pattern="_pass_1.fastq", full.names = TRUE))
  
  # Save sample names of files that will be analyzed
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 2. QC/Filtering
  
  # Save an image to summarize quality profiles of the samples
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
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 3. Dereplication
  
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  names(derepFs) <- sample.names
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 4. Sample inference
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 5. Construct tables
  
  out_path = paste0(base_path,"Output")
  silva_path = paste0(base_path,"Data/taxon")
  
  # Create and save out sequence tables
  seqtab <- makeSequenceTable(dadaFs)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  saveRDS(seqtab.nochim, 
          file = paste0(out_path,"/SeqTables/",name,"_seqtab_nochim.rds"))
  
  # Track reads through pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
  rownames(track) <- sample.names
  saveRDS(track, file=paste0(out_path,"/QC/",name,"_trackedReads.rds"))
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # 
  # 7. Assign taxonomy
 # taxa <- assignTaxonomy(seqtab.nochim, 
 #                        paste0(silva_path,"/silva_nr_v128_train_set.fa.gz"), multithread=TRUE)
 # taxa <- addSpecies(taxa, 
 #                    paste0(silva_path,"/silva_species_assignment_v128.fa.gz"))
 # saveRDS(taxa, file = paste0(out_path, "/Taxa/",name,"_taxa_silva_plus.rds"))
}

# # # # # # # # COMMANDS # # # # # # # # 

folders = c("Helm_DSS", "TMM_AOMDSS_2014", "UCSF_DNR", "Baxter_AOMDSS",
            "TMM_AOMDSS_2016", "TMM_DSS", "UTS_DSS", "UMAA_DSS")

# Completed: Helm_DSS
# dada2("Helm_DSS", 100, 50, 2, 2)
#dada2("TMM_AOMDSS_2014", 0, 0, 2, 2)
#dada2("TMM_AOMDSS_2016", 0, 0, 2, 2)
#dada2("TMM_DSS", 0, 0, 2, 2)
#dada2("UTS_DSS", 100, 100, 2, 2)
#dada2("UMAA_DSS", 0, 0, 2, 2)

# Run for first time: TNBS single-end samples
#dada2_single("UTA_TNBS",25,2)
dada2_single("UCSD_TNBS",0,2)

# Rerun: Helm and UTS; try running Baxter and UCSF
#dada2("Helm_DSS", 0, 0, 2, 2)
dada2("UTS_DSS", 25, 25, 2, 2)
dada2("Baxter_AOMDSS", 0, 0, 2, 2)
dada2("UCSF_DNR", 200, 150, 2, 2)