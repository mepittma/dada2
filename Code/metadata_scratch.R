#!/usr/bin/env Rscript

base_dir = "/Users/student/Documents/PollardRotation/dada2"


# # # # # # RERUN WITH LESS STRICT MERGING # # # # # # # # # #
# Load in the Helm data
name = "Helm_DSS"

helm = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))

# Only keep sample names and response
responses <- helm[,c("SRA_Sample","Treatment")]

# Load in the UTS data
name = "UTS_DSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Treatment")]


# # # # # # # # # # ACTUALLY READY # # # # # # # # # # 
# Load in the TMM_DSS data
name = "TMM_DSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Treatment")]

# Load in the UMAA data
name = "UMAA_DSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Treatment")]

# Load in the TMM_AOMDSS_2014 data
name = "TMM_AOMDSS_2014"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))

# Load in the TMM_AOMDSS_2016 data
name = "TMM_AOMDSS_2016"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))