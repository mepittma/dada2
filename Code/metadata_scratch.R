#!/usr/bin/env Rscript

base_dir = "/Users/student/Documents/PollardRotation/dada2"

# # # # # # # # # # # # # # NEEDS A FILE # # # # #
# Load in the Baxter data
name = "Baxter_AOMDSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Sample_Name")]


# # # # # # # # # # # # # # NEEDS A FILE # # # # #
# Load in the UCSF data
name = "UCSF_DNR"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Sample_Name")]


# # # # # # # # # # # # # # 
# Load in the Helm data
name = "Helm_DSS"
helm = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- helm[,c("SRA_Sample","Treatment")]

# # # # # # # # # # # # # # 
# Load in the UTS data
name = "UTS_DSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Group")]

# # # # # # # # # # # # # # 
# Load in the UTA_TNBS data
name = "UTA_TNBS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Library_Name")]

# # # # # # # # # # # # # # 
# Load in the TMM_DSS data
name = "TMM_DSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Sample_Name")]

# # # # # # # # # # # # # # 
# Load in the UMAA data
name = "UMAA_DSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Sample_Name")]

# # # # # # # # # # # # # # NEEDS MORE INFO # # # # # # # # # # # # # # 
# Load in the UCSD_TNBS
name = "UCSD_TNBS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep stool samples
data <- data[which(data$sample_type == "stool"),]
# Only keep sample names and response
responses <- data[,c("SRA_Sample","Group")]

# # # # # # # # # # # # # # THIS ONE NEEDS CLARIFICATION # # # # # # # # # # # # # # 
# Load in the TMM_AOMDSS_2014 data
name = "TMM_AOMDSS_2014"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
#responses <- data[,c("SRA_Sample","Sample_Name")]

# # # # # # # # # # # # # # THIS ONE NEEDS CLARIFICATION # # # # # # # # # # # # # # 
# Load in the TMM_AOMDSS_2016 data
name = "TMM_AOMDSS_2016"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))