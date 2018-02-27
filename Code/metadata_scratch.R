#!/usr/bin/env Rscript

base_dir = "/Users/student/Documents/PollardRotation/dada2"


# # # # # # # # # # # # # #  BAXTER # # # # # # # # # # # # # # 
name = "Baxter_AOMDSS"
data = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  row.names = 1, sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rename the underscores in the metadata to hyphens (to accomadate necessary filename changes...)
metanames <- gsub("_", "-", row.names(data), fixed=TRUE)
meta <- 
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis


# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD


# Identify other columns that might be useful

# Save out as tab-delimited text file with the format name_processed.txt

# # # # # # # # # # # # # # HELM # # # # # # # # # # # # # # 
# Load in the Helm data
name = "Helm_DSS"
helm = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), sep = '\t', header = TRUE)
# Confirm this worked
QC <- readRDS(paste0(base_dir, "/Output/QC/",name,"_trackedReads.rds"))
# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Only keep sample names and response
responses <- helm[,c("SRA_Sample","Treatment")]
