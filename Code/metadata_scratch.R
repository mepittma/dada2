#!/usr/bin/env Rscript

base_dir = "/Users/student/Documents/PollardRotation/dada2"


# # # # # # # # # # # # # #  BAXTER # # # # # # # # # # # # # # 
name = "Baxter_AOMDSS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  row.names = 1, sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rename the underscores in the metadata to hyphens (to accomadate necessary filename changes...)
metanames <- gsub("_", "-", row.names(data), fixed=TRUE)
row.names(meta) <- metanames
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
meta <- meta[grep("Mus mus", meta$feature),]
meta <- meta[,c("collection_date", "subspecf_gen_lin", "host_subject_id", "inoculum", "disease_stat", "tumors")]
names(meta) <- c("collection_date", "genotype", "Library_Name","inoculum", "disease_stat")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
names(meta)[names(meta) == 'disease_stat'] <- 'response'
meta$response <- gsub('aom/dss', '1', meta$response)
meta$response <- gsub('normal', '0', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)



# # # # # # # # # # # # # # HELM # # # # # # # # # # # # # # 
# Load in the Helm data
name = "Helm_DSS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  row.names = 1, sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
row.names(meta) <- meta$Run
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
meta <- meta[,c("collection_date", "host_genotype", "Library_Name", "Microbiota", "Treatment")]
names(meta) <- c("collection_date", "genotype", "Library_Name", "inoculum", "Treatment")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
names(meta)[names(meta) == 'Treatment'] <- 'response'
meta$response <- gsub('DSS', '1', meta$response)
meta$response <- gsub('none', '0', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)


# # # # # # # # # # # # # # TMM AOM/DSS 2014 - needs response clarification # # # # # # # # # # # # # # 
#name = "TMM_AOMDSS_2014"
#meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
#                  sep = '\t', header = TRUE)

# Load in sequence abundances
#seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
#row.names(meta) <- meta$Run
# Confirm that the two can be merged
#merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
#meta <- meta[,c("collection_date", "Sample_Name", "env_feature")]

# # # # # # # # # # # # # # TMM AOM/DSS 2016 - needs response clarification # # # # # # # # # # # # # # 
#name = "TMM_AOMDSS_2016"
#meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
#                  sep = '\t', header = TRUE)

# Load in sequence abundances
#seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
#row.names(meta) <- meta$Run
# Confirm that the two can be merged
#merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
#meta <- meta[,c("collection_date", "Sample_Name", "env_feature")]

# # # # # # # # # # # # # # TMM DSS # # # # # # # # # # # # # # 
name = "TMM_DSS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
row.names(meta) <- meta$Run
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
meta <- meta[grep("WT_|WT\\+DSS", meta$Sample_Name),]
meta <- meta[,c("collection_date", "Sample_Name", "env_feature")]
names(meta) <- c("collection_date","Library_Name","response")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
meta$response <- gsub('AOM/DSS colon cancer induction', '0', meta$response)
meta$response <- gsub('DSS colon colitis induction', '1', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)


# # # # # # # # # # # # # # UCSD TNBS # # # # # # # # # # # # # # 
name = "UCSD_TNBS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
row.names(meta) <- meta$Run
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
meta <- meta[grep("stool", meta$sample_type),]
meta <- meta[,c("collection_timestamp", "Library_Name", "anonymized_name", "host_subject_id")]
names(meta) <- c("collection_date","Library_Name","inoculum","response")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
meta$response <- gsub('.*A.*', '0', meta$response)
meta$response <- gsub('.*F.*', '1', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)


# # # # # # # # # # # # # # UCSF DNR # # # # # # # # # # # # # # 
name = "UCSF_DNR"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rename samples from underscores to hyphen
metanames <- gsub("_", "-", meta$sample)
row.names(meta) <- metanames
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")
row.names(meta) <- meta$sample
names(meta) <- c("Library_Name","genotype","cohort")
meta$response <- meta$genotype

# Create two separate metadata sets, since the "current" cohort is our out-of-sample data
pilot_meta <- meta[which(meta$cohort == "Pilot"), ]
current_meta <- meta[which(meta$cohort == "Current"),]

# Create response column for the pilot metadata
pilot_meta$response <- gsub('WT', '0', pilot_meta$response)
pilot_meta$response <- gsub('DNR', '1', pilot_meta$response)
current_meta$response <- gsub('WT|DNR', 'unknown', current_meta$response)

# Recombine
meta <- rbind(pilot_meta,current_meta)

# Save out both
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), 
            sep="\t", quote=FALSE)


# # # # # # # # # # # # # # UMAA DSS # # # # # # # # # # # # # # 
name = "UMAA_DSS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
row.names(meta) <- meta$Run
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples not relevant to my hypothesis
meta <- meta[,c("collection_date", "genotype", "mouse_id", "Sample_Name")]
names(meta) <- c("collection_date", "genotype", "Library_Name","response")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
meta$response <- gsub('.*DSS.*', '1', meta$response)
meta$response <- gsub('.*NoAbs.*', '0', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)


# # # # # # # # # # # # # # UTA TNBS # # # # # # # # # # # # # # 
name = "UTA_TNBS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
row.names(meta) <- meta$Run
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples/categories not relevant to my hypothesis
meta <- meta[,c("collection_date", "Sample_Name", "Library_Name")]
names(meta) <- c("collection_date", "Library_Name","response")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
meta$response <- gsub('CONTROL', '0', meta$response)
meta$response <- gsub('TNBS', '1', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)


# # # # # # # # # # # # # # UTS DSS # # # # # # # # # # # # # # 
name = "UTS_DSS"
meta = read.table(paste0(base_dir, "/MetaData/",name,"_metadata.txt"), 
                  sep = '\t', header = TRUE)

# Load in sequence abundances
seqs <- readRDS(paste0(base_dir, "/Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Rownames should be SRR accession
row.names(meta) <- meta$Run
# Confirm that the two can be merged
merged <- merge(seqs, meta, by="row.names")

# Remove samples/categories not relevant to my hypothesis
meta <- meta[,c("Sample_Name", "Group", "Group")]
names(meta) <- c("Library_Name","inoculum","response")

# Make sure to name Response column  "response", with 0 = no IBD and 1 = IBD
meta$response <- gsub('Mock|Sodium tungstate|DSS\\+sodium tungstate', '0', meta$response)
meta$response <- gsub('DSS', '1', meta$response)

# Save out as tab-delimited text file with the format name_processed.txt
write.table(meta, file=paste0(base_dir,"/MetaData/",name,"_processed.txt"), sep="\t", quote=FALSE)

