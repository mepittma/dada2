#!/usr/bin/env Rscript

library(data.table)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # CREATION OF FEATURE MATRICES # # # # # # # # # # # # # 

# Function to combine the datasets for both raw and debatched data, OTUs and taxa
write_filtSeq_mat <- function(list_list){
  
  # Establish index of list to grab
  list_ind = 2
  
  # Establish data.frame to append to
  combined_seq <- data.frame()
  combined_meta <- data.frame()
  combined_tax <- data.frame()
  
  # For each dataset in the list, append the relevant data type
  for (ps_list in list_list){
    
    # Access relevant variables
    ps <- ps_list[list_ind]
    seq <- as.data.frame(otu_table(ps))
    meta <- as.data.frame(sample_data(ps))
    taxa <- as.data.frame(tax_table(ps))
    
    # Append to df
    seq_df <- rbind(combined_seq, seq)
    tax_df <- rbind(combined_tax, taxa)
  }
  
  # Create new phyloseq object from merged dataframes
  
  # Debatch merged ps
  
  # Create taxon feature file
  taxa$C = paste(taxa$Kingdom, taxa$Phylum, taxa$Class, taxa$Order, taxa$Family, 
                 taxa$Genus, sep="_")
  lookup <- as.data.frame(cbind(row.names(taxa), taxa$C))
  OTU2 <- as.data.frame(OTU1)
  names(OTU2) = lookup$V2[match(names(OTU2), lookup$V1)]
  
  # Merge both to relevant metadata columns
  meta <- meta[,c("AOMDSS","DNR","DSS","TNBS","response","study")]
  OTU1 <- merge(OTU1,meta,by="row.names",all.x=TRUE)
  row.names(OTU1) <- OTU1$row.names
  OTU2 <- merge(OTU2,meta,by="row.names",all.x=TRUE)
  
  # Save out datasets
  write.table(OTU1, file = paste0(out_path,name,"_seq.tsv"), sep = "\t", quote = FALSE)
  write.table(OTU2, file = paste0(out_path,name,"_tax.tsv"), sep = "\t", quote = FALSE)
  
}

# # # # # COMMAND # # # # # 

# Define starting variables
feat_path = paste0(base_path, "Output/FeatMat")
list_list = list(Bax_list, Helm_list, TMM_list, UCTN_list, UCIL_list, UMAA_list, 
                 UTA_list, UTS_list, UCSF_list)

# Call the function to create feature matrices

# Desired feature matrices

# 1. Raw, all data (including UCSF test data), both taxon abundances and OTUs
  # NOTE! make sure taxon abundances are being added when a merging event occurs
# 2. 

# Considerations: Only model type

seq_df <- data.frame()
tax_df <- data.frame()

name_list = c("Baxter_AOMDSS","Helm_DSS","TMM_DSS","UCSD_TNBS",
              "UMAA_DSS","UTA_TNBS","UTS_DSS","UCSF_DNR")
for(name in name_list){
  
  # Read in the file
  seq_data <- read.table(paste0(
    "/Users/student/Documents/PollardRotation/scikit/Data/",name,"_seq.tsv"))
  tax_data <- read.table(paste0(
    "/Users/student/Documents/PollardRotation/scikit/Data/",name,"_tax.tsv"))
  
  # add a column indicating study identity
  seq_data$study <- paste(name)
  tax_data$study <- paste(name)
  
  # merge with df, creating 0s where no data for the sequences exists (note: is this ok?)
  seq_df <- rbind(setDT(seq_data), setDT(seq_df), fill=TRUE)
  tax_df <- rbind(setDT(tax_data), setDT(tax_df), fill=TRUE)
  
}

# Read in the test data and merge that too
seq_tdata <- read.table("/Users/student/Documents/PollardRotation/scikit/Data/OOSUCSF_DNR_seq.tsv")
tax_tdata <- read.table("/Users/student/Documents/PollardRotation/scikit/Data/OOSUCSF_DNR_tax.tsv")
seq_tdata$study <- "DNR_test"
tax_tdata$study <- "DNR_test"
### REMOVE THIS SOMEDAY ###
seq_df$Row.names <- NULL
tax_df$Row.names <- NULL

seq_df <- rbind(setDT(seq_tdata), setDT(seq_df), fill=TRUE)
tax_df <- rbind(setDT(tax_tdata), setDT(tax_df), fill=TRUE)

seq_df[is.na(seq_df)] <- 0
tax_df[is.na(tax_df)] <- 0

# Now take out the rows that were part of the test data
test_seq <- seq_df[which(seq_df$study == "DNR_test"),]
test_tax <- tax_df[which(tax_df$study == "DNR_test"),]

tr_seq_df <- seq_df[which(seq_df$study != "DNR_test"),]
tr_tax_df <- tax_df[which(tax_df$study != "DNR_test"),]

# Make one for UCSF_DNR (not test version)
ucsf_seq <- seq_df[which(seq_df$study == "UCSF_DNR"),]
ucsf_tax <- seq_df[which(tax_df$study == "UCSF_DNR"),]

# Remove study column
test_seq$study <- NULL
test_tax$study <- NULL
tr_seq_df$study <- NULL
tr_tax_df$study <- NULL
ucsf_seq$study <- NULL
ucsf_tax$study <- NULL

# Save out the feature files
write.table(tr_seq_df, file = "/Users/student/Documents/PollardRotation/scikit/Data/seq_FF.tsv", sep = "\t", row.names = FALSE)
write.table(tr_tax_df, file = "/Users/student/Documents/PollardRotation/scikit/Data/tax_FF.tsv", sep = "\t", row.names = FALSE)
write.table(test_seq, file = "/Users/student/Documents/PollardRotation/scikit/Data/test_seq_FF.tsv", sep = "\t", row.names = FALSE)
write.table(test_tax, file = "/Users/student/Documents/PollardRotation/scikit/Data/test_tax_FF.tsv", sep = "\t", row.names = FALSE)

write.table(ucsf_seq, file = "/Users/student/Documents/PollardRotation/scikit/Data/UCSF_DNR_seq_FF.tsv", sep = "\t", row.names = FALSE)
write.table(ucsf_tax, file = "/Users/student/Documents/PollardRotation/scikit/Data/UCSF_DNR_tax_FF.tsv.tsv", sep = "\t", row.names = FALSE)
