#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)

#base_path = "/Users/student/Documents/PollardRotation/dada2/"
base_path = "/pollard/home/mpittman/dada2/"
silva_path = paste0(base_path, "Data/taxon")
out_path = paste0(base_path,"Output")

# # # # # # # # # # # # # #  BAXTER # # # # # # # # # # # # # # 
name = "Baxter_AOMDSS"

# Load in sequence abundances
seqs <- readRDS(paste0(base_path, "Output/SeqTables/",name,"_seqtab_nochim.rds"))
# Load in the QC table
QC <- readRDS(paste0(base_path, "Output/QC/", name, "_trackedReads.rds"))

# No taxa. Let's try to run it locally and see what errors we get.
taxa <- assignTaxonomy(seqs, paste0(silva_path,"/silva_nr_v128_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa,paste0(silva_path,"/silva_species_assignment_v128.fa.gz"))
head(taxa)
saveRDS(taxa, file = paste0(out_path, "/Taxa/",name,"_taxa_silva_plus.rds"))

# Load in the RDS file just created
taxa_l <- readRDS(paste0(out_path, "/Taxa/",name,"_taxa_silva_plus.rds"))
