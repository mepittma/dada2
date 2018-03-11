#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(decontam)

#base_path = "/pollard/home/mpittman/dada2/"
base_path = "/Users/student/Documents/PollardRotation/dada2/"
options(tz="America/Los_Angeles")

# # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Load in the seqtab_nochim files and attach metadata

read_seqtab <- function(seq,taxa,meta){ #add to_prune parameter?
  
  ps <- phyloseq(otu_table(seq, taxa_are_rows=FALSE), 
                 sample_data(meta), 
                 tax_table(taxa))
  #if (!missing(to_prune)){
    #ps <- prune_samples(sample_names(ps) != to_prune, ps) # Remove any desired samples
  #}
  
  return(ps)
}


# # # # # # # # # # # # # # # # # # # # # # # # 
#2. Decontamination
deco_ps <- function(ps, taxa, out_plot, out_tab){
  
  # Plot the library sizes - not sure if useful?
  pdf(out_plot)
  df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=response)) + geom_point()
  dev.off()
  
  # Remove contaminants based on negative controls
  sample_data(ps)$is.neg <- sample_data(ps)$response == "0"
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
  to_remove <- contamdf.prev[which(contamdf.prev$contaminant),]
  ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
  
  # Save out filtering info
  contam_tax <- taxa[which(row.names(taxa) %in% row.names(to_remove)),]
  write.table(contam_tax, file = out_tab, sep = "\t")
  
  return(ps.noncontam)
}

# # # # # # # # # # # # # # # # # # # # # # # # 
#3. Visualizations
visualize <- function(ps, var_list, name, suff, base_path){
  
  # Alpha diversity
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_",suff,"_alphaDiversity.pdf")
  pdf(out_file)
  for(var in var_list){
    print(plot_richness(ps, x=var, measures=c("Shannon", "Simpson"), color=var) + theme_bw())
  }
  dev.off()

  # Ordination
  ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_",suff,"_Ordination.pdf")
  pdf(out_file)
  
  for(var in var_list){
    print(plot_ordination(ps, ord.nmds.bray, color = var, title="Bray NMDS"))
  }
  dev.off()
  
  # Taxa Bar plot
  top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_",suff,"_TaxaBar.pdf")
  
  pdf(out_file)
  for(var in var_list){
    plot_bar(ps.top20, x=var, fill="Family") + facet_wrap(~response, scales="free_x")
  }
  dev.off()
  
}

# # # # # # # # # # # # # # # # # # # # # # # # 
#4. Save out feature/response matrices - both taxa and sequence abundances

make_mat <- function(out_path, ps, meta, name){
  
  # Extract abundance matrix from the phyloseq object
  OTU1 = as(otu_table(ps), "matrix")
  if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
  seq_df = as.data.frame(OTU1)
  
  # Get taxon abundance table
  OTU = as.data.frame(as(tax_table(ps), "matrix"))
  OTU$C = paste(OTU$Kingdom, OTU$Phylum, OTU$Class, OTU$Order, OTU$Family, 
                OTU$Genus, OTU$Species, sep="_")
  lookup <- cbind(row.names(OTU), OTU$C)
  lookup <- as.data.frame(lookup)
  OTU2 <- as.data.frame(OTU1)
  match(names(OTU2), lookup$V1)
  names(OTU2) = lookup$V2[match(names(OTU2), lookup$V1)]
  
  # Merge both to relevant metadata columns
  meta <- meta[,c("AOMDSS","DNR","DSS","TNBS","response")]
  OTU1 <- merge(OTU1,meta,by="row.names",all.x=TRUE)
  row.names(OTU1) <- OTU1$row.names
  OTU2 <- merge(OTU2,meta,by="row.names",all.x=TRUE)
  
  # Save out both matrices to scikitlearn folder
  write.table(OTU1, file = paste0(out_path,name,"_seq.tsv"), sep = "\t", quote = FALSE)
  write.table(OTU2, file = paste0(out_path,name,"_tax.tsv"), sep = "\t", quote = FALSE)
  
}


# # # # # # # # # # # # # # # # # # # # # # # # 
#5. Multi-function
phyloML <- function(base_path, out_path, name, var_list){
  
  # Read in necessary data
  seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
  taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
  meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)
  
  # Create phyloseq object, decontaminate it
  ps = read_seqtab(seq,taxa,meta)
  decon = deco_ps(ps,taxa, paste0(base_path, "Output/Decontam/",name,"_LibrarySize.pdf"),
                  paste0(base_path, "Output/Decontam/", name, "_removedTaxa.tsv"))
  
  # Visualize both
  visualize(ps, var_list, name, suff = "raw", base_path)
  visualize(decon, var_list, name, suff = "decon", base_path)
  
  # Save out matrices for ML pipeline
  make_mat(out_path, decon,meta,name)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Commands

base_path = "/Users/student/Documents/PollardRotation/dada2/"
out_path = "/Users/student/Documents/PollardRotation/scikit/Data/"

# Get each individual feature matrix
phyloML(base_path, out_path, "Baxter_AOMDSS", c("collection_date","inoculum","tumors","response"))
phyloML(base_path, out_path, "Helm_DSS", c("collection_date","inoculum","response"))
phyloML(base_path, out_path, "TMM_DSS", c("collection_date","response"))
phyloML(base_path, out_path, "UCSD_TNBS", c("collection_date","inoculum","response"))
phyloML(base_path, out_path, "UMAA_DSS", c("collection_date","response"))
phyloML(base_path, out_path, "UTA_TNBS", c("collection_date","response"))
phyloML(base_path, out_path, "UTS_DSS", c("inoculum","response"))
phyloML(base_path, out_path, "UCSF_DNR", c("response","genotype"))

#### Visualize Current DNR data with Pilot DNR data
# Load in UCSF test data; visualize with other UCSF data; save out test feature matrix

# Read in and combine the two metadata files
name <- "UCSF_DNR"
meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)
test_meta <- read.table(paste0(base_path, "MetaData/", name, "_TEST_processed.txt"), sep = "\t", header = TRUE)
test_meta$DNR <- "unknown"
meta$response<-NULL
c_meta <- rbind(meta, test_meta)

out_path = "/Users/student/Documents/PollardRotation/scikit/Data/OOS"
ps = read_seqtab(seq,taxa,test_meta)
# Remove contaminants based on negative controls
sample_data(ps)$is.neg <- sample_data(ps)$genotype == "WT"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
to_remove <- contamdf.prev[which(contamdf.prev$contaminant),]
decon <- prune_taxa(!contamdf.prev$contaminant, ps)
# Save out filtering info
contam_tax <- taxa[which(row.names(taxa) %in% row.names(to_remove)),]
write.table(contam_tax, 
            file = base_path, "Output/Decontam/", name, "_COMBINED_removedTaxa.tsv", sep = "\t")
visualize(ps, var_list, name, suff = "raw_OOS", base_path)
visualize(decon, var_list, name, suff = "decon_OOS", base_path)
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
seq_df = as.data.frame(OTU1)
# Get taxon abundance table
OTU = as.data.frame(as(tax_table(ps), "matrix"))
OTU$C = paste(OTU$Kingdom, OTU$Phylum, OTU$Class, OTU$Order, OTU$Family, 
              OTU$Genus, OTU$Species, sep="_")
lookup <- cbind(row.names(OTU), OTU$C)
lookup <- as.data.frame(lookup)
OTU2 <- as.data.frame(OTU1)
match(names(OTU2), lookup$V1)
names(OTU2) = lookup$V2[match(names(OTU2), lookup$V1)]
# Save out both matrices to scikitlearn folder - 
# ACTUALLY, no, because these will have wrong # features
# See below where it's added back into conglomerate
write.table(OTU1, file = paste0(out_path,name,"_seq.tsv"), sep = "\t", quote = FALSE)
write.table(OTU2, file = paste0(out_path,name,"_tax.tsv"), sep = "\t", quote = FALSE)

taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
ps = read_seqtab(seq,taxa,c_meta)

# Remove contaminants based on negative controls
sample_data(ps)$is.neg <- sample_data(ps)$DNR == "0"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
to_remove <- contamdf.prev[which(contamdf.prev$contaminant),]
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
# Save out filtering info
contam_tax <- taxa[which(row.names(taxa) %in% row.names(to_remove)),]
write.table(contam_tax, 
            file = base_path, "Output/Decontam/", name, "_COMBINED_removedTaxa.tsv", sep = "\t")

# Visualization
var_list = c("genotype","cohort")
visualize(ps, var_list, name, suff = "raw_COMBINED", base_path)
visualize(ps.noncontam, var_list, name, suff = "decon_COMBINED", base_path)

##################################
### Combine all training feature matrices
library(data.table)
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



# # # # # # # # # # # # # # # # # # 
# Visualize in PhyloSeq

# Reconstruct the model based on columns
seq_df$model <- NA
for (model in c("DNR","DSS","AOMDSS","TNBS")){
  seq_df$model[which(seq_df[[model]] == 1)] <- paste(model)
}

seq_df$model[which(is.na(seq_df$model))] <- "Control"

# Separate metadata
seq_df$DNR <- NULL
seq_df$DSS <- NULL
seq_df$AOMDSS <- NULL
seq_df$TNBS <- NULL
meta <- seq_df[,c("Row.names","response","model")]
row.names(meta) <- meta$Row.names

seq_df$Row.names <- NULL
seq_df$response <- NULL
seq_df$model <- NULL

# Create new phyloseq object, visualize diversity
ps <- phyloseq(otu_table(seq_df, taxa_are_rows=FALSE), 
               sample_data(meta))
var_list = c("response","model")

# Alpha diversity
out_file = paste0(base_path, "Output/PhyloSeq/Combined_alphaDiversity.pdf")
pdf(out_file)
for(var in var_list){
  print(plot_richness(ps, x=var, measures=c("Shannon", "Simpson"), color=var) + theme_bw())
}
dev.off()

# Ordination
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
out_file = paste0(base_path, "Output/PhyloSeq/Combined_Ordination.pdf")
pdf(out_file)

for(var in var_list){
  print(plot_ordination(ps, ord.nmds.bray, color = var, title="Bray NMDS"))
}
dev.off()