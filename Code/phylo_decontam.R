#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(decontam)

#base_path = "/pollard/home/mpittman/dada2/"
base_path = "/Users/student/Documents/PollardRotation/dada2/"
options(tz="America/Los_Angeles")

name = "Baxter_AOMDSS"
seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)
ps = read_seqtab(seq,taxa,meta)
decon = deco_ps(ps,taxa, paste0(base_path, "Output/Decontam/",name,"_LibrarySize.pdf"),
                paste0(base_path, "Output/Decontam/", name, "_removedTaxa.tsv"))
visualize(ps,var_list = c("collection_date","inoculum","response","tumors"),name,base_path)

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
visualize <- function(ps, var_list, name, base_path){
  
  # Alpha diversity
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_alphaDiversity.pdf")
  pdf(out_file)
  for(var in var_list){
    print(plot_richness(ps, x=var, measures=c("Shannon", "Simpson"), color=var) + theme_bw())
  }
  dev.off()

  # Ordination
  ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_Ordination.pdf")
  pdf(out_file)
  
  for(var in var_list){
    print(plot_ordination(ps, ord.nmds.bray, color = var, title="Bray NMDS"))
  }
  dev.off()
  
  # Taxa Bar plot
  top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_TaxaBar.pdf")
  
  for(var in var_list){
    plot_bar(ps.top20, x=var, fill="Family") + facet_wrap(~response, scales="free_x")
  }
  dev.off()
  
}

# # # # # # # # # # # # # # # # # # # # # # # # 
#4. Save out feature/response matrices - both taxa and sequence abundances

make_mat <- function(ps, meta, name){
  out_file = paste0(base_dir,"Output/FeatMat/",name,".tsv")
  
  # Extract abundance matrix from the phyloseq object
  OTU1 = as(otu_table(ps), "matrix")
  if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
  seq_df = as.data.frame(OTU1)
  
  # Repeat with the taxon abundances
  
  # Create dummy variables for IBD model
  
  # For both, append the collection_date and 
  
  
}

# Convert collection_date to number of days
# Create dummy columns for encoding of the model type
# Create matrix
# Save out table

# # # # # # # # # # # # # # # # # # # # # # # # 
#5. Multi-function



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Commands

#For each individual dataset
# Read in the data

# check sample variables
sample_variables(ps)

# Precontam visualization

# Decontam

# Post-decontam visualization

# Save out both sequence and taxa abundance information

