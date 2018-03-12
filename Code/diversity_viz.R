#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(decontam)

#base_path = "/pollard/home/mpittman/dada2/"
base_path = "/Users/student/Documents/PollardRotation/dada2/"
options(tz="America/Los_Angeles")

# Outline of what I want to do
# 1. Read in the dada2 output
# 2. Decontaminate
# 3. Remove batch effects
# 4. Visualize the above three cases
  # o. qplot histogram to show abundance discrepancies among datasets
  # a. PCoA - Bray-Curtis and Unifrac
    # i. with each dataset individually
    # ii. with the combined datasets
  # b. Major taxon abundance
    # i. with each dataset individually
    # ii. with the combined datasets
# 5. Visualize with just decontam/unbatched data
  # a. Tree of life with the combined datasets
  # b. Heatmap
    # i. with each dataset individually, comparing across variables
    # ii. with all datasets, comparing across 
# 6. Feature matrices both with OTUs and taxa

name = "Baxter_AOMDSS"
# Read in necessary data
seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)

# Create phyloseq object, decontaminate it
ps = read_seqtab(seq,taxa,meta)
decon = deco_ps(ps,taxa, paste0(base_path, "Output/Decontam/",name,"_LibrarySize.pdf"),
                paste0(base_path, "Output/Decontam/", name, "_removedTaxa.tsv"))

# # # # # # # # # # # # # # # # # # # # # # # # 
# General functions

# Create heatmap from matrix

# Create dendrogram/heatmap from matrix

# Create taxa-tree from matrix

# Create taxa heatmap from matrix

# PCoA plot
meta_fact <- meta
meta_fact$response <- gsub('0', 'no', meta_fact$response)
meta_fact$response <- gsub('1', 'yes', meta_fact$response)
meta_fact$response <- as.factor(meta_fact$response)
pcoa.plot(as.data.frame(seq), is.OTU=TRUE, meta_fact, factors=c(Diet="inoculum", IBD="response"))

# # # # # # # # # # # # # # # # # # # # # # # # 
# Dataset-specific functions

#1. read in data and create phyloseq object
read_seqtab <- function(seq,taxa,meta){
  
  ps <- phyloseq(otu_table(seq, taxa_are_rows=FALSE), 
                 sample_data(meta), 
                 tax_table(taxa))
  return(ps)
}

#2. Visualizations
viz_diversity <- function(ps, var_list, name, suff, base_path){
  
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
    print(plot_bar(ps.top20, x=var, fill="Family") + facet_wrap(~response, scales="free_x"))
  }
  dev.off()
  
  # 
  
}

viz_heatmaps <- 

#3. Decontamination
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

#4. Batch effect