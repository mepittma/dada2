#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(decontam)

#base_path = "/pollard/home/mpittman/dada2/"
base_path = "/Users/student/Documents/PollardRotation/dada2/"

name = "Baxter_AOMDSS"
# # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Load in the seqtab_nochim files and attach metadata

read_seqtab <- function(name, to_prune=to_prune){
  
  seqtab.nochim <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
  taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
  samdf <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)
  
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                 sample_data(samdf), 
                 tax_table(taxa))
  if (to_prune){
    ps <- prune_samples(sample_names(ps) != to_prune, ps) # Remove any desired samples
  }
}


#2. Visualizations
visualize <- function(ps, out_file){
  pdf(out_file)
  
  # Alpha diversity
  print(plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw())
  
  # Ordination
  ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
  print(plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS"))
  
  # Taxa Bar plot
  top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
  
  dev.off()
}

# check sample variables
sample_variables(ps)

#3. Decontamination
inspect_library_sizes(ps){
  df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
}