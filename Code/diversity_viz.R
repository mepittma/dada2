#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(decontam)
library(sva)
library(RAM)

#base_path = "/pollard/home/mpittman/dada2/"
base_path = "/Users/student/Documents/PollardRotation/dada2/"
options(tz="America/Los_Angeles")

# Create a new column in metadata capturing all the co-variates
#metadata$batch <- ""
#for (var in var_list){
#  metadata$batch <- paste(metadata$batch, var, sep = "_")
#}

cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

# Outline of what I want to do
#c 1. Read in the dada2 output
#c 2. Decontaminate
#c 3. Remove batch effects
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

# read in the dada2 output
name = "Baxter_AOMDSS"
# Read in necessary data
seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)
var_list = c("collection_date","inoculum","response")

seq <- seq[which(row.names(seq) %in% row.names(meta)),]
seq <- seq[,colSums(seq) > 0]

# Create phyloseq object, decontaminate it
ps = read_seqtab(seq,taxa,meta)
decon = deco_ps(ps,taxa, paste0(base_path, "Output/Decontam/",name,"_LibrarySize.pdf"),
                paste0(base_path, "Output/Decontam/", name, "_removedTaxa.tsv"))
debatch = 

# # # # # # # # # # # # # # # # # # # # # # # # 
# General functions

# Create heatmap from matrix

# Create dendrogram/heatmap from matrix

# Create taxa-tree from matrix

# Create taxa heatmap from matrix

# PCoA plot
PCoA <- function(seq, meta, taxa, var_list, suff, base_path){
  # This function takes as input sequence, metadata, and taxonomy information
  
  stringsAsFactors = FALSE
  
  # Get meta matrix into proper form
  meta_fact <- meta
  meta_fact$response <- gsub('0', 'no', meta_fact$response)
  meta_fact$response <- gsub('1', 'yes', meta_fact$response)
  meta_fact$response <- as.factor(meta_fact$response)
  
  # Add taxonomy column to the seq matrix
  taxa<- as.data.frame(taxa)
  taxa$taxonomy <- paste0("k__",taxa$Kingdom,"; ", 
                          "p__",taxa$Phylum, "; ",
                          "c__", taxa$Class, "; ",
                          "o__", taxa$Order, "; ",
                          "f__", taxa$Family, "; ",
                          "g__", taxa$Genus)
  seq_ <- t(seq)
  tax_seq <- as.data.frame(cbind.data.frame(seq_, taxonomy = taxa$taxonomy[match(rownames(seq_), rownames(taxa))]))
  
  # Create and save the pcoa plot
  out_file = paste0(base_path, "Output/PhyloSeq/",name,"_",suff,"_PCoA.pdf")
  pdf(out_file)
  pcoa.plot(tax_seq, is.OTU=TRUE, meta_fact, factors=c(Diet="inoculum", IBD="response"), 
            rank=NULL, sample.labels = FALSE, main = paste0(name, " PCoA plot: ", suff))
  dev.off()
}

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

viz_heatmaps <- function(ps)

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
debatch <- function(abun_mat, meta){
  
  # log-transform the transposed data
  my_data <-log(as.matrix(t(abun_mat)) +1)
  my_data <- my_data[rowSums(my_data) > 0 ,]
  
  # Create a null model and one accounting for my variable of interest
  mod1 = model.matrix(~as.factor(meta$response))
  mod0 = cbind(mod1[,1])
  
  # Calculate batch coefficients
  my_n_sv = num.sv(my_data,mod1,method="leek")
  my_svseq = svaseq(my_data,mod1,mod0,n.sv=my_n_sv-1)
  my_sv<-my_svseq$sv
  
  # Remove batch coefficients
  clean_df <- cleanY(my_data,mod1,my_sv)
  
  # Return abundance matrix in original format
  return(t(clean_df))
}

