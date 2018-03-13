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
source(paste0(base_path,"Code/viz_functions.R"))
options(tz="America/Los_Angeles")

# Outline of what I want to do
#c 1. Read in the dada2 output
#c 2. Decontaminate
#c 3. Remove batch effects
#c 4. Visualize the above three cases with individual datasets + joint, debatched datasets
  #c a. PCoA - Bray-Curtis and Unifrac
  #c b. Major taxon abundance 
# 5. Visualize with just decontam/unbatched data
  # a. Tree of life with the combined datasets
  #c b. Heatmap
    # i. with each dataset individually, comparing across variables
    # ii. with all datasets, comparing across IBD
# 6. Feature matrices both with OTUs and taxa

# read in the dada2 output
name = "Baxter_AOMDSS"

# Read in necessary data
seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), sep = "\t", header = TRUE)
var_list = c("collection_date","inoculum","response")

# Get rid of empty taxonomy, add response factor column
seq <- seq[which(row.names(seq) %in% row.names(meta)),]
seq <- as.data.frame(seq[,colSums(seq) > 0])
meta$response_factor <- meta$response
meta$response_factor <- gsub(1, "yes", meta$response_factor)
meta$response_factor <- gsub(0, "no", meta$response_factor)

# Create phylogenetic tree
fitGR <- get_tree(seq)

# Create phyloseq object, decontaminate it, debatch it
ps = read_seqtab(seq,taxa,meta,fitGR)
##decon = deco_ps(ps,taxa, paste0(base_path, "Output/Decontam/",name,"_LibrarySize.pdf"),
##                paste0(base_path, "Output/Decontam/", name, "_removedTaxa.tsv"))
##debat = debatch(seq, meta)

# Desired visualizations
PCoA(ps, var_list = c("collection_date","inoculum","response_factor"), name, 
     suff="Uncorrected", base_path)
heat_viz(ps, name, "Uncorrected", base_path)
tree_viz(ps, name, var_list = c("collection_date","inoculum","response_factor"),
         suff="Uncorrected", base_path)

# # # # # # # # # # # # # # # # # # # # # # # # 
# General functions

# Create heatmaps from matrix
heat_viz <- function(ps, name, suff, base_path){
  
  out_path = paste0(base_path, "Output/PhyloSeq/Heatmap/",name,"_",suff,"_Heatmaps.pdf")
  pdf(out_path)
    
  theme_set(theme_bw())
  gpt <- subset_taxa(ps, Kingdom=="Bacteria")
  gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:50]), gpt)
  
  # Plot the heatmap
  print(plot_heatmap(gpt, "NMDS", "bray", "response_factor", "Family", 
                     low="#000033", high="#FF3300"))
  
  # Create dendrogram/heatmap from matrix
  top_taxa <- as.data.frame(tax_table(gpt))
  taxa_names(gpt) <- make.names(top_taxa$Family, unique=TRUE)
  heatmap(otu_table(gpt))
  
  dev.off()
  
}



# Create taxa-tree from matrix
tree_viz <- function(ps, name, var_list, suff, base_path){
  
  # Save out to file
  out_path = paste0(base_path, "Output/PhyloSeq/AbundanceTree/",name,"_",suff,"_Trees.pdf")
  pdf(out_path)
  
  # Select top50 taxa
  physeq = prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:50]), ps)
  
  # map color to taxonomic class
  plot_tree(physeq, ladderize="left", color="Class", 
            title = paste0(name, " Phylogenic tree (Class), ", suff ))
  plot_tree(physeq, ladderize="left", color="Phylum",
            title = paste0(name, " Phylogenic tree (Phylum), ", suff ))
  
  for (var in var_list){
    
    # map color to environmental factors
    plot_tree(physeq, ladderize="left", color=var,
              title = paste0(name, " Phylogenic tree by ", var, ", ", suff ))
    
  }
  
}


# PCoA plot
PCoA <- function(ps, var_list, name, suff, base_path){
  
  out_path = paste0(base_path, "Output/PhyloSeq/PCoA/",name,"_",suff,"_PCoA.pdf")
  pdf(out_path)
  
  for (var in var_list){
    
    if (var != "response_factor"){
      ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
      p = plot_ordination(ps, ordu, color= var, 
                            shape="response_factor")
      p = p + geom_point(size=7, alpha=0.75)
      p = p + scale_colour_brewer(type="qual", palette="Set1")
      print(plot(p + ggtitle(paste0("MDS/PCoA on weighted-UniFrac distance, ", name," ", suff))))
      
    }
    
  }
  
  dev.off()
}



# # # # # # # # # # # # # # # # # # # # # # # # 
# Dataset-specific functions

#1. read in data and create phyloseq object


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

