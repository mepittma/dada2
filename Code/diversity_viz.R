#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(decontam)
library(sva)
library(DECIPHER)
library(phangorn)

#base_path = "/pollard/home/mpittman/dada2/"
base_path = "/Users/student/Documents/PollardRotation/dada2/"
source(paste0(base_path,"Code/viz_functions.R"))
options(tz="America/Los_Angeles")

# # # # # # # # # # # FILTERING FUNCTIONS # # # # # # # # # # # # # 

# Decontamination
deco_ps <- function(ps, meta, taxa, out_plot, out_tab){
  
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
  print(dim(taxa))
  contam_tax <- taxa[which(row.names(taxa) %in% row.names(to_remove)),]
  write.table(contam_tax, file = out_tab, sep = "\t")
  
  # Create a new phyloseq object
  new_seq = as.data.frame(otu_table(ps.noncontam))
  new_GF = get_tree(new_seq)
  
  return(read_seqtab(new_seq,taxa,meta,new_GF))

}

# Batch effect
debatch <- function(ps, meta, taxa){
  
  # Extract variables from phyloseq object
  abun_mat <- as.data.frame(otu_table(ps))
  
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
  
  # Create a new phyloseq object
  new_seq = as.data.frame(t(clean_df))
  new_GF = get_tree(new_seq)
  
  return(read_seqtab(seq=new_seq,taxa=taxa,meta=meta,fitGR=new_GF))
}

# # # # # # # # # # # MAIN FUNCTION # # # # # # # # # # # # # 

phylo_viz <- function(name, var_list, base_path){
  
  # Read in data
  seq <- readRDS(paste0(base_path, "Output/SeqTables/", name, "_seqtab_nochim.rds"))
  taxa <- readRDS(paste0(base_path, "Output/Taxa/", name,"_taxa_silva_plus.rds"))
  meta <- read.table(paste0(base_path, "MetaData/", name, "_processed.txt"), 
                     sep = "\t", header = TRUE)
  
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
  decon = deco_ps(ps, meta, taxa, 
                  out_plot = paste0(base_path, "Output/Decontam/",name,"_LibrarySize.pdf"), 
                  out_tab = paste0(base_path, "Output/Decontam/", name, "_removedTaxa.tsv"))
  debat = debatch(ps=decon, meta=meta, taxa=taxa)
  
  # Save RData
  saveRDS(ps, file = paste0(base_path, "Output/RData/", name,"_ps.rds"))
  saveRDS(decon, file = paste0(base_path, "Output/RData/", name,"_decon.rds"))
  saveRDS(debat, file = paste0(base_path, "Output/RData/", name,"_debat.rds"))
  
  # Visualizations for the three datasets
  
  suff_list = c("Uncorrected", "Decontaminated", "Batch-Corrected")
  seq_list <- c(ps, decon, debat)
  
  for (i in 1:3){
    suff = unlist(suff_list[i])
    ps_ = seq_list[[i]]
    
    PCoA(ps_, var_list, name, suff, base_path)
    heat_viz(ps_, name, suff, base_path)
    tree_viz(ps_, name, var_list, suff, base_path)
    
    if(suff != "Batch-Corrected"){
      viz_diversity(ps_, var_list, name, suff, base_path)
    }
    
    # Return the de-batched phyloseq object so we can use it downstream
    return(list(ps, debat))
    
  }
  
  
}

# # # # # # # # # # # VISUALIZATION CALLS # # # # # # # # # # # # # 

Bax_list <- phylo_viz(name = "Baxter_AOMDSS", 
        var_list = c("collection_date","inoculum","response_factor"), base_path)
Helm_list <- phylo_viz("Helm_DSS", c("collection_date","inoculum","response_factor"), base_path)
TMM_list <- phylo_viz("TMM_DSS", c("collection_date","response_factor"), base_path)
UCTN_list <- phylo_viz("UCSD_TNBS", c("collection_date","inoculum","response_factor"), base_path)
UCIL_list <- phylo_viz("UCSD_IL10", c("collection_date","cage","gender","run",
                                      "timepoint","response_factor"), base_path)
UMAA_list <- phylo_viz("UMAA_DSS", c("collection_date","response_factor"), base_path)
UTA_list <- phylo_viz("UTA_TNBS", c("collection_date","response_factor"), base_path)
UTS_list <- phylo_viz("UTS_DSS", c("inoculum","response_factor"), base_path)
UCSF_list <- phylo_viz("UCSF_DNR", c("response_factor","genotype"), base_path)

