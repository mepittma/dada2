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
  write.table(OTU1, file = paste0(out_path,name,"_seq.tsv"))
  write.table(OTU2, file = paste0(out_path,name,"_tax.tsv"))
  
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

phyloML(base_path, out_path, "Baxter_AOMDSS", c("collection_date","inoculum","tumors","response"))
phyloML(base_path, out_path, "Helm_DSS", c("collection_date","inoculum","response"))
phyloML(base_path, out_path, "TMM_DSS", c("collection_date","response"))
phyloML(base_path, out_path, "UCSD_TNBS", c("collection_date","inoculum","response"))
#phyloML(base_path, out_path, "UCSF_DNR", c("response","genotype","cohort")) # why isn't this one working?
phyloML(base_path, out_path, "UMAA_DSS", c("collection_date","response"))
phyloML(base_path, out_path, "UTA_TNBS", c("collection_date","response"))
#phyloML(base_path, out_path, "UTS_DSS", c("inoculum","response")) # why does this one just have EukaryotaNANANA?




