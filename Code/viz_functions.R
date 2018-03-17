# Source this file for miscellaneous functions for diversity_viz.R
# (I'm tired of having to scroll past these...)

# Function to remove batch effects
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


# Create a taxonomy tree using DECIPHER
get_tree <- function(seq){
  
  seqs <- names(seq)
  names(seqs) <- seqs
  
  alignment <- AlignSeqs(DNAStringSet(seqs), iterations=5, refinements=5)
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phang.align)
  
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  
  return(fitGTR)
}

# Read data into a phyloseq object
library(phyloseq)
read_seqtab <- function(seq,taxa,meta,fitGR){
  
  # Create phyloseq object
  ps <- phyloseq(otu_table(seq, taxa_are_rows=FALSE), 
                 sample_data(meta), 
                 tax_table(taxa),
                 phy_tree(fitGR$tree))
  return(ps)
}

viz_diversity <- function(ps, var_list, name, suff, base_path){
  
  # Alpha diversity
  #out_file = paste0(base_path, "Output/PhyloSeq/",name,"_",suff,"_alphaDiversity.pdf")
  #pdf(out_file)
  #for(var in var_list){
  #  print(plot_richness(ps, x=var, measures=c("Shannon", "Simpson"), color=var) + theme_bw())
  #}
  #dev.off()
  
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
  
}

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
    print(plot_tree(physeq, ladderize="left", color=var,
              title = paste0(name, " Phylogenic tree by ", var, ", ", suff )))
    
  }
  dev.off()
  
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