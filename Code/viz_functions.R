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
library(DECIPHER)
library(phangorn)
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
                 phy_tree(fitGTR$tree))
  return(ps)
}