#!/usr/bin/env Rscript

library(dada2)
library(ShortRead)
library(ggplot2)


# # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Get sample files and name

# The argument will be the sample ID of the mouse
#sample = commandArgs(trailingOnly=TRUE)

# The name of the file will be the argument (sample ID)
#R1 <- paste0("/pollard/home/slyalina/work/projects/mouse_ibd_16s/demult/",sample,"-R1.fastq")
#R2 <- paste0("/pollard/home/slyalina/work/projects/mouse_ibd_16s/demult/",sample,"-R2.fastq")

#R1 <- "/Users/student/Documents/PollardRotation/InputData/42_27-R1.fastq"
#R2 <- "/Users/student/Documents/PollardRotation/InputData/42_27-R2.fastq"

sName <- (strsplit(basename(R1), "-"))[[1]][1]


# # # # # # # # # # # # # # # # # # # # # # # # 
# 2. Check sequence quality 
# - have we pre-trimmed the primers? - yes

plotQualityProfile(R1)
plotQualityProfile(R2)

# # # # # # # # # # # # # # # # # # # # # # # # 
# 3. Filtering

# Create a location to contain the filtered data
#filt_path <- "/pollard/home/mpittman/filt_16S"
filt_path <- "/Users/student/Documents/PollardRotation/InputData/filt_16S"

# Just chose tutorial defaults...maybe play with this?
filtF = paste0(filt_path,"/",sName,"-F-filt.fastq.gz")
filtR = paste0(filt_path,"/",sName,"-R-filt.fastq.gz")

out <- filterAndTrim(R1, filtF, 
                     R2, filtR, 
                     truncLen=c(200,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# # # # # # # # # # # # # # # # # # # # # # # # 
# 4. Check that error rates look reasonable

errF <- learnErrors(filtF, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtR, multithread=TRUE)
plotErrors(errR, nominalQ=TRUE)

# # # # # # # # # # # # # # # # # # # # # # # # 
# 5. Dereplication

# In this step we condense the data by collapsing together all reads that encode the same 
# sequence, which significantly reduces later computation times

derepF <- derepFastq(filtF, verbose=TRUE)
derepR <- derepFastq(filtR, verbose=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # 
# 6. Sample Inference

dadaF <- dada(derepF, err=errF, multithread=TRUE)
dadaR <- dada(derepR, err=errR, multithread=TRUE)

# inspect the class object
dadaF[[1]]
dadaR[[1]]

# # # # # # # # # # # # # # # # # # # # # # # # 
# 7. Merge paired reads

# Spurious sequence variants are further reduced by merging overlapping reads. 
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=TRUE)

# Inspect the merger data.frame from the first sample.
# Did most reads merge? If not, maybe I trimmed off too much.
head(mergers[[1]])

# # # # # # # # # # # # # # # # # # # # # # # # 
# 7. Construct sequence table
# Sequences that are much longer or shorter than expected may be the result of non-specific 
# priming, and may be worth removing 
# (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]). 
# This is analogous to “cutting a band” in-silico to get amplicons of the targeted length.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# # # # # # # # # # # # # # # # # # # # # # # # 
# 8. Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# What proportion of sequences was not filtered out as chimeric? 20% of the reads...
sum(seqtab.nochim)/sum(seqtab)


# # # # # # # # # # # # # # # # # # # # # # # # 
# 9. Track reads through the pipeline

track <- cbind(out, sum(getUniques(dadaF)), 
               sum(getUniques(mergers)), rowSums(seqtab), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
head(track)

# # # # # # # # # # # # # # # # # # # # # # # # 
# 10. Assign taxonomy

setwd('/Users/student/Documents/PollardRotation/dada2')

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v128.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# If your reads do not seem to be appropriately assigned, for example lots of your 
# bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, 
# your reads may be in the opposite orientation as the reference database. 
# Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) 
# and see if this fixes the assignments.

# # # # # # # # # # # # # # # # # # # # # # # # 
# 11. Create a table that can be used as input into machine learning program