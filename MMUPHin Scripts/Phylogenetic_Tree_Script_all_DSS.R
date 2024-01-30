#This script was used to generate the phylogenetic tree for all DSS datasets.
#!/usr/bin/Rscript4.2.1

library(phyloseq)
library(dada2)
library(data.table)
library(readr)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(Biostrings)
library(hiReadsProcessor)
library(ShortRead)
library(metagMisc)
library(stringr)
library(microbiome)
library(knitr)
library(DECIPHER)
library(vegan)
library(phangorn)

#Load files
OTU <- readRDS("seqtab_DSS_no_BURELLO_MMUPHin_new.rds")
taxa <- readRDS("tax_DSS_no_BURELLO_MMUPHin_new.rds")
fastafile <- readRDS("fastafile_DSS_no_BURELLO_MMUPHin_new.rds")
metadata <- read_csv(file = "metadata_DSS_no_BURELLO_MMUPHin.csv", col_names=TRUE)

OTU <- as.data.frame(OTU)
metadata<-as.data.frame(metadata)
rownames(metadata) <- metadata[,c(4)]
class(metadata)

rownames(OTU)<-gsub(".fastq.gz","",rownames(OTU)) #remove .fastq.gz extension
rownames(OTU)<-gsub("sample","SRR",rownames(OTU)) #rename the Kazakevych samples

OTUmatrix = as.matrix(OTU) #Convert OTU table to matrix
taxtable <- taxa %>% replace(is.na(.), "Unclassified") #replacing NAs with "Unclassified"
taxamatrix = as.matrix(taxtable) #Convert taxonomy table to matrix


###Run sequence alignment (MSA) using DECIPHER
sequences <- getSequences(fastafile)
names(sequences) <- sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

saveRDS(alignment, file = "alignment_DSS_no_BURELLO_MMUPHin_new.rds")


phang.align <- phyDat(as(alignment, "matrix"), type="DNA") #Change sequence alignment output into a phyDat structure
saveRDS(phang.align, file = "phang.align_DSS_no_BURELLO_MMUPHin_new.rds")

dm <-dist.ml(phang.align) #Create distance matrix
saveRDS(dm, file = "dm_DSS_no_BURELLO_MMUPHin_new.rds")

treeNJ <- NJ(dm) #Perform Neighbour joining
saveRDS(treeNJ, file = "treeNJ_DSS_no_BURELLO_MMUPHin_new.rds")

fit = pml(treeNJ, data=phang.align) #Internal maximum likelihood
saveRDS(fit, file = "fit_DSS_no_BURELLO_MMUPHin_new.rds")

fitGTR <- update(fit, k=4, inv=0.2) #negative edges length changed to 0!
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace=0))

saveRDS(fitGTR, file = "phylogenetic_tree_DSS_no_BURELLO_MMUPHin_new.rds")




