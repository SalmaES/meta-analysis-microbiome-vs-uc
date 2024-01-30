library(data.table)
library(readr)
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
library(dada2)
library(phangorn)
library(speedyseq)
library(microbial)



OTU <- readRDS("OTU_MMUPHin.rds")
taxa <- readRDS("taxtab_MMUPHin.rds")
metadata <- readRDS("meta_MMUPHin.rds")
fastafile <- readRDS("fastafile_MMUPHin.rds")

OTU <- as.data.frame(OTU)
metadata<-as.data.frame(metadata)
rownames(metadata) <- metadata[,c(4)]
class(metadata)

#Check difference between OTU table and medtadata
rownames(OTU) ->OTU_check
rownames(metadata) -> meta_check

#Convert to vectors
OTU_check <- unlist(OTU_check)
meta_check <- unlist(meta_check)

#What's in both files
int <- intersect(OTU_check,meta_check)

#What's different in OTU_check
difs <- setdiff(OTU_check,meta_check)

#What's different in meta_check
difss <- setdiff(meta_check,OTU_check)

rownames(OTU)<-gsub(".fastq.gz","",rownames(OTU)) #remove .fastq.gz extension

#Check difference between OTU table and medtadata
rownames(OTU) ->OTU_check
rownames(metadata) -> meta_check

#Convert to vectors
OTU_check <- unlist(OTU_check)
meta_check <- unlist(meta_check)

#What's in both files
int <- intersect(OTU_check,meta_check)

#What's different in OTU_check
difs <- setdiff(OTU_check,meta_check)

#What's different in meta_check
difss <- setdiff(meta_check,OTU_check)

#Remove the rows with the samples in difss
metadata2 <- metadata[!(row.names(metadata) %in% difss),]
OTU2 <- OTU[!(row.names(OTU) %in% difs),]

#Recheck: Check difference between OTU table and medtadata
rownames(OTU2) ->OTU_check2
rownames(metadata2) -> meta_check2

#Convert to vectors
OTU_check2 <- unlist(OTU_check2)
meta_check2 <- unlist(meta_check2)

#What's in both files
int2 <- intersect(OTU_check2,meta_check2)

#What's different in OTU_check
difs2 <- setdiff(OTU_check2,meta_check2)

#What's different in meta_check
difss2 <- setdiff(meta_check2,OTU_check2)

names(taxa) [1] <- 'Domain'

OTUmatrix2 = as.matrix(OTU2) #Convert OTU table to matrix
taxamatrix = as.matrix(taxa) #Convert taxtable table to matrix


#Create a phyloseq object
ps2 <- phyloseq(otu_table(OTUmatrix2, taxa_are_rows=FALSE), tax_table(taxamatrix), sample_data(metadata2))
ps2
saveRDS(ps2, "phyloseq_object_all_MMUPHin(raw).rds")

ordination.raw <- ordinate(ps2, "PCoA", distance = "bray")
pdf(file="PCoA_raw_all_MMUPHin(Condition).pdf")
plot_ordination(ps2, ordination.raw, color = "Condition", 
                title = "Before Batch Adjustment")  + theme(aspect.ratio=1)
dev.off()


ordination.raw <- ordinate(ps2, "PCoA", distance = "bray")
pdf(file="PCoA_raw_all_MMUPHin(Colitis_Model).pdf")
plot_ordination(ps2, ordination.raw, color = "Colitis_Model", 
                title = "Before Batch Adjustment")  + theme(aspect.ratio=1)
dev.off()


ps2 <- ps2 %>%
  rename_tax_table(Kingdom=Domain) %>%
  rename_with_tax_table(stringr::str_to_title)

raw.permanova <- betatest(ps2, group="Condition", distance = "bray")   
write.csv(raw.permanova, "PERMANOVA.raw.bray.all.MMUPHin(Condition).csv")

raw.permanova <- betatest(ps2, group="Colitis_Model", distance = "bray")   
write.csv(raw.permanova, "PERMANOVA.raw.bray.all.MMUPHin(Colitis_Model).csv")
