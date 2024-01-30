##Data file clean up script

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

rm(list = ls())

#Load files
OTU <- readRDS("seqtab_DSS_no_BURELLO_MMUPHin.rds")
taxa <- readRDS("tax_DSS_MMUPHin.rds")
fastafile <- readRDS("fastafile_MMUPHin.rds")
metadata <- read_csv(file = "metadata_fully_clustering_MMUPHin.csv", col_names=TRUE)

OTU <- as.data.frame(OTU)
metadata<-as.data.frame(metadata)
rownames(metadata) <- metadata[,c(4)]
class(metadata)

rownames(OTU)<-gsub(".fastq.gz","",rownames(OTU)) #remove .fastq.gz extension
rownames(OTU)<-gsub("sample","SRR",rownames(OTU)) #rename the samples so that they're matching across files

OTUmatrix = as.matrix(OTU) #Convert OTU table to matrix
taxtable <- taxa %>% replace(is.na(.), "Unclassified") #replacing NAs with "Unclassified"
taxamatrix = as.matrix(taxtable) #Convert taxonomy table to matrix
metamatrix <- as.matrix(metadata)

#Check for discrepancies between OTU table and medtadata files
rownames(OTUmatrix) -> OTU_check
rownames(metamatrix) -> meta_check

#Convert to vectors
OTU_check <- unlist(OTU_check)
meta_check <- unlist(meta_check)

int <- intersect(OTU_check,meta_check) #What's in both files
difs <- setdiff(OTU_check,meta_check) #What's different in OTU_check
difss <- setdiff(meta_check,OTU_check) #What's different in meta_check


#Remove the rows with the samples in difss, if any
OTU2 <- OTU[!(row.names(OTU) %in% difs),]

OTUmatrix2 = as.matrix(OTU2) #Convert OTU table to matrix


#Check for discrepancies between OTU table and taxonomy table 
colnames(OTUmatrix) -> OTU_check2
rownames(taxamatrix) -> taxa_check

#Convert to vectors
OTU_check2 <- unlist(OTU_check2)
taxa_check <- unlist(taxa_check)


int <- intersect(OTU_check2,taxa_check) #What's in both files
difs <- setdiff(OTU_check2,taxa_check) #What's different in OTU_check
difss <- setdiff(taxa_check,OTU_check2) #What's different in meta_check

saveRDS(OTU2, file = "seqtab_fully_clustering_DSS_MMUPHin.rds")
saveRDS(taxa, "tax_fully_clustering_DSS_MMUPHin.rds")


