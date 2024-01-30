##Preprocessing data script for Contijoch et al., 2019
#######################################################################################################

#!/usr/bin/Rscript3.6.0

library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)


#Specify path where FASTQ files are located
path <- "~/16S_mouse_data/CONTIJOCH/PRJNA413199"

list.files(path)


##Sort files in forward/reverse orders
fnFs <- sort(list.files(path, pattern=".fastq.gz$"))
fnFs

##Extract sample names (assuming filenames have the format: SAMPLENAME_XXX.fastq)
sample.names <- sapply(strsplit(fnFs, "_"), `[`,1)
sample.names



#Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)

pdf(paste0(fnFs,"qualityprofiles_Fs_CONTIJOCH.pdf"), width=10, height=10)
plotQualityProfile(fnFs[1:5])
dev.off()


#Create a new file path for the trimmed data
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))

#Quality filtering and trimming
out <- filterAndTrim(fnFs, filtFs, trimLeft=15, truncLen = 235, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
save(out, file="filterandtrim_CONTIJOCH.Rda") #Saving the output 


#Estimate the error model for DADA2 algorithm using reverse reads
errF <- learnErrors(filtFs, multithread = TRUE)
save(errF, file="errorsF_CONTIJOCH.Rda") #Saving the output 


#Plot error rates for all possible base transitions
pdf("errorplotsFs_CONTIJOCH.pdf", width=10, height=10)
plotErrors(errF, nominalQ=TRUE)
dev.off()

#Dereplicate FASTQ files to speed up computation
derepFs <- derepFastq(filtFs, verbose = TRUE)
save(derepFs, file="dereplicationFs_CONTIJOCH.Rda") #Saving the output 

#Name derep class objects by the sample name
names(derepFs) <- sample.names

#Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
save(dadaFs, file="inferenceFs_CONTIJOCH.Rda") #Saving the output

#Tabulate denoised and merged reads (similar to OTU table)
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab) 

#View the length of all total ribosomal sequence variants (RSVs)
table(nchar(getSequences(seqtab)))
write.table(seqtab, file = "seqtab_CONTIJOCH.txt", sep="\t")

#Chimera checking and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

#Calculate the proportion of non-chimeric RSVs
sum(seqtab.nochim)/sum(seqtab)
write.table(seqtab.nochim, file = "seqtabnoc_CONTIJOCH.txt", sep="\t")

#Read count tracking summary table
getN <- function(x) sum(getUniques(x))

summary_table <- data.frame(row.names=sample.names, Input_reads=out[,1], Filtered_Reads=out[,2], 
                            denoised_F=sapply(dadaFs, getN), 
                            nonchim=rowSums(seqtab.nochim),
                            Percent_Retained_Reads=round(rowSums(seqtab.nochim)/out[,1]*100, 1))

write.table(summary_table, "read-count-tracking_CONTIJOCH.tsv", quote = FALSE, sep = "\t", col.names = NA)


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/TaxonomyClassifiers/RDP/rdp_train_set_18.fa.gz", multithread = TRUE)
unname(head(taxa))


#Assign species
taxa.plus <- addSpecies(taxa, "~/TaxonomyClassifiers/RDP/rdp_species_assignment_18.fa.gz", verbose = TRUE)
colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
write.table(taxa.plus, file = "taxtab_CONTIJOCH.txt", sep="\t")   

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))


# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <-paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_CONTIJOCH.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts_CONTIJOCH.tsv", sep="\t", quote=F, col.names=NA)

