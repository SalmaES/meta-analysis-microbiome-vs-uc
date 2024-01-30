#This script was used to run MMUPhin and the statistical analyses run on the batch adjusted data, with modifications depending on dataset grouping used.

# Cleaning Memory
# --------------------------------------------------
rm(list = ls())

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
library(speedyseq)
library(microbiomeutilities)
library(BiocManager)
library(devtools)
library(DESeq2)
library(microbial)
library(vegan)
library(microbiomeMarker)
library(yingtools2)
library(speedyseq)

rm(list=ls())

setwd("~/Desktop/OneDrive - University of Essex/Year 3/Microbiome Analysis/MMUPHin/Fully Clustering Studies")

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)

#Load files
OTU <- readRDS("seqtab_fully_clustering_DSS_MMUPHin_new.rds")
taxa <- readRDS("tax_fully_clustering_DSS_MMUPHin_new.rds")
metadata <- read_csv("metadata_fully_clustering_MMUPHin(adjusted).csv")
tree <- readRDS("phylogenetic_tree_fully_clustering_DSS_MMUPHin_new.rds")

#Assigning metadata rownames
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[,4]
metadata <- as.matrix(metadata)



###################################### Data file clean up #########################################

rownames(OTU)<-gsub(".fastq.gz","",rownames(OTU)) #remove .fastq.gz extension

#Check difference between OTU table and metadata
rownames(OTU) -> OTU_check
rownames(metadata) -> meta_check

#Convert to vectors
OTU_check <- unlist(OTU_check)
meta_check <- unlist(meta_check)

int <- intersect(OTU_check,meta_check) #What's in both files
difs <- setdiff(OTU_check,meta_check) #What's different in OTU_check
difss <- setdiff(meta_check,OTU_check) #What's different in meta_check

#Check difference between OTU table and taxa table
colnames(OTU) -> OTU_check2
rownames(taxa) -> taxa_check

#Convert to vectors
OTU_check2 <- unlist(OTU_check2)
taxa_check <- unlist(taxa_check)

int <- intersect(OTU_check2,taxa_check) #What's in both files
difs <- setdiff(OTU_check2,taxa_check) #What's different in OTU_check
difss <- setdiff(taxa_check,OTU_check2) #What's different in taxa_check

taxtable <- taxa %>% replace(is.na(.), "Unclassified") #replacing NAs with "Unclassified"
taxamatrix = as.matrix(taxtable) #Convert taxonomy table to matrix
class(taxamatrix) #Check object type (must be a matrix for Phyloseq)
class(OTU) #Check object type (must be a matrix for Phyloseq)
class(metadata) #Check object type (must be a dataframe for Phyloseq) 
metadata <- as.data.frame(metadata)



###################################### Data processing #########################################

#Create a phyloseq object
ps <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE), tax_table(taxamatrix), sample_data(metadata), phy_tree(tree$tree))
ps # (3367 taxa and 101 samples)

#Filter out samples with read counts below 7000 
ps <- prune_samples(microbiome::readcount(ps)>=7000, ps) 
ps # (3367 taxa and 101 samples)


#Adding reference sequences to phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Change taxtable column name from "Domain" to "Kingdom"
ps1 <- ps %>%
  rename_tax_table(Domain=Kingdom) %>%
  rename_with_tax_table(stringr::str_to_title)

ps1 %>% rank_names() #Check the taxtable column names

##Adding the lowest taxonomy for an ASV
ps2 <- add_besthit(ps1, sep=":")
taxa_names(ps2)[1:2]


##Aggregating ASVs at different taxonomic levels
pseq3 <- tax_glom(ps2, taxrank="Species")
###Saving phytree for dowstream analysis (post MMUPHin)
phy <- phy_tree(pseq3)

##Converting phyloseq object to dataframe
pseq3df <- phyloseq_to_df(pseq3, addtax = T, addtot = F, addmaxrank = F,
                          sorting = NULL)





###################################### Data Clean Up for MMUPHin #########################################

OTU_test <- pseq3df %>%
  unite("Taxa", Domain:Phylum:Class:Order:Family:Genus:Species) #combining the taxonomic levels into 1 taxa column
OTU_test %>% 
  mutate(across(where(is.character), str_trim)) -> OTU_table1
OTU_table1 <- OTU_table1[,-c(1)]
OTUtestmatrix <- as.matrix(OTU_table1)
rownames(OTUtestmatrix) <- OTUtestmatrix[,1]
otumatrix <- OTUtestmatrix[,-1]
class(otumatrix) = "numeric"

write.table(otumatrix, file = "otumatrix_fully_clustering_DSS_MMUPHin(adjusted).txt", sep="\t")
write.table(metadata, file = "metadata_fully_clustering_DSS_MMUPHin(adjusted).txt", sep="\t")





###################################### Running MMUPHin #########################################


## Batch effect adjustment with MMUPHin
fit_adjust_batch <- adjust_batch(feature_abd = otumatrix,
                                 batch = "Study",
                                 covariates = "Condition",
                                 data = metadata,
                                 control = list(verbose = FALSE)) $feature_abd_adj


tax <- pseq3df[, -c(1, 9:109)] #creating a new taxonomic table after aggregating ASV for new phyloseq object

rownames(otumatrix) -> rownames
rownames(fit_adjust_batch) <- rownames
rownames(tax) <- rownames
taxmatrix <- as.matrix(tax)

##Creating new phyloseq object with adjusted abundances for data manipulation
ps_mfn <- phyloseq(otu_table(fit_adjust_batch, taxa_are_rows=TRUE), tax_table(taxmatrix), sample_data(metadata))
###Changing taxa_names in phylogenetic tree to extract and merge into new phyloseq object
taxanames1 <- taxa_names(ps_mfn)
taxa_names(pseq3) <- taxanames1
phy <- phy_tree(pseq3)

ps_mfn1 <- phyloseq(otu_table(fit_adjust_batch, taxa_are_rows=TRUE), tax_table(taxmatrix), sample_data(metadata), phy_tree(phy))
ps_mfn1 





############################# Plotting PCoA using Raw and Batch-Adjusted Data (Pre-Filtering) ################################


pdf(file="PCoA_Fully_Clustering_MMUPHin_bray.pdf")
ordination <- ordinate(ps_mfn1, "PCoA", distance = "bray")
plot_ordination(ps_mfn1, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets - MMUPHin (bray)")

pdf(file="PCoA_Fully_Clustering_MMUPHin_unifrac.pdf")
ordination <- ordinate(ps_mfn1, "PCoA", distance = "unifrac")
plot_ordination(ps_mfn1, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets - MMUPHin (unifrac)")

pdf(file="PCoA_Fully_Clustering_MMUPHin_wunifrac.pdf")
ordination <- ordinate(ps_mfn1, "PCoA", distance = "wunifrac")
plot_ordination(ps_mfn1, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets - MMUPHin (wUnifrac)")
dev.off()


pdf(file="PCoA_Fully_Clustering_Raw_bray.pdf")
ordination <- ordinate(ps2, "PCoA", distance = "bray")
plot_ordination(ps2, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets - Raw (bray)")
dev.off()

pdf(file="PCoA_Fully_Clustering_Raw_unifrac.pdf")
ordination <- ordinate(ps2, "PCoA", distance = "unifrac")
plot_ordination(ps2, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets - Raw (Unifrac)")
dev.off()

pdf(file="PCoA_Fully_Clustering_Raw_wunifrac.pdf")
ordination <- ordinate(ps2, "PCoA", distance = "wunifrac")
plot_ordination(ps2, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets - Raw (wUnifrac)")
dev.off()







###################################### PERMANOVA Test (Pre-Filtering) #########################################


D_before <- vegdist(t(otumatrix))
D_after <- vegdist(t(fit_adjust_batch))

set.seed(1)
fit_adonis_before <- adonis2(D_before ~ Study, data = metadata) 
fit_adonis_after <- adonis2(D_after ~ Study, data = metadata) 
print(fit_adonis_before)
print(fit_adonis_after)





###################################### Wilcoxon Rank Sum Test (Pre-Filtering) ######################################### 

rich = estimate_richness(ps_mfn1)
wt1=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Condition)
wt2=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Treatment)
wt3=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Genotype)
wt4=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Intensity_of_Disease)
wt5=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Sex)
wt6=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Extraction_Kit)
wt7=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1)$Sequencing_Kit)


write.csv(wt1[["p.value"]], "Wilcoxon.test.MMUPHin.FC(Condition).csv")
write.csv(wt2[["p.value"]], "Wilcoxon.test.MMUPHin.FC(Treatment).csv")
write.csv(wt3[["p.value"]], "Wilcoxon.test.MMUPHin.FC(Genotype).csv")
write.csv(wt4[["p.value"]], "Wilcoxon.test.MMUPHin.FC(Intensity_of_Disease).csv")
write.csv(wt5[["p.value"]], "Wilcoxon.test.MMUPHin.FC(Sex).csv")
write.csv(wt6[["p.value"]], "Wilcoxon.test.MMUPHin.FC(ExtractionKit).csv")
write.csv(wt7[["p.value"]], "Wilcoxon.test.MMUPHin.FC(SequencingKit).csv")




###################################### Prevalence Filtering ######################################### 

prevdf = apply(X = otu_table(ps_mfn1),
               MARGIN = ifelse(taxa_are_rows(ps_mfn1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_mfn1),
                    tax_table(ps_mfn1))
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps_mfn1, "Phylum"))

pdf(file="Prevalence_Filtering_Plots_MMUPHin.FC.pdf")
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_mfn1),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
dev.off()

prevalenceThreshold = 0.05 * nsamples(ps_mfn1)
nsamples(ps_mfn1)
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps_mfn1.pruned = prune_taxa(keepTaxa, ps_mfn1)
ps_mfn1.pruned 




###################################### Wilcoxon Rank Sum Test (Post-Filtering) ######################################### 

rich = estimate_richness(ps_mfn1.pruned)
wtp1=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Condition)
wtp2=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Treatment)
wtp3=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Genotype)
wtp4=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Intensity_of_Disease)
wtp5=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Sex)
wtp6=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Extraction_Kit)
wtp7=pairwise.wilcox.test(rich$Observed, sample_data(ps_mfn1.pruned)$Sequencing_Kit)


write.csv(wtp1[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(Condition).csv")
write.csv(wtp2[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(Treatment).csv")
write.csv(wtp3[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(Genotype).csv")
write.csv(wtp4[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(Intensity_of_Disease).csv")
write.csv(wtp5[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(Sex).csv")
write.csv(wtp6[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(ExtractionKit).csv")
write.csv(wtp7[["p.value"]], "Wilcoxon.test.MMUPHin.FC.pruned(SequencingKit).csv")




############################# Plotting PCoA Plots using Raw and Batch-Adjusted Data (Post-Filtering) ################################

pdf(file="PCoA_Fully_Clustering_MMUPHin_pruned_bray.pdf")
ordination <- ordinate(ps_mfn1.pruned, "PCoA", distance = "bray")
plot_ordination(ps_mfn1.pruned, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets(P) - MMUPHin (bray)")

pdf(file="PCoA_Fully_Clustering_MMUPHin_pruned_unifrac.pdf")
ordination <- ordinate(ps_mfn1.pruned, "PCoA", distance = "unifrac")
plot_ordination(ps_mfn1.pruned, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets(P) - MMUPHin (unifrac)")

pdf(file="PCoA_Fully_Clustering_MMUPHin_pruned_wunifrac.pdf")
ordination <- ordinate(ps_mfn1.pruned, "PCoA", distance = "wunifrac")
plot_ordination(ps_mfn1.pruned, ordination, color = "Condition", shape = "Study", 
                title = "Fully Clustering Datasets(P) - MMUPHin (wUnifrac)")
dev.off()

#Change taxtable column name from "Domain" to "Kingdom" for downstream analysis
ps_mfn1.pruned1 <- ps_mfn1.pruned %>%
  rename_tax_table(Kingdom=Domain) %>%
  rename_with_tax_table(stringr::str_to_title)



