#This script was used to plot the heatmap for DESeq2 results from DSS datasets group B (figure 3.11).

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
library(ggvenn)
library(RColorBrewer)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)

rm(list=ls())

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)

#Load files
OTU <- readRDS("seqtab_semi_clustering_DSS_MMUPHin_new.rds")
taxa <- readRDS("tax_semi_clustering_DSS_MMUPHin_new.rds")
metadata <- read_csv("metadata_semi_clustering_MMUPHin(adjusted).csv")
tree <- readRDS("phylogenetic_tree_semi_clustering_DSS_MMUPHin_new.rds")

#Assigning metadata rownames
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[,4]
metadata <- as.matrix(metadata)




###################################### Data file clean up #########################################

rownames(OTU)<-gsub(".fastq.gz","",rownames(OTU)) #remove .fastq.gz extension
rownames(OTU)<-gsub("sample","SRR",rownames(OTU)) #rename the Kazakevych samples

#Check difference between OTU table and metadata
rownames(OTU) -> OTU_check
rownames(metadata) -> meta_check

#Convert to vectors
OTU_check <- unlist(OTU_check)
meta_check <- unlist(meta_check)

int <- intersect(OTU_check,meta_check) #What's in both files
difs <- setdiff(OTU_check,meta_check) #What's different in OTU_check
difss <- setdiff(meta_check,OTU_check) #What's different in meta_check

#Remove the rows with the samples in difs
OTU <- OTU[!(row.names(OTU) %in% difs),]
metadata <- metadata[!(row.names(metadata) %in% difss),]

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
ps <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE), 
               tax_table(taxamatrix), sample_data(metadata), phy_tree(tree$tree))

#Filter out samples with read counts below 7000 
ps <- prune_samples(microbiome::readcount(ps)>=7000, ps) 


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


##Adding the lowest taxonomy for an ASV
ps2 <- add_besthit(ps1, sep=":")
taxa_names(ps2)[1:2]


##Aggregating ASVs at different taxonomic levels
pseq3 <- tax_glom(ps, taxrank="Species")
###Saving phytree for dowstream analysis (post MMUPHin)
phy <- phy_tree(pseq3)

metadata2 <-data.frame(sample_data(pseq3))

##Converting phyloseq object to dataframe
pseq3df <- phyloseq_to_df(pseq3, addtax = T, addtot = F, addmaxrank = F,
                          sorting = NULL)




###################################### Data Clean Up for MMUPHin #########################################

OTU_test <- pseq3df %>%
  unite("Taxa", Kingdom:Phylum:Class:Order:Family:Genus:Species) #combining the taxonomic levels into 1 taxa column
OTU_test %>% 
  mutate(across(where(is.character), str_trim)) -> OTU_table1
OTU_table1 <- OTU_table1[,-c(1)]
OTUtestmatrix <- as.matrix(OTU_table1)
rownames(OTUtestmatrix) <- OTUtestmatrix[,1]
otumatrix <- OTUtestmatrix[,-1]
class(otumatrix) = "numeric"

write.table(otumatrix, file = "otumatrix_semi_clustering_DSS_MMUPHin(adjusted).txt", sep="\t")
write.table(metadata2, file = "metadata_semi_clustering_DSS_MMUPHin(adjusted).txt", sep="\t")






###################################### Running MMUPHin #########################################

## Batch effect adjustment with MMUPHin
fit_adjust_batch <- adjust_batch(feature_abd = otumatrix,
                                 batch = "Study",
                                 covariates = "Condition",
                                 data = metadata2,
                                 control = list(verbose = FALSE)) $feature_abd_adj


tax <- pseq3df[, -c(1, 9:1019)] #creating a new taxonomic table after aggregating ASV for new phyloseq object

rownames(otumatrix) -> rownames
rownames(fit_adjust_batch) <- rownames
rownames(tax) <- rownames
taxmatrix <- as.matrix(tax)

##Creating new phyloseq object with adjusted abundances for data manipulation
ps_mfn <- phyloseq(otu_table(fit_adjust_batch, taxa_are_rows=TRUE), 
                   tax_table(taxmatrix), sample_data(metadata2))
###Changing taxa_names in phylogenetic tree to extract and merge into new phyloseq object
taxanames1 <- taxa_names(ps_mfn)
taxa_names(pseq3) <- taxanames1
phy <- phy_tree(pseq3)


ps_mfn1 <- phyloseq(otu_table(fit_adjust_batch, taxa_are_rows=TRUE), 
                    tax_table(taxmatrix), sample_data(metadata2), phy_tree(phy))







###################################### DESeq2 Analysis #########################################

deseq2 <- microbiomeMarker::run_deseq2(ps_mfn,
                                       group = "Condition",
                                       taxa_rank = "none",
                                       sfType = "poscounts",
                                       pvalue_cutoff=0.01, 
                                       transform = "identity", 
                                       p_adjust = "fdr")

deseq2_table <- deseq2@marker_table
write.csv(deseq2_table, "DESeq2_table_semi_clustering_DSS(Condition)_MMUPHin.csv")




###################################### Plotting DESeq2 Heatmap #########################################

#################Functions needed to run this code ##########################
transform_abundances <- function(object,
                                 transform = c("identity", "log10", "log10p")) {
  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  otu <- as(otu_table(object), "matrix")
  
  if (transform == "identity") {
    abd <- otu
  } else if (transform == "log10") {
    abd <- transform_log10(otu)
  } else {
    abd <- transform_log10p(otu)
  }
  
  otu_table(object) <- otu_table(abd, taxa_are_rows = taxa_are_rows(object))
  
  object
}

# the data is transformed using log10(1 + x) if the data contains zeroes
transform_log10 <- function(x) {
  if (min(x) == 0) {
    warning("OTU table contains zeroes. Using log10(1 + x) instead.")
    x_norm <- log10(1 + x)
  } else {
    x_norm <- log10(x)
  }
  
  x_norm
}

# the data is transformed using log10(1 + x)
transform_log10p <- function(x) {
  log10(1 + x)
}


get_features_labels <- function(features, label_level, max_label_len) {
  purrr::map_chr(features, 
                 ~ get_feature_label(.x, label_level, max_label_len))
}

#' get the label of a single feature
#' @noRd
get_feature_label <- function(feature,
                              label_level = 1,
                              max_label_len = 60,
                              sep = "|") {
  if (length(feature) != 1) {
    stop("`feature` muste be a character vector of length 1")
  }
  if (label_level == 0) {
    feature <- feature
  } else {
    feature <- strsplit(feature, split = sep, fixed = TRUE) %>%
      unlist() %>%
      rev()
    feature_level <- length(feature)
    feature <- ifelse(
      label_level > feature_level,
      paste(rev(feature[seq_len(feature_level)]), collapse = sep),
      paste(rev(feature[seq_len(label_level)]), collapse = sep)
    )
  }
  
  feature_len <- nchar(feature)
  if (feature_len > max_label_len) {
    feature_letters <- unlist(strsplit(feature, ""))
    feature <- paste(
      paste(feature_letters[seq_len(max_label_len / 2 - 2)], 
            collapse = ""),
      "..",
      paste(feature_letters[
        (feature_len - max_label_len / 2 + 3):feature_len],
        collapse = ""),
      sep = ""
    )
  }
  # replace "Unknown" label in the species level as "sp."
  feature <- replace_unknown_species(feature)
  
  feature
}

# replace "Unknown" label in the species level as "sp."
replace_unknown_species <- function(feature, sep = "|") {
  species_hased <- grepl("s__", feature, fixed = TRUE)
  if (!species_hased) {
    return(feature)
  }
  
  taxa_lvl <- strsplit(feature, sep, fixed = TRUE)
  n_lvl <- length(taxa_lvl)
  sp <- taxa_lvl[[n_lvl]]
  sp <- gsub("Unknown", "sp.", feature, fixed = TRUE)
  taxa_lvl[[n_lvl]] <- sp
  feature <- paste(taxa_lvl, collapse = sep)
  
  feature
}

############################################################
# HEATMAP STARTS HERE
################################################
# this code was copied from the plot_heatmap function

sampledata <- sample_data(ps_mfn)

deseq2.2 <- transform_abundances(deseq2, transform = "log10")
df <- marker_table(deseq2.2)
df1 <- subset (df, ef_logFC >= 2 | ef_logFC <= -2)
rows <- c("g_Turicibacter", "g_Murimonas", "s_Clostridium.XVIII.cocleatum",
          "s_Helicobacter.typhlonius", "g_Escherichia/Shigella", "g_Marvinbryantia",
          "g_Aestuariispira", "g_Akkermansia", "g_Merdimonas", "g_Clostridium XVIII",
          "f_Peptostreptococcaceae", "g_Mucispirillum", "p_Candidatus_Saccharibacteria",
          "g_Romboutsia", "g_Faecalibaculum", "g_Clostridium XIVb", "g_Paludicola",
          "s_Eggerthella.lenta", "g_Robinsoniella", "g_Anaerobutyricum", "g_Rodentibacter",
          "g_Clostridioides", "s_Lactobacillus.intestinalis",  "g_Fournierella", "g_Bifidobacterium", 
          "g_Streptococcus", "c_Deltaproteobacteria", "g_Anaerostipes", "g_Parasutterella", 
          "g_Eggerthella", "s_Bacteroides.thetaiotamicron", "g_Mailhella", "g_Falsiporyphyromonas",
          "g_Herbinix", "f_Lactobacillaceae", "s_Robinsoniella.peoriensis", "g_Hydrogenoanaerobacterium",
          "g_Ruminococcus", "g_Helicobacter", "s_Anaerostipes.caccae", "g_Caproiciproducens",
          "s_Anaerotruncus.colihominis")



abd <- as(otu_table(deseq2.2), "matrix")
marker_abd <- abd[match(df1$feature, row.names(abd)), ] %>%
  as.data.frame()

labels <- get_features_labels(
  row.names(marker_abd),
  1,
  160
)
row.names(marker_abd) <- labels


annotationData <- data.frame(sampledata$Condition, sampledata$Sample_Accession)
rownames(annotationData) <- annotationData[,2]
samples <- match(rownames(annotationData), colnames(marker_abd))
marker_abd_ordered <- marker_abd[, samples]

abundance_colnames <- colnames(marker_abd_ordered)

annotationData <- annotationData[annotationData$sampledata.Sample_Accession %in% abundance_colnames,] # Remove the rows of annotation data that don't belong in the data set
annotationData <- select(annotationData, 1) # select all the columns by index e.g. 1,2 except for the sample_accession column
colnames(annotationData) <- c('Condition') 
annotation_colours <- list(
  'Condition' = c("Colitogenic" = "#FFC300", "Healthy" = "#1872B0","Colitis" = "#FA7F4A")
)

ht_opt$HEATMAP_LEGEND_PADDING = unit(1,"cm")


colAnnotation <- HeatmapAnnotation(
  df = annotationData,
  which = 'col',
  col = annotation_colours,
  simple_anno_size = unit(0.2, "cm"),
  annotation_name_gp = gpar(
    fontface = "bold", 
    fontsize = 8)
)

max(marker_abd)
min(marker_abd)

col_fun = colorRamp2(c(0,5), c("white", "#C70039"))



hmap <- Heatmap(
  marker_abd_ordered,
  name = "Abundance (log10)",
  cluster_rows = FALSE, #check
  cluster_columns = FALSE, #check
  column_order = order(annotationData$Condition),
  top_annotation = colAnnotation,
  col_fun,
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(direction = "horizontal",
                              labels_gp = gpar(fontsize = 9)
  ),
  height = unit(4, "mm") * nrow(marker_abd_ordered),
  show_column_names = FALSE
) + rowAnnotation(labels=anno_text(rows, which="row"),
                  width = unit(6,"cm"))

pdf(file="DESeq2_semi_clustering_DSS_MMUPHin(Condition+simplelabels+log10).pdf", height=8, width = 15)
draw(hmap, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", legend_grouping = "original")
dev.off()
