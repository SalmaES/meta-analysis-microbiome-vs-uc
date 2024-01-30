##Analysis Script for Kimura et al., 2020
#######################################################################################################
#Load R libraries
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(Biostrings)
library(dplyr)
library("rstatix")
library(stringr)
library(pheatmap)
library(SensoMineR)
library(microbiome)
library(knitr)
library(purrr)
library("ape")
library(microbiomeutilities)
library(BiocManager)
library(devtools)
library(DESeq2)
library(microbial)
library(vegan)
library(microbiomeMarker)
library(yingtools2)
library(speedyseq)
library(ComplexHeatmap)

#Clear workspace
rm(list = ls())

#Assign variables for heatmap colours
colors_teal <- "#66C2A5"
colors_red <- "#FB8072"
colors_blue <- "#80B1D3"
colors_purple <- "#E78AC3"
colors_lime <- "#A6D854"
colors_yellow <- "#FFFFB3"
colors_brown <- "#E5C494"
colors_grey <- "#B3B3B3"
colors_lilac <- "#BEBADA"
colors_orange <- "#FDB462"
colors_pink <- "#FCCDE5"
colors_green <- "#8DD3C7"


########################### LOADING DATASETS ########################### 

OTU <- readRDS("seqtabnoc_KIMURA.rds") #feature counts table
taxa <- read_csv("taxtab_KIMURA.csv") #taxonomy table
sampledata <- read_csv("adjusted_sampledata_KIMURA.csv") #metadata
tree <- readRDS("Phylogenetic_Tree_KIMURA.rds") #phylogenetic tree


########################### DATA WRANGLING ########################### 
##Assigning row names for feature table, taxonomy table and metadata.
##Check for missing data and matching sample names in feature table and metadata
##Convert data to phyloseq friendly format. Feature table and taxonomy table formatted into matrices and metadata to dataframe.

#Assigning sampledata rownames
sampledata <- as.data.frame(sampledata)
rownames(sampledata) <- sampledata[,4]

#Assigning taxonomy table rownames
taxa <- as.data.frame(taxa)
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-1]

#Check for differences between OTU table and sampledata
rownames(OTU) -> OTU_check
rownames(sampledata) -> sample_check

#Convert to vectors
OTU_check <- unlist(OTU_check)
sample_check <- unlist(sample_check)

#What's in both files
int <- intersect(OTU_check,sample_check)

#What's different in OTU_check
difs <- setdiff(OTU_check,sample_check)

#What's different in sample_check
difss <- setdiff(sample_check,OTU_check)

#Check for differences between OTU table and taxa table
colnames(OTU) -> OTU_check2
rownames(taxa) -> taxa_check

#Convert to vectors
OTU_check2 <- unlist(OTU_check2)
taxa_check <- unlist(taxa_check)

#What's in both files
int <- intersect(OTU_check2,taxa_check)

#What's different in OTU_check
difs <- setdiff(OTU_check2,taxa_check)

#What's different in taxa_check
difss <- setdiff(taxa_check,OTU_check2)

taxamatrix = as.matrix(taxa) 
otumatrix = as.matrix(OTU) 
class(taxamatrix) 
class(otumatrix) 
class(sampledata) 

#Create a phyloseq object
ps <- phyloseq(otu_table(otumatrix, taxa_are_rows=FALSE), tax_table(taxamatrix), 
               sample_data(sampledata), phy_tree(tree$tree))
ps
saveRDS(ps, file = "phyloseq_object_raw_KIMURA.rds")

########################### DATA CLEAN UP ########################### 
#Filter out samples with read counts below 7000
ps <- prune_samples(microbiome::readcount(ps)>=7000, ps) 
ps


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
ps2
#plot_bar(ps2, fill="Phylum")
saveRDS(ps2, file = "phyloseq_raw_no_outliers_KIMURA.rds")

#Rarefy the data
set.seed(19920921)
ps.rarefied <- rarefy_even_depth(ps2, rngseed=FALSE, replace=T)
ps.rarefied
saveRDS(ps.rarefied, file = "phyloseq_rarefied_KIMURA.rds")


#Prevalence Filtering
prevdf = apply(X = otu_table(ps.rarefied),
               MARGIN = ifelse(taxa_are_rows(ps.rarefied), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps.rarefied),
                    tax_table(ps.rarefied))
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps.rarefied, "Phylum"))

pdf(file="Prevalence_Filtering_Plots(Rarefied)_KIMURA.pdf")
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps.rarefied),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
dev.off()

prevalenceThreshold = 0.05 * nsamples(ps.rarefied)
nsamples(ps.rarefied)
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps.rarefied.pruned = prune_taxa(keepTaxa, ps.rarefied)
ps.rarefied.pruned
saveRDS(ps.rarefied.pruned, file = "phyloseq_object_rarefied_pruned_KIMURA.rds")

########################### DATA ANALYSIS ########################### 
#Pairwise testing for rarefied and pruned data using Wilcoxon test
rich.pruned = estimate_richness(ps.rarefied.pruned)
wtp1=pairwise.wilcox.test(rich.pruned$Observed, sample_data(ps.rarefied.pruned)$Condition)
wtp2=pairwise.wilcox.test(rich.pruned$Observed, sample_data(ps.rarefied.pruned)$Genotype)
wtp3=pairwise.wilcox.test(rich$Observed, sample_data(ps.rarefied.pruned)$Intensity_of_Disease)

write.csv(wtp1[["p.value"]], "Wilcoxon.test.rarefied.pruned(Condition).KIMURA.csv")
write.csv(wtp2[["p.value"]], "Wilcoxon.test.rarefied.pruned(Genotype).KIMURA.csv")
write.csv(wtp3[["p.value"]], "Wilcoxon.test.rarefied.pruned(Intensity_of_Disease).KIMURA.csv")

#Ordination of rarefied and pruned data
ordination.rarefied.pruned <- ordinate(ps.rarefied.pruned, "PCoA", distance = "bray")
pdf(file="PCoA_Rarified_Pruned_Bray(genotype)_KIMURA.pdf")
plot_ordination(ps.rarefied.pruned, ordination.rarefied.pruned, color = "Genotype", 
                title = "ASVs_Rarefied+Pruned(Bray)") + theme(aspect.ratio=1)
dev.off()


#Change taxtable column name from "Domain" to "Kingdom" for downstream analysis
ps.rarefied.pruned1 <- ps.rarefied.pruned %>%
  rename_tax_table(Kingdom=Domain) %>%
  rename_with_tax_table(stringr::str_to_title)
saveRDS(ps.rarefied.pruned1, file = "phyloseq_rarefied_pruned1_KIMURA.rds")

#PERMANOVA testing (p value for PCoA plots)
beta.rarified.pruned.bray <- betatest(ps.rarefied.pruned1, group="Genotype", distance = "bray") 
write.csv(beta.rarified.pruned.bray, "PERMANOVA.test.bray.rarefied.pruned(genotype).KIMURA.csv")


#Differential abundance testing with DESeq2
deseq2 <- microbiomeMarker::run_deseq2(ps.rarefied.pruned1,
                                       group = "Genotype",
                                       contrast = c("Opg-/-", "WT"),
                                       taxa_rank = "none",
                                       sfType = "poscounts",
                                       pvalue_cutoff=0.01, 
                                       transform = "identity", 
                                       p_adjust = "fdr")

deseq2_table <- deseq2@marker_table
write.csv(deseq2_table, "DESeq2_table(Genotype)_KIMURA.csv")


####DESeq2 Heatmap
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


############################################################ (One Annotation Bar)

# this code was copied from the plot_heatmap function
deseq2.2 <- transform_abundances(deseq2, transform = "log10")
df <- marker_table(deseq2.2)
#df1 <- subset (df, ef_logFC >= 10 | ef_logFC <= -10)


abd <- as(otu_table(deseq2.2), "matrix")
marker_abd <- abd[match(df$feature, row.names(abd)), ,drop=FALSE] %>% ######### drop = false is needed here as the sub query of abd returns just one row, which data frame interprets as a vector
  as.data.frame()

labels <- get_features_labels(
  row.names(marker_abd),
  1,
  60
)
row.names(marker_abd) <- labels

annotationData <- data.frame(sampledata$Genotype, sampledata$Sample_Accession)
abundance_colnames <- colnames(marker_abd)
annotationData <- annotationData[annotationData$sampledata.Sample_Accession %in% abundance_colnames,] # Remove the rows of annotation data that don't belong in the data set
annotationData <- select(annotationData, 1) # select all the columns by index e.g. 1,2 except for the sample_accession column
colnames(annotationData) <- c('Genotype') 
annotation_colours <- list(
  'Genotype' = c("WT" = colors_blue,"Opg-/-" = colors_pink)
)

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
col_fun = colorRamp2(c(0,3.5), c("white", "#C70039"))

hmap <- Heatmap(
  marker_abd,
  name = "Abundance (log10)",
  cluster_rows = FALSE, #check
  cluster_columns = FALSE, #check
  column_order = order(annotationData$Genotype),
  top_annotation = colAnnotation,
  col_fun,
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title_position = "topleft",
    labels_gp = gpar(fontsize = 9)
  ),
  height = unit(4, "mm") * nrow(marker_abd),
  show_column_names = FALSE
)

pdf(file="DESeq2_Heatmap(Genotype+log10)_KIMURA.pdf")
draw(hmap, merge_legend = TRUE, legend_grouping = "original")
dev.off()
