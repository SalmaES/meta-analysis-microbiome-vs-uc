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



ps_mfn <- readRDS("batch_adjusted_phyloseq_object_MMUPHin.rds")

pdf(file="PCoA_all_post_MMUPHin(Condition).pdf")
plot_ordination(ps_mfn, ordinate(ps_mfn, "PCoA", distance = "bray"), color = "Condition")
dev.off()

pdf(file="PCoA_all_post_MMUPHin(Colitis Model).pdf")
plot_ordination(ps_mfn, ordinate(ps_mfn, "PCoA", distance = "bray"), color = "Colitis_Model")
dev.off()

mfn.permanova.condition <- microbial::betatest(ps_mfn, group="Condition", distance = "bray")   
write.csv(mfn.permanova.condition, "PERMANOVA.bray.all.post.MMUPHin(Condition).csv")

mfn.permanova.cm <- betatest(ps_mfn, group="Colitis_Model", distance = "bray")   
write.csv(mfn.permanova.cm, "PERMANOVA.bray.all.post.MMUPHin(Colitis_Model).csv")