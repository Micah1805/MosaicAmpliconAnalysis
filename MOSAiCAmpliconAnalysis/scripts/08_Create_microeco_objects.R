# Script: 08_Create_microeco_objects.R
# Purpose: Create a microeco object used for further analysis (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

####Housekeeping####
#Load packages
library(microeco)

# allow more waiting time to download each package
options(timeout = 1000)
tmp <- c("microeco", "mecoturn", "MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "igraph", "picante", "pheatmap", "rgexf", 
         "ggalluvial", "ggh4x", "rcompanion", "FSA", "gridExtra", "aplot", "NST", "GGally", "ggraph", "networkD3", "poweRlaw", "ggtern", "SRS")

#Check or install
for(x in tmp){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

library("BiocManager")
library("file2meco")
library("MicrobiomeStat")
library("WGCNA") 
library("ggtree") #Bioconductor
library("metagenomeSeq")
library("ALDEx2")
library("ANCOMBC")
library(SpiecEasi)
library(mixedCCA)
library(SPRING)
library(NetCoMi)
library(beemStatic)
library(chorddiag)
library(ggradar)
library(ggnested)
library(ggcor)
library(RJSONIO)
library("biomformat")
library(seqinr)
library(tidyverse)
library(ape)
library(magrittr)

#Fix random number generation to make the results repeatable
set.seed(123)

#Set plotting background
library(ggplot2)
theme_set(theme_bw())

#Set working environment
path="./output"

####Prokaryotes####
#Create microtable object
#Abundance data
normProkAbund <- read.csv("./output/asv_norm_prokaryotes.csv", row.names=1, check.names = FALSE)
otu_table_16S <- normProkAbund
otu_table_16S <- as.data.frame(otu_table_16S)
class(otu_table_16S)

#Taxa data
prok_taxa <- read.csv("./data/taxonomy_prokaryotes.csv", row.names=1)
taxonomy_table_16S <- prok_taxa
taxonomy_table_16S <- as.data.frame(taxonomy_table_16S)
class(taxonomy_table_16S)
taxonomy_table_16S <- taxonomy_table_16S %<>% tidy_taxonomy()

#Environmental data
Env_Prok <- read.csv("./output/environmental_data_prokaryotes.csv", row.names=1)
sample_info_16S <- Env_Prok
sample_info_16S <- as.data.frame(sample_info_16S)
class(sample_info_16S)

dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S)
dataset

#Remove ASVs assigned as “mitochondria” or “chloroplast”
dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dataset

#Make the ASV and sample information consistent across all files
dataset$tidy_dataset()
print(dataset)

#Check the sequence numbers in each sample
dataset$sample_sums() %>% range

#Save all the basic data in microtable object to local files (e.g. feature abundance, metadata, taxonomic table, representative sequences)
dataset$save_table(dirpath = paste0(path, "/Microeco/basic_files"), sep = ",")

#Calculate the taxa abundance at each taxonomic rank
dataset$cal_abund()
class(dataset$taxa_abund)

#Optimal: Show part of the relative abundance at Division level
dataset$taxa_abund$Phylum[1:4, 1:4]
dataset$save_abund(dirpath = paste0(path, "/Microeco/taxa_abund"))

#Merge abundance table into one
dataset$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)
# remove those unclassified
dataset$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)

#Calculate the alpha diversity
dataset$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset$alpha_diversity)
class(dataset$alpha_diversity)

dataset$save_alphadiv(dirpath = paste0(path, "/Microeco/alpha_diversity"))

#Beta diversity
dataset$cal_betadiv(unifrac = FALSE)
class(dataset$beta_diversity)

dataset$save_betadiv(dirpath = paste0(path, "/Microeco/beta_diversity"))

save(dataset,file=paste0(path, "/microeco_prok.Rdata"))


####Eukaryotes####
#Create a microecotable object
#Abundance table
normEukAbund <- read.csv("./output/asv_norm_eukaryotes.csv", row.names=1, check.names = FALSE)
otu_table_18S <- normEukAbund
otu_table_18S <- as.data.frame(otu_table_18S)
class(otu_table_18S)

#Taxa table
euk_taxa <- read.csv("./data/taxonomy_eukaryotes.csv", row.names=1)
taxonomy_table_18S <- euk_taxa
taxonomy_table_18S <- as.data.frame(taxonomy_table_18S) 
class(taxonomy_table_18S)
taxonomy_table_18S <- taxonomy_table_18S %<>% tidy_taxonomy

#Environmental data
Env_Euk <- read.csv("./output/environmental_data_eukaryotes.csv", row.names=1)
sample_info_18S<-Env_Euk
sample_info_18S<-as.data.frame(sample_info_18S)
class(otu_table_18S)

dataset2 <- microtable$new(sample_table = sample_info_18S, otu_table = otu_table_18S, tax_table = taxonomy_table_18S)
dataset2

#Remove ASVs assigned as “mitochondria” or “chloroplast”
dataset2$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dataset2

#Make the ASV and sample information consistent across all files
dataset2$tidy_dataset()
print(dataset2)

#Optimal: Check the sequence numbers in each sample
dataset2$sample_sums() %>% range

#Save all the basic data in microtable object to local files (e.g. feature abundance, metadata, taxonomic table, representative sequences)
dataset2$save_table(dirpath = paste0(path, "/Microeco/basic_files", sep = ","))

#Calculate the taxa abundance at each taxonomic rank
dataset2$cal_abund()
class(dataset2$taxa_abund)

#Optimal: Show part of the relative abundance at Division level
dataset2$taxa_abund$Division[1:4, 1:4]

dataset2$save_abund(dirpath = paste0(path, "/Microeco/taxa_abund"))

#Merge abundance table into one
dataset2$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)

# remove those unclassified
dataset2$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)

#Calculate the alpha diversity
dataset2$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset2$alpha_diversity)
class(dataset2$alpha_diversity)

dataset2$save_alphadiv(dirpath = paste0(path, "/Microeco/alpha_diversity"))

#Beta diversity
dataset2$cal_betadiv(unifrac = FALSE)
class(dataset2$beta_diversity)

dataset2$save_betadiv(dirpath = paste0(path, "/Microeco/beta_diversity"))
save(dataset2,file=paste0(path, "/microeco_euk.Rdata"))


#####Create a microeco object for Prokaryotes and Eukaryotes together####
#Create microtable object
normProkEukAbund <- read.csv("./output/asv_norm_prokaryotes_eukaryotes.csv", row.names=1, check.names = FALSE)
Taxa_Prok_Euk <- read.csv("./output/taxonomy_prokaryotes_eukaryotes.csv", row.names=1)
Env_Prok_Euk<- read.csv("./output/environmental_data_prokaryotes_eukaryotes.csv", row.names=1)

#ASV table
otu_table<-normProkEukAbund
otu_table<-as.data.frame(otu_table)
class(otu_table)

#Taxonomy table
taxonomy_table<-Taxa_Prok_Euk
taxonomy_table<-as.data.frame(taxonomy_table)
class(taxonomy_table)
taxonomy_table<-taxonomy_table %<>% tidy_taxonomy()

#Sample table
sample_info<-Env_Prok_Euk
sample_info<-as.data.frame(sample_info)
class(otu_table)

dataset3 <- microtable$new(sample_table = sample_info,otu_table = otu_table,tax_table = taxonomy_table)
dataset3

#Remove ASVs assigned as “mitochondria” or “chloroplast”
dataset3$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dataset3

#Make the ASV and sample information consistent across all files
dataset3$tidy_dataset()
print(dataset3)

#Check the sequence numbers in each sample
dataset3$sample_sums() %>% range

#Save all the basic data in microtable object to local files (e.g. feature abundance, metadata, taxonomic table, representative sequences
dataset3$save_table(dirpath = paste0(path, "/Microeco/basic_files"), sep = ",")

#Calculate the taxa abundance at each taxonomic rank
dataset3$cal_abund()
class(dataset3$taxa_abund)

#Optimal: Show part of the relative abundance at Subdivision_Phylum level
dataset3$taxa_abund$Subdivision_Phylum[1:4, 1:4]

dataset3$save_abund(dirpath = paste0(path, "/Microeco/taxa_abund"))

#Merge abundance table into one
dataset3$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)

#Remove those unclassified
dataset3$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)

#Calculate the alpha diversity
dataset3$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset3$alpha_diversity)
class(dataset3$alpha_diversity)

dataset3$save_alphadiv(dirpath = paste0(path, "/Microeco/alpha_diversity"))

#Beta diversity
dataset3$cal_betadiv(unifrac = FALSE)
class(dataset3$beta_diversity)

dataset3$save_betadiv(dirpath = paste0(path, "/Microeco/beta_diversity"))

save(dataset3,file=paste0(path, "./microeco_prok_euk.Rdata"))



