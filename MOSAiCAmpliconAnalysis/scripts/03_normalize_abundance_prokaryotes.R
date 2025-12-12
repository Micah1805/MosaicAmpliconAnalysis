# Script: 01_normalize_abundance_prokaryotes.R
# Purpose: Normalizes abundance values for prokaryotic ASVs.
# Input: Metadata_16S.csv, asv_table_prokaryotes.csv
# Main output: asv_categorised_prokaryotes.csv,
#              asv_abundance_prokaryotes.csv, asv_norm_prokaryotes.csv
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

#Load required libraries
library(stringr)
library(plyr)
library(tidyverse)
library(readxl)

#Define working environment
path="./output"

####Create Abundance tables separated by environment (Water and ice)####

#Load data
#Metadata
Meta_data <- read.csv("./data/Metadata_16S.csv", sep=";")
#str(Meta_data)

#Raw abundance data
abundMatRaw_prok <- read.csv("./output/asv_table_prokaryotes.csv", row.names=1)

##Prepare Metadata##
#Define sample names (extracted from read name to be in line with Fastq-files used for ASV table)
sampNames <- Meta_data$`forward_read_file_name`
sample <- sub("_S.*",'',sampNames)
sample

sample2 <- sub("-corr.*",'',sample)
sample2

#Subset metadata and add sample names
t_Met <- Meta_data %>% dplyr::select("collection date",
                                   "geographic location (latitude)",
                                   "geographic location (longitude)",
                                   "geographic location (depth)",
                                   sample_description) 
t_Met$ID <- sample2
t_Met <- na.omit(t_Met)
t_Met$ID <- str_replace_all(t_Met$ID,"-",".")
t_Met$ID <- str_replace_all(t_Met$ID,"_",".")
t_Met <- as.data.frame(t_Met)
class(t_Met)

#Extract Leg number from sample desription
t_Met$sample_description <- sub("PS122/",'',t_Met$sample_description)
t_Met$sample_description <- sub("*_(\\D*\\d+)-(\\D*\\d+)",'',t_Met$sample_description)
colnames(t_Met)[5] <- "Leg"

##Prepare abundance data##
t_Abund <- as.data.frame(t(abundMatRaw_prok))
t_Abund <- tibble::rownames_to_column(t_Abund, "ID")

##Merge metadata and abundance tables
Abundance_merged <- dplyr::inner_join(t_Met, t_Abund, by ='ID')
#Outcome includes only 107 samples as one sample name was included twice (two fasta-files with same sample name, but different station numbers) 
Abundance_merged <- Abundance_merged%>% relocate("ID")
#save(Abundance_merged,file=paste0(path,"/Abundance_Prok_water_CTD.Rdata"))

####Create depth categories####
#Make sample depth numeric
Abundance_merged$`geographic location (depth)` <- as.numeric(Abundance_merged$`geographic location (depth)`)

#Create rounded depth values
Abundance_merged$Depth_rounded <- round(Abundance_merged$`geographic location (depth)`)
Abundance_merged <- Abundance_merged %>% relocate("Depth_rounded")
Abundance_merged$`geographic location (depth)` <- as.character(Abundance_merged$`geographic location (depth)`)
Abundance_merged <- Abundance_merged[order(Abundance_merged$Depth_rounded),]


#Sort samples into different depth ranges (Categories)
TopBottom = gsub( " .*$", "", Abundance_merged$Depth_rounded)
TopBottom[which(TopBottom == "1"| 
                  TopBottom == "2"| 
                  TopBottom == "5"| 
                  TopBottom == "10")]="Surface"

TopBottom[which(TopBottom == "12"|
                  TopBottom == "15"|
                  TopBottom == "16"|
                  TopBottom == "18"|
                  TopBottom == "20"|
                  TopBottom == "22"| 
                  TopBottom == "25"| 
                  TopBottom == "30"| 
                  TopBottom == "31"| 
                  TopBottom == "35"| 
                  TopBottom == "36"| 
                  TopBottom == "40"| 
                  TopBottom == "49")]="Chlmax"

TopBottom[which(TopBottom == "50"|
                  TopBottom == "51"|
                  TopBottom == "52"|
                  TopBottom == "52")]="50m"

TopBottom[which(TopBottom == "100"|
                  TopBottom == "101")]="100m"

#Assign categories
Abundance_merged$Category = paste0("Water","_",TopBottom)
Water_sort_Prok <- Abundance_merged %>% relocate("Category") 
Water_sort_Prok <- Water_sort_Prok[order(Water_sort_Prok$Depth_rounded),]
Water_sort_Prok <- Water_sort_Prok[,-2]
Water_sort_Prok <- Water_sort_Prok %>% relocate("ID")


#Create a new column Date_Category
Water_sort_Prok$Date_Category <- paste(Water_sort_Prok$"collection date",Water_sort_Prok$Category)
Water_sort_Prok$Date_Category <- str_replace_all(Water_sort_Prok$Date_Category,"-","_")
Water_sort_Prok$Date_Category <- str_replace_all(Water_sort_Prok$Date_Category," ","_")
Water_sort_Prok<-Water_sort_Prok %>% relocate("Date_Category")
#save(Water_sort_Prok,file=paste0(path,"/Water_sort_Prok.Rdata"))
write.csv(Water_sort_Prok,file=paste0(path,"/asv_categorised_prokaryotes.csv"), row.names=TRUE)

####Transposed categorized abundance tables (to get abundance tables only containing the ASVs)####
class(Water_sort_Prok)
Water_sort_Prok <- as.data.frame(Water_sort_Prok)

class(Water_sort_Prok)

#Remove duplicated samples (by calculating the mean)
#After averaging only 91 samples are left
Water_sort_Prok_m <- Water_sort_Prok %>% group_by(Date_Category,Category) %>% dplyr::summarise(across(c(7:5602),mean)) %>% distinct 
Water_sort_Prok_m <- Water_sort_Prok_m[order(Water_sort_Prok_m$Category,decreasing=TRUE),]
Water_sort_Prok_m <- Water_sort_Prok_m %>% remove_rownames %>% column_to_rownames(var="Date_Category")
Water_sort_Prok_m <- Water_sort_Prok_m[,-1]
Prok_abund_water <- as.data.frame(t(Water_sort_Prok_m))
#save(Prok_abund_water,file=paste0(path,"/Abundance_Prok.Rdata"))
write.csv(Prok_abund_water,file=paste0(path,"/asv_abundance_prokaryotes.csv"), row.names=TRUE)

####Normalize data####
#Hellinger transformation and normalization abundance data
normProkAbund <- apply(Prok_abund_water,2,function(x)sqrt(x/sum(x)))
#save(normProkAbund,file=paste0(path,"/asv_norm_prokaryotes.Rdata"))
write.csv(normProkAbund,file=paste0(path,"/asv_norm_prokaryotes.csv"), row.names=TRUE)



