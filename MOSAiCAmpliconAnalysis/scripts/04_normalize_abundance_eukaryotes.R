# Script: 04_normalize_abundance_eukaryotes.R
# Purpose: Normalizes abundance values for eukaryotic ASVs.
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

#Load required libraries
library(stringr)
library(plyr)
library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)

#Optional: Define working environment
path="./output"

#####Prepare tables for microeco####
#Load data
#Metadata
Meta_Euk <- read.csv("./data/Metadata_18S.csv", sep=";") 

#Raw abundances
abundMatRaw_euk <- read.csv("./output/asv_table_eukaryotes.csv", row.names=1)

#Prepare Metadata
#Define sample names ( extracted from read name to be in line with Fastq-files used  for ASV table)
sampNames <- Meta_Euk$`forward_read_file_name`
sample <- sub("_S.*",'',sampNames)
sample

sample2 <- sub("-corr.*",'',sample)
sample2

#Subset metadata and add sample names
Euk_Date<-Meta_Euk %>% dplyr::select("collection.date..start.",
                                     "geographic.location..longitude.",
                                     "geographic.location..latitude.",
                                     "geographic.location..depth.",
                                     sample_description,"environmental.package") 
Euk_Date$ID <- sample2
Euk_Date <- na.omit(Euk_Date)
Euk_Date$ID <- str_replace_all(Euk_Date$ID,"-",".")
Euk_Date$ID <- str_replace_all(Euk_Date$ID,"_",".")
Euk_Date <- as.data.frame(Euk_Date)
class(Euk_Date)

#Extract Leg number from sample desription
Euk_Date$sample_description <- sub("PS122/",'',Euk_Date$sample_description)
Euk_Date$sample_description <- sub("*_(\\D*\\d+)-(\\D*\\d+)",'',Euk_Date$sample_description)
colnames(Euk_Date)[5]<-"Leg"

#Prepare abundance data
Euk <- as.data.frame(t(abundMatRaw_euk))
Euk <- tibble::rownames_to_column(Euk, "ID")
#Euk$ID <- str_replace_all(Euk_Date$ID,"-",".")
#Euk$ID <- str_replace_all(Euk_Date$ID,"_",".")
class(Euk)

##Merge metadata and abundance tables
Abundance_merged_Euk <- dplyr::inner_join(Euk_Date, Euk, by ='ID')
#Outcome includes only 107 samples as one sample name was icluded twice (two Fasta-files with same sample name, but different station numbers) 
Abundance_merged_Euk <- Abundance_merged_Euk%>% relocate("ID")
#save(Abundance_merged_Euk,file=paste0(path,"/Abundance_Euk_water_CTD.Rdata"))

####Create depth categories####
#Make sample depth numeric
Abundance_merged_Euk$geographic.location..depth. <- as.numeric(Abundance_merged_Euk$geographic.location..depth.)

#Create rounded depth values
Abundance_merged_Euk$Depth_rounded <- round(Abundance_merged_Euk$geographic.location..depth.)
Abundance_merged_Euk <- Abundance_merged_Euk %>% relocate("Depth_rounded")
Abundance_merged_Euk$geographic.location..depth. <- as.character(Abundance_merged_Euk$geographic.location..depth.)
Abundance_merged_Euk <- Abundance_merged_Euk[order(Abundance_merged_Euk$Depth_rounded),]


#Sort samples into different depth ranges (Categories)
TopBottom = gsub( " .*$", "", Abundance_merged_Euk$Depth_rounded)
TopBottom[which(TopBottom == "1"| 
                  TopBottom == "2"| 
                  TopBottom == "5"| 
                  TopBottom == "surface"|
                  TopBottom == "10")]="Surface"


TopBottom[which(TopBottom == "chlmax"|
                  TopBottom == "12"|
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
                  TopBottom == "52")]=">50m"

TopBottom[which(TopBottom == "100"|
                  TopBottom == "101")]=">75m"


#Assign categories
Abundance_merged_Euk$Category = paste0("Water","_",TopBottom)
Water_sort_Euk <- Abundance_merged_Euk %>% relocate("Category") 
Water_sort_Euk <- Water_sort_Euk[order(Water_sort_Euk$Depth_rounded),]
Water_sort_Euk <- Water_sort_Euk[,-2]
Water_sort_Euk <- Water_sort_Euk %>% relocate("ID")


#Create a new column Date_Category
Water_sort_Euk$Date_Category <- paste(Water_sort_Euk$collection.date..start.,Water_sort_Euk$Category)
Water_sort_Euk$Date_Category <- str_replace_all(Water_sort_Euk$Date_Category,"-","_")
Water_sort_Euk$Date_Category <- str_replace_all(Water_sort_Euk$Date_Category," ","_")
Water_sort_Euk <- Water_sort_Euk %>% relocate("Date_Category")

#save for later work
#save(Water_sort_Euk,file=paste0(path,"/Water_sort_Euk.Rdata"))
write.csv(Water_sort_Euk,file=paste0(path,"/asv_categorised_eukaryotes.csv"), row.names=TRUE)

####Transposed categorized abundance tables (to get abundance tables only containing the ASVs)####
#If package plyr is loaded dplyr my not work accurate in next step, therefore detach it
detach("package:plyr", unload = TRUE)

#Remove duplicated samples (by calculating the mean)
#After averaging only 91 samples are left
Euk_water_m <- Water_sort_Euk %>% group_by(Date_Category,Category) %>% dplyr::summarise(across(c(8:6136),mean)) %>% distinct 
Euk_water_m <- Euk_water_m[order(Euk_water_m$Category,decreasing=TRUE),]
Euk_water_m <- Euk_water_m %>% remove_rownames %>% column_to_rownames(var="Date_Category")
Euk_water_m <- Euk_water_m[,-1]
Euk_abund_water<-as.data.frame(t(Euk_water_m))
#save(Euk_abund_water,file=paste0(path,"/Abundance_Euk.Rdata"))
write.csv(Euk_abund_water,file=paste0(path,"/asv_abundance_eukaryotes.csv"), row.names=TRUE)

####Normlize data####
#Hellinger transformation and normalization abundance data
normEukAbund <- apply(Euk_abund_water,2,function(x)sqrt(x/sum(x)))
#save(normEukAbund,file=paste0(path,"/asv_norm_eukaryotes.Rdata"))
write.csv(normEukAbund,file=paste0(path,"/asv_norm_eukaryotes.csv"), row.names=TRUE)

