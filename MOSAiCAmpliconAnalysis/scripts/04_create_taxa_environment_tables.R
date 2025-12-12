# Script: 04_create_taxa_environment_tables.R
# Purpose: Creates taxonomic and environmental metadata tables.
# Input: taxonomy_eukaryotes.csv, 
#        asv_categorised_prokaryotes.csv, asv_categorised_eukaryotes.csv
# Main output: taxonomy_eukaryotes_clean.csv,
#         environmental_data_prokaryotes.csv, environmental_data_eukaryotes.csv,
#         environmental_data_prokaryotes_eukaryotes.csv
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

####Housekeeping####
#Load required libraries
library(stringr)
library(plyr)
library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)

####Taxa tables####
#Clean taxa tables eukaryotes
#Load data
euk_taxa <- read.csv("./output/taxonomy_eukaryotes.csv", row.names=1)

#Rename euk_taxa columns
colnames(euk_taxa)[2] ="Supergroup"
colnames(euk_taxa)[3] ="Division"
colnames(euk_taxa)[4] ="Subdivision"
colnames(euk_taxa)[5] ="Class"
colnames(euk_taxa)[6] ="Order"
colnames(euk_taxa)[7] ="Family"
colnames(euk_taxa)[8] ="Genus"
colnames(euk_taxa)[9] ="Species"

#save(euk_taxa,file=paste0(path,"/Taxa_Euk.Rdata"))
write.csv(euk_taxa,file="./output/taxonomy_eukaryotes_clean.csv", row.names=TRUE)


####Environmental data####
#Load data
Water_sort_Prok <- read.csv("./output/asv_categorised_prokaryotes.csv", row.names=1)
Water_sort_Prok$`geographic location (depth)`<- as.numeric(Water_sort_Prok$`geographic location (depth)`)
Water_sort_Prok$Leg <- as.numeric(Water_sort_Prok$Leg)

Water_sort_Euk <- read.csv("./output/asv_categorised_eukaryotes.csv", row.names=1) 
Water_sort_Euk$`geographic location (depth)` <- as.numeric(Water_sort_Euk$`geographic location (depth)`)
Water_sort_Euk$Leg <- as.numeric(Water_sort_Euk$Leg)

#Transform tables
#If package plyr is loaded dplyr my not work accurate in next step, therefore detach it
detach("package:plyr", unload = TRUE)

#Create environmental table for prokaryotes
Env_Prok <- Water_sort_Prok %>% group_by(Date_Category,Category) %>% dplyr::summarise(across(c(3:5602),mean)) %>% distinct 
Env_Prok <- Env_Prok[order(Env_Prok$Category,decreasing=TRUE),]
Env_Prok$Date_Category <- str_replace(Env_Prok$Date_Category,"_","-")
Env_Prok$Date_Category <- str_replace(Env_Prok$Date_Category,"_","-")
Env_Prok$Date_Category <- str_replace(Env_Prok$Date_Category,"_"," ")
Env_Prok[c("Date", "Category")] <- str_split_fixed(Env_Prok$Date_Category, " ", 2) 
Env_Prok$Date_Category <- str_replace(Env_Prok$Date_Category," ","_")
Env_Prok <- Env_Prok %>% relocate("Leg") %>% relocate(`geographic location (depth)`)
Env_Prok <- Env_Prok %>% relocate(geographic location (longitude)") %>% relocate("geographic location (latitude)")
Env_Prok <- Env_Prok %>% relocate("Date") %>% relocate("Date_Category")
Env_Prok <- Env_Prok[1:7]
Env_Prok$Date_Category <- str_replace_all(Env_Prok$Date_Category,"-","_")
Env_Prok <- Env_Prok %>% remove_rownames %>% column_to_rownames(var="Date_Category")

#save(Env_Prok,file="./output/environmental_data_prokaryotes.Rdata")
write.csv(Env_Prok,file="./output/environmental_data_prokaryotes.csv", row.names=TRUE)

#Create environmental table for eukaryotes
Env_Euk <- Water_sort_Euk %>% group_by(Date_Category,Category) %>% dplyr::summarise(across(c(3:6136),mean)) %>% distinct 
Env_Euk <- Env_Euk[order(Env_Euk$Category,decreasing=TRUE),]
Env_Euk$Date_Category <- str_replace(Env_Euk$Date_Category,"_","-")
Env_Euk$Date_Category <- str_replace(Env_Euk$Date_Category,"_","-")
Env_Euk$Date_Category <- str_replace(Env_Euk$Date_Category,"_"," ")
Env_Euk[c("Date", "Category")] <- str_split_fixed(Env_Euk$Date_Category, " ", 2) 
Env_Euk$Date_Category<-str_replace(Env_Euk$Date_Category," ","_")
Env_Euk <- Env_Euk %>% relocate("Leg") %>% relocate("geographic location (depth)")
Env_Euk <- Env_Euk %>% relocate("geographic location (longitude)") %>% relocate("geographic location (latitude)")
Env_Euk <- Env_Euk %>% relocate("Date") %>% relocate("Date_Category")
Env_Euk <- Env_Euk[1:7]
Env_Euk$Date_Category <- str_replace_all(Env_Euk$Date_Category,"-","_")
Env_Euk <- Env_Euk %>% remove_rownames %>% column_to_rownames(var="Date_Category")

#save(Env_Euk,file=paste0(path,"/Env_data_Euk.Rdata"))
write.csv(Env_Euk,file="./output/environmental_data_eukaryotes.csv", row.names=TRUE)

#Build one table for prokaryotes and eukaryotes
Env_Euk <- tibble::rownames_to_column(Env_Euk, "Date_Category")
Env_Prok <- tibble::rownames_to_column(Env_Prok, "Date_Category")

Env_Euk_Prok <- full_join(Env_Prok,Env_Euk,by="Date_Category")
Env_Euk_Prok <- Env_Euk_Prok %>% remove_rownames %>% column_to_rownames(var="Date_Category")
Env_Euk_Prok <- Env_Euk_Prok[,1:6]
colnames(Env_Euk_Prok)[1] ="Date"
colnames(Env_Euk_Prok)[2] ="Lat"
colnames(Env_Euk_Prok)[3] ="Long"
colnames(Env_Euk_Prok)[4] ="Depth"
colnames(Env_Euk_Prok)[5] ="Leg"
colnames(Env_Euk_Prok)[6] ="Category"

#save(Env_Euk_Prok,file=paste0(path,"/Env_Euk_Prok.Rdata"))
write.csv(Env_Euk_Prok,file="./output/environmental_data_prokaryotes_eukaryotes.csv", row.names=TRUE)


