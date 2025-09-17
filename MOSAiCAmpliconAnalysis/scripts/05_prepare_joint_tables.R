# Script: 05_prepare_joint_tables.R
# Purpose: Combines prokaryotic and eukaryotic tables into unified formats.
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

#Load required packages
library(dplyr)
library(tidyr)

####Prepare abundance table####
#Load data
Prok_abund_water <- read.csv("./output/asv_abundance_prokaryotes.csv", row.names=1,check.names = FALSE)
Euk_abund_water <- read.csv("./output/asv_abundance_eukaryotes.csv", row.names=1,check.names = FALSE)

Abundance_Prok_Euk <- bind_rows(as.data.frame(Euk_abund_water),as.data.frame(Prok_abund_water))
#Remove NA (if some were produced)
Abundance_Prok_Euk_no_NA <- Abundance_Prok_Euk %>% select_if(~ !any(is.na(.)))

#save(Abundance_Prok_Euk,file=paste0("./output","/Abundance_Prok_Euk_NA.Rdata"))
#save(Abundance_Prok_Euk_no_NA,file=paste0("./output","/Abundance_Prok_Euk_no_NA.Rdata"))
write.csv(Abundance_Prok_Euk,file="./output/asv_abundance_prokaryotes_eukaryotes.csv", row.names=TRUE)
write.csv(Abundance_Prok_Euk_no_NA,file="./output/asv_abundance_prokaryotes_eukaryotes_no_NAs.csv", row.names=TRUE)

#Normalize abundances
#Hellinger transformation and normalization abundance data
normProkEukAbund <- apply(Abundance_Prok_Euk_no_NA,2,function(x)sqrt(x/sum(x)))

#save(normProkEukAbund,file=paste0("./output","/normProkEukAbund.Rdata"))
write.csv(normProkEukAbund,file="./output/asv_norm_prokaryotes_eukaryotes.csv", row.names=TRUE)


####Merge taxa tables prokaryotes and eukaryotes####
#Load data
prok_taxa <- read.csv("./output/taxonomy_prokaryotes.csv", row.names=1)
euk_taxa <- read.csv("./output/taxonomy_eukaryotes.csv", row.names=1)

#Remove "Supergroup" and rename "Devision" column from Eukaryotes
euk_taxa <- euk_taxa[,-2]
colnames(euk_taxa)[3] ="Subdivision_Phylum"
euk_taxa <- euk_taxa[,-2]

#Rename Phylum column from Prokaryotes
colnames(prok_taxa)[2] ="Subdivision_Phylum"

Taxa_Prok_Euk <- rbind(as.data.frame(prok_taxa),as.data.frame(euk_taxa))

#save(Taxa_Prok_Euk,file="./output/taxonomy_prokaryotes_eukaryotes.Rdata")
write.csv(Taxa_Prok_Euk,file="./output/taxonomy_prokaryotes_eukaryotes.csv", row.names=TRUE)
