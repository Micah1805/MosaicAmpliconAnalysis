################################################################################
##############Prepare tables for prokaryotes and eukaryotes together############
################################################################################

#Load required packages
library(dplyr)
library(tidyr)

####Prepare abundance table####
#Load data
load("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/Abundance_Euk.Rdata")
load("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/Abundance_Prok.Rdata")

Abundance_Prok_Euk<-bind_rows(as.data.frame(Euk_abund_water),as.data.frame(Prok_abund_water))
#Remove NA (if some were produced)
Abundance_Prok_Euk_no_NA<-Abundance_Prok_Euk %>% select_if(~ !any(is.na(.)))

save(Abundance_Prok_Euk,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Abundance_Prok_Euk_NA.Rdata"))
save(Abundance_Prok_Euk_no_NA,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Abundance_Prok_Euk_no_NA.Rdata"))

#Normlize abundances
#Hellinger transformation and normalization abundance data
normProkEukAbund<-apply(Abundance_Prok_Euk_no_NA,2,function(x)sqrt(x/sum(x)))

save(normProkEukAbund,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/normProkEukAbund.Rdata"))


####Merge taxa tables prokaryotes and eukaryotes####
#Load data
load("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/taxa_euk.Rdata")
load("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/taxa_prok.Rdata")

#Remove "Supergroup" and rename "Devision" column from Eukaryotes
euk_taxa<-euk_taxa[,-2]
colnames(euk_taxa)[3] ="Subdivision_Phylum"
euk_taxa<-euk_taxa[,-2]

#Rename Phylum column from Prokaryotes
colnames(prok_taxa)[2] ="Subdivision_Phylum"

Taxa_Prok_Euk<-rbind(as.data.frame(prok_taxa),as.data.frame(euk_taxa))

save(Taxa_Prok_Euk,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Taxa_Prok_Euk.Rdata"))
