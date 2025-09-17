################################################################################
####Phaseplot III-Calculate frequency table for Prokaryotes and Eukaryotes######
################################################################################

#Housekeeping
library(ggplot2)
library(tidyverse)
library(stringr)
library(xlsx)

#Create a table including Kingdom, ASV and frequency
taxa <- tibble::rownames_to_column(taxa, "ASV")
yy <- tibble::rownames_to_column(yy, "ASV")
Frequ <- full_join(taxa,yy,by="ASV") 
Frequ<- Frequ %>% relocate("freq") %>% relocate("Kingdom") %>% relocate("ASV")
Frequ<-Frequ[1:3]
Frequ<-Frequ[order(Frequ$freq),]

#Rename all Bacteria and Archaea into Prokaryotes (to have it uniform with other data presented)
Frequ$Kingdom <- str_replace(Frequ$Kingdom, "Bacteria", "Prokaryotes")
Frequ$Kingdom <- str_replace(Frequ$Kingdom, "Archaea", "Prokaryotes")

#barplot(Frequ)

#str(Frequ)

#Count number of Archaea, Bacteria, Eukaryotes occuring at the specific pecentages
Frequ_sort<-Frequ %>% count( Kingdom,freq,sort = TRUE)

#Calculate percentage of Archaea, Bacteria, Eukaryotes occuring at the specific pecentages
Frequ_perc<-Frequ %>% count( Kingdom,freq) %>%
  group_by(Kingdom) %>%
  mutate(percent = n/sum(n))
Frequ_perc<-Frequ_perc[order(Frequ_perc$freq),]

#str(Frequ_perc)

#Rename columns
colnames(Frequ_perc)[2] <- 'Frequency'
colnames(Frequ_perc)[4] <- 'Percentage'

#Transfor decimal percentage into percentual percentage 
Frequ_perc$Percentage_perc<-Frequ_perc$Percentage*100


#Graphical summary
ggplot(Frequ_perc,aes(x=Kingdom, y = Percentage_perc,group=Frequency,fill = Frequency)) + 
  geom_bar(stat="identity",position="dodge", colour="black")+
  scale_y_continuous(breaks=seq(0,100,5))+
  labs(y="Number of ASVs [%]")

#Save data
#save(Frequ_Prok_Euk,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/Frequ_Prok_Euk.Rdata"))
save(Frequ_perc,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/Frequ_perc_100m.Rdata"))

#Export frequency table
write.xlsx(as.data.frame(Frequ_perc), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Frequ_100m.xlsx",
           col.names = TRUE, row.names = TRUE)




















Frequ_Prok_Euk<-cbind(taxa$Kingdom,taxa$ASV,Prok_Euk_data$phaseFreqInfo$freq)
Frequ_Prok_Euk<-as.data.frame(Frequ_Prok_Euk)

Frequ_Prok_Euk$V1 <- str_replace(Frequ_Prok_Euk$V1, "Bacteria", "Prokaryotes")
Frequ_Prok_Euk$V1 <- str_replace(Frequ_Prok_Euk$V1, "Archaea", "Prokaryotes")

Frequ_Prok_Euk<-Frequ_Prok_Euk[order(Frequ_Prok_Euk$V3),]
#barplot(Frequ)

str(Frequ_Prok_Euk)
Frequ_Prok_Euk$V3<-as.numeric(Frequ_Prok_Euk$V3)



#Frequ_Prok_Euk<-Frequ_Prok_Euk %>% count( V1,V3,sort = TRUE)


Frequ_Prok_Euk<-Frequ_Prok_Euk %>% count( V1,V3) %>%
  group_by(V1) %>%
  mutate(percent = n/sum(n))
Frequ_Prok_Euk<-Frequ_Prok_Euk[order(Frequ_Prok_Euk$V3),]




colnames(Frequ_Prok_Euk)[1] <- 'Kingdom'
colnames(Frequ_Prok_Euk)[2] <- 'Frequency'
colnames(Frequ_Prok_Euk)[4] <- 'Percentage'

Frequ_Prok_Euk$Percentage_perc<-Frequ_Prok_Euk$Percentage*100


ggplot(Frequ_Prok_Euk,aes(x=Kingdom, y = Percentage_perc,group=Frequency,fill = Frequency)) + 
  geom_bar(stat="identity",position="dodge", colour="black")+
  scale_y_continuous(breaks=seq(0,100,5))+
  labs(y="Number of ASVs [%]")


save(Frequ_Prok_Euk,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/Frequ_Prok_Euk.Rdata"))
save(Frequ_perc,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/Frequ_perc.Rdata"))

write.xlsx(as.data.frame(Frequ_Prok_Euk), file="Frequ_Prok_Euk.xlsx",
           col.names = TRUE, row.names = TRUE)


class(Frequ_Prok_Euk)
