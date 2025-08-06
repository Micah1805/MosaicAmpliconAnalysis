################################################################################
#####################Prepare environmental data for RDA#########################
################################################################################

#Load packages
library(stringr)
library(plyr)
library(tidyverse)
library(readxl)

#Load environmental data
Env_Meta <- read_excel("~/AWI/MOSAiC/Sequencing_Prok_Euk/Env-Data/Enviroment_Meta_Nutrients_chlorophyll.xlsx",
                       col_types = c("numeric","date","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","text","text","numeric","numeric","text","text","numeric","numeric","numeric","text","text","text"))
View(Env_Meta)
str(Env_Meta)

####Subdivide environmental data by depth####
Env_Meta<- Env_Meta[order(Env_Meta$`Depth [m]`),]

Surface<-Env_Meta[1:1076,]
Surface$Category<-c("Water_Surface")
Surface<-Surface[-1]
Surface$Date<-format(as.POSIXct(Surface$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%Y-%m-%d')
Surface$Time<-format(as.POSIXct(Surface$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%H:%M:%S')
Surface$Date_Category<-paste(Surface$Date,Surface$Category)
Surface$Date_Category<-str_replace_all(Surface$Date_Category,"-","_")
Surface$Date_Category<-str_replace_all(Surface$Date_Category," ","_")
Surface<-Surface %>% relocate("Date_Category","Date")
Surface<- Surface %>%  group_by(Date_Category,Date) %>% summarise_if(is.numeric, mean)

Chlmax<-Env_Meta[1077:2672,]
Chlmax$Category<-c("Water_Chlmax")
Chlmax<-Chlmax[-1]
Chlmax$Date<-format(as.POSIXct(Chlmax$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%Y-%m-%d')
Chlmax$Time<-format(as.POSIXct(Chlmax$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%H:%M:%S')
Chlmax$Date_Category<-paste(Chlmax$Date,Chlmax$Category)
Chlmax$Date_Category<-str_replace_all(Chlmax$Date_Category,"-","_")
Chlmax$Date_Category<-str_replace_all(Chlmax$Date_Category," ","_")
Chlmax<-Chlmax %>% relocate("Date_Category","Date")
Chlmax<- Chlmax %>%  group_by(Date_Category,Date) %>% summarise_if(is.numeric, mean)

Fünfzig<-Env_Meta[2673:3104,]
Fünfzig$Category<-c("Water_>50m")
Fünfzig<-Fünfzig[-1]
Fünfzig$Date<-format(as.POSIXct(Fünfzig$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%Y-%m-%d')
Fünfzig$Time<-format(as.POSIXct(Fünfzig$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%H:%M:%S')
Fünfzig$Date_Category<-paste(Fünfzig$Date,Fünfzig$Category)
Fünfzig$Date_Category<-str_replace_all(Fünfzig$Date_Category,"-","_")
Fünfzig$Date_Category<-str_replace_all(Fünfzig$Date_Category," ","_")
Fünfzig<-Fünfzig %>% relocate("Date_Category","Date")
Fünfzig<- Fünfzig %>%  group_by(Date_Category,Date) %>% summarise_if(is.numeric, mean)

Hundert<-Env_Meta[3105:3555,]
Hundert$Category<-c("Water_>75m")
Hundert<-Hundert[-1]
Hundert$Date<-format(as.POSIXct(Hundert$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%Y-%m-%d')
Hundert$Time<-format(as.POSIXct(Hundert$`Date/Time`,format='%m/%d/%Y %H:%M:%S'),format='%H:%M:%S')
Hundert$Date_Category<-paste(Hundert$Date,Hundert$Category)
Hundert$Date_Category<-str_replace_all(Hundert$Date_Category,"-","_")
Hundert$Date_Category<-str_replace_all(Hundert$Date_Category," ","_")
Hundert<-Hundert %>% relocate("Date_Category","Date")
Hundert<- Hundert %>%  group_by(Date_Category,Date) %>% summarise_if(is.numeric, mean)

Env_tot<-rbind(Surface,Chlmax,Fünfzig,Hundert)
Env_tot<-Env_tot %>% remove_rownames %>% column_to_rownames(var="Date_Category")
Env_tot[is.na(Env_tot)] <- 0
save(Env_tot,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Env_tot.Rdata"))

#Calculate mean values per collection date
Env_mean<-Env_tot %>% relocate("Date")
#Env_Prok_mean<-Env_Prok_mean[-1]
Env_mean<-Env_mean %>% 
  group_by(Date) %>%
  summarise_at(vars(c(1:17)), list(name = mean))
Env_mean<-Env_mean %>% remove_rownames %>% column_to_rownames(var="Date")

save(Env_mean,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Env_mean.Rdata"))


####Add Light to the environmental data####

#Load data and modify it
Light <- read_excel("~/AWI/MOSAiC/Sequencing_Prok_Euk/Env-Data/Light.xlsx")
colnames(Light)[1]<-"Date"
Light$Category<-c("Water_Surface") 
Light$ID<-paste(Light$Date,Light$Category)

Light$ID<-str_replace_all(Light$ID,"-","_")
Light$ID<-str_replace_all(Light$ID," ","_")
Light<-Light %>% relocate("ID")

#Modify environmental data
str(Env_tot)
Env_tot$Date<-as.Date(Env_tot$Date)
Env_tot <- tibble::rownames_to_column(Env_tot, "ID")
colnames(Env_tot)[6] <- "Depth"

#Merge both data sets
Env_tot_Light<-full_join(Env_tot,Light,by=c("Date","Depth","ID"))
Env_tot_Light<-Env_tot_Light[1:20]
Env_tot_Light<-Env_tot_Light[order(Env_tot_Light$Depth),]

#Remove duplicated samles by averaging values
Env_tot_Light_m<- Env_tot_Light %>% group_by(ID) %>% summarise(across(c(1:19),mean)) %>% distinct 

Env_tot_Light_m<-Env_tot_Light_m %>% remove_rownames %>% column_to_rownames(var="ID")

Env_tot_Light_m$`Transmitted 350 - 920nm`[is.na(Env_tot_Light_m$`Transmitted 350 - 920nm`)] <- 0

save(Env_tot_Light_m,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Env_tot_Light.Rdata"))




##################
Env_tot$Date<-as.Date(Env_tot$Date)
Env_tot<-tibble::rownames_to_column(Env_tot, "ID")
Env_tot_Light_2<-full_join(Env_tot,Light,by=c("ID"))

Env_tot_Light_2<-Env_tot_Light_2[1:21]
Env_tot_Light_2<-Env_tot_Light_2[1:733,]
Env_tot_Light_2<-Env_tot_Light_2[order(Env_tot_Light_2$Depth),]
Env_tot_Light_2<- Env_tot_Light_2 %>% relocate("Date.y")
Env_tot_Light_2<- Env_tot_Light_2[-1]

Env_tot_Light_m_1<- Env_tot_Light_2 %>% group_by(ID) %>% summarise(across(c(1:19),mean)) %>% distinct 

Env_tot_Light_m_1<-Env_tot_Light_m_1 %>% remove_rownames %>% column_to_rownames(var="ID")

Env_tot_Light_m_2<-Env_tot_Light_m_1
Env_tot_Light_m_2[is.na(Env_tot_Light_m_1)] <- 0

save(Env_tot_Light_m_2,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Env_tot_Light_2.Rdata"))

save(Env_tot_Light_m_1,file=paste0("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep","/Env_tot_Light_2_NA.Rdata"))

