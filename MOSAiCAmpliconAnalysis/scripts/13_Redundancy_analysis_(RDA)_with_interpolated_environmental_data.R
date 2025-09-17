# Script: 13_Redundancy_analysis_(RDA)_with_interpolated_environmental_data.R
# Purpose: Create RDA plots (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

####Housekeeping####
#Load packages
library(zoo)
library(ggplot2)
library(microeco)
library(dplyr)
library(magrittr)
library(stringr)

####Interpolate environmental data####
#Load environmental data
environmental_data_light_RDA <- read.csv("./output/environmental_data_light_RDA.csv", row.names=1)

#Interpolate
Env_inter <- na.approx(environmental_data_light_RDA[9:19],rule=2)
Inter <- cbind(environmental_data_light_RDA[1:8],as.data.frame(Env_inter))
df <- as.data.frame(Env_inter)
Inter[is.na(Inter)] <- 0

#Change labels of columns in the environmental data table
str(Inter)
Inter$Date <- Inter$Date %>% as.Date
colnames(Inter)[1] <- "Time"
colnames(Inter)[5] <- "Depth [m]"
colnames(Inter)[7] <- "Salinity"
colnames(Inter)[8] <- "Temperature [°C]"
colnames(Inter)[9] <- "PO4 [µmol/l]"
colnames(Inter)[10] <- "NO2 [µmol/l]"
colnames(Inter)[11] <- "N3N2 [µmol/l]"
colnames(Inter)[12] <- "Si [µmol/l]"
colnames(Inter)[13] <- "NH4 [µmol/l]"
colnames(Inter)[17] <- "Chl a [µg/l]"
colnames(Inter)[19] <- "Light (350-920 nm (transmission))"
#str(Inter)

#Plot interpolated data
ggplot(as.data.frame(Inter),aes(Time,`Light (350-920 nm (transmission))`))+
  geom_point()+
  geom_line()+
  geom_point(data=environmental_data_light_RDA,aes(Date.x,Transmitted.350...920nm),color="red")


ggplot(as.data.frame(Inter),aes(Time,`NO2 [µmol/l]`))+
  geom_point()+
  geom_line()+
  geom_point(data=environmental_data_light_RDA,aes(Date.x,NO2_mean),color="red")

ggplot(as.data.frame(Inter),aes(Time,`PO4 [µmol/l]`))+
  geom_point()+
  geom_line()+
  geom_point(data=environmental_data_light_RDA,aes(Date.x,PO4_mean),color="red")


ggplot(as.data.frame(Inter),aes(Time,`Si [µmol/l]`))+
  geom_point()+
  geom_line()+
  geom_point(data=environmental_data_light_RDA,aes(Date.x,Si_mean),color="red")


ggplot(as.data.frame(Inter),aes(Time,`NH4 [µmol/l]`))+
  geom_point()+
  geom_line()+
  geom_point(data=environmental_data_light_RDA,aes(Date.x,NH4_mean),color="red")

#Move to right directory
path="./output"


####Prokaryotes####
load("./output/microeco_prok.Rdata")
dataset$sample_table[c("Water", "Water_depth")] <- str_split_fixed(dataset$sample_table$Category, "_", 2) 
dataset$sample_table$Water_depth %<>% factor(., levels = c("Surface","Chlmax",">50m",">75m"))

#Add the environmental data
t1 <- trans_env$new(dataset = dataset, add_data = Inter[c(1,3:5,7:13,17,19)])

#Calculate RDA using taxonomic level
t1$cal_ordination(method = "RDA", taxa_level = "Class")

#Select 10 features and adjust the arrow length
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
t1$res_ordination_trans$df_sites$Water_depth = factor(t1$res_ordination_trans$df_sites$Water_depth, c("Surface","Chlmax",">50m",">75m"))

#Plot RDA
RDA_prok <- t1$plot_ordination(plot_color = "Water_depth")

ggsave(RDA_prok, 
       filename = "./figures/Figure_6_a.jpeg",
       device = "jpeg",
       width = 15,
       height = 15)


####Eukaryotes####
load("./output/microeco_euk.Rdata")
dataset2$sample_table[c("Water", "Water_depth")] <- str_split_fixed(dataset2$sample_table$Category, "_", 2) 

dataset2$sample_table$Water_depth %<>% factor(., levels = c("Surface","Chlmax",">50m",">75m"))

#Define order of the sample lables in the legend
dataset2$sample_table$Water_depth %<>% factor(., levels = c("Surface","Chlmax",">50m",">75m"))

#Add the environmental data
t2 <- trans_env$new(dataset = dataset2, add_data = Inter[c(1,3:5,7:13,17,19)])

#Calculate RDA using taxonomic level
t2$cal_ordination(method = "RDA", taxa_level = "Class")

#Select 10 features and adjust the arrow length
t2$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
t2$res_ordination_trans$df_sites$Water_depth = factor(t2$res_ordination_trans$df_sites$Water_depth, c("Surface","Chlmax",">50m",">75m"))

#Plot RDA
RDA_euk <- t2$plot_ordination(plot_color = "Water_depth")

ggsave(RDA_euk, 
       filename = "./figures/Figure_6_b.jpeg",
       device = "jpeg",
       width = 15,
       height = 15)


####Prokaryotes und eukaryotes####
#Load dataset
load("./output/microeco_prok_euk.Rdata")

#Remove :nucl from Eukaryota:nucl
dataset3$tax_table$Kingdom<-as.character(stringr::str_remove(dataset3$tax_table$Kingdom, ":nucl"))

#Define order of the sample lables in the legend
dataset3$sample_table$Category %<>% factor(., levels = c("Water_Surface","Water_Chlmax","Water_>50m","Water_>75m"))

#Add the environmental data
t3 <- trans_env$new(dataset = dataset3, add_data = Inter[c(1,3:5,7:13,17,19)])

#Calculate RDA using taxonomic level
t3$cal_ordination(method = "RDA", taxa_level = "Class")

#Select 10 features and adjust the arrow length
t3$trans_ordination(show_taxa = 30, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

#Assign each Class to its Kingdom
t3$res_ordination_trans$df_arrows_spe <- tibble::rownames_to_column(t3$res_ordination_trans$df_arrows_spe, "ID")
Prok_class <- t3$res_ordination_trans$df_arrows_spe[grepl("prok_asv*",t3$res_ordination_trans$df_arrows_spe$ID),]
Prok_class$Kingdom<-c("Prokaryotes")
Euk_class <- t3$res_ordination_trans$df_arrows_spe[grepl("euk_asv*",t3$res_ordination_trans$df_arrows_spe$ID),]
Euk_class$Kingdom<-c("Eukaryotes")
Classes <- rbind(Prok_class,Euk_class)

#Replace df_arrows_spe by own variable
t3$res_ordination_trans$df_arrows_spe = Classes

#Plot RDA
RDA_prok_euk <- t3$plot_ordination(plot_color = "Category", taxa_arrow_color=ifelse(t3$res_ordination_trans$df_arrows_spe$Kingdom == "Prokaryotes","red","blue"),taxa_text_color=ifelse(t3$res_ordination_trans$df_arrows_spe$Kingdom == "Prokaryotes","red","blue"),
                   env_text_size = 6,
                   taxa_text_size = 5)+
                   theme(axis.title.x = element_text(size=15),
                   axis.title.y = element_text(size=15),
                   legend.text = element_text(size=15),
                   legend.title = element_text(size=15))

ggsave(RDA_prok_euk, 
       filename = "./figures/Figure_6_c.jpeg",
       device = "jpeg",
       width = 25,
       height = 25)
