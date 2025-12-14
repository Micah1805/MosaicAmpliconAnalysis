# Script: 8_Change_in_species_composition_removed_taxa.R
# Purpose: Create bar plots showing the relative abundance without most abundant phyla (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: microeco_prok.Rdata, microeco_euk.Rdata
# Main output: SF_2_a.jpeg, SF_3_b.jpeg, SF_2_c.jpeg, SF_3_d.jpeg
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

####Housekeeping####
#Load packages
library(microeco)
library(dplyr)
library(ggplot2)
library(magrittr)
library(paletteer)
library(grid)
library(RColorBrewer)
library(stringr)
library(tidyr)

#Set background plots
theme_set(theme_gray())

####Prokaryotes####
#Load data
load("./output/microeco_prok.Rdata")
dataset$sample_table$Date <- as.Date(dataset$sample_table$Date)

dataset$sample_table[c("Water", "Category")] <- str_split_fixed(dataset$sample_table$Category, "_", 2) 
dataset$sample_table$Category<-str_replace(dataset$sample_table$Category,"_","-")

#Create date lables for the x-axis
dataset$sample_table[c("Year", "Date_Month")] <- str_split_fixed(dataset$sample_table$Date, "-", 2)
dataset$sample_table[c("Month", "Day")] <- str_split_fixed(dataset$sample_table$Date_Month, "-", 2)
dataset$sample_table$Month <- as.numeric(dataset$sample_table$Month)
dataset$sample_table$Month_name <- month.abb[c(dataset$sample_table$Month)]
dataset$sample_table$Date_Month <- factor(dataset$sample_table$Date_Month, 
                                          levels =  unique(dataset$sample_table$Date_Month))
dataset$sample_table$Day <- factor(dataset$sample_table$Day, 
                                          levels =  unique(dataset$sample_table$Day))
dataset$sample_table <- dataset$sample_table[order(dataset$sample_table$Date,decreasing=FALSE),]
dataset$sample_table <- dataset$sample_table %>%
  mutate(
    M = ifelse(lag(Month_name) != Month_name |
                 is.na(lag(Month_name)), Month_name, ""),
    Y = ifelse(lag(Year) != Year |
                 is.na(lag(Year)), Year, "")
  )
dataset$sample_table$Day <- factor(dataset$sample_table$Day,
                                   levels =  unique(dataset$sample_table$Day))

dataset$sample_table <- dataset$sample_table %>%
  unite("Day_Month", Day, Month_name, sep = " ")
dataset$sample_table$Day_Month <- factor(dataset$sample_table$Day_Month,
                                   levels =  unique(dataset$sample_table$Day_Month))

#Remove "p__" in front of Phyla names
dataset$tax_table$Phylum <- gsub("p__","",dataset$tax_table$Phylum) 


#Exclude Proteobacteria and Bacteroidota
No_proteo<-clone(dataset)
No_proteo$tax_table<-subset(No_proteo$tax_table,Phylum!="Proteobacteria")
No_proteo$tax_table<-subset(No_proteo$tax_table,Phylum!="Bacteroidota")

#Remove ASVs with the taxonomic assignments “mitochondria” or “chloroplast”
No_proteo$filter_pollution(taxa = c("mitochondria", "chloroplast"))
No_proteo

#Make the ASV and sample information consistent across all files
No_proteo$tidy_dataset()
print(No_proteo)

#Check the sequence numbers in each sample
No_proteo$sample_sums() %>% range

#Calculate the taxa abundance at each taxonomic rank
No_proteo$cal_abund()

#Show part of the relative abundance at Phylum level
No_proteo$taxa_abund$Phylum[1:5, 1:5]

#Calculate the alpha diversity
No_proteo$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(No_proteo$alpha_diversity)

#Beta diversity
No_proteo$cal_betadiv(unifrac = FALSE)

#Plot species composition at class level
test1 <- trans_abund$new(dataset = No_proteo, taxrank = "Class", ntaxa = 10, high_level = "Phylum")

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax","50m","100m"))
test1$data_abund$Day_Month <- factor(test1$data_abund$Day_Month, levels = c("14 Nov", "21 Nov", "29 Nov", "26 Dec",
                                                                            "02 Jan", "09 Jan", "16 Jan", "23 Jan",
                                                                            "30 Jan", "06 Feb", "20 Feb", "05 Mar",
                                                                            "07 Mar", "12 Mar", "27 Mar", "04 Apr",
                                                                            "09 Apr", "23 Apr", "07 May", "16 May",
                                                                            "16 Jun", "27 Jun", "01 Jul", "09 Jul",
                                                                            "16 Jul", "23 Jul", "30 Jul", "27 Aug",
                                                                            "03 Sep", "11 Sep", "17 Sep"))


#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1","2"= "Leg 2","3"="Leg 3","4"="Leg 4","5"="Leg 5"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1","Leg 2","Leg 3","Leg 4","Leg 5"))

#Plot abundance
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Day_Month",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=24,face="bold"),
        axis.text.x = element_text(size = 24,angle = 90, margin = margin(t = 5), vjust = 0.5),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, 'cm'),
        legend. = element_text(size = 22),
        legend.title =  element_text(size=24))

#Add color code to the panels 
g <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g)

ggsave(g, 
       filename = "./figures/SF_2_a.jpeg",
       device = "jpeg",
       width = 12,
       height = 10)

#Plot species composition at genus level
test1 <- trans_abund$new(dataset = No_proteo, taxrank = "Genus", ntaxa = 30, high_level = "Phylum")

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax","50m","100m"))
test1$data_abund$Day_Month <- factor(test1$data_abund$Day_Month, levels = c("14 Nov", "21 Nov", "29 Nov", "26 Dec",
                                                                            "02 Jan", "09 Jan", "16 Jan", "23 Jan",
                                                                            "30 Jan", "06 Feb", "20 Feb", "05 Mar",
                                                                            "07 Mar", "12 Mar", "27 Mar", "04 Apr",
                                                                            "09 Apr", "23 Apr", "07 May", "16 May",
                                                                            "16 Jun", "27 Jun", "01 Jul", "09 Jul",
                                                                            "16 Jul", "23 Jul", "30 Jul", "27 Aug",
                                                                            "03 Sep", "11 Sep", "17 Sep"))


#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1","2"= "Leg 2","3"="Leg 3","4"="Leg 4","5"="Leg 5"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1","Leg 2","Leg 3","Leg 4","Leg 5"))

#Plot abundance
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=24,face="bold"),
        axis.text.x = element_text(size = 24,angle = 90, margin = margin(t = 5), vjust = 0.5),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, 'cm'),
        legend. = element_text(size = 22),
        legend.title =  element_text(size=24))

#Add color code to the panels 
g_genus <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g_genus$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_genus$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_genus$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_genus)

ggsave(g_genus, 
       filename = "./figures/SF_3_b.jpeg",
       device = "jpeg",
       width = 14,
       height = 10)


####Eukaryotes####
load("./output/microeco_euk.Rdata")
dataset2$sample_table$Date <- as.Date(dataset2$sample_table$Date)

dataset2$sample_table[c("Water", "Category")] <- str_split_fixed(dataset2$sample_table$Category, "_", 2) 
dataset2$sample_table$Category<-str_replace(dataset2$sample_table$Category,"_","-")

#Create date lables for the x-axis
dataset2$sample_table[c("Year", "Date_Month")] <- str_split_fixed(dataset2$sample_table$Date, "-", 2)
dataset2$sample_table[c("Month", "Day")] <- str_split_fixed(dataset2$sample_table$Date_Month, "-", 2)
dataset2$sample_table$Month <- as.numeric(dataset2$sample_table$Month)
dataset2$sample_table$Month_name <- month.abb[c(dataset2$sample_table$Month)]
dataset2$sample_table$Date_Month <- factor(dataset2$sample_table$Date_Month, 
                                           levels =  unique(dataset2$sample_table$Date_Month))
dataset2$sample_table$Day <- factor(dataset2$sample_table$Day, 
                                    levels =  unique(dataset2$sample_table$Day))
dataset2$sample_table <- dataset2$sample_table[order(dataset2$sample_table$Date,decreasing=FALSE),]
dataset2$sample_table <- dataset2$sample_table %>%
  mutate(
    M = ifelse(lag(Month_name) != Month_name |
                 is.na(lag(Month_name)), Month_name, ""),
    Y = ifelse(lag(Year) != Year |
                 is.na(lag(Year)), Year, "")
  )
dataset2$sample_table$Day <- factor(dataset2$sample_table$Day,
                                    levels =  unique(dataset2$sample_table$Day))

dataset2$sample_table <- dataset2$sample_table %>%
  unite("Day_Month", Day, Month_name, sep = " ")
dataset2$sample_table$Day_Month <- factor(dataset2$sample_table$Day_Month,
                                          levels =  unique(dataset2$sample_table$Day_Month))

#Remove "s__" in front of subdivision names
dataset2$tax_table$Subdivision <- gsub("s__","",dataset2$tax_table$Subdivision)


col<-colorRampPalette( color=c(brewer.pal(10,"Dark2"),"Darkred","Darkblue"))(10)

#Exclude Dinoflagellata
No_Dino<-clone(dataset2)

No_Dino$tax_table<-subset(No_Dino$tax_table,Subdivision!="Dinoflagellata")

#Remove ASVs with the taxonomic assignments “mitochondria” or “chloroplast”
No_Dino$filter_pollution(taxa = c("mitochondria", "chloroplast"))
No_Dino

#Make the ASV and sample information consistent across all files
No_Dino$tidy_dataset()
print(No_Dino)

#Check the sequence numbers in each sample
No_Dino$sample_sums() %>% range

#Calculate the taxa abundance at each taxonomic rank
No_Dino$cal_abund()

#Show part of the relative abundance at Phylum level
No_Dino$taxa_abund$Division[1:5, 1:5]

#Calculate the alpha diversity
No_Dino$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(No_Dino$alpha_diversity)

#Beta diversity
No_Dino$cal_betadiv(unifrac = FALSE)

#Plot species composition
test1 <- trans_abund$new(dataset = No_Dino, taxrank = "Class", ntaxa = 10, high_level = "Subdivision")

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax","50m","100m"))
test1$data_abund$Day_Month <- factor(test1$data_abund$Day_Month, levels = c("14 Nov", "21 Nov", "29 Nov", "26 Dec",
                                                                            "02 Jan", "09 Jan", "16 Jan", "23 Jan",
                                                                            "30 Jan", "06 Feb", "20 Feb", "05 Mar",
                                                                            "07 Mar", "12 Mar", "27 Mar", "04 Apr",
                                                                            "09 Apr", "23 Apr", "07 May", "16 May",
                                                                            "16 Jun", "27 Jun", "01 Jul", "09 Jul",
                                                                            "16 Jul", "23 Jul", "30 Jul", "27 Aug",
                                                                            "03 Sep", "11 Sep", "17 Sep"))


#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1","2"= "Leg 2","3"="Leg 3","4"="Leg 4","5"="Leg 5"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1","Leg 2","Leg 3","Leg 4","Leg 5"))

#Plot relative abundance at class level
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Day_Month",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  #theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=24,face="bold"),
        axis.text.x = element_text(size = 24,angle = 90, margin = margin(t = 5), vjust = 0.5),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, 'cm'),
        legend. = element_text(size = 22),
        legend.title =  element_text(size=24))

#Add color code to the panels 
g_euk <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g_euk$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_euk$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_euk$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_euk)

ggsave(g_euk, 
       filename = "./figures/SF_2_c.jpeg",
       device = "jpeg",
       width = 12,
       height = 10)

#Plot species composition at  genus level
test1 <- trans_abund$new(dataset = No_Dino, taxrank = "Genus", ntaxa = 30, high_level = "Subdivision")

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax","50m","100m"))
test1$data_abund$Day_Month <- factor(test1$data_abund$Day_Month, levels = c("14 Nov", "21 Nov", "29 Nov", "26 Dec",
                                                                            "02 Jan", "09 Jan", "16 Jan", "23 Jan",
                                                                            "30 Jan", "06 Feb", "20 Feb", "05 Mar",
                                                                            "07 Mar", "12 Mar", "27 Mar", "04 Apr",
                                                                            "09 Apr", "23 Apr", "07 May", "16 May",
                                                                            "16 Jun", "27 Jun", "01 Jul", "09 Jul",
                                                                            "16 Jul", "23 Jul", "30 Jul", "27 Aug",
                                                                            "03 Sep", "11 Sep", "17 Sep"))


#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1","2"= "Leg 2","3"="Leg 3","4"="Leg 4","5"="Leg 5"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1","Leg 2","Leg 3","Leg 4","Leg 5"))

#Plot relative abundance at genus level
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Day_Month",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  #theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=24,face="bold"),
        axis.text.x = element_text(size = 24,angle = 90, margin = margin(t = 5), vjust = 0.5),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, 'cm'),
        legend. = element_text(size = 22),
        legend.title =  element_text(size=24))

#Add color code to the panels 
g_genus_euk <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g_genus_euk$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_genus_euk$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_genus_euk$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g)

ggsave(g_genus, 
       filename = "./figures/SF_3_d.jpeg",
       device = "jpeg",
       width = 14,
       height = 10)

