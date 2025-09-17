# Script: 12_Change_in_species_composition_removed_taxa.R
# Purpose: Create bar plots showing the relative abundance without most abundant phyla (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
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

#Set background plots
theme_set(theme_gray())

####Prokaryotes####
#Load data
load("./output/microeco_prok.Rdata")
dataset$sample_table[c("Water", "Category")] <- str_split_fixed(dataset$sample_table$Category, "_", 2) 
 
#Exclude Proteobacteria and Bacteroidota
No_proteo<-clone(dataset)
No_proteo$tax_table<-subset(No_proteo$tax_table,Phylum!="p__Proteobacteria")
No_proteo$tax_table<-subset(No_proteo$tax_table,Phylum!="p__Bacteroidota")

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
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot abundance
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=15,face="bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        legend.title =  element_text(size=15))

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
       filename = "./figures/SF_5_a.jpeg",
       device = "jpeg",
       width = 12,
       height = 10)

#Plot species composition at genus level
test1 <- trans_abund$new(dataset = No_proteo, taxrank = "Genus", ntaxa = 30, high_level = "Phylum")

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot abundance
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=15,face="bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        legend.title =  element_text(size=15))

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
       filename = "./figures/SF_6_b.jpeg",
       device = "jpeg",
       width = 14,
       height = 10)


####Eukaryotes####
load("./output/microeco_euk.Rdata")
dataset2$sample_table[c("Water", "Category")] <- str_split_fixed(dataset2$sample_table$Category, "_", 2) 
col<-colorRampPalette( color=c(brewer.pal(10,"Dark2"),"Darkred","Darkblue"))(10)

#Exclude Dinoflagellata
No_Dino<-clone(dataset2)

No_Dino$tax_table<-subset(No_Dino$tax_table,Subdivision!="s__Dinoflagellata")

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
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot relative abundance at class level
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  #theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=15,face="bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        legend.title =  element_text(size=15))

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
       filename = "./figures/SF_5_c.jpeg",
       device = "jpeg",
       width = 12,
       height = 10)

#Plot relative abundance at genus level
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  #theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="white",size=15,face="bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        legend.title =  element_text(size=15))

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
       filename = "./figures/SF_6_d.jpeg",
       device = "jpeg",
       width = 14,
       height = 10)
