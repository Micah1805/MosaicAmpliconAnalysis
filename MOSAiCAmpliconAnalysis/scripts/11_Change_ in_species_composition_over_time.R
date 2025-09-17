# Script: 11_Changes_in_species_composition_over_time.R
# Purpose: Create bar plots showing the relative abundance (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: microeco_prok.Rdata, microeco_euk.Rdata
# Main output: Figure_5.jpeg, SF_5_b.jpeg, SF_6_a.jpeg, Figure_5.jpeg, SF_5_d.jpeg, SF_6_c.jpeg
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
library(xlsx)

#Set background plots
theme_set(theme_gray())

####Prokaryotes####
#Load data
load("./output/microeco_prok.Rdata")
dataset$sample_table[c("Water", "Category")] <- 
  str_split_fixed(dataset$sample_table$Category, "_", 2) 

#Plot relative abundances at class level
#Create color scheme
col<-colorRampPalette( color=c(brewer.pal(10,"Dark2"),"Darkred","Darkblue"))(10)

test1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 10, high_level = "Phylum") 

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot abundance
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="White"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

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
       filename = "./figures/Figure_5.jpeg",
       device = "jpeg",
       width = 12,
       height = 8)

#Optional: Prepare dataset (to show class level instead of phylum for Bacteroidota)
#Filter out phylum Bacteroidota 
Bact_moved <- dataset$tax_table %>% filter(grepl("p__Bacteroidota",Phylum))

Bact_moved$Species_2 <- ""
Bact_moved<-Bact_moved[,-2]
colnames(Bact_moved)[2] <- "Phylum"
colnames(Bact_moved)[3] <- "Class"
colnames(Bact_moved)[4] <- "Order"
colnames(Bact_moved)[5] <- "Family"
colnames(Bact_moved)[6] <- "Genus"
colnames(Bact_moved)[7] <- "Species"
Bact_moved <- as.data.frame(Bact_moved)

Bact_moved_all <- dataset$tax_table %>% filter(!grepl("p__Bacteroidota",Phylum))
Bact_moved_all <- rbind(Bact_moved,Bact_moved_all)

dataset_moved <- clone(dataset)
dataset_moved$tax_table<-Bact_moved_all

#Prepare transformed microtable object and calculate new abundances
#Make the ASV and sample information consistent across all files 
dataset_moved$tidy_dataset()
print(dataset_moved)

#Check the sequence numbers in each sample
dataset_moved$sample_sums() %>% range

#Calculate the taxa abundance at each taxonomic rank
dataset_moved$cal_abund()
dataset_moved$taxa_abund$Phylum[1:5, 1:5]

#Calculate the alpha diversity
dataset_moved$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset_moved$alpha_diversity)

#Beta diversity
dataset_moved$cal_betadiv(unifrac = FALSE)

##Plot relative abundances at order level
#Create color scheme
col<-colorRampPalette( color=c(brewer.pal(20,"Dark2"),"Darkred","Darkblue"))(20)

test1 <- trans_abund$new(dataset = dataset_moved, taxrank = "Order", ntaxa = 20, high_level = "Phylum") 

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))


#Plot abundance
p<-test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="White"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#Add color code to the panels 
g_order <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g_order$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_order$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_order$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_order)

ggsave(g_order, 
       filename = "./figures/SF_5_b.jpeg",
       device = "jpeg",
       width = 12,
       height = 8)

#Plot relative abundances at genus level
#Create color scheme
col <- colorRampPalette( color=c(brewer.pal(20,"Dark2"),"Darkred","Darkblue"))(20)

test1 <- trans_abund$new(dataset = dataset_moved, taxrank = "Genus", ntaxa = 20, high_level = "Phylum") #Adjust parameters to show different taxonomic levels or increase or decrease the number of taxa shown

#Define order labels legend
test1$data_abund$Category = factor(test1$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test1$data_abund <- test1$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test1$data_abund$Leg = factor(test1$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot abundance
p <- test1$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="White"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#Add color code to the panels 
g_genus <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g_genus$layout$name))
col_strip <- c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_genus$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_genus$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_genus)

ggsave(g_genus, 
       filename = "./figures/SF_6_a.jpeg",
       device = "jpeg",
       width = 12,
       height = 8)


####Eukaryotes####
#Load data
load("./output/microeco_euk.Rdata")
dataset2$sample_table[c("Water", "Category")] <- str_split_fixed(dataset2$sample_table$Category, "_", 2) 

#Plot relative abundance at class level
#Create color scheme
col<-colorRampPalette( color=c(brewer.pal(10,"Dark2"),"Darkred","Darkblue"))(10)

test2 <- trans_abund$new(dataset = dataset2, taxrank = "Class", ntaxa = 10, high_level = "Subdivision")

#Define order labels legend
test2$data_abund$Category = factor(test2$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test2$data_abund <- test2$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test2$data_abund$Leg = factor(test2$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot abundance
p<-test2$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="White"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

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
       filename = "./figures/Figure_4.jpeg",
       device = "jpeg",
       width = 12,
       height = 8)

#Optional: Prepare dataset (to show class level instead of subdivision for Gyrista)
#Filter out Gyrista (subdivision)
Euk_moved <- dataset2$tax_table %>% filter(grepl("s__Gyrista",Subdivision))

Euk_moved$Species_2 <- ""
Euk_moved<-Euk_moved[,-2]
colnames(Euk_moved)[2] <- "Supergroup"
colnames(Euk_moved)[3] <- "Division"
colnames(Euk_moved)[4] <- "Subdivision"
colnames(Euk_moved)[5] <- "Class"
colnames(Euk_moved)[6] <- "Order"
colnames(Euk_moved)[7] <- "Family"
colnames(Euk_moved)[8] <- "Genus"
colnames(Euk_moved)[9] <- "Species"
Euk_moved <- as.data.frame(Euk_moved)

Euk_moved_all <- dataset2$tax_table %>% filter(!grepl("s__Gyrista",Subdivision))
Euk_moved_all <- rbind(Euk_moved,Euk_moved_all)

dataset_moved_euk <- clone(dataset2)
dataset_moved_euk$tax_table<-Euk_moved_all

#Make the ASV and sample information consistent across all files
dataset_moved_euk$tidy_dataset()
print(dataset_moved_euk)

#Check the sequence numbers in each sample
dataset_moved_euk$sample_sums() %>% range

#Calculate the taxa abundance at each taxonomic rank
dataset_moved_euk$cal_abund()
class(dataset_moved_euk$taxa_abund)

#Show part of the relative abundance at Division level
dataset_moved_euk$taxa_abund$Division[1:5, 1:5]

#Calculate the alpha diversity
dataset_moved_euk$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset_moved_euk$alpha_diversity)
class(dataset_moved_euk$alpha_diversity)

#Beta diversity
dataset_moved_euk$cal_betadiv(unifrac = FALSE)
class(dataset_moved_euk$beta_diversity)

dataset_moved_euk$save_betadiv(dirpath = "beta_diversity")

#Plot relative abundance at order level
#Create color scheme
col <- colorRampPalette(color=c(brewer.pal(20,"Dark2"),"Darkred","Darkblue"))(20)

test2 <- trans_abund$new(dataset = dataset_moved_euk, taxrank = "Order", ntaxa = 20, high_level = "Subdivision")

#Define order labels legend
test2$data_abund$Category = factor(test2$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test2$data_abund <- test2$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test2$data_abund$Leg = factor(test2$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))

#Plot abundance
p<-test2$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="White"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#Add color code to the panels 
g_order_euk <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g_order_euk$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_order_euk$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_order_euk$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_order_euk)

ggsave(g_order_euk, 
       filename = "./figures/SF_5_d.jpeg",
       device = "jpeg",
       width = 12,
       height = 8)

#Plot relative abundance at class level
#Create color scheme
col<-colorRampPalette( color=c(brewer.pal(20,"Dark2"),"Darkred","Darkblue"))(20)

test2 <- trans_abund$new(dataset = dataset_moved_euk, taxrank = "Genus", ntaxa = 20, high_level = "Subdivision")

#Define order labels legend
test2$data_abund$Category = factor(test2$data_abund$Category, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
test2$data_abund <- test2$data_abund %>%
  mutate(Leg= recode(Leg, "1" = "Leg 1: \n 20.09.2019- \n 13.12.2019","2"= "Leg 2: \n 13.12.2019- \n 24.02.2020","3"="Leg 3: \n 24.02.2020- \n 04.06.2020","4"="Leg 4: \n 04.06.2020- \n 12.08.2020","5"="Leg 5: \n 12.08.2020- \n 12.10.2020"))
test2$data_abund$Leg = factor(test2$data_abund$Leg, c("Leg 1: \n 20.09.2019- \n 13.12.2019","Leg 2: \n 13.12.2019- \n 24.02.2020","Leg 3: \n 24.02.2020- \n 04.06.2020","Leg 4: \n 04.06.2020- \n 12.08.2020","Leg 5: \n 12.08.2020- \n 12.10.2020"))


#Plot abundance
p<-test2$plot_bar(ggnested = TRUE,facet=c("Category","Leg"),others_color = "grey70",x_axis_name = "Date",color_values = col,  xtext_angle = 90,xtext_size = 8)+
  facet_grid(Category~Leg,scales = "free_x") +
  theme(strip.text = element_text(size=10,face="bold"))+
  ylim(0,100)+
  theme(strip.text = element_text(color ="White"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#Add color code to the panels 
g_genus_euk <- ggplot_gtable(ggplot_build(p))
striptr <- which(grepl('strip-', g$layout$name))
col_strip<-c("darkred","orange","darkgreen","darkblue","purple","grey","grey","grey","grey")
k <- 1
for (i in striptr) {
  j <- which(grepl('rect', g_genus_euk$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_genus_euk$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_genus_euk)

ggsave(g_genus_euk, 
       filename = "./figures/SF_6_c.jpeg",
       device = "jpeg",
       width = 12,
       height = 8)
