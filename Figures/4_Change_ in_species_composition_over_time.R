################################################################################
##################Changes in species composition over time######################
################################################################################

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
load("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/dataset.Rdata")
dataset$sample_table[c("Water", "Category")] <- 
  str_split_fixed(dataset$sample_table$Category, "_", 2) 

##Optional: Prepare dataset (to show class level instead of phylum for Bacteroidota)##
#Filter out phylum Bacteroidota 
Bact_moved <- dataset$tax_table %>% filter(grepl("p__Bacteroidota",Phylum))


Bact_moved$Species_2<-""
Bact_moved<-Bact_moved[,-2]
colnames(Bact_moved)[2] <-"Phylum"
colnames(Bact_moved)[3] <-"Class"
colnames(Bact_moved)[4] <-"Order"
colnames(Bact_moved)[5] <-"Family"
colnames(Bact_moved)[6] <-"Genus"
colnames(Bact_moved)[7] <-"Species"
Bact_moved<-as.data.frame(Bact_moved)

Bact_moved_all <- dataset$tax_table %>% filter(!grepl("p__Bacteroidota",Phylum))
Bact_moved_all <- rbind(Bact_moved,Bact_moved_all)

#write.xlsx(as.data.frame(Bact_moved_all), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Bact_moved_all.xlsx",
#          col.names = TRUE, row.names = TRUE)

dataset$tax_table<-Bact_moved_all

#Prepare transformed microecoobject and calculate new abundances
#Make the OTU and sample information consistent across all files in the dataset object
dataset$tidy_dataset()
print(dataset)

#Check the sequence numbers in each sample
dataset$sample_sums() %>% range

#Save all the basic data in microtable object to local files, including feature abundance, metadata, taxonomic table, phylogenetic tree and representative sequences

#Calculate the taxa abundance at each taxonomic rank
dataset$cal_abund()
class(dataset$taxa_abund)

dataset$taxa_abund$Division[1:5, 1:5]

dataset$save_abund(dirpath = "taxa_abund")

#Calculate the alpha diversity
dataset$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset$alpha_diversity)
class(dataset$alpha_diversity)

#Beta diversity
dataset$cal_betadiv(unifrac = FALSE)
class(dataset$beta_diversity)



##Plot relative abundances##
#Create color scheme
col<-colorRampPalette( color=c(brewer.pal(20,"Dark2"),"Darkred","Darkblue"))(10)

test1 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 20, high_level = "Phylum")

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


####Eukaryotes####

#Load data
load("~/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/dataset2.Rdata")
dataset2$sample_table[c("Water", "Category")] <- str_split_fixed(dataset2$sample_table$Category, "_", 2) 

##Optional: Prepare dataset (to show class level instead of subdivision for Gyrista)## 
#Filter out Gyrista (subdivision)
Euk_moved <- dataset2$tax_table %>% filter(grepl("s__Gyrista",Subdivision))

Euk_moved$Species_2<-""
Euk_moved<-Euk_moved[,-2]
colnames(Euk_moved)[2] <-"Supergroup"
colnames(Euk_moved)[3] <-"Division"
colnames(Euk_moved)[4] <-"Subdivision"
colnames(Euk_moved)[5] <-"Class"
colnames(Euk_moved)[6] <-"Order"
colnames(Euk_moved)[7] <-"Family"
colnames(Euk_moved)[8] <-"Genus"
colnames(Euk_moved)[9] <-"Species"
Euk_moved<-as.data.frame(Euk_moved)

Euk_moved_all <- dataset2$tax_table %>% filter(!grepl("s__Gyrista",Subdivision))
Euk_moved_all <- rbind(Euk_moved,Euk_moved_all)

#write.xlsx(as.data.frame(Bact_moved_all), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Bact_moved_all.xlsx",
#          col.names = TRUE, row.names = TRUE)

dataset2$tax_table<-Euk_moved_all

#Make the OTU and sample information consistent across all files in the dataset object
dataset2$tidy_dataset()
print(dataset2)

#Check the sequence numbers in each sample
dataset2$sample_sums() %>% range

#Calculate the taxa abundance at each taxonomic rank
dataset2$cal_abund()
class(dataset2$taxa_abund)

#Show part of the relative abundance at Division level
dataset2$taxa_abund$Division[1:5, 1:5]

#Calculate the alpha diversity
dataset2$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(dataset2$alpha_diversity)
class(dataset2$alpha_diversity)

#Beta diversity
dataset2$cal_betadiv(unifrac = FALSE)
class(dataset2$beta_diversity)

dataset2$save_betadiv(dirpath = "beta_diversity")



##Plot data##
#Create color scheme
col<-colorRampPalette( color=c(brewer.pal(20,"Dark2"),"Darkred","Darkblue"))(10)

test2 <- trans_abund$new(dataset = dataset2, taxrank = "Order", ntaxa = 20, high_level = "Subdivision")

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
