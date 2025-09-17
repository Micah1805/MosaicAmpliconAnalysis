# Script: 10_Create_a_heatmap_showing_the_beta_diversity.R
# Purpose: Create a heatmap showing the beta diversity
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

####Housekeeping####
#Load packages
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gplots)

####Prokaryotes####
#Load data
load("./output/microeco_prok.Rdata")

#Create a data.frame containing the beta-diversity (Bray-Curtis dissimilarity)
Beta <- as.data.frame(dataset$beta_diversity$bray)

Beta<-tibble::rownames_to_column(Beta, "Samples")
Beta$Samples <-str_replace(Beta$Samples,"_","-")
Beta$Samples <-str_replace(Beta$Samples,"_","-")
Beta$Samples <-str_replace(Beta$Samples,"_"," ")
Beta[c("Date", "Category")] <- str_split_fixed(Beta$Samples, " ", 2) 
Beta<-Beta %>% relocate("Date") %>% relocate("Category")
Beta_m <- Beta %>% group_by(Category) %>% dplyr :: summarise(across(c(4:92),mean)) %>% distinct 
Beta_m<-as.data.frame(Beta_m) %>% remove_rownames %>% column_to_rownames(var="Category")
#str(Beta_m)

Beta_t<-as.data.frame(t(Beta_m))
Beta_t<-tibble::rownames_to_column(Beta_t, "Samples")
Beta_t$Samples <-str_replace(Beta_t$Samples,"_","-")
Beta_t$Samples <-str_replace(Beta_t$Samples,"_","-")
Beta_t$Samples <-str_replace(Beta_t$Samples,"_"," ")
Beta_t[c("Date", "Category")] <- str_split_fixed(Beta_t$Samples, " ", 2) 
Beta_t<-Beta_t %>% relocate("Date") %>% relocate("Category")
Beta_m_t <- Beta_t %>% group_by(Category) %>% dplyr :: summarise(across(c(3:6),mean)) %>% distinct 
Beta_m_t<-as.data.frame(Beta_m_t) %>% remove_rownames %>% column_to_rownames(var="Category")

Beta_heat<-as.matrix(round(t(Beta_m_t),2))

#Optional: Create a heat map (show beta-diversity per category)
heatmap(Beta_heat,
        col=heat.colors(3))
legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))

#Subdivide the diversity table into the legs of the expedition
Beta_sort<-Beta[order(Beta$Date),]
Beta_sort[,"Leg"]<-"NA"
Beta_sort$Leg[1:11][Beta_sort$Leg[1:11] == 'NA'] <- '1'
Beta_sort$Leg[12:35][Beta_sort$Leg[12:35] == 'NA'] <- '2'
Beta_sort$Leg[36:55][Beta_sort$Leg[36:55] == 'NA'] <- '3'
Beta_sort$Leg[56:76][Beta_sort$Leg[56:76] == 'NA'] <- '4'
Beta_sort$Leg[77:91][Beta_sort$Leg[77:91] == 'NA'] <- '5'
Beta_sort<-Beta_sort %>% relocate("Leg") %>% relocate("Category")
Beta_sort$Category_Leg<-paste(Beta_sort$Category,Beta_sort$Leg)
Beta_sort<-Beta_sort %>% relocate("Category_Leg")
Beta_sort$Category_Leg<-str_replace_all(Beta_sort$Category_Leg," ","_")
Beta_sort_m <- Beta_sort %>% group_by(Category_Leg) %>% dplyr :: summarise(across(c(5:95),mean)) %>% distinct 
Beta_sort_m<-as.data.frame(Beta_sort_m) %>% remove_rownames %>% column_to_rownames(var="Category_Leg")

Beta_sort_t<-as.data.frame(t(Beta_sort_m))
Beta_sort_t<-tibble::rownames_to_column(Beta_sort_t, "Samples")
Beta_sort_t$Samples <-str_replace(Beta_sort_t$Samples,"_","-")
Beta_sort_t$Samples <-str_replace(Beta_sort_t$Samples,"_","-")
Beta_sort_t$Samples <-str_replace(Beta_sort_t$Samples,"_"," ")
Beta_sort_t[c("Date", "Category")] <- str_split_fixed(Beta_sort_t$Samples, " ", 2) 
Beta_sort_t<-Beta_sort_t %>% relocate("Date") %>% relocate("Category")
Beta_sort_t<-Beta_sort_t[order(Beta_sort_t$Date),]
Beta_sort_t[,"Leg"]<-"NA"
Beta_sort_t$Leg[1:11][Beta_sort_t$Leg[1:11] == 'NA'] <- '1'
Beta_sort_t$Leg[12:36][Beta_sort_t$Leg[12:36] == 'NA'] <- '2'
Beta_sort_t$Leg[37:55][Beta_sort_t$Leg[37:55] == 'NA'] <- '3'
Beta_sort_t$Leg[56:76][Beta_sort_t$Leg[56:76] == 'NA'] <- '4'
Beta_sort_t$Leg[77:91][Beta_sort_t$Leg[77:91] == 'NA'] <- '5'
Beta_sort_t<-Beta_sort_t %>% relocate("Leg") %>% relocate("Category")
Beta_sort_t$Category_Leg<-paste(Beta_sort_t$Category,Beta_sort_t$Leg)
Beta_sort_t<-Beta_sort_t %>% relocate("Category_Leg")
Beta_sort_t$Category_Leg<-str_replace_all(Beta_sort_t$Category_Leg," ","_")
Beta_sort_t_m <- Beta_sort_t %>% group_by(Category_Leg) %>% dplyr :: summarise(across(c(5:24),mean)) %>% distinct 
Beta_sort_t_m<-as.data.frame(Beta_sort_t_m) %>% remove_rownames %>% column_to_rownames(var="Category_Leg")

Beta_heat_Leg<-as.matrix(round(t(Beta_sort_t_m),2))
Beta_heat_Leg<-tibble::rownames_to_column(as.data.frame(Beta_heat_Leg), "Samples")
Beta_heat_Leg[c("Water", "Depth_Leg")] <- str_split_fixed(Beta_heat_Leg$Samples, "_", 2) 
Beta_heat_Leg<-Beta_heat_Leg %>% relocate("Depth_Leg","Water","Samples")
Beta_heat_Leg<-as.data.frame(Beta_heat_Leg) %>% remove_rownames %>% column_to_rownames(var="Depth_Leg")
Beta_heat_Leg<-Beta_heat_Leg[3:22]
Beta_heat_Leg<-t(Beta_heat_Leg)
Beta_heat_Leg<-tibble::rownames_to_column(as.data.frame(Beta_heat_Leg), "Samples")
Beta_heat_Leg[c("Water", "Depth_Leg")] <- str_split_fixed(Beta_heat_Leg$Samples, "_", 2) 
Beta_heat_Leg<-Beta_heat_Leg %>% relocate("Depth_Leg","Water","Samples")
Beta_heat_Leg<-as.data.frame(Beta_heat_Leg) %>% remove_rownames %>% column_to_rownames(var="Depth_Leg")
Beta_heat_Leg<-Beta_heat_Leg[3:22]
Beta_heat_Leg<-t(Beta_heat_Leg)
Beta_heat_Leg<-as.matrix(Beta_heat_Leg)

#Plot beta-diversity per leg
heatmap(Beta_heat_Leg,
        col=heat.colors(3))
legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))

#Add a color code to the heatmap
pdf(file = "./figures/SF_4_a.pdf", width = 12)
my_palette <- colorRampPalette(c("red","yellow","green"))(n=299)
col_breaks=c(seq(-1,0,length=100),
             seq(0,0.8,length=100),
             seq(0.81,1,length=100))

heatmap.2(Beta_heat_Leg,
          cellnote=Beta_heat_Leg,
          notecol="black",
          density.info="none",
          trace="none",
          margins = c(9,10),
          col=my_palette,
          dendrogram="both")
dev.off()


####Eukaryotes####
load("./output/microeco_euk.Rdata")

#Create a data.frame containing the beta diversity (Bray-Curtis dissimilarity)
Beta <- as.data.frame(dataset2$beta_diversity$bray)

Beta<-tibble::rownames_to_column(Beta, "Samples")
Beta$Samples <- str_replace(Beta$Samples,"_","-")
Beta$Samples <- str_replace(Beta$Samples,"_","-")
Beta$Samples <- str_replace(Beta$Samples,"_"," ")
Beta[c("Date", "Category")] <- str_split_fixed(Beta$Samples, " ", 2) 
Beta <- Beta %>% relocate("Date") %>% relocate("Category")
Beta_m <- Beta %>% group_by(Category) %>% dplyr ::summarise(across(c(4:92),mean)) %>% distinct 
Beta_m <- as.data.frame(Beta_m) %>% remove_rownames %>% column_to_rownames(var="Category")
#str(Beta_m)

Beta_t <- as.data.frame(t(Beta_m))
Beta_t <- tibble::rownames_to_column(Beta_t, "Samples")
Beta_t$Samples <- str_replace(Beta_t$Samples,"_","-")
Beta_t$Samples <- str_replace(Beta_t$Samples,"_","-")
Beta_t$Samples <- str_replace(Beta_t$Samples,"_"," ")
Beta_t[c("Date", "Category")] <- str_split_fixed(Beta_t$Samples, " ", 2) 
Beta_t <- Beta_t %>% relocate("Date") %>% relocate("Category")
Beta_m_t <- Beta_t %>% group_by(Category) %>% dplyr :: summarise(across(c(3:6),mean)) %>% distinct 
Beta_m_t <- as.data.frame(Beta_m_t) %>% remove_rownames %>% column_to_rownames(var="Category")

Beta_heat <- as.matrix(round(t(Beta_m_t),2))

#Create a heat map (show beta-diversity per category)
heatmap(Beta_heat,
        #Rowv = NA, 
        #Colv = NA,
        col=heat.colors(3))
legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))

#Subdivide the diversity table into the legs of the expedition
Beta_sort<-Beta[order(Beta$Date),]
Beta_sort[,"Leg"]<-"NA"
Beta_sort$Leg[1:11][Beta_sort$Leg[1:11] == 'NA'] <- '1'
Beta_sort$Leg[12:35][Beta_sort$Leg[12:35] == 'NA'] <- '2'
Beta_sort$Leg[36:55][Beta_sort$Leg[36:55] == 'NA'] <- '3'
Beta_sort$Leg[56:76][Beta_sort$Leg[56:76] == 'NA'] <- '4'
Beta_sort$Leg[77:91][Beta_sort$Leg[77:91] == 'NA'] <- '5'
Beta_sort<-Beta_sort %>% relocate("Leg") %>% relocate("Category")
Beta_sort$Category_Leg<-paste(Beta_sort$Category,Beta_sort$Leg)
Beta_sort<-Beta_sort %>% relocate("Category_Leg")
Beta_sort$Category_Leg<-str_replace_all(Beta_sort$Category_Leg," ","_")
Beta_sort_m <- Beta_sort %>% group_by(Category_Leg) %>% dplyr :: summarise(across(c(5:95),mean)) %>% distinct 
Beta_sort_m<-as.data.frame(Beta_sort_m) %>% remove_rownames %>% column_to_rownames(var="Category_Leg")

Beta_sort_t<-as.data.frame(t(Beta_sort_m))
Beta_sort_t<-tibble::rownames_to_column(Beta_sort_t, "Samples")
Beta_sort_t$Samples <-str_replace(Beta_sort_t$Samples,"_","-")
Beta_sort_t$Samples <-str_replace(Beta_sort_t$Samples,"_","-")
Beta_sort_t$Samples <-str_replace(Beta_sort_t$Samples,"_"," ")
Beta_sort_t[c("Date", "Category")] <- str_split_fixed(Beta_sort_t$Samples, " ", 2) 
Beta_sort_t<-Beta_sort_t %>% relocate("Date") %>% relocate("Category")
Beta_sort_t<-Beta_sort_t[order(Beta_sort_t$Date),]
Beta_sort_t[,"Leg"]<-"NA"
Beta_sort_t$Leg[1:11][Beta_sort_t$Leg[1:11] == 'NA'] <- '1'
Beta_sort_t$Leg[12:36][Beta_sort_t$Leg[12:36] == 'NA'] <- '2'
Beta_sort_t$Leg[37:55][Beta_sort_t$Leg[37:55] == 'NA'] <- '3'
Beta_sort_t$Leg[56:76][Beta_sort_t$Leg[56:76] == 'NA'] <- '4'
Beta_sort_t$Leg[77:91][Beta_sort_t$Leg[77:91] == 'NA'] <- '5'
Beta_sort_t<-Beta_sort_t %>% relocate("Leg") %>% relocate("Category")
Beta_sort_t$Category_Leg<-paste(Beta_sort_t$Category,Beta_sort_t$Leg)
Beta_sort_t<-Beta_sort_t %>% relocate("Category_Leg")
Beta_sort_t$Category_Leg<-str_replace_all(Beta_sort_t$Category_Leg," ","_")
Beta_sort_t_m <- Beta_sort_t %>% group_by(Category_Leg) %>% dplyr :: summarise(across(c(5:24),mean)) %>% distinct 
Beta_sort_t_m<-as.data.frame(Beta_sort_t_m) %>% remove_rownames %>% column_to_rownames(var="Category_Leg")

Beta_heat_Leg<-as.matrix(round(t(Beta_sort_t_m),2))
Beta_heat_Leg<-tibble::rownames_to_column(as.data.frame(Beta_heat_Leg), "Samples")
Beta_heat_Leg[c("Water", "Depth_Leg")] <- str_split_fixed(Beta_heat_Leg$Samples, "_", 2) 
Beta_heat_Leg<-Beta_heat_Leg %>% relocate("Depth_Leg","Water","Samples")
Beta_heat_Leg<-as.data.frame(Beta_heat_Leg) %>% remove_rownames %>% column_to_rownames(var="Depth_Leg")
Beta_heat_Leg<-Beta_heat_Leg[3:22]
Beta_heat_Leg<-t(Beta_heat_Leg)
Beta_heat_Leg<-tibble::rownames_to_column(as.data.frame(Beta_heat_Leg), "Samples")
Beta_heat_Leg[c("Water", "Depth_Leg")] <- str_split_fixed(Beta_heat_Leg$Samples, "_", 2) 
Beta_heat_Leg<-Beta_heat_Leg %>% relocate("Depth_Leg","Water","Samples")
Beta_heat_Leg<-as.data.frame(Beta_heat_Leg) %>% remove_rownames %>% column_to_rownames(var="Depth_Leg")
Beta_heat_Leg<-Beta_heat_Leg[3:22]
Beta_heat_Leg<-t(Beta_heat_Leg)
Beta_heat_Leg<-as.matrix(Beta_heat_Leg)

#Add a color code to the heatmap
pdf(file = "./figures/SF_4_b.pdf", width = 12)
my_palette <- colorRampPalette(c("red","yellow","green"))(n=299)
col_breaks=c(seq(-1,0,length=100),
             seq(0,0.8,length=100),
             seq(0.81,1,length=100))

Beta_div_heatmap <- heatmap.2(Beta_heat_Leg,
          cellnote=Beta_heat_Leg,
          notecol="black",
          density.info="none",
          trace="none",
          margins = c(8,10),
          col=my_palette,
          dendrogram="both")
dev.off()

