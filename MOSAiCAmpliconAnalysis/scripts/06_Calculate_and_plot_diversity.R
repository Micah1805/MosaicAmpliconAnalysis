# Script: 06_Calculate_and_plot_diversity.R
# Purpose: Create plots showing alpha diversity and  PCoA (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: microeco_prok.Rdata, microeco_euk.Rdata
# Main output: ST_1.csv, ST_2.csv, Figure_3.jpeg 
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

####Housekeeping####
#Load packages
library(microeco)
library(ggplot2)
library(tidyverse)
library(grid)
library(xlsx)
library(cowplot)
library(ggpubr)

#Set background
theme_set(theme_bw())

set.seed(123)

####Prokaryotes####
load("./output/microeco_prok.Rdata")
dataset$sample_table[c("Water", "Water_depths")] <- str_split_fixed(dataset$sample_table$Category, "_", 2) 
dataset$sample_table$Water_depths<-str_replace(dataset$sample_table$Water_depths,"_","-")

##Alpha diversity
t1 <- trans_alpha$new(dataset = dataset, group = "Water_depths")
head(t1$data_stat)

#Calculate difference between the groups
#ANOVA
t1$cal_diff(method = "anova")
head(t1$res_diff)

t1$cal_diff(method = "anova")

#Set order labels legend
t1$data_alpha$`Water_depths` = factor(t1$data_alpha$`Water_depths` , c("Surface","Chlmax","50m","100m"))

#Plot alpha-diversity
alpha_div_proks <- t1$plot_alpha(measure = "Shannon", add_sig_text_size = 6, add = "jitter", xtext_size = 15)+
  scale_y_continuous(breaks = seq(0, 7, 0.5))

#Beta-diversity
t1 <- trans_beta$new(dataset = dataset, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`, c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "1" = "Leg 1","2"= "Leg 2","3"="Leg 3","4"="Leg 4","5"="Leg 5"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 1","Leg 2","Leg 3","Leg 4","Leg 5"))

#Plot the PCoA result with confidence ellipse
#PCoA_proks <- t1$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point","ellipse"), ellipse_chull_fill = FALSE)+
#  theme(axis.title = element_text(size=15,face="bold"),
#        axis.text = element_text(size=15,face = "bold"),
#        legend.text = element_text(size=15),
#        legend.title = element_text(size=15))+
#  guides(fill=guide_legend(title="Water_depths"))


#Separated by leg
#p<-t1$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"), ellipse_chull_fill = FALSE)+
#  facet_wrap(~Leg,ncol=1)+
#  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
#  theme(strip.text = element_text(color ="White",size = 12),
#        axis.title = element_text(size=12),
#        axis.text = element_text(size=12),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        legend.position="none")

#Add color code to the panels 
#g_prok <- ggplot_gtable(ggplot_build(p))
#stript <- which(grepl('strip-t', g_prok$layout$name))
#col_strip<-c("purple","darkblue","darkgreen","orange","darkred")
#k <- 1
#for (i in stript) {
#  j <- which(grepl('rect', g_prok$grobs[[i]]$grobs[[1]]$childrenOrder))
#  g_prok$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
#  k <- k+1
#}
#grid.draw(g_prok)


##PERMANOVA
#Separated by legs
t1$cal_manova(manova_all = FALSE,group = "Water_depths", by_group = "Leg")
PERMANOVA_Leg_Prok<-t1$res_manova

#Export results
#write.csv(as.data.frame(PERMANOVA_Prok), file="./output/PERMANOVA_proks_no_legs.csv", row.names = TRUE)
write.csv(as.data.frame(PERMANOVA_Leg_Prok), file="./output/ST_1.csv", row.names = TRUE)


####Eukaryotes####
#Load data
load("./output/microeco_euk.Rdata")
dataset2$sample_table[c("Water", "Water_depths")] <- str_split_fixed(dataset2$sample_table$Category, "_", 2) 
dataset2$sample_table$Water_depths<-str_replace(dataset2$sample_table$Water_depths,"_","-")

##Alpha diversity
#Create trans_alpha object
t2 <- trans_alpha$new(dataset = dataset2, group = "Water_depths")
head(t2$data_stat)

#Calculate difference between the groups
#ANOVA
t2$cal_diff(method = "anova")
# return t1$res_diff
head(t2$res_diff)

#Set order lables legend
t2$data_alpha$`Water_depths` = factor(t2$data_alpha$`Water_depths`, c("Surface","Chlmax","50m","100m"))


#Plot alpha-diversity
alpha_div_euks <- t2$plot_alpha(measure = "Shannon", add_sig_text_size = 6, add = "jitter")+
  scale_y_continuous(breaks = seq(0, 7, 0.5))+
  coord_cartesian(ylim = c(4,7))

##Beta-diversity
t2 <- trans_beta$new(dataset = dataset2, group = "Water_depths", measure = "bray")

#Calculate PCoA
t2$cal_ordination(method = "PCoA")
class(t2$res_ordination)

#Define order labels legend
t2$res_ordination$scores$`Water_depths` = factor(t2$res_ordination$scores$`Water_depths`, c("Surface","Chlmax","50m","100m"))

#Rename Legs
t2$res_ordination$scores <- t2$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "1" = "Leg 1","2"= "Leg 2","3"="Leg 3","4"="Leg 4","5"="Leg 5"))
t2$res_ordination$scores$Leg = factor(t2$res_ordination$scores$Leg, c("Leg 1","Leg 2","Leg 3","Leg 40","Leg 5"))

#Plot the PCoA result with confidence ellipse
#PCoA_euks <- t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point","ellipse"),, ellipse_chull_fill = FALSE)+
#  theme(axis.title = element_text(size=15,face="bold"),
#        axis.text = element_text(size=15,face = "bold"),
#        legend.text = element_text(size=15),
#        legend.title = element_text(size=15))

#Separated by leg
#p2<-t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"), ellipse_chull_fill = FALSE)+
#  facet_wrap(~Leg,ncol=1)+
#  scale_y_continuous(breaks = seq(-0.4, 0.5, 0.1))+
#  theme(strip.text = element_text(color ="White"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        legend.text = element_text(size = 12),
#        legend.title = element_text())

#p3<-get_legend(p2)

#Separated by leg
#p2<-t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"),ellipse_chull_fill = FALSE)+
#  facet_wrap(~Leg,ncol=1)+
#  scale_y_continuous(breaks = seq(-0.4, 0.5, 0.1))+
#  theme(strip.text = element_text(color ="White",size = 12),
#        axis.title = element_text(size=12),
#        axis.text = element_text(size=12),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        legend.position = "none")

#Add color code to the panels 
#g_euk <- ggplot_gtable(ggplot_build(p2))
#stript <- which(grepl('strip-t', g_euk$layout$name))
#col_strip<-c("purple","darkblue","darkgreen","orange","darkred")
#k <- 1
#for (i in stript) {
#  j <- which(grepl('rect', g_euk$grobs[[i]]$grobs[[1]]$childrenOrder))
#  g_euk$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
#  k <- k+1
#}
#grid.draw(g_euk)

#Plot PCoAs of prokaryotes and eukaryotes together
#library(ggpubr)
#g_prok_euks <- ggarrange(g_prok,g_euk,p3,
#          ncol=3)

##PERMANOVA
#Separated by legs
t2$cal_manova(manova_all = FALSE,group = "Water_depths", by_group = "Leg")
PERMANOVA_Leg_Euk<-t2$res_manova

#Export results
#write.csv(as.data.frame(PERMANOVA_Euk), file="./output/PERMANOVA_euks_no_leg.csv", row.names = TRUE)
write.csv(as.data.frame(PERMANOVA_Leg_Euk), file="./output/ST_2.csv", row.names = TRUE)






####Calculate PCoA for each leg individually#####
#Prokaryotes
#Leg 1
Leg_1_Prok <- clone(dataset)
Leg_1_Prok$sample_table <- subset(dataset$sample_table,Leg == 1)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_1_Prok$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_1_Prok

#Make the OTU and sample information consistent across all files in the dataset object
Leg_1_Prok$tidy_dataset()
print(Leg_1_Prok)

#Check the sequence numbers in each sample
Leg_1_Prok$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
# use default parameters
Leg_1_Prok$cal_abund()

# return dataset$taxa_abund
class(Leg_1_Prok$taxa_abund)

####Calculate the alpha diversity####
Leg_1_Prok$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)

####Beta diversity#### --> needs phyl_tree
Leg_1_Prok$cal_betadiv(unifrac = FALSE)
class(Leg_1_Prok$beta_diversity)


#Leg 2
Leg_2_Prok <- clone(dataset)
Leg_2_Prok$sample_table <- subset(dataset$sample_table,Leg == 2)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_2_Prok$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_2_Prok

#Make the OTU and sample information consistent across all files in the dataset object
Leg_2_Prok$tidy_dataset()
print(Leg_2_Prok)

#Check the sequence numbers in each sample
Leg_2_Prok$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
# use default parameters
Leg_2_Prok$cal_abund()

# return dataset$taxa_abund
class(Leg_2_Prok$taxa_abund)

####Calculate the alpha diversity####
Leg_2_Prok$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)

# return dataset$alpha_diversity
class(Leg_2_Prok$alpha_diversity)

####Beta diversity#### --> needs phyl_tree
Leg_2_Prok$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_2_Prok$beta_diversity)


#Leg 3
Leg_3_Prok <- clone(dataset)
Leg_3_Prok$sample_table <- subset(dataset$sample_table,Leg == 3)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_3_Prok$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_3_Prok

#Make the OTU and sample information consistent across all files in the dataset object
Leg_3_Prok$tidy_dataset()
print(Leg_3_Prok)

#Check the sequence numbers in each sample
Leg_3_Prok$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
Leg_3_Prok$cal_abund()


####Calculate the alpha diversity####
Leg_3_Prok$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(Leg_3_Prok$alpha_diversity)

# return dataset$alpha_diversity
class(Leg_3_Prok$alpha_diversity)

####Beta diversity#### --> needs phyl_tree
Leg_3_Prok$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_3_Prok$beta_diversity)


#Leg 4
Leg_4_Prok <- clone(dataset)
Leg_4_Prok$sample_table <- subset(dataset$sample_table,Leg == 4)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_4_Prok$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_4_Prok

#Make the OTU and sample information consistent across all files in the dataset object
Leg_4_Prok$tidy_dataset()
print(Leg_4_Prok)

#Check the sequence numbers in each sample
Leg_4_Prok$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
# use default parameters
Leg_4_Prok$cal_abund()

####Calculate the alpha diversity####
Leg_4_Prok$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)
plot(Leg_4_Prok$alpha_diversity)

# return dataset$alpha_diversity
class(Leg_4_Prok$alpha_diversity)

####Beta diversity#### --> needs phyl_tree
Leg_4_Prok$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_4_Prok$beta_diversity)


#Leg 5
Leg_5_Prok <- clone(dataset)
Leg_5_Prok$sample_table <- subset(dataset$sample_table,Leg == 5)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_5_Prok$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_5_Prok

#Make the OTU and sample information consistent across all files in the dataset object
Leg_5_Prok$tidy_dataset()
print(Leg_5_Prok)

#Check the sequence numbers in each sample
Leg_5_Prok$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
# use default parameters
Leg_5_Prok$cal_abund()

####Calculate the alpha diversity####
Leg_5_Prok$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)

####Beta diversity#### --> needs phyl_tree
Leg_5_Prok$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_5_Prok$beta_diversity)


#Plot Legs
#Leg 1
#dataset$sample_table$Category<-str_replace(dataset$sample_table$Category,"_","-")
t1 <- trans_beta$new(dataset = Leg_1_Prok, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "1" = "Leg 1"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 1"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  #scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Prok_1_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Prok_1_Leg$layout$name))
col_strip<-c("darkred")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Prok_1_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Prok_1_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Prok_1_Leg)


#Leg 2
t1 <- trans_beta$new(dataset = Leg_2_Prok, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "2"= "Leg 2"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 2"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  #scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Prok_2_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Prok_2_Leg$layout$name))
col_strip<-c("orange")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Prok_2_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Prok_2_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Prok_2_Leg)


#Leg 3
t1 <- trans_beta$new(dataset = Leg_3_Prok, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg,"3"="Leg 3"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 3"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  #scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Prok_3_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Prok_3_Leg$layout$name))
col_strip<-c("darkgreen")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Prok_3_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Prok_3_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Prok_3_Leg)


#Leg 4
t1 <- trans_beta$new(dataset = Leg_4_Prok, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg,"4"="Leg 4"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 4"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  #scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Prok_4_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Prok_4_Leg$layout$name))
col_strip<-c("darkblue")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Prok_4_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Prok_4_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Prok_4_Leg)


#Leg 5
t1 <- trans_beta$new(dataset = Leg_5_Prok, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg,"5"="Leg 5"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 5"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  #scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Prok_5_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Prok_5_Leg$layout$name))
col_strip<-c("purple")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Prok_5_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Prok_5_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Prok_5_Leg)



#Eukaryotes
#Leg 1
Leg_1_Euk <- clone(dataset2)
Leg_1_Euk$sample_table <- subset(dataset$sample_table,Leg == 1)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_1_Euk$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_1_Euk

#Make the OTU and sample information consistent across all files in the dataset object
Leg_1_Euk$tidy_dataset()
print(Leg_1_Euk)

#Check the sequence numbers in each sample
Leg_1_Euk$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
Leg_1_Euk$cal_abund()


####Calculate the alpha diversity####
Leg_1_Euk$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)


####Beta diversity#### --> needs phyl_tree
Leg_1_Euk$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_1_Euk$beta_diversity)


#Leg 2
Leg_2_Euk <- clone(dataset2)
Leg_2_Euk$sample_table <- subset(dataset$sample_table,Leg == 2)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_2_Euk$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_2_Euk

#Make the OTU and sample information consistent across all files in the dataset object
Leg_2_Euk$tidy_dataset()
print(Leg_2_Euk)

#Check the sequence numbers in each sample
Leg_2_Euk$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
Leg_2_Euk$cal_abund()

####Calculate the alpha diversity####
Leg_2_Euk$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)

####Beta diversity#### --> needs phyl_tree
Leg_2_Euk$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_2_Euk$beta_diversity)


#Leg 3
Leg_3_Euk <- clone(dataset2)
Leg_3_Euk$sample_table <- subset(dataset$sample_table,Leg == 3)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_3_Euk$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_3_Euk

#Make the OTU and sample information consistent across all files in the dataset object
Leg_3_Euk$tidy_dataset()
print(Leg_3_Euk)

#Check the sequence numbers in each sample
Leg_3_Euk$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
Leg_3_Euk$cal_abund()

# return dataset$taxa_abund
class(Leg_3_Euk$taxa_abund)

# show part of the relative abundance at Phylum level
Leg_3_Euk$taxa_abund$Division_Phylum[1:5, 1:5]


####Calculate the alpha diversity####
Leg_3_Euk$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)

####Beta diversity#### --> needs phyl_tree
Leg_3_Euk$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_3_Euk$beta_diversity)


#Leg 4
Leg_4_Euk <- clone(dataset2)
Leg_4_Euk$sample_table <- subset(dataset$sample_table,Leg == 4)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_4_Euk$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_4_Euk

#Make the OTU and sample information consistent across all files in the dataset object
Leg_4_Euk$tidy_dataset()
print(Leg_4_Euk)

#Check the sequence numbers in each sample
Leg_4_Euk$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
Leg_4_Euk$cal_abund()

# return dataset$taxa_abund
class(Leg_4_Euk$taxa_abund)

####Calculate the alpha diversity####
Leg_4_Euk$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)

####Beta diversity#### --> needs phyl_tree
Leg_4_Euk$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_4_Euk$beta_diversity)


#Leg 5
Leg_5_Euk <- clone(dataset2)
Leg_5_Euk$sample_table <- subset(dataset$sample_table,Leg == 5)

#Remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”
Leg_5_Euk$filter_pollution(taxa = c("mitochondria", "chloroplast"))
Leg_5_Euk

#Make the OTU and sample information consistent across all files in the dataset object
Leg_5_Euk$tidy_dataset()
print(Leg_5_Euk)

#Check the sequence numbers in each sample
Leg_5_Euk$sample_sums() %>% range

####Calculate the taxa abundance at each taxonomic rank####
Leg_5_Euk$cal_abund()

# return dataset$taxa_abund
class(Leg_5_Euk$taxa_abund)


####Calculate the alpha diversity####
Leg_5_Euk$cal_alphadiv(measures= c("Shannon","Simpson","Coverage","InvSimpson","PD"), PD = FALSE)



####Beta diversity#### --> needs phyl_tree
Leg_5_Euk$cal_betadiv(unifrac = FALSE)
#return dataset$beta_diversity
class(Leg_5_Euk$beta_diversity)


#Plot Legs
#Leg 1
t1 <- trans_beta$new(dataset = Leg_1_Euk, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "1" = "Leg 1"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 1"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(position = "right") +  # Move y-axis to the right
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        axis.title.y.right = element_text(angle = 90, hjust = 0.5), # Rotate Y-axis title to horizontal
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Euk_1_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Euk_1_Leg$layout$name))
col_strip<-c("darkred")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Euk_1_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Euk_1_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Euk_1_Leg)


#Leg 2
t1 <- trans_beta$new(dataset = Leg_2_Euk, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "2"= "Leg 2"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 2"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(position = "right") +  # Move y-axis to the right
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        axis.title.y.right = element_text(angle = 90, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Euk_2_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Euk_2_Leg$layout$name))
col_strip<-c("orange")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Euk_2_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Euk_2_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Euk_2_Leg)


#Leg 3
t1 <- trans_beta$new(dataset = Leg_3_Euk, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg,"3"="Leg 3"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 3"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(position = "right") +  # Move y-axis to the right
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        axis.title.y.right = element_text(angle = 90, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Euk_3_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Euk_3_Leg$layout$name))
col_strip<-c("darkgreen")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Euk_3_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Euk_3_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Euk_3_Leg)

#Leg 4
t1 <- trans_beta$new(dataset = Leg_4_Euk, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg,"4"="Leg 4"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 4"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(position = "right") +  # Move y-axis to the right
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        axis.title.y.right = element_text(angle = 90, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Euk_4_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Euk_4_Leg$layout$name))
col_strip<-c("darkblue")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Euk_4_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Euk_4_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Euk_4_Leg)


#Leg 5
t1 <- trans_beta$new(dataset = Leg_5_Euk, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`,  c("Surface","Chlmax","50m","100m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg,"5"="Leg 5"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 5"))

#Separated by leg
p<-t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(position = "right") +  # Move y-axis to the right
  #scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White",size = 24),
        axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        axis.title.y.right = element_text(angle = 90, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
Euk_5_Leg <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', Euk_5_Leg$layout$name))
col_strip<-c("purple")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', Euk_5_Leg$grobs[[i]]$grobs[[1]]$childrenOrder))
  Euk_5_Leg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(Euk_5_Leg)


p2 <- t1$plot_ordination(point_size=5, plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(position = "right") +  # Move y-axis to the right
  # scale_color_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  #scale_shape_discrete(name = "Depth", labels = c("Surface", "Chlmax", ">50m",">75m"))+
  theme(strip.text = element_text(color ="White"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.position="bottom")+
        #legend.key.size = unit(1, 'cm'), #change legend key size
        #legend.key.height = unit(10, 'cm'), #change legend key height
        #legend.key.width = unit(10, 'cm'))+ #change legend key width
labs(color="Water depths",
     shape="Water depths")

p3 <- get_legend(p2)


#Merge plots
#Plot PCoAs together
Prok <- ggarrange(Prok_1_Leg,Prok_2_Leg,Prok_3_Leg,Prok_4_Leg,Prok_5_Leg,
                  #labels = c("A","B"),
                  nrow=5)
Euks <- ggarrange(Euk_1_Leg,Euk_2_Leg,Euk_3_Leg,Euk_4_Leg,Euk_5_Leg,nrow=5)
Figure_3 <- ggarrange(Prok, Euks, p3,
          nrow= 1, ncol = 3, widths = c(1,1,0.75))

ggsave(Figure_3, 
       filename = "./figures/Figure_3.jpeg",
       device = "jpeg",
       width = 10,
       height = 10,
       unit = "cm")
