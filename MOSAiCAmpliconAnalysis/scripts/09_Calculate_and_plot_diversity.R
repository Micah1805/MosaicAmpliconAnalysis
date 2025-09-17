# Script: 09_Calculate_and_plot_diversity.R
# Purpose: Create plots showing alpha diversity and  PCoA (steps modified from the Microeco-Tutorial (https://chiliubio.github.io/microeco_tutorial/))
# Input: microeco_prok.Rdata, microeco_euk.Rdata
# Main output: SF_3_a.jpeg, SF_2_a.jpeg, ST_1.csv, SF_3_b.jpeg, SF_2_b.jpeg, Figure_3.jpeg, ST_2.csv
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

#Set background
theme_set(theme_bw())

####Prokaryotes####
load("./output/microeco_prok.Rdata")
dataset$sample_table[c("Water", "Water_depths")] <- str_split_fixed(dataset$sample_table$Category, "_", 2) 

##Alpha diversity
t1 <- trans_alpha$new(dataset = dataset, group = "Water_depths")
head(t1$data_stat)

#Calculate difference between the groups
#ANOVA
t1$cal_diff(method = "anova")
head(t1$res_diff)

t1$cal_diff(method = "anova")

#Set order labels legend
t1$data_alpha$`Water_depths` = factor(t1$data_alpha$`Water_depths` , c("Surface","Chlmax",">50m",">75m"))

#Plot alpha-diversity
alpha_div_proks <- t1$plot_alpha(measure = "Shannon", add_sig_text_size = 6, add = "jitter", xtext_size = 15)+
  scale_y_continuous(breaks = seq(0, 7, 0.5))

ggsave(alpha_div_proks, 
       filename = "./figures/SF_3_a.jpeg",
       device = "jpeg",
       width = 10,
       height = 10)

#Beta-diversity
t1 <- trans_beta$new(dataset = dataset, group = "Water_depths", measure = "bray")

#Calculate PCoA
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)

#Set order labels legend
t1$res_ordination$scores$`Water_depths` = factor(t1$res_ordination$scores$`Water_depths`, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
t1$res_ordination$scores <- t1$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "1" = "Leg 1: 20.09.2019-13.12.2019","2"= "Leg 2: 13.12.2019-24.02.2020","3"="Leg 3: 24.02.2020-04.06.2020","4"="Leg 4: 04.06.2020-12.08.2020","5"="Leg 5: 12.08.2020-12.10.2020"))
t1$res_ordination$scores$Leg = factor(t1$res_ordination$scores$Leg, c("Leg 1: 20.09.2019-13.12.2019","Leg 2: 13.12.2019-24.02.2020","Leg 3: 24.02.2020-04.06.2020","Leg 4: 04.06.2020-12.08.2020","Leg 5: 12.08.2020-12.10.2020"))

#Plot the PCoA result with confidence ellipse
PCoA_proks <- t1$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point","ellipse"), ellipse_chull_fill = FALSE)+
  theme(axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=15,face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))+
  guides(fill=guide_legend(title="Water_depths"))

ggsave(PCoA_proks, 
       filename = "./figures/SF_2_a.jpeg",
       device = "jpeg",
       width = 10,
       height = 8)

#Separated by leg
p<-t1$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1))+
  theme(strip.text = element_text(color ="White",size = 12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Add color code to the panels 
g_prok <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g_prok$layout$name))
col_strip<-c("purple","darkblue","darkgreen","orange","darkred")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g_prok$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_prok$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_prok)

ggsave(g_prok, 
       filename = "./figures/PCoA_legs_proks.jpeg",
       device = "jpeg",
       width = 4,
       height = 10)

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
t2$data_alpha$`Water_depths` = factor(t2$data_alpha$`Water_depths`, c("Surface","Chlmax",">50m",">75m"))


#Plot alpha-diversity
alpha_div_euks <- t2$plot_alpha(measure = "Shannon", add_sig_text_size = 6, add = "jitter")+
  scale_y_continuous(breaks = seq(0, 7, 0.5))+
  coord_cartesian(ylim = c(4,7))

ggsave(alpha_div_euks, 
       filename = "./figures/SF_3_b.jpeg",
       device = "jpeg",
       width = 10,
       height = 10)

##Beta-diversity
t2 <- trans_beta$new(dataset = dataset2, group = "Water_depths", measure = "bray")

#Calculate PCoA
t2$cal_ordination(method = "PCoA")
class(t2$res_ordination)

#Define order labels legend
t2$res_ordination$scores$`Water_depths` = factor(t2$res_ordination$scores$`Water_depths`, c("Surface","Chlmax",">50m",">75m"))

#Rename Legs
t2$res_ordination$scores <- t2$res_ordination$scores %>%
  # Rename 
  mutate(Leg= recode(Leg, "1" = "Leg 1: 20.09.2019-13.12.2019","2"= "Leg 2: 13.12.2019-24.02.2020","3"="Leg 3: 24.02.2020-04.06.2020","4"="Leg 4: 04.06.2020-12.08.2020","5"="Leg 5: 12.08.2020-12.10.2020"))
t2$res_ordination$scores$Leg = factor(t2$res_ordination$scores$Leg, c("Leg 1: 20.09.2019-13.12.2019","Leg 2: 13.12.2019-24.02.2020","Leg 3: 24.02.2020-04.06.2020","Leg 4: 04.06.2020-12.08.2020","Leg 5: 12.08.2020-12.10.2020"))

#Plot the PCoA result with confidence ellipse
PCoA_euks <- t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point","ellipse"),, ellipse_chull_fill = FALSE)+
  theme(axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=15,face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

ggsave(PCoA_euks, 
       filename = "./figures/SF_2_b.jpeg",
       device = "jpeg",
       width = 10,
       height = 8)

#Separated by leg
p2<-t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(breaks = seq(-0.4, 0.5, 0.1))+
  theme(strip.text = element_text(color ="White"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text())

p3<-get_legend(p2)

#Separated by leg
p2<-t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"),ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(breaks = seq(-0.4, 0.5, 0.1))+
  theme(strip.text = element_text(color ="White",size = 12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

#Add color code to the panels 
g_euk <- ggplot_gtable(ggplot_build(p2))
stript <- which(grepl('strip-t', g_euk$layout$name))
col_strip<-c("purple","darkblue","darkgreen","orange","darkred")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g_euk$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_euk$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col_strip[k]
  k <- k+1
}
grid.draw(g_euk)

ggsave(g_euk, 
       filename = "./figures/PCoA_legs_euks.jpeg",
       device = "jpeg",
       width = 4,
       height = 10)

#Plot PCoAs of prokaryotes and eukaryotes together
library(ggpubr)
g_prok_euks <- ggarrange(g_prok,g_euk,p3,
          ncol=3)

ggsave(g_prok_euks, 
       filename = "./figures/Figure_3.jpeg",
       device = "jpeg",
       width = 10,
       height = 10)

##PERMANOVA
#Separated by legs
t2$cal_manova(manova_all = FALSE,group = "Water_depths", by_group = "Leg")
PERMANOVA_Leg_Euk<-t2$res_manova

#Export results
#write.csv(as.data.frame(PERMANOVA_Euk), file="./output/PERMANOVA_euks_no_leg.csv", row.names = TRUE)
write.csv(as.data.frame(PERMANOVA_Leg_Euk), file="./output/ST_2.csv", row.names = TRUE)
