################################################################################
####################Calculate and plot diversity################################
################################################################################

#Load packages
library(microeco)
library(ggplot2)
library(tidyverse)
library(grid)
library(xlsx)

#Set background
theme_set(theme_bw())

####Prokaryotes####
load("C:/Users/jahuem001/OneDrive/Dokumente/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/dataset.Rdata")
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

#Plot alpha-diversityt
t1$plot_alpha(measure = "Shannon", y_increase = 0.3)
t1$plot_alpha(measure = "Simpson", y_increase = 0.1)

t1$plot_alpha(measure = "Shannon", add_sig_text_size = 6, add = "jitter", xtext_size = 15)+
  scale_y_continuous(breaks = seq(0, 7, 0.5))

#Diversity in different legs
t1 <- trans_alpha$new(dataset = dataset, group = "Water_depths",by_group = "Leg")
anova<-t1$data_stat

#Calculate difference between the groups
#ANOVA
t1$cal_diff(method = "anova")
anova<-t1$res_diff

#Set order labels legend
t1$data_alpha$`Water_depths` = factor(t1$data_alpha$`Water_depths` , c("Surface","Chlmax",">50m",">75m"))

#Plot alpha-diversity
t1$plot_alpha(measure = "Shannon")

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
t1$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point","ellipse"), ellipse_chull_fill = FALSE)+
  theme(axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=15,face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))+
  guides(fill=guide_legend(title="Water_depths"))


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


#Plot and compare group distances
t1 <- trans_beta$new(dataset = dataset, group = "Water_depths", measure = "bray")
#Calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
#Perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")

t1$plot_group_distance(add = "mean")

# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "wilcox")
t1$plot_group_distance(plot_type = "ggviolin", add = "mean_se")
t1$plot_group_distance(add = "mean")

##PERMANOVA
#For all
t1$cal_manova(manova_all = FALSE)
PERMANOVA_Prok<-t1$res_manova

#Separated by legs
t1$cal_manova(manova_all = FALSE,group = "Water_depths", by_group = "Leg")
PERMANOVA_Leg_Prok<-t1$res_manova

#Export results
write.xlsx(as.data.frame(PERMANOVA_Prok), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/PERMANOVA_Prok.xlsx",
           col.names = TRUE, row.names = TRUE)
write.xlsx(as.data.frame(PERMANOVA_Leg_Prok), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/PERMANOVA_Leg_Prok.xlsx",
           col.names = TRUE, row.names = TRUE)

##ANOSIM
#For all
t1$cal_anosim(paired=TRUE)
ANOSIM_Prok<-t1$res_anosim

#Separated by legs
t1$cal_anosim(paired=TRUE,group = "Water_depths", by_group = "Leg")
ANOSIM_Leg_Prok<-t1$res_anosim

#Export ANOSIM results
write.xlsx(as.data.frame(ANOSIM_Prok), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/ANOSIM_Prok.xlsx",
           col.names = TRUE, row.names = TRUE)
write.xlsx(as.data.frame(ANOSIM_Leg_Prok), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/ANOSIM_Leg_Prok.xlsx",
           col.names = TRUE, row.names = TRUE)

####Eukaryotes####

#Load data
load("C:/Users/jahuem001/OneDrive/Dokumente/AWI/MOSAiC/Sequencing_Prok_Euk/Rdata/Tables_microeco/Prok_Euk_sep/dataset2.Rdata")
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


t2$cal_diff(method = "anova")

#Set order lables legend
t2$data_alpha$`Water_depths` = factor(t2$data_alpha$`Water_depths`, c("Surface","Chlmax",">50m",">75m"))


#Plot alpha-diversity
t2$plot_alpha(measure = "Shannon", y_increase = 0.3)
t2$plot_alpha(measure = "Shannon", y_increase = 0.1)

t2$plot_alpha(measure = "Shannon", add_sig_text_size = 6, add = "jitter")+
  scale_y_continuous(breaks = seq(0, 7, 0.5))+
  coord_cartesian(ylim = c(4,7))

#Diversity in different legs
t2 <- trans_alpha$new(dataset = dataset2, group = "Water_depths",by_group = "Leg")
anova<-t2$data_stat


#Calculate difference between the groups
#ANOVA
t2$cal_diff(method = "anova")
anova2<-t2$res_diff

#Set order labels legend
t2$data_alpha$`Water_depths` = factor(t2$data_alpha$`Water_depths` , c("Surface","Chlmax",">50m",">75m"))

#Plot alpha-diversity
t2$plot_alpha(measure = "Shannon")



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
t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point","ellipse"),, ellipse_chull_fill = FALSE)+
  theme(axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=15,face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

#Separated by leg
p2<-t2$plot_ordination(plot_color = "Water_depths", plot_shape = "Water_depths", plot_type = c("point", "chull"), ellipse_chull_fill = FALSE)+
  facet_wrap(~Leg,ncol=1)+
  scale_y_continuous(breaks = seq(-0.4, 0.5, 0.1))++
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


#Plot PCoAs of prokaryotes and eukaryotes together
library(ggpubr)
ggarrange(g_prok,g_euk,p3,
          ncol=3)


##PERMANOVA
#For all
t2$cal_manova(manova_all = FALSE)
PERMANOVA_Euk<-t2$res_manova

#Separated by legs
t2$cal_manova(manova_all = FALSE,group = "Water_depths", by_group = "Leg")
PERMANOVA_Leg_Euk<-t2$res_manova

#Export results
write.xlsx(as.data.frame(PERMANOVA_Euk), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/PERMANOVA_Euk.xlsx",
           col.names = TRUE, row.names = TRUE)
write.xlsx(as.data.frame(PERMANOVA_Leg_Euk), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/PERMANOVA_Leg_Euk.xlsx",
           col.names = TRUE, row.names = TRUE)

##ANOSIM
#For all
t2$cal_anosim(paired=TRUE)
ANOSIM_Euk<-t2$res_anosim

#Separated by legs
t2$cal_anosim(paired=TRUE,group = "Water_depths", by_group = "Leg")
ANOSIM_Leg_Euk<-t2$res_anosim

#Export ANOSIM results
write.xlsx(as.data.frame(ANOSIM_Euk), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/ANOSIM_Euk.xlsx",
           col.names = TRUE, row.names = TRUE)
write.xlsx(as.data.frame(ANOSIM_Leg_Euk), file="~/AWI/MOSAiC/Sequencing_Prok_Euk/Graphics/Final/10.03.2025/ANOSIM_Leg_Euk.xlsx",
           col.names = TRUE, row.names = TRUE)
