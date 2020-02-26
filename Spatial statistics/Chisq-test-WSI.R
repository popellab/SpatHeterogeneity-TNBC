#################################################################################
# This script is used to perform Chi-square test for the whole sample (Fig. S5) #
#################################################################################

library(R.matlab)
library(ggplot2)
library(tidyverse)
library(grid)
library(sp)
library(plyr)
library(reshape2)

# for formatting grids



# combine

Profile_collect <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(Profile_collect) <- c('Section', 'Area', 'Real','Binomial', 'Density', 'Marker', 'Case')

# Profile path
Profile_prefix <- 'Final_diction_'
for(case in seq(1,6)){
  case_path <- switch (case,'./Case 0/', './Case 1/', './Case 2/', './Case 3/', './Case 4/', './Case 5/')
  case_name <- switch (case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  stripe_name <- switch (case,'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  
  for(marker in seq(1,5)){
    marker_name <- switch(marker, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
    # get profile
    Profile <- readRDS(paste(case_path, 'Final_diction/', Profile_prefix, marker_name, '.rds', sep = ''))
    theo_num <- (sum(Profile$Number)*Profile$Area)/sum(Profile$Area)
    
    wTheo_sub <- data.frame(cbind(Profile$Section, Profile$Area, Profile$Number, theo_num, Profile$Density,marker_name, stripe_name))
    colnames(wTheo_sub) <- c('Section', 'Area', 'Real','Binomial', 'Density', 'Marker', 'Case')
    
    wTheo_sub$Real <- as.numeric(as.character(wTheo_sub$Real))
    wTheo_sub$Binomial <- as.numeric(as.character(wTheo_sub$Binomial))
    
    csDat <- data.frame(t(wTheo_sub[,3:4]))
    cs_test <- chisq.test(csDat)
    print(cs_test$p.value < 0.0001)
    Profile_collect <- rbind(Profile_collect, wTheo_sub)
    
  }
}

Profile_collect$Real <- as.numeric(as.character(Profile_collect$Real))
Profile_collect$Binomial <- as.numeric(as.character(Profile_collect$Binomial))
Profile_collect$Area <- as.numeric(as.character(Profile_collect$Area))
Profile_collect$Section <- as.numeric(as.character(Profile_collect$Section))

Profile_melt <- melt(Profile_collect, variable.name = 'type', value.names = 'number', id.vars = c('Section','Marker','Case','Area', 'Density'))
Profile_melt$type <- as.character(Profile_melt$type)
csDat <- data.frame((Profile_collect[Profile_collect$Marker == 'CD3' & Profile_collect$Case == 'Case 4',3:4]))
csDat_Area <- data.frame((Profile_collect[Profile_collect$Marker == 'CD3' & Profile_collect$Case == 'Case 4',2]))/sum(data.frame((Profile_collect[Profile_collect$Marker == 'CD3' & Profile_collect$Case == 'Case 4',2])))

chisq.test(csDat$Real, p = csDat_Area[,1])

chisq.test(csDat)
kruskal.test(csDat)

# filter 
Profile_melt <- Profile_melt[Profile_melt$Section > 0 & Profile_melt$Area > 1.1,]
jpeg('./Chisq_test.jpeg', units="in", width=10, height=12, res=600)
ggplot(Profile_melt, aes(type, value)) + 
  geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
  theme_bw() +
  #stat_compare_means(method = "chisq.test") +
  facet_grid(Case ~ Marker, scales = 'free') +
  geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = TRUE, test = function(a, b){
    list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 6, vjust=1.2) +
  background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
  xlab('Category') +
  ylab('Cell counts(Frequencies)')+
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
        strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
  #coord_flip() +
  panel_border()
dev.off()




ggpaired(Dat_ROI, x = "criteria", y = "value",
         color = "criteria", line.color = "gray", line.size = 0.6, point.size = 5)+
  theme_bw() + 
  
  xlab('Methods') +
  ylab('DSC Scores') +
  ################################
{
  CD3 <- Profile_melt[Profile_melt$Mark == 'CD3' & Profile_melt$Case == 'Case 1A',]
  CD4 <- Profile_melt[Profile_melt$Mark == 'CD4' & Profile_melt$Case == 'Case 1A',]
  CD8 <- Profile_melt[Profile_melt$Mark == 'CD8' & Profile_melt$Case == 'Case 1A',]
  CD20 <- Profile_melt[Profile_melt$Mark == 'CD20' & Profile_melt$Case == 'Case 1A',]
  FoxP3 <- Profile_melt[Profile_melt$Mark == 'FoxP3' & Profile_melt$Case == 'Case 1A',]
  
  p1 <- 
    ggplot(CD3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  p2 <- 
    ggplot(CD4, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p3 <- 
    ggplot(CD8, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  
  p4 <- 
    ggplot(CD20, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p5 <- 
    ggplot(FoxP3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab('Case 1A')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank()) +
    #coord_flip() +
    panel_border()
  Case1A <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Profile_melt[Profile_melt$Mark == 'CD3' & Profile_melt$Case == 'Case 1B',]
  CD4 <- Profile_melt[Profile_melt$Mark == 'CD4' & Profile_melt$Case == 'Case 1B',]
  CD8 <- Profile_melt[Profile_melt$Mark == 'CD8' & Profile_melt$Case == 'Case 1B',]
  CD20 <- Profile_melt[Profile_melt$Mark == 'CD20' & Profile_melt$Case == 'Case 1B',]
  FoxP3 <- Profile_melt[Profile_melt$Mark == 'FoxP3' & Profile_melt$Case == 'Case 1B',]
  
  
  p1 <- 
    ggplot(CD3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  p2 <- 
    ggplot(CD4, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p3 <- 
    ggplot(CD8, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  
  p4 <- 
    ggplot(CD20, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p5 <- 
    ggplot(FoxP3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab('Case 1B')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank()) +
    #coord_flip() +
    panel_border()
  Case1B <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Profile_melt[Profile_melt$Mark == 'CD3' & Profile_melt$Case == 'Case 2',]
  CD4 <- Profile_melt[Profile_melt$Mark == 'CD4' & Profile_melt$Case == 'Case 2',]
  CD8 <- Profile_melt[Profile_melt$Mark == 'CD8' & Profile_melt$Case == 'Case 2',]
  CD20 <- Profile_melt[Profile_melt$Mark == 'CD20' & Profile_melt$Case == 'Case 2',]
  FoxP3 <- Profile_melt[Profile_melt$Mark == 'FoxP3' & Profile_melt$Case == 'Case 2',]
  
  
  p1 <- 
    ggplot(CD3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  p2 <- 
    ggplot(CD4, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p3 <- 
    ggplot(CD8, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  
  p4 <- 
    ggplot(CD20, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p5 <- 
    ggplot(FoxP3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab('Case 2')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank()) +
    #coord_flip() +
    panel_border()
  Case2<-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Profile_melt[Profile_melt$Mark == 'CD3' & Profile_melt$Case == 'Case 3',]
  CD4 <- Profile_melt[Profile_melt$Mark == 'CD4' & Profile_melt$Case == 'Case 3',]
  CD8 <- Profile_melt[Profile_melt$Mark == 'CD8' & Profile_melt$Case == 'Case 3',]
  CD20 <- Profile_melt[Profile_melt$Mark == 'CD20' & Profile_melt$Case == 'Case 3',]
  FoxP3 <- Profile_melt[Profile_melt$Mark == 'FoxP3' & Profile_melt$Case == 'Case 3',]
  
  
  p1 <- 
    ggplot(CD3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  p2 <- 
    ggplot(CD4, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p3 <- 
    ggplot(CD8, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  
  p4 <- 
    ggplot(CD20, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p5 <- 
    ggplot(FoxP3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab('Case 3')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank()) +
    #coord_flip() +
    panel_border()
  Case3 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Profile_melt[Profile_melt$Mark == 'CD3' & Profile_melt$Case == 'Case 4',]
  CD4 <- Profile_melt[Profile_melt$Mark == 'CD4' & Profile_melt$Case == 'Case 4',]
  CD8 <- Profile_melt[Profile_melt$Mark == 'CD8' & Profile_melt$Case == 'Case 4',]
  CD20 <- Profile_melt[Profile_melt$Mark == 'CD20' & Profile_melt$Case == 'Case 4',]
  FoxP3 <- Profile_melt[Profile_melt$Mark == 'FoxP3' & Profile_melt$Case == 'Case 4',]
  
  
  p1 <- 
    ggplot(CD3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  p2 <- 
    ggplot(CD4, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p3 <- 
    ggplot(CD8, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  
  p4 <- 
    ggplot(CD20, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p5 <- 
    ggplot(FoxP3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab(NULL) +
    ylab('Case 4')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank()) +
    #coord_flip() +
    panel_border()
  Case4 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Profile_melt[Profile_melt$Mark == 'CD3' & Profile_melt$Case == 'Case 5',]
  CD4 <- Profile_melt[Profile_melt$Mark == 'CD4' & Profile_melt$Case == 'Case 5',]
  CD8 <- Profile_melt[Profile_melt$Mark == 'CD8' & Profile_melt$Case == 'Case 5',]
  CD20 <- Profile_melt[Profile_melt$Mark == 'CD20' & Profile_melt$Case == 'Case 5',]
  FoxP3 <- Profile_melt[Profile_melt$Mark == 'FoxP3' & Profile_melt$Case == 'Case 5',]
  
  
  p1 <- 
    ggplot(CD3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab('CD3') +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  p2 <- 
    ggplot(CD4, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab('CD4') +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p3 <- 
    ggplot(CD8, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab('CD8') +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  
  p4 <- 
    ggplot(CD20, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab('CD20') +
    ylab(NULL)+
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25)) +
    #coord_flip() +
    panel_border()
  
  p5 <- 
    ggplot(FoxP3, aes(type, value)) + 
    geom_boxplot(aes(fill = factor(type), alpha = 0.5), outlier.shape = NA,show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE,aes(colour = factor(type)), alpha = 0.5) +
    theme_bw() +
    geom_signif(comparisons = list(c("Real", "Binomial")), map_signif_level = FALSE, test = function(a, b){
      list(p.value = chisq.test(data.frame(cbind(a, b)))$p.value)}, textsize = 4, vjust=0.2) +
    background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
    xlab('FoxP3') +
    ylab('Case 5')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 30), 
          strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 25), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank()) +
    #coord_flip() +
    panel_border()
  Case5 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
jpeg('Chisq-test-1225.jpeg', units="in", width=15, height=20, res=300)
plot_grid(Case1A, Case1B, Case2, Case3,Case4, Case5, labels = c('','','',''),  ncol = 1, align = 'hv')+
  ylab('Density')
dev.off()