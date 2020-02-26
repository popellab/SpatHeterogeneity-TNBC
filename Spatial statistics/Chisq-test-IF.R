#################################################################################
# This script is used to perform Chi-square test whitn invasive front (Fig. S6) #
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

Profile_collect <- readRDS('./Invasive_densityProfile_combined.rds')
colnames(Profile_collect) <-  c('Section', 'Area', 'Real','Density','Marker', 'Case')

Profile_wTheo <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(Profile_wTheo) <- c('Section', 'Area', 'Real', 'Density', 'Marker', 'Case', 'Binomial')

for(case in seq(1,6)){
  case_name <- switch (case,'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  for(marker in seq(1,5)){
    marker_name <- switch(marker, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
    Profile_sub <- Profile_collect[Profile_collect$Case == case_name & Profile_collect$Marker == marker_name,]
    total_area <- sum(Profile_sub$Area)
    total_num <- sum(Profile_sub$Real)
    
    Profile_sub$Binomial <- (total_num*Profile_sub$Area)/total_area
    
    Profile_wTheo <- rbind(Profile_wTheo, Profile_sub)
  }
}





Profile_melt_IF <- melt(Profile_wTheo, variable.name = 'type', value.names = 'number', id.vars = c('Section','Marker','Case','Area', 'Density'))
Profile_melt_IF$type <- as.character(Profile_melt_IF$type)
csDat <- data.frame((Profile_collect[Profile_collect$Marker == 'CD3' & Profile_collect$Case == 'Case 1A',3:4]))
chisq.test(csDat)

# filter 
jpeg('./Chisq_test_IF.jpeg', units="in", width=10, height=12, res=600)
ggplot(Profile_melt_IF, aes(type, value)) + 
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
  ##################################################
{
  CD3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD3' & Profile_melt_IF$Case == 'Case 1A',]
  CD4 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD4' & Profile_melt_IF$Case == 'Case 1A',]
  CD8 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD8' & Profile_melt_IF$Case == 'Case 1A',]
  CD20 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD20' & Profile_melt_IF$Case == 'Case 1A',]
  FoxP3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'FoxP3' & Profile_melt_IF$Case == 'Case 1A',]
  
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
  CD3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD3' & Profile_melt_IF$Case == 'Case 1B',]
  CD4 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD4' & Profile_melt_IF$Case == 'Case 1B',]
  CD8 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD8' & Profile_melt_IF$Case == 'Case 1B',]
  CD20 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD20' & Profile_melt_IF$Case == 'Case 1B',]
  FoxP3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'FoxP3' & Profile_melt_IF$Case == 'Case 1B',]
  
  
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
  CD3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD3' & Profile_melt_IF$Case == 'Case 2',]
  CD4 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD4' & Profile_melt_IF$Case == 'Case 2',]
  CD8 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD8' & Profile_melt_IF$Case == 'Case 2',]
  CD20 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD20' & Profile_melt_IF$Case == 'Case 2',]
  FoxP3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'FoxP3' & Profile_melt_IF$Case == 'Case 2',]
  
  
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
  CD3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD3' & Profile_melt_IF$Case == 'Case 3',]
  CD4 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD4' & Profile_melt_IF$Case == 'Case 3',]
  CD8 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD8' & Profile_melt_IF$Case == 'Case 3',]
  CD20 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD20' & Profile_melt_IF$Case == 'Case 3',]
  FoxP3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'FoxP3' & Profile_melt_IF$Case == 'Case 3',]
  
  
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
  plot(p1)
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
  CD3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD3' & Profile_melt_IF$Case == 'Case 4',]
  CD4 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD4' & Profile_melt_IF$Case == 'Case 4',]
  CD8 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD8' & Profile_melt_IF$Case == 'Case 4',]
  CD20 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD20' & Profile_melt_IF$Case == 'Case 4',]
  FoxP3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'FoxP3' & Profile_melt_IF$Case == 'Case 4',]
  
  
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
  CD3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD3' & Profile_melt_IF$Case == 'Case 5',]
  CD4 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD4' & Profile_melt_IF$Case == 'Case 5',]
  CD8 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD8' & Profile_melt_IF$Case == 'Case 5',]
  CD20 <- Profile_melt_IF[Profile_melt_IF$Mark == 'CD20' & Profile_melt_IF$Case == 'Case 5',]
  FoxP3 <- Profile_melt_IF[Profile_melt_IF$Mark == 'FoxP3' & Profile_melt_IF$Case == 'Case 5',]
  
  
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
jpeg('Chisq-test-IF-1225.jpeg', units="in", width=15, height=20, res=300)
plot_grid(Case1A, Case1B, Case2, Case3,Case4, Case5, labels = c('','','',''),  ncol = 1, align = 'hv')+
  ylab('Density')
dev.off()