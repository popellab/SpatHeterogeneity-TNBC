  #######################################################################
##### This script is used for analyzing intra-tumoral heterogeenity####
#######################################################################

# annotation.R script should be run first, to obtain polygon information

library(spatstat)
library(R.matlab)
library(ggplot2)
library(scales) # to access break formatting functions
library(dplyr)
library(kulife)
library(rgeos)
library(gridExtra)
library(concaveman)
library(cowplot)
library(epade)
library(pracma) # to calculate AUC
library(RANN)
library(ggpubr)
source('MiFunction.R')


# Find outliers
library(tsoutliers)
library(devtools)
p_th = 5E-2
# Define global variable
#sub_path <- './Case_2/'
Var_suffix <- '_Points.mat'
outname = "fittedResult_"


##########################################################################
############ First-order property and Spatial model fitting ##############
##########################################################################

########
# The definitions of parameters are defined in the manuscript
#######
p_all <- data.frame(matrix(nrow=0,ncol=15))
header <- c('index',"xstart","xend", "ystart","yend","n","intensity",'notCSR?','fit','kappa','sigma2','mu','location', 'Marker', 'Case')

colnames(p_all) <- header

for(case in seq(1,6,1)){
  
  case_name <- switch(case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  
  # Change case here
  sub_path <- switch(case,'./Case_0/', './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  # Loop to calculate statistics
  for(label_idx in seq(1,5,1)){
    # For case 2
    Current_lab <- switch (label_idx,'CD3','CD4','CD8','CD20','FoxP3')
    
    
    Coords_path <- paste(sub_path, Current_lab, '/',Current_lab, Var_suffix, sep = '')
    marker_Coords <- as.data.frame(readMat(Coords_path))/125
    
    #######################
    if(case == 1){
      InvasiveFile <- readRDS(paste(sub_path, 'Invasive_1a.rds', sep = ''))
      
      ##########
      
      TumorFile <- readRDS(paste(sub_path, 'Tumor_1a.rds', sep = ''))
      NormalFile <- readRDS(paste(sub_path, 'Normal_1a.rds', sep = ''))
      }
    if(case == 2){
      
      InvasiveFile <- readRDS(paste(sub_path, 'Invasive_1b.rds', sep = ''))
      
      TumorFile <- readRDS(paste(sub_path, 'Tumor_1b.rds', sep = ''))
      
      NormalFile <- readRDS(paste(sub_path, 'Normal_1b.rds', sep = ''))
      }
    if(case > 2){
      
      InvasiveFilePath <- paste(sub_path, 'Invasive.rds', sep = '')
      InvasiveFile <- readRDS(InvasiveFilePath)
      
      TumorFilePath <- paste(sub_path, 'Tumor.rds', sep = '')
      TumorFile <- readRDS(TumorFilePath)
      
      NormalFile <- readRDS(paste(sub_path, 'Normal.rds', sep = ''))
      
    }
    # Characterize all points
    

    # input file-related path
    # output file-related path
    # Scanning the whole space
    window <- 0.15
    
    # Index for each window
    index <- 0
    # Create files for store fit results for current case 
    for(y_origin in seq(0, 25-window, window)){
        for(x_origin in seq(0,35-window, window)){
          window_coords <- data.frame(cbind(c(x_origin, x_origin, x_origin+window, x_origin+window), c(y_origin, y_origin + window, y_origin + window, y_origin)))
          
          # Characterize the window
          if(isTRUE(ncol(NormalFile) == 3)){
            normal_polygon_num <- NormalFile[nrow(NormalFile),3]
            normal_count <- 0
            for(num in seq(1, normal_polygon_num)){
              subNormalFile <- NormalFile[NormalFile$label == num, ]
              normal <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], subNormalFile[,1], subNormalFile[,2]))
              subnormal_count <- nrow(data.frame(normal[normal[,1] >0,]))
              normal_count <- normal_count + subnormal_count
            }
          }
          if(isTRUE(ncol(NormalFile) == 2)){
            normal_count <- 0
            normal <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], NormalFile[,1], NormalFile[,2]))
            normal_count <- nrow(data.frame(normal[normal[,1] >0,]))
          }
          
          
          invasive <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], InvasiveFile[,1], InvasiveFile[,2]))
          invasive_count <- nrow(data.frame(invasive[invasive[,1] >0,]))
          
          tumor <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], TumorFile[,1], TumorFile[,2]))
          tumor_count <- nrow(data.frame(tumor[tumor[,1] >0,]))
          
          
          if(isTRUE(invasive_count != 0)){
            loc <- 'IF'
          } 
          if(isTRUE(invasive_count == 0 & tumor_count != 0)){
            loc <- 'CT'
          }
          if(isTRUE(invasive_count == 0 & tumor_count == 0 & normal_count != 0)){
            loc <- 'N'
          }
          if(isTRUE(invasive_count == 0 & tumor_count == 0 & normal_count == 0)){
            loc <- 'BG'
          }
          if(isTRUE(loc != 'BG')){
            SubRegion <- marker_Coords[marker_Coords[,1] < (x_origin+window)  & marker_Coords[,1] >= x_origin & marker_Coords[,2] < (y_origin + window) & marker_Coords[,2] >= y_origin,]
            SubRegion <- SubRegion[complete.cases(SubRegion), ]
            fitParameter <- subRegionFit(SubRegion, x_origin, y_origin, window, alpha = p_th)
            if(isTRUE(nrow(SubRegion) != 0)){
              index <- index + 1
              number <- nrow(SubRegion)
              density <- number/(window * window)
              p_sub <- data.frame(cbind(index,x_origin,x_origin+window,y_origin,y_origin+window, number, density,fitParameter, loc,Current_lab, case_name ))
              names(p_sub) <- header
              p_all <- rbind(p_all, p_sub)
            } 
         
          }
        }
      }  
    # write results to file
    print('finish')

  }
}

write.csv(p_all, './waterfall_data_fit.csv') # run after the above loop finishing


####################
###### Plotting ####
####################

waterfall <- read.csv('./waterfall_data_fit.csv')
waterfall<- waterfall[waterfall$fit == 1,]

Density_dat <- waterfall[,c(8, 14, 15,16)]
Density_dat$Criteira <- 'Density'
colnames(Density_dat) <- c('value', 'location', 'Marker', 'Case', 'Criteria')
Number_dat <- waterfall[,c(13, 14, 15,16)]
Number_dat$Criteira <- 'Number'
colnames(Number_dat) <- c('value', 'location', 'Marker', 'Case', 'Criteria')

Distance_dat <- waterfall[,c(12, 14, 15,16)]
Distance_dat$sigma2 <- sqrt(Distance_dat$sigma2*pi/2)
Distance_dat$Criteira <- 'Distance'
colnames(Distance_dat) <- c('value', 'location', 'Marker', 'Case', 'Criteria')

Area_dat <- waterfall[,c(12, 14, 15,16)]
Area_dat$sigma2 <- 2*pi*Area_dat$sigma2*log(20)
Area_dat$Criteira <- 'Area'
colnames(Area_dat) <- c('value', 'location', 'Marker', 'Case', 'Criteria')

Stat_comb <- rbind(Density_dat, Number_dat,Distance_dat,  Area_dat)

{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Density',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Density',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Density',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Density',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Density',]
  
  p1 <- ggplot(data = CD3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(expression(paste("Density, ", mm^{-2}))) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- ggplot(data = CD4, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('CD4') +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- ggplot(data = CD8, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('CD8') +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- ggplot(data = CD20, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('CD20') +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- ggplot(data = FoxP3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('FoxP3') +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  Density <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Number',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Number',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Number',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Number',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Number',]
  
  p1 <- 
    ggplot(data = CD3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab('Cell/cluster') +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10, paired = TRUE) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- ggplot(data = CD4, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +  
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- ggplot(data = CD8, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- ggplot(data = CD20, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- ggplot(data = FoxP3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  Count <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Distance',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Distance',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Distance',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Distance',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Distance',]
  
  p1 <- 
    ggplot(data = CD3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(expression(paste("Mean radius, ", mm))) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10, paired = TRUE) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- ggplot(data = CD4, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +  
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- ggplot(data = CD8, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- ggplot(data = CD20, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- ggplot(data = FoxP3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  Radius <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Area',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Area',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Area',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Area',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Area',]
  
  p1 <- 
    ggplot(data = CD3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(expression(paste("Cluster area, ", mm^{2}))) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10, paired = TRUE) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- ggplot(data = CD4, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +  
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- ggplot(data = CD8, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- ggplot(data = CD20, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- ggplot(data = FoxP3, aes(location, value, color = factor(location), fill = factor(location))) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE) +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    stat_compare_means(aes(group = factor(location)), method = 't.test' ,label = "p.signif",ref.group = 'N', show.legend = FALSE, size = 10) +
    
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  Area <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}


jpeg('./Spatial_fit_stat.jpeg', units="in", width=14, height=15, res=300)
plot_grid(Density, Count, Radius, Area,labels = c('', '','',''),  ncol = 1, align = 'hv') +
  xlab('Locations')
dev.off()

#######################
###### Figs.5D, S8 ####
#######################

# for Case 1A
waterfall <- read.csv('./waterfall_data_fit.csv')
waterfall<- waterfall[waterfall$fit == 1,]
waterfall<- waterfall[waterfall$sigma2 != 0,]

waterfall_1A<- waterfall[waterfall$Case == 'Case 0' & waterfall$Marker == 'FoxP3',]

QCoD_density <- (quantile(waterfall_1A$intensity, 0.75) - quantile(waterfall_1A$intensity, 0.25))/ (quantile(waterfall_1A$intensity, 0.75) + quantile(waterfall_1A$intensity, 0.25))
print(QCoD_density)



# cell/ cluster
QCoD_mu <- (quantile(waterfall_1A$mu, 0.75) - quantile(waterfall_1A$mu, 0.25))/ (quantile(waterfall_1A$mu, 0.75) + quantile(waterfall_1A$mu, 0.25))
print(QCoD_mu)

mean_dist <- sqrt(waterfall_1A$sigma2*pi/2)

QCoD_sigma <- (quantile(mean_dist, 0.75) - quantile(mean_dist, 0.25))/ (quantile(mean_dist, 0.75) + quantile(mean_dist, 0.25))
print(QCoD_sigma)


QCoD_dist <- (quantile(waterfall_1A$sigma2*2*log(20)*pi, 0.75) - quantile(waterfall_1A$sigma2*2*log(20)*pi, 0.25))/ (quantile(waterfall_1A$sigma2*2*log(20)*pi, 0.75) + quantile(waterfall_1A$sigma2*2*log(20)*pi, 0.25))
print(QCoD_dist)

waterfall <- waterfall[,c(2,8,14,15,16)]



waterfall_new <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(waterfall_new) <- c('index','intensity','location','Marker','Case')
for(case in seq(1,6,1)){
  
  case_name <- switch(case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')

  # Loop to calculate statistics
  for(label_idx in seq(1,5,1)){
    # For case 2
    Current_lab <- switch (label_idx,'CD3','CD4','CD8','CD20','FoxP3')
    waterfall_sub <- waterfall[waterfall$Marker == Current_lab & waterfall$Case == case_name,]
    waterfall_sub <- waterfall_sub[order(-waterfall_sub$intensity),]
    waterfall_sub$index <- seq(1, nrow(waterfall_sub))
    waterfall_new <- rbind(waterfall_new, waterfall_sub)
  }
   
  
}

write.csv(waterfall_new, './waterfall_data_new.csv', row.names = FALSE)

waterfall_new <- read.csv('./waterfall_data_new.csv')


# Density plot A
{
  CD3 <- waterfall_new[waterfall_new$Marker == 'CD3' & waterfall_new$Case == 'Case 0',]
  CD4 <- waterfall_new[waterfall_new$Marker == 'CD4' & waterfall_new$Case == 'Case 0',]
  CD8 <- waterfall_new[waterfall_new$Marker == 'CD8' & waterfall_new$Case == 'Case 0',]
  CD20 <- waterfall_new[waterfall_new$Marker == 'CD20' & waterfall_new$Case == 'Case 0',]
  FoxP3 <- waterfall_new[waterfall_new$Marker == 'FoxP3' & waterfall_new$Case == 'Case 0',]
  
  CD3$Case <- 'Case 1A'
  CD4$Case <- 'Case 1A'
  CD8$Case <- 'Case 1A'
  CD20$Case <- 'Case 1A'
  FoxP3$Case <- 'Case 1A'
  
  mean = round(mean(CD3$intensity), digit = 2)
  CV = round(sd(CD3$intensity)/mean, digit = 2)
  
  
  jpeg('./CD3_Case1A.jpeg', units="in", width= 45, height=40, res= 600)

  p1 <- 
    ggplot(CD3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    scale_y_continuous(breaks = c(0, 3000, 6000, 9000),labels = c(0, 3000, 6000, 9000)) +
    
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p1)
  dev.off()
  
  mean = round(mean(CD4$intensity), digit = 2)
  CV = round(sd(CD4$intensity)/mean, digit = 2)
  jpeg('./CD4_Case1A.jpeg', units="in", width= 45, height=40, res= 600)
  
  p2 <- 
    ggplot(CD4, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0, 1000, 2000),labels = c(0, 1000, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p2)
  dev.off()
  
  
  
  mean = round(mean(CD8$intensity), digit = 2)
  CV = round(sd(CD8$intensity)/mean, digit = 2)
  
  jpeg('./CD8_Case1A.jpeg', units="in", width= 45, height=40, res= 600)
  
  p3 <- 
    ggplot(CD8, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p3)
  dev.off()
  
  mean = round(mean(CD20$intensity), digit = 2)
  CV = round(sd(CD20$intensity)/mean, digit = 2)
  
  jpeg('./CD20_Case1A.jpeg', units="in", width= 45, height=40, res= 600)
  
  p4 <- 
    ggplot(CD20, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0, 1000, 2000),labels = c(0, 1000, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p4)
  dev.off()
  
  mean = round(mean(FoxP3$intensity), digit = 2)
  CV = round(sd(FoxP3$intensity)/mean, digit = 2)
  
  jpeg('./FoxP3_Case1A.jpeg', units="in", width= 45, height=40, res= 600)
  
  p5 <- 
    ggplot(FoxP3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p5)
  dev.off()
  

  
  
  Case1A <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

{
  CD3 <- waterfall_new[waterfall_new$Marker == 'CD3' & waterfall_new$Case == 'Case 1',]
  CD4 <- waterfall_new[waterfall_new$Marker == 'CD4' & waterfall_new$Case == 'Case 1',]
  CD8 <- waterfall_new[waterfall_new$Marker == 'CD8' & waterfall_new$Case == 'Case 1',]
  CD20 <- waterfall_new[waterfall_new$Marker == 'CD20' & waterfall_new$Case == 'Case 1',]
  FoxP3 <- waterfall_new[waterfall_new$Marker == 'FoxP3' & waterfall_new$Case == 'Case 1',]
  
  CD3$Case <- 'Case 1B'
  CD4$Case <- 'Case 1B'
  CD8$Case <- 'Case 1B'
  CD20$Case <- 'Case 1B'
  FoxP3$Case <- 'Case 1B'
  
  mean = round(mean(CD3$intensity), digit = 2)
  CV = round(sd(CD3$intensity)/mean, digit = 2)
  
  jpeg('./CD3_Case1B.jpeg', units="in", width= 45, height=40, res= 600)
  p1 <- 
    ggplot(CD3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab('CD3')+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0,500, 1000,1500, 2000),labels = c(0,500, 1000,1500, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p1)
  dev.off()
  
  mean = round(mean(CD4$intensity), digit = 2)
  CV = round(sd(CD4$intensity)/mean, digit = 2)
  
  jpeg('./CD4_Case1B.jpeg', units="in", width= 45, height=40, res= 600)
  p2 <- 
    ggplot(CD4, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab('CD4')+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0,500, 1000,1500, 2000),labels = c(0,500, 1000,1500, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p2)
  dev.off()
  
  

  mean = round(mean(CD8$intensity), digit = 2)
  CV = round(sd(CD8$intensity)/mean, digit = 2)
  jpeg('./CD8_Case1B.jpeg', units="in", width= 45, height=40, res= 600)
  
  p3 <- 
    ggplot(CD8, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab('CD8')+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0,500, 1000,1500, 2000),labels = c(0,500, 1000,1500, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p3)
  dev.off()
  
  
  mean = round(mean(CD20$intensity), digit = 2)
  CV = round(sd(CD20$intensity)/mean, digit = 2)
  
  jpeg('./CD20_Case1B.jpeg', units="in", width= 45, height=40, res= 600)
  p4 <- 
    ggplot(CD20, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab('CD20')+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0,500, 1000,1500, 2000),labels = c(0,500, 1000,1500, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p4)
  dev.off()
  
  
  mean = round(mean(FoxP3$intensity), digit = 2)
  CV = round(sd(FoxP3$intensity)/mean, digit = 2)

  jpeg('./FoxP3_Case1B.jpeg', units="in", width= 45, height=40, res= 600)
  p5 <- 
    ggplot(FoxP3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab('FoxP3')+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0,500, 1000,1500, 2000),labels = c(0,500, 1000,1500, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  plot(p5)
  dev.off()
  
  
  Case1B <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# Case 2
{
  CD3 <- waterfall_new[waterfall_new$Marker == 'CD3' & waterfall_new$Case == 'Case 2',]
  CD4 <- waterfall_new[waterfall_new$Marker == 'CD4' & waterfall_new$Case == 'Case 2',]
  CD8 <- waterfall_new[waterfall_new$Marker == 'CD8' & waterfall_new$Case == 'Case 2',]
  CD20 <- waterfall_new[waterfall_new$Marker == 'CD20' & waterfall_new$Case == 'Case 2',]
  FoxP3 <- waterfall_new[waterfall_new$Marker == 'FoxP3' & waterfall_new$Case == 'Case 2',]
  
  mean = round(mean(CD3$intensity), digit = 2)
  CV = round(sd(CD3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case2_CD3.jpeg', units="in", width= 45, height=40, res= 600)
  
  #p1 <- 
    ggplot(CD3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    #scale_x_continuous(breaks = c(0,500, 1000,1500, 2000),labels = c(0,500, 1000,1500, 2000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  mean = round(mean(CD4$intensity), digit = 2)
  CV = round(sd(CD4$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case2_CD4.jpeg', units="in", width= 45, height=40, res= 600)
  
  ggplot(CD4, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  mean = round(mean(CD8$intensity), digit = 2)
  CV = round(sd(CD8$intensity)/mean, digit = 2)
  jpeg('./waterfall_Case2_CD8.jpeg', units="in", width= 45, height=40, res= 600)
  
  ggplot(CD8, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2500, 5000,7500),labels = c(0,2500, 5000,7500))
    
  dev.off()
  
  mean = round(mean(CD20$intensity), digit = 2)
  CV = round(sd(CD20$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case2_CD20.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(CD20, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2500, 5000,7500),labels = c(0,2500, 5000,7500))
    
  dev.off()
  
  mean = round(mean(FoxP3$intensity), digit = 2)
  CV = round(sd(FoxP3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case2_FoxP3.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(FoxP3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_x_continuous(breaks = c(0,2000, 4000),labels = c(0,2500, 4000))
  
dev.off()
  
  #Case2 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
  dev.off()
}
# Case 3
{
  CD3 <- waterfall_new[waterfall_new$Marker == 'CD3' & waterfall_new$Case == 'Case 3',]
  CD4 <- waterfall_new[waterfall_new$Marker == 'CD4' & waterfall_new$Case == 'Case 3',]
  CD8 <- waterfall_new[waterfall_new$Marker == 'CD8' & waterfall_new$Case == 'Case 3',]
  CD20 <- waterfall_new[waterfall_new$Marker == 'CD20' & waterfall_new$Case == 'Case 3',]
  FoxP3 <- waterfall_new[waterfall_new$Marker == 'FoxP3' & waterfall_new$Case == 'Case 3',]
  
  mean = round(mean(CD3$intensity), digit = 2)
  CV = round(sd(CD3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case3_CD3.jpeg', units="in", width= 45, height=40, res= 600)
    ggplot(CD3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +

    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  mean = round(mean(CD4$intensity), digit = 2)
  CV = round(sd(CD4$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case3_CD4.jpeg', units="in", width= 45, height=40, res= 600)
  
    ggplot(CD4, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +

    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  
  mean = round(mean(CD8$intensity), digit = 2)
  CV = round(sd(CD8$intensity)/mean, digit = 2)
  jpeg('./waterfall_Case3_CD8.jpeg', units="in", width= 45, height=40, res= 600)
  
    ggplot(CD8, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
  
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  
  mean = round(mean(CD20$intensity), digit = 2)
  CV = round(sd(CD20$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case3_CD20.jpeg', units="in", width= 45, height=40, res= 600)
  
    ggplot(CD20, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
  
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)+
    scale_x_continuous(breaks = c(0,1000, 2000),labels = c(0,1000, 2000)) +
    scale_y_continuous(breaks = c(0,1000, 2000, 3000, 4000),labels = c(0,1000, 2000, 3000, 4000))
    
  dev.off()
    
  
  mean = round(mean(FoxP3$intensity), digit = 2)
  CV = round(sd(FoxP3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case3_FoxP3.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(FoxP3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
  
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_x_continuous(breaks = c(0,2000, 4000),labels = c(0,2500, 4000))
  
  dev.off()
  
  Case3 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}
# Case 4
{
  CD3 <- waterfall_new[waterfall_new$Marker == 'CD3' & waterfall_new$Case == 'Case 4',]
  CD4 <- waterfall_new[waterfall_new$Marker == 'CD4' & waterfall_new$Case == 'Case 4',]
  CD8 <- waterfall_new[waterfall_new$Marker == 'CD8' & waterfall_new$Case == 'Case 4',]
  CD20 <- waterfall_new[waterfall_new$Marker == 'CD20' & waterfall_new$Case == 'Case 4',]
  FoxP3 <- waterfall_new[waterfall_new$Marker == 'FoxP3' & waterfall_new$Case == 'Case 4',]
  
  mean = round(mean(CD3$intensity), digit = 2)
  CV = round(sd(CD3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case4_CD3.jpeg', units="in", width= 45, height=40, res= 600)
  
  #p1 <- 
  ggplot(CD3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_y_continuous(breaks = c(0, 2000, 4000, 6000, 7500),labels = c(0, 2000,4000, 6000, 7500)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  mean = round(mean(CD4$intensity), digit = 2)
  CV = round(sd(CD4$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case4_CD4.jpeg', units="in", width= 45, height=40, res= 600)
  
  ggplot(CD4, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  mean = round(mean(CD8$intensity), digit = 2)
  CV = round(sd(CD8$intensity)/mean, digit = 2)
  jpeg('./waterfall_Case4_CD8.jpeg', units="in", width= 45, height=40, res= 600)
  
  ggplot(CD8, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2000, 4000,6000),labels =c(0,2000, 4000,6000)) +
    scale_x_continuous(breaks = c(0,2000, 4000),labels =c(0,2000, 4000))
  
  
  dev.off()
  
  mean = round(mean(CD20$intensity), digit = 2)
  CV = round(sd(CD20$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case4_CD20.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(CD20, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2000, 4000,6000),labels = c(0,2000, 4000,6000))+
    scale_x_continuous(breaks = c(0,1000, 2000),labels =c(0,1000, 2000))
  
  dev.off()
  
  mean = round(mean(FoxP3$intensity), digit = 2)
  CV = round(sd(FoxP3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case4_FoxP3.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(FoxP3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,1000,2000,3000),labels = c(0,1000, 2000, 3000))+
    scale_x_continuous(breaks = c(0,2000,4000),labels = c(0, 2000,4000))
  
  dev.off()
  
  #Case2 <-
  plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
  dev.off()
}
# Case 5

{
  CD3 <- waterfall_new[waterfall_new$Marker == 'CD3' & waterfall_new$Case == 'Case 5',]
  CD4 <- waterfall_new[waterfall_new$Marker == 'CD4' & waterfall_new$Case == 'Case 5',]
  CD8 <- waterfall_new[waterfall_new$Marker == 'CD8' & waterfall_new$Case == 'Case 5',]
  CD20 <- waterfall_new[waterfall_new$Marker == 'CD20' & waterfall_new$Case == 'Case 5',]
  FoxP3 <- waterfall_new[waterfall_new$Marker == 'FoxP3' & waterfall_new$Case == 'Case 5',]
  
  mean = round(mean(CD3$intensity), digit = 2)
  CV = round(sd(CD3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case5_CD3.jpeg', units="in", width= 45, height=40, res= 600)
  
  #p1 <- 
  ggplot(CD3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_y_continuous(breaks = c(0, 2000, 4000, 6000, 8000, 10000),labels = c(0, 2000,4000, 6000, 8000, 10000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE)
  dev.off()
  
  mean = round(mean(CD4$intensity), digit = 2)
  CV = round(sd(CD4$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case5_CD4.jpeg', units="in", width= 45, height=40, res= 600)
  
  ggplot(CD4, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    # annotate(geom="text", x=1700, y= 8500, label= paste('mean = ',mean, sep = ''), size = 5) +
    #annotate(geom="text", x=1700, y= 7600, label= paste('CoV = ',CV, sep = ''), size = 5) +
    scale_x_continuous(breaks = c(0,2000, 4000),labels =c(0,2000, 4000)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge")
  dev.off()
  
  mean = round(mean(CD8$intensity), digit = 2)
  CV = round(sd(CD8$intensity)/mean, digit = 2)
  jpeg('./waterfall_Case5_CD8.jpeg', units="in", width= 45, height=40, res= 600)
  
  ggplot(CD8, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2000, 4000,6000, 8000),labels =c(0,2000, 4000,6000, 8000)) +
    scale_x_continuous(breaks = c(0,2000, 4000),labels =c(0,2000, 4000))
  
  
  dev.off()
  
  mean = round(mean(CD20$intensity), digit = 2)
  CV = round(sd(CD20$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case5_CD20.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(CD20, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2000, 4000,6000, 8000, 10000, 12000),labels =c(0,2000, 4000,6000, 8000, 10000, 12000)) +
    scale_x_continuous(breaks = c(0,1000, 2000),labels =c(0,1000, 2000))
  
  dev.off()
  
  mean = round(mean(FoxP3$intensity), digit = 2)
  CV = round(sd(FoxP3$intensity)/mean, digit = 2)
  
  jpeg('./waterfall_Case5_FoxP3.jpeg', units="in", width= 45, height=40, res= 600)
  ggplot(FoxP3, aes(x = index, y = intensity, color = factor(location))) + 
    theme_bw(base_rect_size = 15) +
    scale_color_manual(name = 'Regions',values = c("#F2F2C6", "#F26969",'#02f527')) +
    xlab(NULL)+
    ylab(NULL) +
    theme(axis.text = element_text(size = 250), panel.grid = element_line(size = 8), axis.title = element_text(size = 280)) +
    geom_bar(stat = 'identity', fill = 'white', width  = 0.1, position = "dodge", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,2000, 4000,6000, 8000, 10000),labels =c(0,2000, 4000,6000, 8000, 10000)) +
    scale_x_continuous(breaks = c(0,1000,2000, 3000),labels =c(0,1000,2000, 3000))
  
  dev.off()
  
  #Case2 <-
  plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
  dev.off()
}
plot(Case1B)
jpeg('./waterfall_2345.jpeg', units="in", width= 500, height=450, res= 600)
plot_grid(Case2, Case3,Case4,Case5, labels = c('','','',''),  ncol = 1, align = 'hv')+
  ylab('Density')
dev.off()

#################################################################################
###### Data for Figs.5F, S10, S11, visualization is done by software Blender ####
#################################################################################

allMarker_case_all <- data.frame(matrix(nrow = 0, ncol  = 10))
names(allMarker_case_all) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location','marker', 'case')
for(case in seq(1,5,1)){
  allMarker <- data.frame(matrix(nrow = 0, ncol  = 9))
  names(allMarker) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location','marker')
  # Change case here
  sub_path <- switch(case, './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  case_name <- switch(case, 'Case 1' , 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  
  # Loop to calculate statistics
  for(label_idx in seq(1,5,1)){
    

    
    #######################
    if(case == 1){
      Current_lab <- switch (label_idx,'CD3','CD4','CD8','CD20','FoxP3')
      Coords_path <- paste(sub_path, Current_lab, Var_suffix, sep = '')
      marker_Coords <- as.data.frame(readMat(Coords_path))/125
      
      InvasiveFile_1 <- readRDS(paste(sub_path, 'Invasive_1.rds', sep = ''))
      
      InvasiveFile_2 <- readRDS(paste(sub_path, 'Invasive_2.rds', sep = ''))
      ##########
      
      TumorFile_1 <- readRDS(paste(sub_path, 'Tumor_1.rds', sep = ''))
      
      
      TumorFile_2 <- readRDS(paste(sub_path, 'Tumor_2.rds', sep = ''))
      ########
      InvasiveFile <- rbind(InvasiveFile_1, InvasiveFile_2)
      TumorFile <- rbind(TumorFile_1, TumorFile_2)
      NormalFile <- readRDS(paste(sub_path, 'Normal.rds', sep = ''))
    }
    
    if(case > 1){
      Current_lab <- switch (label_idx,'CD3','CD4','CD8','CD20','FoxP3')
      Coords_path <- paste(sub_path, Current_lab, '/',Current_lab, Var_suffix, sep = '')
      marker_Coords <- as.data.frame(readMat(Coords_path))/125
      
      InvasiveFile <- readRDS(paste(sub_path, 'Invasive.rds', sep = ''))
      
      TumorFile <- readRDS(paste(sub_path, 'Tumor.rds', sep = ''))
      
      NormalFile <- readRDS(paste(sub_path, 'Normal.rds', sep = ''))
    }
    
    header = c('index',"xstart","xend", "ystart","yend","n","intensity",'location')
    # Scanning the whole space
    window <- 0.15
    
    # Index for each window
    index <- 0
    x_origin <- 0
    y_origin <- 0
    p_all <- data.frame(matrix(nrow=0,ncol=8))
    colnames(p_all) = header
    # Create files for store fit results for current case 
    
    for(y_origin in seq(0, 25-window, window)){
      for(x_origin in seq(0,35-window, window)){
        window_coords <- data.frame(cbind(c(x_origin, x_origin, x_origin+window, x_origin+window), c(y_origin, y_origin + window, y_origin + window, y_origin)))
        
        # Characterize the window
        if(isTRUE(ncol(NormalFile) == 3)){
          normal_polygon_num <- NormalFile[nrow(NormalFile),3]
          normal_count <- 0
          for(num in seq(1, normal_polygon_num)){
            subNormalFile <- NormalFile[NormalFile$label == num, ]
            normal <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], subNormalFile[,1], subNormalFile[,2]))
            subnormal_count <- nrow(data.frame(normal[normal[,1] >0,]))
            normal_count <- normal_count + subnormal_count
          }
        }
        if(isTRUE(ncol(NormalFile) == 2)){
          normal_count <- 0
          normal <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], NormalFile[,1], NormalFile[,2]))
          normal_count <- nrow(data.frame(normal[normal[,1] >0,]))
        }

 
        invasive <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], InvasiveFile[,1], InvasiveFile[,2]))
        invasive_count <- nrow(data.frame(invasive[invasive[,1] >0,]))
    
        tumor <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], TumorFile[,1], TumorFile[,2]))
        tumor_count <- nrow(data.frame(tumor[tumor[,1] >0,]))

        
        if(isTRUE(invasive_count != 0)){
          loc <- 'IF'
        } 
        if(isTRUE(invasive_count == 0 & tumor_count != 0)){
          loc <- 'CT'
        }
        if(isTRUE(invasive_count == 0 & tumor_count == 0 & normal_count != 0)){
          loc <- 'N'
        }
        if(isTRUE(invasive_count == 0 & tumor_count == 0 & normal_count == 0)){
          loc <- 'BG'
        }
        if(isTRUE(loc != 'BG')){
          SubRegion <- marker_Coords[marker_Coords[,1] < (x_origin+window)  & marker_Coords[,1] >= x_origin & marker_Coords[,2] < (y_origin + window) & marker_Coords[,2] >= y_origin,]
          SubRegion <- SubRegion[complete.cases(SubRegion), ]
          if(isTRUE(nrow(SubRegion) != 0)){
            number <- nrow(SubRegion)
            density <- number/(window * window)
          } 
          if(isTRUE(nrow(SubRegion) == 0)){
            number <- 0
            density <- 0
          }
        
          p_sub <- data.frame(cbind(index,x_origin,x_origin+window,y_origin,y_origin+window, number, density, loc))
          names(p_sub) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location')
          p_all <- rbind(p_all, p_sub)
        }
      }
    }
    
    p_all$marker <- Current_lab
    names(p_all) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location','marker')
    names(allMarker) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location','marker')
    allMarker <- rbind(allMarker, p_all)
  }
  allMarker$case <- case_name
  #names(allMarker_case) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location','marker', 'case')
  names(allMarker_case_all) <- c('index',"xstart","xend", "ystart","yend","n","intensity",'location','marker', 'case')
  allMarker_case_all <- rbind(allMarker_case_all, allMarker)
}

write.csv(allMarker_case_all, './3D_plot_data.csv',row.names=FALSE)
allMarker_case_all <- read.csv('./3D_plot_data.csv')

