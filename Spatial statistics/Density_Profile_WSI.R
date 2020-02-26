##################################################################################
#### This script is used to generate density profile with SE (Figs. 5C, S7) ######
##################################################################################

library(ggplot2)
library(R.matlab)
library(cowplot) # panel_border function

# Calculate SE for each marker across each slide
rawFile <- read.csv('./3D_plot_data_with0.csv')

SE <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(SE) <- c('Marker', 'Case', 'SE', 'SD')
for(case in seq(1,6,1)){
  case_name <- switch (case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  for(marker in seq(1, 5)){
    marker_name <- switch(marker, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
    for(loc in seq(1,3)){
      loc_name <- switch(loc, 'N', 'IF', 'CT')
      statFile <- rawFile[rawFile$case == case_name & rawFile$marker == marker_name & rawFile$location == loc_name,]

      num_window <- nrow(statFile)
      sd <- sd(statFile$intensity)
      std.error <- sd/sqrt(num_window)
      SE_sub <- cbind(marker_name, case_name, loc_name, std.error, sd)
      colnames(SE_sub) <- c('Marker', 'Case','Location', 'SE', 'SD')
      SE <- rbind(SE, SE_sub)
    }
  }
}

# Calculate confidence interval
z_value <- 1.96 # 95% confidence interval
window_size <- 0.15

Profile_collect <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(Profile_collect) <- c('Section', 'Area', 'Number', 'Density', 'CI_low', 'CI_high', 'Marker', 'Case')
# Profile path
Profile_prefix <- 'Final_diction_'
for(case in seq(1,6)){
  case_path <- switch (case,'./Case_0/', './Case_1/', './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  case_name <- switch (case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  stripe_name <- switch (case,'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  for(marker in seq(1,5)){
    marker_name <- switch(marker, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
    # get profile
    Profile <- readRDS(paste(case_path, 'Final_diction/', Profile_prefix, marker_name, '.rds', sep = ''))
    # get SE value for current case
    SE_cur <- SE[SE$Marker == marker_name & SE$Case == case_name,]
    SE_cur$SD <- as.numeric(as.character(SE_cur$SD))
    
    SD_N <- SE_cur[SE_cur$Location == 'N',]$SD
    SD_IF <- SE_cur[SE_cur$Location == 'IF',]$SD
    SD_CT <- SE_cur[SE_cur$Location == 'CT',]$SD
    
    # divide the Profile according to tissue
    
    # normal
    Profile_N <- Profile[Profile$Section < -0.5,]
    
    # invasive
    Profile_IF <- Profile[Profile$Section >= -0.5 & Profile$Section <= 0.5,]
    
    # tumor
    Profile_CT <- Profile[Profile$Section > 0.5,]
    # calculate CI for N 
    
    CI_low <- Profile_N$Density - z_value*SD_N/(sqrt(Profile_N$Area/(window_size*window_size)))
    
    CI_high <- Profile_N$Density + z_value*SD_N/(sqrt(Profile_N$Area/(window_size*window_size)))
    
    CI_N <- cbind(CI_low, CI_high)
    # calculate CI for IF 
    
    CI_low <- Profile_IF$Density - z_value*SD_IF/(sqrt(Profile_IF$Area/(window_size*window_size)))
    CI_high <- Profile_IF$Density + z_value*SD_IF/(sqrt(Profile_IF$Area/(window_size*window_size)))
    CI_IF <- cbind(CI_low, CI_high)
    
    # calculate CI for CT 
    
    CI_low <- Profile_CT$Density - z_value*SD_CT/(sqrt(Profile_CT$Area/(window_size*window_size)))
    CI_high <- Profile_CT$Density + z_value*SD_CT/(sqrt(Profile_CT$Area/(window_size*window_size)))
    CI_CT <- cbind(CI_low, CI_high)
    
    CI<- rbind(CI_N, CI_IF, CI_CT)
    
    Profile_CI <- cbind(Profile, CI,marker_name, stripe_name)

      
    colnames(Profile_CI) <- c('Section', 'Area', 'Number', 'Density', 'CI_low', 'CI_high', 'Marker', 'Case')
    
    Profile_collect <- rbind(Profile_collect, Profile_CI) 
  }
}

Profile_collect[Profile_collect$CI_low < 0,]$CI_low <- 0

Discard_sub <- data.frame(matrix(nrow = 0, ncol = 8))
for(pts in seq(1, nrow(Profile_collect))){
  subProfile <- Profile_collect[pts,]
  if(isTRUE(subProfile$Density*1.8 < subProfile$CI_high | subProfile$Density*0.2 > subProfile$CI_low)){
    Discard_sub <- rbind(Discard_sub,subProfile)
  }
}

########################
# Plotting for each case
########################

# Density plot A
{
  CD3 <- Profile_collect[Profile_collect$Mark == 'CD3' & Profile_collect$Case == 'Case 1A',]
  CD4 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 1A',]
  CD8 <- Profile_collect[Profile_collect$Mark == 'CD8' & Profile_collect$Case == 'Case 1A',]
  CD20 <- Profile_collect[Profile_collect$Mark == 'CD20' & Profile_collect$Case == 'Case 1A',]
  FoxP3 <- Profile_collect[Profile_collect$Mark == 'FoxP3' & Profile_collect$Case == 'Case 1A',]
  
  Discard_CD3 <- Discard_sub[Discard_sub$Mark == 'CD3' & Discard_sub$Case == 'Case 1A',]
  Discard_CD4 <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 1A',]
  Discard_CD8 <- Discard_sub[Discard_sub$Mark == 'CD8' & Discard_sub$Case == 'Case 1A',]
  Discard_CD20 <- Discard_sub[Discard_sub$Mark == 'CD20' & Discard_sub$Case == 'Case 1A',]
  Discard_FoxP3 <- Discard_sub[Discard_sub$Mark == 'FoxP3' & Discard_sub$Case == 'Case 1A',]

  p1 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL) +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD4, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD4, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL) +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD8, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD8, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD8, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL) +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD20, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD20, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD20, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL) +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
 p5 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = FoxP3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_FoxP3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = FoxP3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab('Case 1A') +
    xlab(NULL) +
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank())
  
  Case1A <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# Density plot B
{
  CD3 <- Profile_collect[Profile_collect$Mark == 'CD3' & Profile_collect$Case == 'Case 1B',]
  CD4 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 1B',]
  CD8 <- Profile_collect[Profile_collect$Mark == 'CD8' & Profile_collect$Case == 'Case 1B',]
  CD20 <- Profile_collect[Profile_collect$Mark == 'CD20' & Profile_collect$Case == 'Case 1B',]
  FoxP3 <- Profile_collect[Profile_collect$Mark == 'FoxP3' & Profile_collect$Case == 'Case 1B',]
  Discard_CD3 <- Discard_sub[Discard_sub$Mark == 'CD3' & Discard_sub$Case == 'Case 1B',]
  Discard_CD4 <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 1B',]
  Discard_CD8 <- Discard_sub[Discard_sub$Mark == 'CD8' & Discard_sub$Case == 'Case 1B',]
  Discard_CD20 <- Discard_sub[Discard_sub$Mark == 'CD20' & Discard_sub$Case == 'Case 1B',]
  Discard_FoxP3 <- Discard_sub[Discard_sub$Mark == 'FoxP3' & Discard_sub$Case == 'Case 1B',]
  p1 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    
    ylab(NULL) +
    xlab('CD3') +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD4, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    
    geom_line(data = CD4, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab('CD4') +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD8, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD8, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD8, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab('CD8') +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD20, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD20, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD20, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab('CD20') +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = FoxP3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_FoxP3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = FoxP3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab('Case 1B') +
    xlab('FoxP3') +
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank())
  
  Case1B <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}


# Density plot 2
{
  CD3 <- Profile_collect[Profile_collect$Mark == 'CD3' & Profile_collect$Case == 'Case 2',]
  CD4 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 2',]
  CD8 <- Profile_collect[Profile_collect$Mark == 'CD8' & Profile_collect$Case == 'Case 2',]
  CD20 <- Profile_collect[Profile_collect$Mark == 'CD20' & Profile_collect$Case == 'Case 2',]
  FoxP3 <- Profile_collect[Profile_collect$Mark == 'FoxP3' & Profile_collect$Case == 'Case 2',]
  Discard_CD3 <- Discard_sub[Discard_sub$Mark == 'CD3' & Discard_sub$Case == 'Case 2',]
  Discard_CD4 <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 2',]
  Discard_CD8 <- Discard_sub[Discard_sub$Mark == 'CD8' & Discard_sub$Case == 'Case 2',]
  Discard_CD20 <- Discard_sub[Discard_sub$Mark == 'CD20' & Discard_sub$Case == 'Case 2',]
  Discard_FoxP3 <- Discard_sub[Discard_sub$Mark == 'FoxP3' & Discard_sub$Case == 'Case 2',]
  p1 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD4, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    
    geom_line(data = CD4, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD8, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD8, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD8, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD20, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD20, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD20, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = FoxP3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_FoxP3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = FoxP3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab('Case 2') +
    xlab(NULL)+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank())
  
  Case2 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# Density plot 3
{
  CD3 <- Profile_collect[Profile_collect$Mark == 'CD3' & Profile_collect$Case == 'Case 3',]
  CD4 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 3',]
  CD8 <- Profile_collect[Profile_collect$Mark == 'CD8' & Profile_collect$Case == 'Case 3',]
  CD20 <- Profile_collect[Profile_collect$Mark == 'CD20' & Profile_collect$Case == 'Case 3',]
  FoxP3 <- Profile_collect[Profile_collect$Mark == 'FoxP3' & Profile_collect$Case == 'Case 3',]
  Discard_CD3 <- Discard_sub[Discard_sub$Mark == 'CD3' & Discard_sub$Case == 'Case 3',]
  Discard_CD4 <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 3',]
  Discard_CD8 <- Discard_sub[Discard_sub$Mark == 'CD8' & Discard_sub$Case == 'Case 3',]
  Discard_CD20 <- Discard_sub[Discard_sub$Mark == 'CD20' & Discard_sub$Case == 'Case 3',]
  Discard_FoxP3 <- Discard_sub[Discard_sub$Mark == 'FoxP3' & Discard_sub$Case == 'Case 3',]
  p1 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD4, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    
    geom_line(data = CD4, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD8, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD8, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD8, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD20, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD20, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD20, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = FoxP3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_FoxP3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = FoxP3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab('Case 3') +
    xlab(NULL)+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank())
  
  Case3 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# Density plot 4

{
  CD3 <- Profile_collect[Profile_collect$Mark == 'CD3' & Profile_collect$Case == 'Case 4',]
  CD4 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 4',]
  CD8 <- Profile_collect[Profile_collect$Mark == 'CD8' & Profile_collect$Case == 'Case 4',]
  CD20 <- Profile_collect[Profile_collect$Mark == 'CD20' & Profile_collect$Case == 'Case 4',]
  FoxP3 <- Profile_collect[Profile_collect$Mark == 'FoxP3' & Profile_collect$Case == 'Case 4',]
  Discard_CD3 <- Discard_sub[Discard_sub$Mark == 'CD3' & Discard_sub$Case == 'Case 4',]
  Discard_CD4 <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 4',]
  Discard_CD8 <- Discard_sub[Discard_sub$Mark == 'CD8' & Discard_sub$Case == 'Case 4',]
  Discard_CD20 <- Discard_sub[Discard_sub$Mark == 'CD20' & Discard_sub$Case == 'Case 4',]
  Discard_FoxP3 <- Discard_sub[Discard_sub$Mark == 'FoxP3' & Discard_sub$Case == 'Case 4',]
  p1 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD4, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    
    geom_line(data = CD4, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD8, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD8, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD8, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD20, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD20, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD20, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab(NULL)+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = FoxP3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_FoxP3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = FoxP3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab('Case 4') +
    xlab(NULL)+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank())
  
  Case4 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# Density plot 5
{
  CD3 <- Profile_collect[Profile_collect$Mark == 'CD3' & Profile_collect$Case == 'Case 5',]
  CD4 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 5',]
  CD8 <- Profile_collect[Profile_collect$Mark == 'CD8' & Profile_collect$Case == 'Case 5',]
  CD20 <- Profile_collect[Profile_collect$Mark == 'CD20' & Profile_collect$Case == 'Case 5',]
  FoxP3 <- Profile_collect[Profile_collect$Mark == 'FoxP3' & Profile_collect$Case == 'Case 5',]
  Discard_CD3 <- Discard_sub[Discard_sub$Mark == 'CD3' & Discard_sub$Case == 'Case 5',]
  Discard_CD4 <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 5',]
  Discard_CD8 <- Discard_sub[Discard_sub$Mark == 'CD8' & Discard_sub$Case == 'Case 5',]
  Discard_CD20 <- Discard_sub[Discard_sub$Mark == 'CD20' & Discard_sub$Case == 'Case 5',]
  Discard_FoxP3 <- Discard_sub[Discard_sub$Mark == 'FoxP3' & Discard_sub$Case == 'Case 5',]
  p1 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    
    ylab(NULL) +
    xlab('CD3')+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p2 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD4, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    
    geom_line(data = CD4, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab('CD4')+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p3 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD8, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD8, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD8, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab('CD8')+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  
  p4 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = CD20, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_CD20, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = CD20, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab(NULL) +
    xlab('CD20')+
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
  
  p5 <- 
    ggplot() + 
    theme_bw(base_rect_size = 1) +
    geom_ribbon(data = FoxP3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
    geom_point(data = Discard_FoxP3, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
    geom_line(data = FoxP3, aes(Section,Density),  show.legend = FALSE, size = 1) +
    ylab('Case 5') +
    xlab('FoxP3')+
    scale_y_continuous(position = 'right', sec.axis = dup_axis()) + 
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color = "blue") +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5), axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), axis.title.y.left = element_blank())
  
  Case5 <-
    plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}


jpeg('./Density_profile_2345.jpeg', units="in", width=15, height=10, res=300)
plot_grid(Case2, Case3,Case4, Case5, labels = c('','','',''),  ncol = 1, align = 'hv')+
  ylab('Density')
dev.off()





