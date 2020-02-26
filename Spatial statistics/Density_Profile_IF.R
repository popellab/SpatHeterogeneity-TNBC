#######################################################################
##### This script examine the heterogeneiry along Invasive front ######
#######################################################################

# import packages
library(R.matlab)
library(jpeg)
library(ggplot2)
library(ggforce)
library(sp)
library(raster)
library(cowplot)
library(concaveman)
library(rgeos)
library(SDraw)


#####################################################
# Generate Voronoi tesselation along Invasive front #
#####################################################


###################################
# Generate sample points along IF #
###################################

# Read spatial patterns
ref_line <- readRDS('./Case_0/reference_line.rds') # Case_0 the upper left one
ref_line <- SpatialLines(list(Lines(Line(ref_line), ID="Case0_ref")))
numOfPoints  <-  lineLength(ref_line) / 0.2

sample_pts <- data.frame(spsample(ref_line, n = numOfPoints, type = "regular"))

invasive_front <- readRDS('./Case_0/Invasive_1a.rds')



#  Get pixels in invasive front
sample_idx <- point.in.polygon(sample_pts[,1],sample_pts[,2], invasive_front[,1], invasive_front[,2])

# Count numbers
sample_pts <- cbind(sample_pts, data.frame(sample_idx))
sample_pts <- sample_pts[sample_pts[,3] != 0,-3]

# reorganize, no need for case 3, 4, 5
#sample_pts_sec <- rbind(sample_pts[seq(179,183),],sample_pts[seq(1,178),]) # case 2
#
sample_pts_sec <- rbind(sample_pts[seq(67, 82),],sample_pts[seq(1,66),]) # case 1
#sample_pts_sec <- rbind(sample_pts[seq(65, 91),],sample_pts[seq(1,64),]) # case 0

sample_pts_sec <- sample_pts
sample_pts_sec$label <- seq(1, nrow(sample_pts_sec))



dat <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dat) <- c('x','y')
for(i in seq(0.1, 28.072, 0.1)){
  for(j in seq(0.1, 20.824,0.1)){
    dat_sub <- data.frame(cbind(i+0.05,j+0.05))
    colnames(dat_sub) <- c('x','y')
    dat <- rbind(dat, dat_sub)
  }
}

dot_idx <- point.in.polygon(dat[,1],dat[,2], invasive_front[,1], invasive_front[,2])
dat_pts <- cbind(dat, data.frame(dot_idx))
dat_pts <- dat_pts[dat_pts[,3] != 0,-3]

##################################
# Plot sampled points, (Fig. 4D) #
##################################

jpeg('./Case1A_invasive.jpeg', units="in", width=4, height=3, res=300)
ggplot() + 
  scale_y_reverse() +
  theme_bw(base_rect_size = 2) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_point(data = dat_pts, aes(x = dat_pts[,1], y = dat_pts[,2]), color = 'red', size = 1, alpha = 0.8) +
  geom_point(data = sample_pts_sec, aes(x = sample_pts_sec[,1], y = sample_pts_sec[,2])) +
  geom_polygon(data = invasive_front, aes(x = invasive_front[,1], y = invasive_front[,2]), alpha = 0.5) +
  coord_equal()
dev.off()

#########################################
# Generate Voronoi tesselation along IF #
#########################################

Pixel_coords <- read.csv('./Case_5/Pixel_coords.csv')/125
Invasive_pixels <- point.in.polygon(Pixel_coords[,1],Pixel_coords[,2], invasive_front[,1], invasive_front[,2])
Pixel_coords <- cbind(Pixel_coords, data.frame(Invasive_pixels))
Pixel_coords <- Pixel_coords[Pixel_coords[,3] != 0,-3]



Assigned_df <- data.frame(matrix(nrow = 0, ncol = 3))
names(Assigned_df) <- c('x','y','label')


for(pts_id in seq(1, nrow(Pixel_coords))){
  Pts <- Pixel_coords[pts_id, ]
  Distance_profile <- data.frame(pointDistance(Pts, sample_pts_sec[,1:2], lonlat = FALSE))# calculate the distance of all pixel point within the invasive front to the current sample point
  sample_pts_wLocation <- cbind(sample_pts_sec, Distance_profile)
  names(sample_pts_wLocation) <- c('x','y','label', 'distance')
  Optimized_label <- which.min(sample_pts_wLocation$distance)
  
  # assign the current pixel coordinate with correct label
  Assign_coords <- cbind(Pts, Optimized_label)
  Assigned_df <- rbind(Assigned_df, Assign_coords)
  #print(pts_id)
}
names(Assigned_df) <- c('x','y','label')

# save data to file
Assigned_df$label <- as.numeric(Assigned_df$label)

saveRDS(Assigned_df, './Case_5/segment_invasive_CD3.rds')

Assigned_df <- readRDS('./Case_0/segment_invasive.rds')/125

#########################
# Plot partition result #
#########################

jpeg('./Case1A_section_IF.jpeg', units="in", width=4, height=3, res=300)
ggplot() + 
  scale_y_reverse()+
  theme_bw(base_rect_size = 2) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_point(data = Assigned_df, aes(x = Assigned_df[,1], y = Assigned_df[,2], color = factor(Assigned_df$label)), size = 1, alpha = 0.8, show.legend = FALSE) +
  #xlim(10,25) +#0, 15
  #ylim(22,9) + #15, 2
  coord_equal()
dev.off()

###################################################################
######### Calculate SE for reigon densities within IF  ############
###################################################################

#####################################
# Calculate SE for IF for each case #
#####################################

rawFile <- read.csv('./3D_plot_data_with0.csv')

SE <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(SE) <- c('Marker', 'Case', 'SE', 'SD')
for(case in seq(1,6,1)){
  case_name <- switch (case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  Case_name <- switch (case,'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  for(marker in seq(1, 5)){
    marker_name <- switch(marker, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
    for(loc in seq(1,3)){
      loc_name <- switch(loc, 'N', 'IF', 'CT')
      statFile <- rawFile[rawFile$case == case_name & rawFile$marker == marker_name & rawFile$location == loc_name,]
      
      num_window <- nrow(statFile)
      sd <- sd(statFile$intensity)
      std.error <- sd/sqrt(num_window)
      SE_sub <- cbind(marker_name, Case_name, loc_name, std.error, sd)
      colnames(SE_sub) <- c('Marker', 'Case','Location', 'SE', 'SD')
      SE <- rbind(SE, SE_sub)
    }
  }
}

###################################################
# Calculate density and associated CI based on SE #
###################################################

Density_profile <- data.frame(matrix(nrow = 0, ncol = 6))
names(Density_profile) <- c('section_num', 'Area', 'Number','density', 'marker', 'case')
shape_descriptor <- data.frame(matrix(nrow = 0, ncol = 3))
names(shape_descriptor) <- c('x','y','label')

# Read cell points
for(label_idx in seq(1, 5, 1)){
  marker <- switch(label_idx, 'CD3/', 'CD4/', 'CD8/', '/CD20/', 'FoxP3/')
  marker_name <- switch(label_idx, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
  
  for(idx in seq(1,6,1)){
    

    
    Case <- switch(idx, './Case_0/', './Case_1/', './Case_2/', './Case_3/', './Case_4/', './Case_5/')
    Cell_pts <- data.frame(readMat(paste(Case,  marker, marker_name, '_Points.mat', sep = '')))/125
    
    ######### for naming
    case_name <- switch(idx, 'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
    
    # create data frames
    Assigned_df<- readRDS(paste(Case,'segment_invasive.rds',sep = '' ))
    Assigned_df <- Assigned_df[order(Assigned_df$label),] 

    names(Assigned_df) <- c('x','y','label')

    # Extract polygon for each section

    for(label_num in seq(1, Assigned_df[nrow(Assigned_df),3])){
  
      sec_pts <- Assigned_df[Assigned_df$label == label_num,]
      ctr_pts <- as.data.frame(concaveman(as.matrix(sec_pts[,1:2]), concavity = 2))
      Cell_in_Sec <- point.in.polygon(Cell_pts[,1],Cell_pts[,2], ctr_pts[,1], ctr_pts[,2])
      Cell_in_Sec <- cbind(Cell_pts, Cell_in_Sec) 
  

      # Calculate numbers of cell within polygon
      Number <- nrow(Cell_in_Sec[Cell_in_Sec[,3] != 0,])
  
      # Calculate section area
      Area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(ctr_pts)),1))))
  
      # Calculate density
      Density <- Number/Area
  
      ############################
      ### prepare data for plot###
      ############################
      density_profile <- cbind(label_num, Area, Number, Density, marker_name, case_name)
      names(density_profile) <- c('section_num','Area' , 'Number' ,'density', 'marker', 'case')
      Density_profile <- rbind(Density_profile, density_profile)

      # shape
  
      #shape <- cbind(ctr_pts, label_num)
      #names(shape) <- c('x', 'y', 'label')
      #shape_descriptor <- rbind(shape_descriptor, shape)
      }
  }
}

############### This section calculate CI
Density_profile$Density <- as.numeric(as.character(Density_profile$Density))
Density_profile$label_num <- as.numeric(as.character(Density_profile$label_num))
Density_profile$Area <- as.numeric(as.character(Density_profile$Area))
Density_profile$Number <- as.numeric(as.character(Density_profile$Number))

saveRDS(Density_profile,'./Invasive_densityProfile_combined.rds')

###############
Density_profile <- readRDS('./Invasive_densityProfile_combined.rds')

z_value <- 1.96 # 95% confidence interval
window_size <- 0.15
#Density_profile <- readRDS('./Invasive_densityProfile_combined.rds')


Profile_collect_sub <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(Profile_collect_sub) <- c('Section', 'Area','Number' ,'Density','Marker', 'Case', 'CI_low', 'CI_high')

Profile_collect <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(Profile_collect) <- c('Section', 'Area', 'Number','Density','Marker', 'Case', 'CI_low', 'CI_high')
# Profile path
Profile_prefix <- 'Final_diction_'
for(case in seq(1,6)){
  case_path <- switch (case,'./Case_0/', './Case_1/', './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  case_name <- switch (case,'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  for(marker in seq(1,5)){
    marker_name <- switch(marker, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
    # get profile
    Profile <- Density_profile[Density_profile$case_name == case_name & Density_profile$marker_name == marker_name,]
    # get SE value for current case
    SE_cur <- SE[SE$Marker == marker_name & SE$Case == case_name,]
    SE_cur$SD <- as.numeric(as.character(SE_cur$SD))
    
    SD_IF <- SE_cur[SE_cur$Location == 'IF',]$SD
    # divide the Profile according to tissue
    
    # calculate CI for IF 
    
    CI_low <- Profile$Density - z_value*SD_IF/(sqrt(Profile$Area/(window_size*window_size)))
    CI_high <- Profile$Density + z_value*SD_IF/(sqrt(Profile$Area/(window_size*window_size)))
    CI_IF <- cbind(CI_low, CI_high)
    
    Profile_collect_sub <- cbind(Profile, CI_IF) 
    Profile_collect <- rbind(Profile_collect, Profile_collect_sub)
    
  }
  #Profile_collect_sub <- cbind(Profile, CI_collect) 
}

colnames(Profile_collect) <- c('Section', 'Area', 'Number','Density','Marker', 'Case', 'CI_low', 'CI_high')


Profile_collect[Profile_collect$CI_low < 0,]$CI_low <- 0

Discard_sub <- data.frame(matrix(nrow = 0, ncol = 9))
for(pts in seq(1, nrow(Profile_collect))){
  subProfile <- Profile_collect[pts,]
  if(isTRUE(subProfile$Density*1.8 < subProfile$CI_high | subProfile$Density*0.2 > subProfile$CI_low)){
    Discard_sub <- rbind(Discard_sub,subProfile)
  }
}


#saveRDS(Profile_collect,  './Invasive_densityProfile_combined.rds')

Profile_collect <- readRDS('./Invasive_densityProfile_combined.rds')

CD4_Case1A <- Profile_collect[Profile_collect$Marker == 'CD4' & Profile_collect$Case == 'Case 1A',]
CD4_Case1B <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 1B',]

CD4_Case2 <- Profile_collect[Profile_collect$Marker == 'CD4' & Profile_collect$Case == 'Case 2',]
CD4_Case3 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 3',]

CD4_Case4 <- Profile_collect[Profile_collect$Marker == 'CD4' & Profile_collect$Case == 'Case 4',]
CD4_Case5 <- Profile_collect[Profile_collect$Mark == 'CD4' & Profile_collect$Case == 'Case 5',]



CD4_Case1A_Discard <- Discard_sub[Discard_sub$Marker == 'CD4' & Discard_sub$Case == 'Case 1A',]
CD4_Case1B_Discard <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 1B',]

CD4_Case2_Discard <- Discard_sub[Discard_sub$Marker == 'CD4' & Discard_sub$Case == 'Case 2',]
CD4_Case3_Discard <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 3',]

CD4_Case4_Discard <- Discard_sub[Discard_sub$Marker == 'CD4' & Discard_sub$Case == 'Case 4',]
CD4_Case5_Discard <- Discard_sub[Discard_sub$Mark == 'CD4' & Discard_sub$Case == 'Case 5',]

p1 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_ribbon(data = CD4_Case1A, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
  geom_point(data = CD4_Case1A_Discard, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
  geom_line(data = CD4_Case1A, aes(Section,Density),  show.legend = FALSE, size = 0.5) +
  ylab(NULL) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

p2 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_ribbon(data = CD4_Case1B, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
  geom_point(data = CD4_Case1B_Discard, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
  geom_line(data = CD4_Case1B, aes(Section,Density),  show.legend = FALSE, size = 0.5) +
  ylab(NULL) +
  xlab('x,mm') +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

p3 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_ribbon(data = CD4_Case2, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
  geom_point(data = CD4_Case2_Discard, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
  geom_line(data = CD4_Case2, aes(Section,Density),  show.legend = FALSE, size = 0.5) +
  ylab(NULL) +
  xlab(NULL) +  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

p4 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_ribbon(data = CD4_Case3, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
  geom_point(data = CD4_Case3_Discard, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
  geom_line(data = CD4_Case3, aes(Section,Density),  show.legend = FALSE, size = 0.5) +
  ylab(NULL) +
  xlab(NULL) +  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

p5 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_ribbon(data = CD4_Case4, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
  geom_point(data = CD4_Case4_Discard, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
  geom_line(data = CD4_Case4, aes(Section,Density),  show.legend = FALSE, size = 0.5) +
  ylab(NULL) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

p6 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_ribbon(data = CD4_Case5, aes(Section, ymin = CI_low,ymax = CI_high),  fill = "grey70", show.legend = FALSE, size = 20) +
  geom_point(data = CD4_Case5_Discard, aes(Section, Density),  color = "red", show.legend = FALSE, size = 5) +
  geom_line(data = CD4_Case5, aes(Section,Density),  show.legend = FALSE, size = 0.5) +
  ylab(NULL) +
  xlab('x,mm') +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))


##########
# Fig 5E #
##########

jpeg('./IF_Case1AB.jpg', units="in", width=9, height=12, res=300)
plot_grid(p1, p2, labels = c('', ''), label_size = 30, ncol = 1)
dev.off()

##########
# Fig S9 #
##########

jpeg('./IF_Case2345.jpg', units="in", width=9, height=8, res=300)
   plot_grid(p3, p4, p5, p6, labels = c('', '','',''), label_size = 30, ncol = 1)
dev.off()



#jpeg( './Case_3/IF_section.jpeg',units="in", width=12, height= 8.015, res=1200)
#ggplot() + geom_polygon(data = Assigned_df, aes(x = x, y =y, color =factor(label)), show.legend = FALSE) +
#  theme_bw() + 
  #geom_point(data = Cell_pts, aes(x = Cell_pts[,1], y = Cell_pts[,2])) +
#  theme(axis.text = element_text(size = 32), axis.title = element_text(size = 30)) +
#  xlab('x, mm') +
#  ylab('y, mm') +
#  xlim(0,28.088) +
#  ylim(25, 0)
#dev.off()



# Calculat QCoD (Numbers on plot)
# Case 1a (Case 0)
CoV_collect <- data.frame(matrix(nrow = 0, ncol = 2))
names(CoV_collect) <- c('num', 'CoV')
for(case_id in seq(0,5,1)){
  filepath <- paste('./Case_',case_id, '/','Invasive_densityProfile.rds', sep = '')
  Density_profile <- readRDS(filepath)
  sd <- sd(Density_profile$Density)
  mean <- mean(Density_profile$Density)
  CoV <- sd/mean # intra-invasive frontal heterogeneity
  CoV <- cbind(case_id, CoV)
  names(CoV) <- c('num', 'CoV')
  CoV_collect <- rbind(CoV_collect, CoV)
}

sd <- sd(CoV_collect$CoV)
mean <- mean(CoV_collect$CoV)
CoV_inter <- sd/mean # inter-invasive frontal heterogeneity
print(paste('Intra-invasive frontal heterogeneity:', min(CoV_collect$CoV), '-', max(CoV_collect$CoV)))

print(paste('Inter-invasive frontal heterogeneity :', CoV_inter))

