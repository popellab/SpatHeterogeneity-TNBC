#############################################################################################
#### This script calculate the area for central tumor, invasive front, and normal tissue ####
#############################################################################################

library(R.matlab)
library(jpeg)
library(ggpubr)
library(alphahull)
library(gginnards)
library(ggplot2)
library(spatstat)
library(tidyverse)
library(grid)
library(sp)
library(rgeos)
library(sf)
library(plyr)
library(dplyr)
library(concaveman)
source('MiFunction.R')

####################################
# Case 1############################
####################################


imgage <- readJPEG('./Case1_CD4.jpeg')


mean_matrix <- data.frame(rowMeans(imgage, dims = 2))
mean_matrix <- t(mean_matrix)

indice_x <- row(mean_matrix)[which(mean_matrix <= 0.7)]
indice_y <- col(mean_matrix)[which(mean_matrix  <= 0.7)]

Index <- data.frame(cbind(indice_x,indice_y))/125
ggplot() + geom_point(data = Index, aes(x = Index[,1], y = Index[,2])) 

Tumor_region_1 <- readRDS('./Case 0/Tumor_1.rds')
Invasive_region_1 <- readRDS('./Case 0/Invasive_1.rds')
Invasive_outer_1 <- readRDS('./Case 0/Invasive_outer_1.rds')

Tumor_region_2 <- readRDS('./Case 0/Tumor_2.rds')
Invasive_region_2 <- readRDS('./Case 0/Invasive_2.rds')
Invasive_outer_2 <- readRDS('./Case 0/Invasive_outer_2.rds')


Tumor_region <- rbind(Tumor_region_1, Tumor_region_2)
Invasive_region <- rbind(Invasive_region_1, Invasive_region_2)
#################
Normal_dist <- as.matrix(point.in.polygon(Index[,1], Index[,2], Invasive_region[,1], Invasive_region[,2])) # invasive pts
normal_idx <- which(Normal_dist != 0)
sub1 <- Index[-normal_idx,]

###############
Normal_dist <- as.matrix(point.in.polygon(sub1[,1], sub1[,2], Tumor_region[,1], Tumor_region[,2])) # tumor pts

normal_idx <- which(Normal_dist != 0) 

Normal_complete <- sub1[-normal_idx,]


dat <- as.matrix(Normal_complete)
case4vis <- largeVis(t(dat), dim=2, K = 4, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)
clusters <- largeVis::hdbscan(case4vis, K=4, minPts = 50, verbose = FALSE, threads = 1)

cluster_pts <- as.data.frame(clusters$clusters) 
# Characterize clusters based on their location
#Cluster_dist <- clusterCharacterize(Cluster_dats, TumorFile, InvasiveFile)

Pts <- cbind(Normal_complete,cluster_pts)
Pts <- Pts[is.na(Pts[,3]) == FALSE, ]

#ggplot() + geom_point(data = Pts, aes(x = Pts[,1], y = Pts[,2]))


Pts <- Pts[!(Pts[,2] > 13.63 & Pts[,1] < 12.7143),]

Pts <- Pts[!(Pts[,1]> 15 & Pts[,2] < 3.211),]

Pts <- Pts[!(Pts[,1] > 25 & Pts[,2] < 11.392),]

# Extract first polygon:
Poly1 <- Pts[(Pts[,1] < 15 & Pts[,2] < 15),]
Poly1 <- Pts[!(Pts[,1] > 9.344 & Pts[,2] > 6.164),]
Poly1 <- Poly1[(Poly1[,1] < 15),]
Poly1 <- Poly1[!(Poly1[,1] > 13.744 & Poly1[,2] > 5.279),]
Poly1 <- Poly1[!(Poly1[,1] > 10.8805 & Poly1[,1] < 11.666 & Poly1[,2] > 5.5027),]

names(Poly1) <- c('x', 'y')

Poly1 <- rbind(Poly1[,1:2], Invasive_outer_2)
Poly1_path <- as.data.frame(concaveman(as.matrix(Poly1[,1:2]), concavity = 2))
saveRDS(Poly1_path, './Case 0/normal_1a.rds')


# Extract Second polygon:
Poly2 <- Pts[(Pts[,1] > 11 & Pts[,2] > 10),]
Poly2 <- Poly2[!(Poly2[,2] > 20.3),]
Poly2 <- Poly2[!(Poly2[,1] < 13.543 & Poly2[,2] > 19.7527),]
Poly2 <- Poly2[!(Poly2[,1] < 21.8779 & Poly2[,1] > 20 & Poly2[,2] > 13.175),]
Poly2 <- Poly2[!(Poly2[,1] > 15 &Poly2[,2] > 15),]

names(Poly2) <- c('x', 'y')
Poly2 <- rbind(Poly2[,1:2], Invasive_outer_1)
Poly2_path <- as.data.frame(concaveman(as.matrix(Poly2[,1:2]), concavity = 2))
saveRDS(Poly2_path,'./Case 1/normal_1b.rds')
ggplot() + #geom_point(data = Poly2, aes(x = Poly2[,1], y = Poly2[,2])) +
  geom_polygon(data = Poly2_path, aes(x = Poly2_path[,1], y = Poly2_path[,2]))+
  
  xlim(0,30) +
  ylim(22,0)

# Extract the third polygon:
Poly3 <- Normal_complete[(Normal_complete[,1] > 8 & Normal_complete[,2] < 12.5 & Normal_complete[,1] < 26),]
Poly3 <- Poly3[,1:2]
ggplot() + geom_point(data = Poly3, aes(x = Poly3[,1], y = Poly3[,2]))

Poly3 <- Poly3[!(Poly3[,2] < 5.3727 & Poly3[,1] < 14.6569),]
Poly3 <- Poly3[!(Poly3[,2] < 6.21118 & Poly3[,1] < 13.52818),]
Poly3 <- Poly3[!(Poly3[,2] > 10.4307 & Poly3[,1] > 14.8),]
Poly3 <- Poly3[!(Poly3[,2] > 10 & Poly3[,1] > 20),]
Poly3 <- Poly3[!(Poly3[,2] < 4.2331 & Poly3[,1] < 18.1875),]
Poly3 <- Poly3[!(Poly3[,2] < 4.1718 & Poly3[,1] > 20.2496),]
Poly3 <- Poly3[!(Poly3[,2] < 3.75),]
Poly3 <- Poly3[!(Poly3[,1] > 24.09),]
Poly3 <- Poly3[!(Poly3[,1] > 10.5296 & Poly3[,1] < 14.405 & Poly3[,2] > 9.6307),]
Poly3 <- Poly3[!(Poly3[,1] > 10.5296 & Poly3[,1] < 14.405 & Poly3[,2] > 9.6307),]

Poly3 <- Poly3[!(Poly3[,1] < 9.346),]
Poly3 <- Poly3[!(Poly3[,2] > 11.319),]
Poly3 <- Poly3[!(Poly3[,2] > 7.77 & Poly3[,1] > 22.1875),]
Poly3 <- Poly3[!(Poly3[,2] > 7.77 & Poly3[,1] > 22.1875),]
Poly3 <- Poly3[!(Poly3[,1] > 23.4668 & Poly3[,2] < 5.2469),]
Poly3_path <- as.data.frame(concaveman(as.matrix(Poly3[,1:2]), concavity = 2))


#Poly3_path <- as.data.frame(concaveman(as.matrix(Poly3[,1:2]), concavity = 2))

ggplot() + geom_point(data = Poly3, aes(x = Poly3[,1], y = Poly3[,2])) +
  geom_path(data = Poly3_path, aes(x = Poly3_path[,1], y = Poly3_path[,2]))+
  xlim(0,30) +
  ylim(22,0)


# Combine polygons
Poly1_path$label <- 1
Poly2_path$label <- 2
Poly3_path$label <- 3
Normal_path <- rbind(Poly1_path,Poly2_path,Poly3_path)
saveRDS(Normal_path, './Case 1/Normal_case1.rds')

###################
# Calculate areas #
###################

# Stroma area
Area_1 <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Poly1_path[,1:2])),1))))
Area_2 <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Poly2_path[,1:2])),1))))
Area_3 <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Poly3_path[,1:2])),1))))
Normal_area <- Area_1 + Area_2 + Area_3
# Tumor area
Tumor_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Tumor_region)),1))))
# Invasive area
Invasive_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Invasive_region)),1))))

for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  28.072, -20.824, 0) +
  xlab('x,mm')+
  ylab('y,mm')

# visualization
jpeg('./Case1AB.jpeg', units="in", width=9, height=6.675, res=300)
move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_polygon(data = Poly1_path, aes(x = Poly1_path[,1], y = Poly1_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly2_path, aes(x = Poly2_path[,1], y = Poly2_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly3_path, aes(x = Poly3_path[,1], y = Poly3_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3, show.legend = FALSE) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  theme_bw(base_rect_size = 1) +
  ylim(22,0) + xlim(0,30) +
  xlab(NULL) +
  ylab(NULL) +
  theme(text = element_text(size = 22))+
  labs(fill = "Annotations") +
  #theme(legend.text = element_text(size = 16)) +
  #theme(legend.title = element_text(size = 16)) +
  #ggtitle('Case 1-CD4') +
  #theme(plot.title = element_text(size = 20)) +
  coord_equal() 
dev.off()


####################################
# Case 2############################
####################################
imgage <- readJPEG('./Case2_CD4.jpeg')


mean_matrix <- data.frame(rowMeans(imgage, dims = 2))
mean_matrix <- t(mean_matrix)

indice_x <- row(mean_matrix)[which(mean_matrix <= 0.4)]
indice_y <- col(mean_matrix)[which(mean_matrix  <= 0.4)]

Index <- data.frame(cbind(indice_x,indice_y))/125

###################
# Calculate areas #
###################

Tumor_region <- readRDS('./Case 2/Tumor.rds')
Invasive_region <- readRDS('./Case 2/Invasive.rds')
Invasive_outer <- readRDS('./Case 2/Invasive_outer.rds')

Tumor_area <- SpatialPolygons(list(Polygons(list(Polygon(Tumor_region)),1)))
Tumor_area <- gArea(Tumor_area)

Invasive_area <- SpatialPolygons(list(Polygons(list(Polygon(Invasive_region)),1)))

Invasive_area <- gArea(Invasive_area)
#################
Normal_dist <- as.matrix(point.in.polygon(Index[,1], Index[,2], Invasive_region[,1], Invasive_region[,2])) # invasive pts
normal_idx <- which(Normal_dist != 0)
sub1 <- Index[-normal_idx,]

###############
Normal_dist <- as.matrix(point.in.polygon(sub1[,1], sub1[,2], Tumor_region[,1], Tumor_region[,2])) # tumor pts

normal_idx <- which(Normal_dist != 0) 

Normal_complete <- sub1[-normal_idx,]


# Define polygon for normal region 1
Normal_poly1 <- Normal_complete[Normal_complete$indice_x > 3.2, ]
Normal_poly1 <- Normal_poly1[!(Normal_poly1$indice_x >9 & Normal_poly1$indice_x <25 & Normal_poly1$indice_y > 7.5),]
Normal_poly1 <- Normal_poly1[!(Normal_poly1$indice_y > 14),]
Normal_poly1 <- Normal_poly1[!(Normal_poly1$indice_y > 10 & Normal_poly1$indice_x <4),]
Normal_poly1 <- Normal_poly1[!(Normal_poly1$indice_y > 11.5 & Normal_poly1$indice_x <10),]
names(Normal_poly1) <- c('x','y')
Normal_poly1 <- rbind(Normal_poly1, Invasive_outer)
Normal_poly1 <- Normal_poly1[!(Normal_poly1$y > 8.6 & Normal_poly1$x >6.3 & Normal_poly1$x < 10),]

Normal_path1 <- as.data.frame(concaveman(as.matrix(Normal_poly1), concavity = 2))

ggplot() + geom_point(data = Normal_poly1, aes(x = Normal_poly1[,1], y = Normal_poly1[,2])) +
  geom_path(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2]))


# Define polygon for normal region: 2
Normal_poly2 <- Normal_complete[(Normal_complete$indice_x >13.59 & Normal_complete$indice_x < 22.458 & Normal_complete$indice_y > 9.067),]
Normal_poly2 <- rbind(Normal_poly2,c(16.023,10.4))
Normal_poly2 <- rbind(Normal_poly2,c(17.23,10.057))

Normal_path2 <- as.data.frame(concaveman(as.matrix(Normal_poly2), concavity = 2))

ggplot() + geom_point(data = Normal_poly2, aes(x = Normal_poly2[,1], y = Normal_poly2[,2])) +
  geom_path(data = Normal_path2, aes(x = Normal_path2[,1], y = Normal_path2[,2]))

# Define polygon for normal region: 3
Normal_poly3 <- Normal_complete[(Normal_complete$indice_x < 5 & Normal_complete$indice_y > 12.5),]
Normal_path3 <- as.data.frame(concaveman(as.matrix(Normal_poly3), concavity = 2))

ggplot() + geom_point(data = Normal_path, aes(x = Normal_poly1[,1], y = Normal_poly1[,2])) +
  xlim(0, 30) +
  ylim(20,0)

####### Combine polygon
Normal_path1$label <- 1
Normal_path2$label <- 2
Normal_path3$label <- 3

Normal_path <- rbind(Normal_path1, Normal_path2, Normal_path3)



saveRDS(Normal_path, 'Normal_Case2.rds')
##### Calculate areas

Area1 <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Normal_path1[,1:2])),1))))
Area2 <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Normal_path2[,1:2])),1))))
Area3 <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Normal_path3[,1:2])),1))))

Area <- Area1 + Area2 + Area3
imgage <- readJPEG('./Case2_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  31.872, -18.024, 0) +
  xlab('x,mm')+
  ylab('y,mm')
jpeg('./Case2_Regions.jpeg', units="in", width=15.9, height=9, res=300)
move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_polygon(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2], fill = 'Stroma',group = label), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3,show.legend = FALSE) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  theme_bw() +
  ylim(20,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 60))+
  theme(axis.title = element_text(size = 62)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(plot.title = element_text(size = 40)) +
  coord_equal() 
  #ggtitle('Case 2') 
dev.off()

####################################
# Case 3############################
####################################


imgage <- readJPEG('./Case3_CD4.jpeg')


mean_matrix <- data.frame(rowMeans(imgage, dims = 2))
mean_matrix <- t(mean_matrix)

indice_x <- row(mean_matrix)[which(mean_matrix <= 0.4)]
indice_y <- col(mean_matrix)[which(mean_matrix  <= 0.4)]

Index <- data.frame(cbind(indice_x,indice_y))/125
ggplot() + geom_point(data = Index, aes(x = Index[,1], y = Index[,2])) 

Tumor_region <- readRDS('./Case 3/Tumor.rds')
Invasive_region <- readRDS('./Case 3/Invasive.rds')
Invasive_outer <- readRDS('./Case 3/Invasive_outer.rds')



#################
Normal_dist <- as.matrix(point.in.polygon(Index[,1], Index[,2], Invasive_region[,1], Invasive_region[,2])) # invasive pts
normal_idx <- which(Normal_dist != 0)
sub1 <- Index[-normal_idx,]

###############
Normal_dist <- as.matrix(point.in.polygon(sub1[,1], sub1[,2], Tumor_region[,1], Tumor_region[,2])) # tumor pts

normal_idx <- which(Normal_dist != 0) 

Normal_complete <- sub1[-normal_idx,]
Normal_complete <- Normal_complete[!(Normal_complete[,1] > 10),]
Normal_complete <- Normal_complete[Normal_complete[,1] > 2.29,]
Normal_complete <- Normal_complete[Normal_complete[,1] < 6,]
Normal_complete <- Normal_complete[Normal_complete[,2] < 14.8848,]

names(Normal_complete) <- c('x','y')

Normal_path <- rbind(Normal_complete, Invasive_outer)
Normal_path <- Normal_path[Normal_path[,1] > 2.29,]
Normal_path <- Normal_path[Normal_path[,2] > 7.5,]
Normal_path <- rbind(Normal_path, c(2.5, 10), c(2.38, 11.74))
Normal_path <- rbind(Normal_path, c(2.39, 11.019), c(3.43, 8.28),c(2.238, 13.73),c(2.201, 13.92),c(3.805, 7.903),c(3.059, 8.8557), c(3.2462,8.4807),c(3.3955, 8.25), c(2.4253, 10.067))
Normal_path <- rbind(Normal_path, c(3.4328, 8.1057), c(2.5, 9.8943),c(4.067,7.644), c(3.9179,7.817))

Normal_path <- as.data.frame(concaveman(as.matrix(Normal_path), concavity = 2))
saveRDS(Normal_path, './Case 3/Normal_case3.rds')

###################
# Calculate areas #
###################

# Stroma area
Area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Normal_path[,1:2])),1))))
# Tumor area
Tumor_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Tumor_region)),1))))
# Invasive area
Invasive_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Invasive_region)),1))))

ggplot() + geom_path(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2]))

imgage <- readJPEG('./Case3_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  28.088, -18.76, 0) +
  xlab('x,mm')+
  ylab('y,mm')
contour_Case3 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  #eom_path(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2])) +
  geom_polygon(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2], fill = 'Stroma'), alpha = 0.2) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  theme_bw() +
  ylim(20,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 26))+
  theme(axis.title = element_text(size = 20)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  ggtitle('Case 3-CD4') + 
  theme(plot.title = element_text(size  = 20))

plot(contour_Case3) 


####################################
# Case 4############################
####################################


imgage <- readJPEG('./Case4_CD4.jpeg')


mean_matrix <- data.frame(rowMeans(imgage, dims = 2))
mean_matrix <- t(mean_matrix)

indice_x <- row(mean_matrix)[which(mean_matrix <= 0.7)]
indice_y <- col(mean_matrix)[which(mean_matrix  <= 0.7)]

Index <- data.frame(cbind(indice_x,indice_y))/125
ggplot() + geom_point(data = Index, aes(x = Index[,1], y = Index[,2])) 

Tumor_region <- readRDS('./Case 4/Tumor.rds')
Invasive_region <- readRDS('./Case 4/Invasive.rds')
Invasive_outer <- readRDS('./Case 4/Invasive_outer.rds')



#################
Normal_dist <- as.matrix(point.in.polygon(Index[,1], Index[,2], Invasive_region[,1], Invasive_region[,2])) # invasive pts
normal_idx <- which(Normal_dist != 0)
sub1 <- Index[-normal_idx,]

###############
Normal_dist <- as.matrix(point.in.polygon(sub1[,1], sub1[,2], Tumor_region[,1], Tumor_region[,2])) # tumor pts

normal_idx <- which(Normal_dist != 0) 

Normal_complete <- sub1[-normal_idx,]

dat <- as.matrix(Normal_complete)
case4vis <- largeVis(t(dat), dim=2, K = 4, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)
clusters <- hdbscan(case4vis, K=4, minPts = 60, verbose = FALSE, threads = 1)

cluster_pts <- as.data.frame(clusters$clusters) 
# Characterize clusters based on their location
#Cluster_dist <- clusterCharacterize(Cluster_dats, TumorFile, InvasiveFile)

Pts <- cbind(Normal_complete,cluster_pts)
Pts <- Pts[is.na(Pts[,3]) == FALSE, ]
Pts <- Pts[Pts[,2]< 22.8,]

Pts <- Pts[!(Pts[,1]< 5 & Pts[,2]> 10.98),]
Pts <- Pts[!(Pts[,1]> 20 & Pts[,2]> 21.922),]
Pts <- Pts[!(Pts[,2] < 1.025),]
Pts <- Pts[!(Pts[,1] > 22.674 & Pts[,2] < 2.5),]
Pts <- Pts[!(Pts[,1] > 28.568),]
Pts <- rbind(Pts, c(20.9339, 21.473,1), c(20.7729, 21.586,1),c(19.935, 21.6997,1))
Pts <- rbind(Pts, c(22.4798, 20.948,1), c(22.09339, 21.1331,1),c(21.8035, 21.1614,1))
Pts <- rbind(Pts, c(17.3638, 21.728,1), c(16.7793, 21.983,1),c(8.18, 22.5212,1), c(19.388,21.6997,1 ))
Pts <- rbind(Pts, c(1.6103, 6.60,1), c(1.6425, 7.0538,1))

ggplot() + geom_point(data = Normal_complete, aes(x = Normal_complete[,1], y = Normal_complete[,2]))

Normal_path <- as.data.frame(concaveman(as.matrix(Pts[,1:2]), concavity = 2))
names(Normal_path) <- c('x', 'y')

#Normal_path <- Normal_path[-nrow(Normal_path),]
Normal_path <- rbind(Normal_path[-seq(1,72,1),], Normal_path[seq(1,72,1),])
Normal_path = Normal_path[order(nrow(Normal_path):1),] #invert row order

Normal_path <- rbind(Normal_path, Invasive_outer,Normal_path[1,])

ggplot() + geom_point(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2]))+
  geom_path(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2]))
saveRDS(Normal_path, './Case 4/Normal_case4.rds')

###################
# Calculate areas #
###################

# Stroma area
Area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Normal_path[,1:2])),1))))
# Tumor area
Tumor_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Tumor_region)),1))))
# Invasive area
Invasive_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Invasive_region)),1))))

imgage <- readJPEG('./Case4_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  31.904, -22.904, 0) +
  xlab('x,mm')+
  ylab('y,mm')

contour_Case4 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_polygon(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2], fill = 'Stroma'), alpha = 0.2) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  theme_bw() +
  ylim(24,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 26))+
  theme(axis.title = element_text(size = 20)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  ggtitle('Case 4-CD4') +
  theme(plot.title = element_text(size = 16))

plot(contour_Case4) 


####################################
# Case 5############################
####################################


imgage <- readJPEG('./Case5_CD4.jpeg')


mean_matrix <- data.frame(rowMeans(imgage, dims = 2))
mean_matrix <- t(mean_matrix)

indice_x <- row(mean_matrix)[which(mean_matrix <= 0.7)]
indice_y <- col(mean_matrix)[which(mean_matrix  <= 0.7)]

Index <- data.frame(cbind(indice_x,indice_y))/125
ggplot() + geom_point(data = Index, aes(x = Index[,1], y = Index[,2])) 

Tumor_region <- readRDS('./Case 5/Tumor.rds')
Invasive_region <- readRDS('./Case 5/Invasive.rds')
Invasive_outer <- readRDS('./Case 5/Invasive_outer.rds')



#################
Normal_dist <- as.matrix(point.in.polygon(Index[,1], Index[,2], Invasive_region[,1], Invasive_region[,2])) # invasive pts
normal_idx <- which(Normal_dist != 0)
sub1 <- Index[-normal_idx,]

###############
Normal_dist <- as.matrix(point.in.polygon(sub1[,1], sub1[,2], Tumor_region[,1], Tumor_region[,2])) # tumor pts

normal_idx <- which(Normal_dist != 0) 

Normal_complete <- sub1[-normal_idx,]

ggplot() + geom_point(data = Normal_complete, aes(x = Normal_complete[,1], y = Normal_complete[,2]))

dat <- as.matrix(Normal_complete)
case4vis <- largeVis(t(dat), dim=2, K = 4, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)
clusters <- hdbscan(case4vis, K=4, minPts = 60, verbose = FALSE, threads = 1)
gplot(clusters, dat)
cluster_pts <- as.data.frame(clusters$clusters) 
# Characterize clusters based on their location
#Cluster_dist <- clusterCharacterize(Cluster_dats, TumorFile, InvasiveFile)

Pts <- cbind(Normal_complete,cluster_pts)
Pts <- Pts[is.na(Pts[,3]) == FALSE, ]

#ggplot() + geom_point(data = Pts, aes(x = Pts[,1], y = Pts[,2]))


Pts <- Pts[Pts[,2]< 15,]

Pts <- Pts[(Pts[,1]< 22.5),]
Pts <- Pts[!(Pts[,1]> 18.5 & Pts[,2] < 2.5),]

Pts <- Pts[!(Pts[,2] < 1.025),]
Pts <- Pts[!(Pts[,1] > 20 & Pts[,2] < 4.75),]

Pts <- Pts[!(Pts[,1] > 21.73 & Pts[,2] <8.1624),]
Pts <- Pts[!(Pts[,1] < 4 & Pts[,2] < 5.38),]
names(Pts) <- c('x', 'y')
Pts <- rbind(Pts, c(2.75, 7.25,1))
Normal_path <- rbind(Pts[,1:2], Invasive_outer)

Normal_path <- as.data.frame(concaveman(as.matrix(Normal_path[,1:2]), concavity = 2))
Normal_path <- Normal_path[!(Normal_path[,1]>13.435 & Normal_path[,1] < 13.83755 & Normal_path[,2] > 5),]
Normal_path <- Normal_path[!(Normal_path[,1]>12.5 & Normal_path[,1] < 15 & Normal_path[,2] < 7.5 & Normal_path[,2] > 5),]
saveRDS(Normal_path, './Case 5/Normal_Case5.rds')

ggplot() + geom_point(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2]))+
  geom_path(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2]))

###################
# Calculate areas #
###################

# Stroma area
Area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Normal_path[,1:2])),1))))
# Tumor area
Tumor_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Tumor_region)),1))))
# Invasive area
Invasive_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Invasive_region)),1))))

for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  25.232, -21.424, 0) +
  xlab('x,mm')+
  ylab('y,mm')

jpeg('./Case 5/Annotations.jpeg', units="in", width=12, height=8, res=300)
contour_Case5 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_polygon(data = Normal_path, aes(x = Normal_path[,1], y = Normal_path[,2], fill = 'Stroma'), alpha = 0.2) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  theme_bw() +
  ylim(22,0) + xlim(0,26) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 24))+
  theme(axis.title = element_text(size = 20)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(plot.title = element_text(size = 20)) +
  ggtitle('Case 5-CD4') 
dev.off()
plot(contour_Case5) 
