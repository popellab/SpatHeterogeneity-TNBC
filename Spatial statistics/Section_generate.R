####################################################################################
# This script is used to generate equidistant sections for Case 1 (Example script) #
# The generation of each section requires fine-tuning #
####################################################################################

library(R.matlab)
library(jpeg)
library(alphahull)
library(ggplot2)
library(ggforce)
library(spatstat)
library(tidyverse)
library(grid)
library(sp)
library(rgeos)
library(sf)
library(plyr)
library(concaveman)
library(largeVis)
source('MiFunction.R')
# for formatting grids

# The pixel-coordinate file is required as input
Indice1 <- as.data.frame(readRDS('./Case 1/Indice_case1_1.rds'))/125 
Indice2 <- as.data.frame(readRDS('./Case 1/Indice_case1_2.rds'))/125

Invasive_1 <- as.data.frame(readRDS('./Case 1/Invasive_1.rds'))
Invasive_2 <- as.data.frame(readRDS('./Case 1/Invasive_2.rds'))

# Case 1B
sp_tumor1 <- SpatialPoints(Indice1)
sp_Invasive1 <- SpatialPolygons(list(Polygons(list(Polygon(Invasive_1)),1)))
sp_Invasive2 <- SpatialPolygons(list(Polygons(list(Polygon(Invasive_2)),1)))


Tumor_dist1 <- gDistance(sp_Invasive1, sp_tumor1, byid = TRUE)


# Generate sections
Hallmark <- TRUE
Area_profile <- data.frame(matrix(nrow = 0, ncol = 2))
Poly_profile <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_profile) <- c('x','y','section')
n <- 1
step_size <- 0.15
low_boundary <- 0
while(isTRUE(Hallmark)){
  Sec <- which(low_boundary < Tumor_dist1 & Tumor_dist1 <= low_boundary + step_size)
  Sec_Coords <- Indice1[Sec,]
  if(nrow(Sec_Coords) != 0){
    Sec_poly <- as.data.frame(concaveman(as.matrix(Sec_Coords[,1:2]), concavity = 2))
    Sec_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Sec_poly)),1))))
    Area_sub <- c(n, Sec_area)
    names(Area_sub) <- c('Section', 'Area')
    Area_profile <- rbind(Area_profile, Area_sub)
    
    # Poly sub
    Poly_sub <- cbind(Sec_poly, n)
    Poly_profile <- rbind(Poly_profile, Poly_sub)
    low_boundary <- low_boundary + step_size
    n <- n + 1
  } else{
    
    Hallmark <- FALSE
    
  }
}

saveRDS(Poly_profile, './Case 1/Tumor_poly1.rds')

# Case 1A
sp_tumor2 <- SpatialPoints(Indice2)
sp_Invasive2 <- SpatialPolygons(list(Polygons(list(Polygon(Invasive_2)),1)))


Tumor_dist2 <- gDistance(sp_Invasive2, sp_tumor2, byid = TRUE)


# Generate 
Hallmark <- TRUE
Area_profile <- data.frame(matrix(nrow = 0, ncol = 2))
Poly_profile2 <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_profile2) <- c('x','y','section')
n <- 1
step_size <- 0.15
low_boundary <- 0
while(isTRUE(Hallmark)){
  Sec <- which(low_boundary < Tumor_dist2 & Tumor_dist2 <= low_boundary + step_size)
  Sec_Coords <- Indice2[Sec,]
  if(nrow(Sec_Coords) != 0){
    Sec_poly <- as.data.frame(concaveman(as.matrix(Sec_Coords[,1:2]), concavity = 2))
    Sec_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Sec_poly)),1))))
    Area_sub <- c(n, Sec_area)
    names(Area_sub) <- c('Section', 'Area')
    Area_profile <- rbind(Area_profile, Area_sub)
    
    # Poly sub
    Poly_sub <- cbind(Sec_poly, n)
    Poly_profile2 <- rbind(Poly_profile2, Poly_sub)
    low_boundary <- low_boundary + step_size
    n <- n + 1
  } else{
    
    Hallmark <- FALSE
    
  }
}

saveRDS(Poly_profile2, './Case 1/Tumor_poly2.rds')


#Combine Case 1A and 1B
Poly_profile$n <- Poly_profile$n + 0.1
Combine_poly <- rbind(Poly_profile, Poly_profile2)
Combine_poly <- Combine_poly[order(Combine_poly$n),]

ggplot() + #geom_path(data = Poly_profile, aes(x = Poly_profile[,1], y = Poly_profile[,2], group = n, fill = factor(n)), alpha = 0.6) +
  
  geom_path(data = Combine_poly, aes(x = Combine_poly[,1], y = Combine_poly[,2], group = n, fill = factor(n)), alpha = 0.6) +
  geom_vline(xintercept= 13.160, linetype="dashed", color = "red") +
  geom_hline(yintercept= 12.832, linetype="dashed", color = "blue") +  
  xlim(0, 30) +
  ylim(25, 0)

saveRDS(Combine_poly, './Case 1/tumor_stratified_regions.rds')

##############################
# Extract invasive pixels ####
##############################
Tumor_poly1 <- readRDS('./Case 1/Tumor_poly1.rds')

ggplot() + geom_path(data = Combine_poly, aes(x = Combine_poly[,1], y = Combine_poly[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = tumor, aes(x = tumor[,1], y = tumor[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = Invasive_1, aes(x = Invasive_1[,1], y = Invasive_1[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = tumor_poly1, aes(x = tumor_poly1[,1], y = tumor_poly1[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  
  ylim(20, 0) +
  xlim(0, 30)


# Case 1A
tumor_poly1 <- readRDS('./Case 1/Tumor_poly1.rds')

tumor_poly1 <- as.data.frame(concaveman(as.matrix(tumor_poly1), concavity = 2))

sp_tumor_poly1 <- SpatialPolygons(list(Polygons(list(Polygon(tumor_poly1[,1:2])),1)))


Pixel_coords <- read.csv('./Case 1/Pixel_coords.csv')/125
Invasive_pixels1 <- point.in.polygon(Pixel_coords[,1],Pixel_coords[,2], Invasive_1[,1], Invasive_1[,2])

Pixel_coords <- cbind(Pixel_coords, data.frame(Invasive_pixels1))
Pixel_coords <- Pixel_coords[Pixel_coords[,3] != 0,-3]
sp_invasive_pixels1 <- SpatialPoints(Pixel_coords)

Invasive_dist <- gDistance(sp_tumor_poly1, sp_invasive_pixels1,byid = TRUE)
saveRDS(Invasive_dist, './Case 1/invasive_dist.rds')

# Generate 
Hallmark <- TRUE
Area_profile <- data.frame(matrix(nrow = 0, ncol = 2))
Poly_profile <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_profile) <- c('x','y','section')
n <- 1
step_size <- 0.15
low_boundary <- 0
while(isTRUE(Hallmark)){
  Sec <- which(low_boundary <Invasive_dist & Invasive_dist <= low_boundary + step_size)
  Sec_Coords <- Pixel_coords[Sec,]
  if(nrow(Sec_Coords) != 0){
    Sec_poly <- as.data.frame(concaveman(as.matrix(Sec_Coords[,1:2]), concavity = 2))
    Sec_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Sec_poly)),1))))
    Area_sub <- c(n, Sec_area)
    names(Area_sub) <- c('Section', 'Area')
    Area_profile <- rbind(Area_profile, Area_sub)
    
    # Poly sub
    Poly_sub <- cbind(Sec_poly, n)
    Poly_profile <- rbind(Poly_profile, Poly_sub)
    low_boundary <- low_boundary + step_size
    n <- n + 1
  } else{
    
    Hallmark <- FALSE
    
  }
}

# Case 1B
tumor_poly2 <- readRDS('./Case 1/Tumor_poly2.rds')

tumor_poly2 <- as.data.frame(concaveman(as.matrix(tumor_poly2), concavity = 2))

sp_tumor_poly2 <- SpatialPolygons(list(Polygons(list(Polygon(tumor_poly2[,1:2])),1)))


Pixel_coords <- read.csv('./Case 1/Pixel_coords.csv')/125
Invasive_pixels2 <- point.in.polygon(Pixel_coords[,1],Pixel_coords[,2], Invasive_2[,1], Invasive_2[,2])

Pixel_coords <- cbind(Pixel_coords, data.frame(Invasive_pixels2))
Pixel_coords <- Pixel_coords[Pixel_coords[,3] != 0,-3]
sp_invasive_pixels2 <- SpatialPoints(Pixel_coords)

Invasive_dist <- gDistance(sp_tumor_poly2, sp_invasive_pixels2, byid = TRUE)
#saveRDS(Invasive_dist, './Case 5/invasive_dist.rds')

# Generate 
Hallmark <- TRUE
Area_profile <- data.frame(matrix(nrow = 0, ncol = 2))
Poly_profile2 <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_profile2) <- c('x','y','section')
n <- 1
step_size <- 0.15
low_boundary <- 0
while(isTRUE(Hallmark)){
  Sec <- which(low_boundary <Invasive_dist & Invasive_dist <= low_boundary + step_size)
  Sec_Coords <- Pixel_coords[Sec,]
  if(nrow(Sec_Coords) != 0){
    Sec_poly <- as.data.frame(concaveman(as.matrix(Sec_Coords[,1:2]), concavity = 2))
    Sec_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Sec_poly)),1))))
    Area_sub <- c(n, Sec_area)
    names(Area_sub) <- c('Section', 'Area')
    Area_profile <- rbind(Area_profile, Area_sub)
    
    # Poly sub
    Poly_sub <- cbind(Sec_poly, n)
    Poly_profile2 <- rbind(Poly_profile2, Poly_sub)
    low_boundary <- low_boundary + step_size
    n <- n + 1
  } else{
    
    Hallmark <- FALSE
    
  }
}

# adjust the shape
Combine_poly <- Poly_profile[Poly_profile[,3] <= 7,]
Combine_poly2 <- Poly_profile2[Poly_profile2[,3] <= 7,]
Combine_poly2$n <- Combine_poly2$n + 0.1
CombinePoly <- rbind(Combine_poly, Combine_poly2)
ggplot() + geom_polygon(data = CombinePoly, aes(x = CombinePoly[,1], y = CombinePoly[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = Invasive_1, aes(x = Invasive_1[,1], y = Invasive_1[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = tumor_poly1, aes(x = tumor_poly1[,1], y = tumor_poly1[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  ylim(20, 0) +
  xlim(0, 30)

saveRDS(CombinePoly, './Case 1/invasive_stratified_regions.rds')

##############################
### Extract normal pixels ####
##############################

# Case 1A
Pixel_coords <- read.csv('./Case 1/Pixel_coords.csv')/125

Normal_region <- readRDS('./Case 1/Normal_case1.rds')

Normal_poly1 <- Normal_region[Normal_region[,3] == 1,]

sp_normal <- SpatialPolygons(list(Polygons(list(Polygon(Normal_poly1[,1:2])),1)))

Normal_pixels <- point.in.polygon(Pixel_coords[,1],Pixel_coords[,2], Normal_poly1[,1], Normal_poly1[,2])
Normal_pixels <- cbind(Pixel_coords, data.frame(Normal_pixels))
Normal_pixels <- Normal_pixels[Normal_pixels[,3] != 0,-3]
sp_normal_pixels <- SpatialPoints(Normal_pixels)

Normal_dist <- gDistance(sp_Invasive2, sp_normal_pixels,byid = TRUE)
# Generate 
Hallmark <- TRUE
Area_profile <- data.frame(matrix(nrow = 0, ncol = 2))
Poly_profile <- data.frame(matrix(nrow = 0, ncol = 3))
Poly_coords <- data.frame(matrix(nrow = 0, ncol = 2))

names(Poly_profile) <- c('x','y','section')
names(Poly_coords) <- c('x','y')

n <- 1
step_size <- 0.15
low_boundary <- 0
while(isTRUE(Hallmark)){
  Sec <- which(low_boundary <Normal_dist & Normal_dist <= low_boundary + step_size)
  Sec_Coords <- Normal_pixels[Sec,-3]
  names(Sec_Coords) <- c('x','y')
  
  if(nrow(Sec_Coords) != 0){
    Poly_coords_sub <- cbind(Sec_Coords, n)
    Poly_coords <- rbind(Poly_coords, Poly_coords_sub)
    
    Sec_poly <- as.data.frame(concaveman(as.matrix(Sec_Coords[,1:2]), concavity = 2))
    Sec_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Sec_poly)),1))))
    Area_sub <- c(n, Sec_area)
    names(Area_sub) <- c('Section', 'Area')
    Area_profile <- rbind(Area_profile, Area_sub)
    
    # Poly sub
    Poly_sub <- cbind(Sec_poly, n)
    
    Poly_profile <- rbind(Poly_profile, Poly_sub)
    
    low_boundary <- low_boundary + step_size
    n <- n + 1
  } else{
    
    Hallmark <- FALSE
    
  }
}
Poly_collect <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_collect) <- c('x', 'y', 'n')
for(poly in seq(4,7,1)){
  Polygon <- Poly_profile[Poly_profile[,3] == poly,]
  Polygon1 <-  Polygon[(Polygon[,1]  < 7),]
  Polygon2 <-  Polygon[(Polygon[,1]  > 9),]
  if(nrow(Polygon1) != 0){
  Polygon1$n <- poly + 0.1
  names(Polygon1) <- c('x', 'y', 'n')
  }
  if(nrow(Polygon2) != 0){
  Polygon2$n <- poly + 0.2
  names(Polygon2) <- c('x', 'y', 'n')
  }
  Poly_collect <- rbind(Poly_collect, Polygon1, Polygon2)
  
}


##### adjust the shape
Poly8 <- Poly_profile[Poly_profile[,3] == 8,]
Poly8.1 <-  Poly8[(Poly8[,2]  > 10),]
Poly8.1$n <- 8.1
Poly8.2 <-  Poly8[(Poly8[,2]  > 5 & Poly8[,2] < 7.5 & Poly8[,1] < 5),]
Poly8.2$n <- 8.2

Poly8.3 <-  Poly8[(Poly8[,1] > 9),]
Poly8.3$n <- 8.3
Poly8 <-rbind(Poly8.1, Poly8.2, Poly8.3)

# for poly 9-29
Poly_collect1 <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_collect1) <- c('x', 'y', 'n')
for(poly in seq(9,29,1)){
  Polygon <- Poly_profile[Poly_profile[,3] == poly,]
  Polygon1 <-  Polygon[(Polygon[,1]  < 2.4),]
  Polygon2 <-  Polygon[(Polygon[,1]  > 10),]
  if(nrow(Polygon1) != 0){
    Polygon1$n <- poly + 0.1
    names(Polygon1) <- c('x', 'y', 'n')
  }
  if(nrow(Polygon2) != 0){
    Polygon2$n <- poly + 0.2
    names(Polygon2) <- c('x', 'y', 'n')
  }
  Poly_collect1 <- rbind(Poly_collect1, Polygon1, Polygon2)
  
}
names(Poly_profile) <- c('x', 'y', 'n')
names(Poly8) <- c('x', 'y', 'n')
names(Poly_collect1) <- c('x', 'y', 'n')

Normal_poly1 <- rbind(Poly_profile[Poly_profile[,3] <= 3,], Poly_collect, Poly8, Poly_collect1)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = Normal_poly1, aes(x = Normal_poly1[,1], y = Normal_poly1[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  
  ylim(20, 0) +
  xlim(0, 30)


#Case 1B 
Pixel_coords <- read.csv('./Case 1/Pixel_coords.csv')/125

Normal_region <- readRDS('./Case 1/Normal_case1.rds')

Normal_poly2 <- Normal_region[Normal_region[,3] == 2,]


sp_normal <- SpatialPolygons(list(Polygons(list(Polygon(Normal_poly2[,1:2])),1)))

Normal_pixels <- point.in.polygon(Pixel_coords[,1],Pixel_coords[,2], Normal_poly2[,1], Normal_poly2[,2])
Normal_pixels <- cbind(Pixel_coords, data.frame(Normal_pixels))
Normal_pixels <- Normal_pixels[Normal_pixels[,3] != 0,-3]
sp_normal_pixels <- SpatialPoints(Normal_pixels)

Normal_dist <- gDistance(sp_Invasive1, sp_normal_pixels,byid = TRUE)
# Generate 
Hallmark <- TRUE
Area_profile <- data.frame(matrix(nrow = 0, ncol = 2))
Poly_profile <- data.frame(matrix(nrow = 0, ncol = 3))
Poly_coords <- data.frame(matrix(nrow = 0, ncol = 2))

names(Poly_profile) <- c('x','y','section')
names(Poly_coords) <- c('x','y')

n <- 1
step_size <- 0.15
low_boundary <- 0
while(isTRUE(Hallmark)){
  Sec <- which(low_boundary <Normal_dist & Normal_dist <= low_boundary + step_size)
  Sec_Coords <- Normal_pixels[Sec,-3]
  names(Sec_Coords) <- c('x','y')
  
  if(nrow(Sec_Coords) != 0){
    Poly_coords_sub <- cbind(Sec_Coords, n)
    Poly_coords <- rbind(Poly_coords, Poly_coords_sub)
    
    Sec_poly <- as.data.frame(concaveman(as.matrix(Sec_Coords[,1:2]), concavity = 2))
    Sec_area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Sec_poly)),1))))
    Area_sub <- c(n, Sec_area)
    names(Area_sub) <- c('Section', 'Area')
    Area_profile <- rbind(Area_profile, Area_sub)
    
    # Poly sub
    Poly_sub <- cbind(Sec_poly, n)
    
    Poly_profile <- rbind(Poly_profile, Poly_sub)
    
    low_boundary <- low_boundary + step_size
    n <- n + 1
  } else{
    
    Hallmark <- FALSE
    
  }
}


Poly_collect <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_collect) <- c('x', 'y', 'n')
for(poly in seq(1,8,1)){
  Polygon <- Poly_profile[Poly_profile[,3] == poly,]
  Polygon1 <-  Polygon[(Polygon[,2]  < 15),]
  Polygon2 <-  Polygon[(Polygon[,2]  > 15),]
  if(nrow(Polygon1) != 0){
    Polygon1$n <- poly + 0.1
    names(Polygon1) <- c('x', 'y', 'n')
  }
  if(nrow(Polygon2) != 0){
    Polygon2$n <- poly + 0.2
    names(Polygon2) <- c('x', 'y', 'n')
  }
  Poly_collect <- rbind(Poly_collect, Polygon1, Polygon2)
  
}

Poly_collect1 <- data.frame(matrix(nrow = 0, ncol = 3))
names(Poly_collect1) <- c('x', 'y', 'n')
for(poly in seq(9,31,1)){
  Polygon <- Poly_profile[Poly_profile[,3] == poly,]
  Polygon1 <-  Polygon[(Polygon[,1]  < 18),]
  Polygon2 <-  Polygon[(Polygon[,1]  > 19),]
  if(nrow(Polygon1) != 0){
    Polygon1$n <- poly + 0.1
    names(Polygon1) <- c('x', 'y', 'n')
  }
  if(nrow(Polygon2) != 0){
    Polygon2$n <- poly + 0.2
    names(Polygon2) <- c('x', 'y', 'n')
  }
  Poly_collect1 <- rbind(Poly_collect1, Polygon1, Polygon2)
  
}
Polygon <- Poly_profile[Poly_profile[,3] == 3,-3]
Polygon <- Polygon[Polygon[,2] < 15,]

Polygon1 <- Poly_profile[Poly_profile[,3] == 2,]
Polygon1 <- Polygon1[Polygon1[,2] < 15,]
Polygon1 <- Polygon1[-seq(2341, 2447,1),-3]
Polygon1 <- Polygon1[-seq(1, 1100,1),]
Polygon1 <- Polygon1[order(nrow(Polygon1):1),] #invert row order
Polygon3 <- rbind(Polygon1, Polygon, Polygon1[1,])
Polygon3$n <- 3
names(Polygon3) <- c('x', 'y', 'n')

Normal_poly2 <- rbind(Poly_collect[Poly_collect[,3] < 3,], Polygon3, Poly_collect[Poly_collect[,3] >=4,], Poly_collect1[Poly_collect1[,3]<= 25,])


saveRDS(Normal_poly1, './Case 1/Normal_poly1_final.rds')#
saveRDS(Normal_poly2, './Case 1/Normal_poly2_final.rds')#

#############################

Normal_poly1 <- readRDS('./Case 1/Normal_poly1_final.rds')
Normal_poly2 <-  readRDS('./Case 1/Normal_poly2_final.rds')

# Modify section info
Normal_poly2_mod <- Normal_poly2[Normal_poly2[,3] >=18,]
Normal_poly2_mod$n <- Normal_poly2_mod$n - 0.1

Normal_poly2_mod3 <- Normal_poly2[Normal_poly2[,3] == 3,]
Normal_poly2_mod3$n <- Normal_poly2_mod3$n + 0.1

Normal_poly2 <- rbind(Normal_poly2[Normal_poly2[,3] < 3,], Normal_poly2_mod3,Normal_poly2[Normal_poly2[,3] >=4  & Normal_poly2[,3] <18,], Normal_poly2_mod)
  
Normal_2 <- Normal_poly1[Normal_poly1[,3] <= 2,]
Normal_2$n <- Normal_2$n + 0.3

Normal_3 <- Normal_poly1[Normal_poly1[,3] == 3,]
Normal_3$n <- Normal_3$n + 0.2

Normal_48 <- Normal_poly1[Normal_poly1[,3] >= 4 & Normal_poly1[,3] < 9,]
Normal_48$n <- Normal_48$n + 0.1

Normal_913 <- Normal_poly1[Normal_poly1[,3] >= 9 & Normal_poly1[,3] < 14,]
Normal_913$n <- Normal_913$n + 0.2

Normal_1417 <- Normal_poly1[Normal_poly1[,3] >= 14 & Normal_poly1[,3] < 18,]
Normal_1417$n <- Normal_1417$n + 0.1

Normal_1824 <- Normal_poly1[Normal_poly1[,3] >= 18 & Normal_poly1[,3] < 25,]

Normal_2529 <- Normal_poly1[Normal_poly1[,3] >= 25 & Normal_poly1[,3] < 30,]
Normal_2529$n <- Normal_2529$n - 0.1

# Combine to poly 1
Normal_poly1 <- rbind(Normal_2, Normal_3, Normal_38, Normal_913, Normal_1417, Normal_1824, Normal_2529)

Normal_poly <- rbind(Normal_poly1, Normal_poly2)
Normal_poly <- Normal_poly[order(Normal_poly$n),] #invert row order

saveRDS(Normal_poly, './Case 1/normal_stratified.rds')#

############################################
##### Calculate density ####################
############################################


Invasive_stratified <- readRDS('./Case 1/invasive_stratified.rds')

Tumor_stratified <- readRDS('./Case 1/tumor_stratified.rds')

Normal_stratified <- readRDS('./Case 1/normal_stratified.rds')#
names(Normal_stratified) <- c('V1', 'V2', 'n')


Normal_stratified$n <- as.numeric(as.character(Normal_stratified$n))
# Case 1A

Invasive_stratified$n <- Invasive_stratified$n + 24
Tumor_stratified$n <- Tumor_stratified$n + 31

final_poly <- rbind(Normal_stratified, Invasive_stratified, Tumor_stratified)
final_poly$n <- as.numeric(final_poly$n)

jpeg( './Case 1/Section_Case1.jpeg',units="in", width=12, height=8, res=1200)

ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = final_poly, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_point(data=CD8_Pts, aes(CD8_Pts[,1], CD8_Pts[,2])) +
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw() +
  ylim(22, 0) +
  xlim(0, 30) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text = element_text(size = 26)) 
dev.off()
 # Calculate cell number and density
CD8_Pts <- data.frame(readMat('./Case 1/FoxP3_Points.mat'))/125
Profile_int <- data.frame(matrix(nrow = 0, ncol = 4))
names(Profile_int) <- c('Section', 'Area', 'Number', 'Density')


idx <- 1
Hallmark <- TRUE
while (Hallmark == TRUE) {
  Current_poly <- final_poly[ (final_poly[,3] < idx + 1 & final_poly[,3] >= idx), ]
  if(nrow(Current_poly) != 0){
    dist <- point.in.polygon(CD8_Pts[,1], CD8_Pts[,2], Current_poly[,1], Current_poly[,2])
    Current_poly_distProfile <- cbind(CD8_Pts, dist)
    if(round(Current_poly[1,3]) == Current_poly[1,3]){ 
      Current_poly_distProfile <- Current_poly_distProfile[Current_poly_distProfile[,3] != 0,]
      Area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(Current_poly[,1:2])),1))))
      Density <- nrow(Current_poly_distProfile)/Area
      Sub_profile <- cbind(idx, Area, nrow(Current_poly_distProfile), Density)
      names(Sub_profile) <- c('Section', 'Area', 'Number', 'Density')
      
    } else {
      Area <- 0
      Number <- 0
      sub_count <- 1
      subMark <- TRUE
      while (subMark == TRUE) {
        sub_idx <- as.character(paste(idx, '.', sub_count, sep = ''))
        print(sub_idx)
        sub_current_poly <- Current_poly[ (Current_poly[,3] == sub_idx), ]
        if(nrow(sub_current_poly) != 0){
          subDist <- point.in.polygon(CD8_Pts[,1], CD8_Pts[,2], sub_current_poly[,1], sub_current_poly[,2])
          sub_distProfile <- cbind(CD8_Pts, subDist)
          sub_Profile <- sub_distProfile[sub_distProfile[,3] != 0,]
          subArea <- gArea(SpatialPolygons(list(Polygons(list(Polygon(sub_current_poly[,1:2])),1))))
          subCount <- nrow(sub_Profile)
          Area <- Area + subArea
          Number <- Number + subCount
          
          sub_count <- sub_count + 1
        } else{
          subMark <- FALSE
        }
        Density <- Number/Area
        Sub_profile <- cbind(idx, Area, Number, Density)
        names(Sub_profile) <- c('Section', 'Area', 'Number', 'Density')
        
      }
      
    }
    Sub_profile <- data.frame(Sub_profile)
    names(Sub_profile) <- c('Section', 'Area', 'Number', 'Density')
    names(Profile_int) <- c('Section', 'Area', 'Number', 'Density')
    
    Profile_int <- rbind(Profile_int, Sub_profile)
    
    idx <- idx + 1
  } else{
    Hallmark = FALSE
  }
}


#########################
# Ordering all sections #
#########################

Normal_diction  <- Profile_int[Profile_int[,1] < 25, ]
Invasive_diction  <- Profile_int[Profile_int[,1] >=25 & Profile_int[,1] < 32, ]
tumor_diction  <- Profile_int[Profile_int[,1] >= 32, ]


Normal_diction <- Normal_diction[order(nrow(Normal_diction):1),] #invert row order
Normal_diction$Section <- -0.15*(Normal_diction$Section) - 0.425
Invasive_diction <- Invasive_diction[order(nrow(Invasive_diction):1),] #invert row order

Invasive_diction$Section <- c(-0.45, -0.325, -0.175, -0.025, 0.125, 0.275, 0.425)

tumor_diction$Section <- seq(0.5, 3.8, 0.15) + 0.075

Final_diction <- rbind(Normal_diction,Invasive_diction,tumor_diction)
#saveRDS(Final_diction, './Case 1/Final_diction_FoxP3.rds')

