######################################################################################
# This script is used to obtain the coordinates for invasive front and central tumor #
######################################################################################

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
library(concaveman)
source('MiFunction.R')
# for formatting grids
library(ggpubr)
# for matrix smoothing
library(oce)


##########################
# Case 1 Invasive front###
##########################

# Contour 1
Color_coords <- read.csv('./Case 1/Pixel_coords.csv',sep = ',',header = T)
imgage <- readJPEG('./Case1_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  28.072, -20.824, 0) +
  xlab('x,mm')+
  ylab('y,mm')


CD3_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD3_1.mat'))
CD4_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD4_1.mat'))
CD8_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD8_1.mat'))
CD20_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD20_1.mat'))
FoxP3_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_FoxP3_1.mat'))
List <- infer_invasiveFrontColror(CD3_Contour_1,CD4_Contour_1,CD8_Contour_1,CD20_Contour_1,FoxP3_Contour_1, Color_coords,2.5, 2603)

# outer

Invasive_outer <- List$outter
Invasive_outer <- Invasive_outer[!(Invasive_outer[,1] > 14.6 & Invasive_outer[,2] > 16.225),]
Invasive_outer <- Invasive_outer[!((Invasive_outer[,2] < 16.225 & Invasive_outer[,2] > 13.2605) & Invasive_outer[,1] > 18.114),]
ggplot() + geom_point(data = Invasive_outer, aes(x = Invasive_outer[,1], y = Invasive_outer[,2]))

# Inner
Invasive_inner <- List$tumor
Invasive_inner <- Invasive_inner[!(Invasive_inner[,1] > 14.6 & Invasive_inner[,2] > 16.225),]

Invasive_inner <- Invasive_inner[!( (Invasive_inner[,2] < 16.225 & Invasive_inner[,2] > 13.21875) & Invasive_inner[,1] > 18.114),]
ggplot() + geom_point(data = Invasive_inner, aes(x = Invasive_inner[,1], y = Invasive_inner[,2]))

Invasive_inner <- Invasive_inner[-1,]

Invasive_inner <- rbind(Invasive_inner[-seq(1,60,1),],Invasive_inner[seq(1,60,1),])

Invasive_outer <- Invasive_outer[-1,]

Invasive_outer <- rbind(Invasive_outer[-seq(1,72,1),],Invasive_outer[seq(1,72,1),])
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order

Band <- rbind(Invasive_outer, Invasive_inner, Invasive_outer[1,])
ggplot() + geom_path(data = Band, aes(x = Band[,1], y = Band[,2]))

# Rest part
rest <- List$outter
rest <- rest[(rest[,1] > 14.6),] #408/729*20, 768/910*25

rest <- rest[(rest[,2] > 16.225),] #408/729*20, 768/910*25
rest <- rest[!(rest[,2] < 10 & rest[,1] < 8),] #408/729*20, 768/910*25

names(rest) <- c('x','y')
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Final_tumor <- rbind(Invasive_inner, rest,Invasive_inner[1,] )
ggplot() + geom_path(data = Final_tumor, aes(x = Final_tumor[,1], y = Final_tumor[,2]))



# Contour 2
CD3_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD3_2.mat'))
CD4_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD4_2.mat'))
CD8_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD8_2.mat'))
CD20_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD20_2.mat'))
FoxP3_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_FoxP3_2.mat'))
List2 <- infer_invasiveFrontColror(CD3_Contour_2,CD4_Contour_2,CD8_Contour_2,CD20_Contour_2,FoxP3_Contour_2, Color_coords,2.5, 2603, 0.5)

# outer

Invasive_outer <- List2$outter
Invasive_outer <- Invasive_outer[!(Invasive_outer[,1] > 3.94 & Invasive_outer[,2] > 8.0996),]
Invasive_outer <- Invasive_outer[!((Invasive_outer[,2] < 8.0996 & Invasive_outer[,2] > 5.919) & Invasive_outer[,1] > 5.166),]

# Inner
Invasive_inner <- List2$tumor
Invasive_inner <- Invasive_inner[!(Invasive_inner[,1] > 4.244 & Invasive_inner[,2] > 7.788),]

Invasive_inner <- Invasive_inner[!( (Invasive_inner[,2] < 7.788 & Invasive_inner[,2] > 5.66978) & Invasive_inner[,1] > 5.166),]

Invasive_inner <- Invasive_inner[-1,]

Invasive_inner <- rbind(Invasive_inner[-seq(1,176,1),],Invasive_inner[seq(1,176,1),])

Invasive_outer <- Invasive_outer[-1,]

Invasive_outer <- rbind(Invasive_outer[-seq(1,179,1),],Invasive_outer[seq(1,179,1),])
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order

Band_2 <- rbind(Invasive_outer, Invasive_inner, Invasive_outer[1,])

# Rest part
rest <- List2$outter
rest <- rest[(rest[,1] > 3.94),] #408/729*20, 768/910*25

rest <- rest[(rest[,2] > 5.919),] #408/729*20, 768/910*25
rest <- rest[!(rest[,2] < 10 & rest[,1] < 8),] #408/729*20, 768/910*25

names(rest) <- c('x','y')

Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Final_tumor_2 <- rbind(Invasive_inner, rest,Invasive_inner[1,] )






contour_Case1 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2])) +
  geom_path(data = List$outter, aes(x = List$outter[,1], y = List$outter[,2])) +
  geom_polygon(data = List$invasiveCoords, aes(x = List$invasiveCoords[,1], y = List$invasiveCoords[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$tumor, aes(x = List$tumor[,1], y =List$tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  #geom_path(data = Band, aes(x = Band[,1], y = Band[,2])) +
  #geom_polygon(data = Band, aes(x = Band[,1], y = Band[,2], fill = 'Invasive front'),alpha = 0.2) +
  #geom_polygon(data = Final_tumor, aes(x = Final_tumor[,1], y = Final_tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  # Contour2
  geom_path(data = List2$tumor, aes(x = List2$tumor[,1], y = List2$tumor[,2])) +
  geom_path(data = List2$outter, aes(x = List2$outter[,1], y = List2$outter[,2])) +
  geom_polygon(data = List2$invasiveCoords, aes(x = List2$invasiveCoords[,1], y = List2$invasiveCoords[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List2$tumor, aes(x = List2$tumor[,1], y =List2$tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  #geom_path(data = Band_2, aes(x = Band_2[,1], y = Band_2[,2])) +
  #geom_polygon(data = Band_2, aes(x = Band_2[,1], y = Band_2[,2], fill = 'Invasive front'),alpha = 0.2) +
  #geom_polygon(data = Final_tumor_2, aes(x = Final_tumor_2[,1], y = Final_tumor_2[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(22,0) + xlim(0,30) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  labs(tag = 'Case1-CD4')
plot(contour_Case1)




##########################
# Case 2 Invasive front###
##########################
# Color method 
Color_coords <- read.csv('./Case 2/Pixel_coords.csv',sep = ',',header = T)

CD3_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD3.mat'))
CD4_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD4.mat'))
CD8_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD8.mat'))
CD20_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD20.mat'))
FoxP3_Contour <- as.data.frame(readMat('./Case 2/Ctr_FoxP3.mat'))



####
List <- infer_invasiveFrontColror(CD3_Contour,CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour, Color_coords,2.5, 2253,0.5)
# outer

Invasive_outer <- List$outter
Invasive_outer <- Invasive_outer[!(Invasive_outer[,1] < 28.3 & Invasive_outer[,2] > 9.7),]
Invasive_outer <- Invasive_outer[!( (Invasive_outer[,1] > 6.614 & Invasive_outer[,2] > 8.814) & (Invasive_outer[,1] < 28.3 & Invasive_outer[,2] > 8.814)),]
Invasive_outer <- rbind(Invasive_outer,Invasive_outer[2,])
Invasive_outer <- Invasive_outer[-c(1,2),]

# Inner
Invasive_inner <- List$tumor
Invasive_inner <- Invasive_inner[!(Invasive_inner[,1] < 28.3 & Invasive_inner[,2] > 9.7),]
 
Invasive_inner <- Invasive_inner[!( (Invasive_inner[,1] > 6.806 & Invasive_inner[,2] > 8.814) & (Invasive_inner[,1] < 28.3 & Invasive_inner[,2] > 8.814)),]
Invasive_inner <- rbind(Invasive_inner,Invasive_inner[2:11,])
Invasive_inner <- Invasive_inner[-seq(1,11,1),]
# Combine
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Band <- rbind(Invasive_outer, Invasive_inner)
Band <- rbind(Band, Band[1,])
# Fine lower boundary
Lower <- List$refence
Lower <- Lower[(Lower[,1] < 28.3 & Lower[,1] > 6.614),]

Lower <- Lower[(Lower[,2] > 8.814),]
names(Lower) <- c('x','y')

Final_tumor <- rbind(Lower, Invasive_inner)
# save data to file
#saveRDS(Final_tumor, file = "Case2_tumor.rds")
#saveRDS(Band, file = "Case2_Invasive.rds")
imgage <- readJPEG('./Case2_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  31.872, -18.024, 0) +
  xlab('x,mm')+
  ylab('y,mm')

contour_Case2 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  #geom_path(data = List$outter, aes(x = List$tumor[,1], y = List$tumor[,2])) +
  #geom_path(data = List$outter, aes(x = List$outter[,1], y = List$outter[,2])) +
  #geom_polygon(data = List$invasiveCoords, aes(x = List$invasiveCoords[,1], y = List$invasiveCoords[,2], fill = 'Invasive front'),alpha = 0.2) +
  #geom_polygon(data = List$tumor, aes(x = List$tumor[,1], y =List$tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  geom_path(data = Band, aes(x = Band[,1], y = Band[,2])) +
  geom_polygon(data = Band, aes(x = Band[,1], y = Band[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = Final_tumor, aes(x = Final_tumor[,1], y =Final_tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(20,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  labs(tag = 'Case 2-CD4')

plot(contour_Case2) 

###########################
####Case 3#################
###########################

Color_coords <- read.csv('./Case 3/Pixel_coords.csv',sep = ',',header = T)

CD3_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD3.mat'))
CD4_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD4.mat'))
CD8_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD8.mat'))
CD20_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD20.mat'))
FoxP3_Contour <- as.data.frame(readMat('./Case 3/Ctr_FoxP3.mat'))

List <- infer_invasiveFrontColror(CD3_Contour,CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour, Color_coords,2.5, 2345)
####
# Inner

Invasive_outer <- List$outter
Invasive_outer <- Invasive_outer[(Invasive_outer[,1] < 6.486),]
Invasive_outer <- Invasive_outer[ (Invasive_outer[,2] < 15 & Invasive_outer[,2] > 5.92),]

# Inner
Invasive_inner <- List$tumor
Invasive_inner <- Invasive_inner[(Invasive_inner[,1] < 6.72),]

Invasive_inner <- Invasive_inner[( Invasive_inner[,2] > 6.28 & Invasive_inner[,2] < 15.3),]

# Combine
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Band <- rbind(Invasive_outer, Invasive_inner)
Band <- rbind(Band, Band[1,])
# Fine lower boundary
rest <- List$refence
rest <- rest[!(rest[,1] < 6.667 & (rest[,2] < 15.1 & rest[,2] > 5.96)),]

names(rest) <- c('x','y')
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order

Final_tumor <- rbind(rest[seq(1,92,1),], Invasive_inner, rest[-seq(1,92,1),])

#####
imgage <- readJPEG('./Case3_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  28.088, -18.76, 0) +
  xlab('x,mm')+
  ylab('y,mm')

contour_Case3 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = Band, aes(x = Band[,1], y = Band[,2])) +
  geom_polygon(data = Band, aes(x = Band[,1], y = Band[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = Final_tumor, aes(x = Final_tumor[,1], y = Final_tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  #geom_path(data = List$outter, aes(x = List$outter[,1], y = List$outter[,2])) +
  #geom_path(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2])) +
  #geom_polygon(data = List$invasiveCoords, aes(x = List$invasiveCoords[,1], y = List$invasiveCoords[,2], fill = 'Invasive front'),alpha = 0.2) +
  #geom_polygon(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(20,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  labs(tag = 'Case3-CD4')
plot(contour_Case3) 

###########################
####Case 4#################
###########################
Color_coords <- read.csv('./Case 4/Pixel_coords.csv',sep = ',',header = T)

CD3_Contour <- as.data.frame(readMat('./Case 4/CD3_Contour.mat'))
CD4_Contour <- as.data.frame(readMat('./Case 4/CD4_Contour.mat'))
CD8_Contour <- as.data.frame(readMat('./Case 4/CD8_Contour.mat'))
CD20_Contour <- as.data.frame(readMat('./Case 4/CD20_Contour.mat'))
FoxP3_Contour <- as.data.frame(readMat('./Case 4/FoxP3_Contour.mat'))

List <- infer_invasiveFrontColror(CD3_Contour,CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour, Color_coords,2.5, 2863)
####
# outter

Invasive_outer <- List$outter
Invasive_outer <- Invasive_outer[!(Invasive_outer[,2] < 20.9 & Invasive_outer[,2] > 11.19 & Invasive_outer[,1] < 5.55),] #408/729*20, 768/910*25

# inner
Invasive_inner <- List$tumor

Invasive_inner <- Invasive_inner[!(Invasive_inner[,2] < 21.098 & Invasive_inner[,2] > 11.19 & Invasive_inner[,1] < 5.55),] #408/729*20, 768/910*25
Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Band <- rbind(Invasive_outer, Invasive_inner)
Band <- rbind(Band, Band[1,])
# Fine lower boundary
rest <- List$outter
rest <- rest[(rest[,1] < 5.55 & (rest[,2] < 21.098 & rest[,2] > 11.19)),]

names(rest) <- c('x','y')
rest_sec1 <- rest[seq(387,425,1),]
rest_sec2 <- rest[-seq(387,425,1),]

Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Final_tumor <- rbind(Invasive_inner, rest_sec1,rest_sec2, Invasive_inner[1,])



imgage <- readJPEG('./Case4_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  31.904, -22.904, 0) +
  xlab('x,mm')+
  ylab('y,mm')


contour_Case4 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$outter, aes(x = List$outter[,1], y = List$outter[,2])) +
  geom_path(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2])) +
  geom_polygon(data = List$invasiveCoords, aes(x = List$invasiveCoords[,1], y = List$invasiveCoords[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  #geom_path(data = Band, aes(x = Band[,1], y = Band[,2])) +
  #geom_polygon(data = Band, aes(x = Band[,1], y = Band[,2], fill = 'Invasive front'),alpha = 0.2) +
  #geom_polygon(data = Final_tumor, aes(x = Final_tumor[,1], y = Final_tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(24,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  labs(tag = 'Case4-CD4')
plot(contour_Case4) 

###########################
####Case 5#################
###########################
Color_coords <- read.csv('./Case 5/Pixel_coords.csv',sep = ',',header = T)

CD3_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD3.mat'))
CD4_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD4.mat'))
CD8_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD8.mat'))
CD20_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD20.mat'))
FoxP3_Contour <- as.data.frame(readMat('./Case 5/Ctr_FoxP3.mat'))

List <- infer_invasiveFrontColror(CD3_Contour,CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour, Color_coords,2.5, 2678)
###########################3
Invasive_outer <- List$outter
Invasive_outer <- Invasive_outer[(Invasive_outer[,2] < 13.8986),] #408/729*20, 768/910*25
Invasive_outer <- Invasive_outer[!(Invasive_outer[,2] > 10.14 & Invasive_outer[,1] < 10),] #408/729*20, 768/910*25
ggplot() + geom_path(data = Invasive_outer, aes(x = Invasive_outer[,1], y = Invasive_outer[,2]))

# inner
Invasive_inner <- List$tumor

Invasive_inner <- Invasive_inner[(Invasive_inner[,2] < 13.6363),] #408/729*20, 768/910*25
Invasive_inner <- Invasive_inner[!(Invasive_inner[,2] > 10.59 & Invasive_inner[,1] < 10),] #408/729*20, 768/910*25
Invasive_inner <- Invasive_inner[-1,]
Invasive_inner <- rbind(Invasive_inner[-seq(1,35,1),],Invasive_inner[seq(1,35,1),])

Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order

Band <- rbind(Invasive_outer, Invasive_inner)
Band <- rbind(Band, Band[1,])

# Fine lower boundary
rest <- List$outter
rest <- rest[(rest[,2] > 10.14),] #408/729*20, 768/910*25

rest <- rest[!(rest[,2] < 13.8968 & rest[,1] > 10),] #408/729*20, 768/910*25

names(rest) <- c('x','y')
rest_sec1 <- rest[seq(1,832,1),]
rest_sec2 <- rest[-seq(1,832,1),]

Invasive_inner = Invasive_inner[order(nrow(Invasive_inner):1),] #invert row order
Final_tumor <- rbind(Invasive_inner, rest_sec2,rest_sec1,Invasive_inner[1,] )

imgage <- readJPEG('./Case5_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  25.232, -21.424, 0) +
  xlab('x,mm')+
  ylab('y,mm')

contour_Case5 <- move_layers(for_plot, idx = 1L, position = "bottom") +
  geom_path(data = List$outter, aes(x = List$outter[,1], y = List$outter[,2])) +
  geom_path(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2])) +
  geom_polygon(data = List$invasiveCoords, aes(x = List$invasiveCoords[,1], y = List$invasiveCoords[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$tumor, aes(x = List$tumor[,1], y = List$tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  #geom_path(data = Band, aes(x = Band[,1], y = Band[,2])) +
  #geom_polygon(data = Band, aes(x = Band[,1], y = Band[,2], fill = 'Invasive front'),alpha = 0.2) +
  #geom_polygon(data = Final_tumor, aes(x = Final_tumor[,1], y = Final_tumor[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(22,0) + xlim(0,26) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  labs(tag = 'Case5-CD4')
plot(contour_Case5) 















######
Invasive_inner = List$tumor[order(nrow(List$tumor):1),] #invert row order
# Create real-invasive front
#Invasive_pts <- filter(List$outter, List$invasiveCoords[,1] < 6.8 | List$invasiveCoords[,2] < 8.7 | List$invasiveCoords[,1] >= 28.3)
Invasive_out <- List$outter[!(List$outter[,1] < 6.9 & List$outter[,2] > 9.8),]
Invasive_out <- Invasive_out[!(Invasive_out[,1] < 28.3 & Invasive_out[,2] > 8.814 & Invasive_out[,1] >6.87),]

Invasive_inner <- Invasive_inner[!(Invasive_inner[,1] < 6.9 & Invasive_inner[,2] > 9.8),]
Invasive_inner <- Invasive_inner[!(Invasive_inner[,1] < 28.3 & Invasive_inner[,2] > 8.814 & Invasive_inner[,1] >6.87),]
Invasive_out = Invasive_out[order(nrow(Invasive_out):1),] #invert row order

bind <- rbind(Invasive_out,Invasive_inner )

ggplot() + #geom_point(data = Invasive_inner, aes(x = Invasive_inner[,1], y = Invasive_inner[,2])) +
  geom_point(data = Invasive_out, aes(x = Invasive_out[,1], y = Invasive_out[,2]))
#geom_path(data = bind, aes(x = bind[,1], y = bind[,2]))
imgage <- readJPEG('./Case2_CD4.jpeg')
