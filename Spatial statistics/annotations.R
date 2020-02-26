# This file is to visualize annotations
library(R.matlab)
library(jpeg)
library(ggpubr)
library(alphahull)
library(gginnards)
library(ggplot2)
library(magick)
library(tidyverse)
library(grid)
library(sp)
library(rgeos)
library(sf)
library(concaveman)
source('MiFunction.R')
# for formatting grids
library(ggpubr)
# read background image
imgage <- readJPEG('./Case4_CD4.jpeg')
bg <- ggplot()
for_plot <- bg +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0, 31.904, -22.904, 0) +
  xlab('x,mm')+
  ylab('y,mm')
plot(for_plot)

##########################
# Case 1 Invasive front###
##########################

# Contour 1
imgage <- readJPEG('./Case1_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  28.072, -20.824, 0) +
  xlab('x,mm')+
  ylab('y,mm')

CD3_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD3_1.mat'))/125
CD4_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD4_1.mat'))/125
CD8_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD8_1.mat'))/125
CD20_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_CD20_1.mat'))/125
FoxP3_Contour_1 <- as.data.frame(readMat('./Case 1/Ctr_FoxP3_1.mat'))/125
List <- infer_invasiveFront(CD3_Contour_1, CD4_Contour_1,CD8_Contour_1,CD20_Contour_1,FoxP3_Contour_1)

# Contour 2
CD3_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD3_2.mat'))/125
CD4_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD4_2.mat'))/125
CD8_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD8_2.mat'))/125
CD20_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_CD20_2.mat'))/125
FoxP3_Contour_2 <- as.data.frame(readMat('./Case 1/Ctr_FoxP3_2.mat'))/125
List2 <- infer_invasiveFront(CD3_Contour_2, CD4_Contour_2,CD8_Contour_2,CD20_Contour_2,FoxP3_Contour_2)



contour_Case1 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$`1`, aes(x = List$`1`[,1], y = List$`1`[,2]),alpha = 0.3) +
  geom_path(data = List$`2`, aes(x = List$`2`[,1], y = List$`2`[,2])) +
  geom_path(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2])) +
  geom_polygon(data = List$`4`, aes(x = List$`4`[,1], y = List$`4`[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2], fill = 'Tumor'),alpha = 0.3) +
  # Contour2
  geom_path(data = List2$`1`, aes(x = List2$`1`[,1], y = List2$`1`[,2]),alpha = 0.3) +
  geom_path(data = List2$`2`, aes(x = List2$`2`[,1], y = List2$`2`[,2])) +
  geom_path(data = List2$`3`, aes(x = List2$`3`[,1], y = List2$`3`[,2])) +
  geom_polygon(data = List2$`4`, aes(x = List2$`4`[,1], y = List2$`4`[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List2$`3`, aes(x = List2$`3`[,1], y = List2$`3`[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(22,0) + xlim(0,30) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))
plot(contour_Case1)

##########################
# Case 2 Invasive front###
##########################

# Contour 1
imgage <- readJPEG('./Case2_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  31.872, -18.024, 0) +
  xlab('x,mm')+
  ylab('y,mm')

CD3_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD3.mat'))/125
CD4_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD4.mat'))/125
CD8_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD8.mat'))/125
CD20_Contour <- as.data.frame(readMat('./Case 2/Ctr_CD20.mat'))/125
FoxP3_Contour <- as.data.frame(readMat('./Case 2/Ctr_FoxP3.mat'))/125
List <- infer_invasiveFront(CD3_Contour, CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour)

contour_Case2 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$`1`, aes(x = List$`1`[,1], y = List$`1`[,2]),alpha = 0.3) +
  geom_path(data = List$`2`, aes(x = List$`2`[,1], y = List$`2`[,2])) +
  geom_path(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2])) +
  geom_polygon(data = List$`4`, aes(x = List$`4`[,1], y = List$`4`[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(20,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))
plot(contour_Case2)


##########################
# Case 3 Invasive front###
##########################

# Contour 1
imgage <- readJPEG('./Case3_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  28.088, -18.76, 0) +
  xlab('x,mm')+
  ylab('y,mm')

CD3_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD3.mat'))/125
CD4_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD4.mat'))/125
CD8_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD8.mat'))/125
CD20_Contour <- as.data.frame(readMat('./Case 3/Ctr_CD20.mat'))/125
FoxP3_Contour <- as.data.frame(readMat('./Case 3/Ctr_FoxP3.mat'))/125
List <- infer_invasiveFront(CD3_Contour, CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour)

contour_Case3 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$`1`, aes(x = List$`1`[,1], y = List$`1`[,2]),alpha = 0.3) +
  geom_path(data = List$`2`, aes(x = List$`2`[,1], y = List$`2`[,2])) +
  geom_path(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2])) +
  geom_polygon(data = List$`4`, aes(x = List$`4`[,1], y = List$`4`[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2], fill = 'Tumor'),alpha = 0.3) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(20,0) + xlim(0,30) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))
plot(contour_Case3)
##########################
# Case 4 Invasive front###
##########################
CD3_Contour <- as.data.frame(readMat('./Case 4/CD3_Contour.mat'))/125
CD4_Contour <- as.data.frame(readMat('./Case 4/CD4_Contour.mat'))/125
CD8_Contour <- as.data.frame(readMat('./Case 4/CD8_Contour.mat'))/125
CD20_Contour <- as.data.frame(readMat('./Case 4/CD20_Contour.mat'))/125
FoxP3_Contour <- as.data.frame(readMat('./Case 4/FoxP3_Contour.mat'))/125
# Return data order: inter_Coords,inter_Coords_large, inter_Coords_small,Band_contour
List <- infer_invasiveFront(CD3_Contour, CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour)
  
contour_Case4 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$`1`, aes(x = List$`1`[,1], y = List$`1`[,2]),alpha = 0.3) +
  geom_path(data = List$`2`, aes(x = List$`2`[,1], y = List$`2`[,2])) +
  geom_path(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2])) +
  geom_polygon(data = List$`4`, aes(x = List$`4`[,1], y = List$`4`[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2], fill = 'Tumor'),alpha = 0.2) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(24,0) + xlim(0,32) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))
plot(contour_Case4)

##########################
# Case 5 Invasive front###
##########################
imgage <- readJPEG('./Case5_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  25.232, -21.424, 0) +
  xlab('x,mm')+
  ylab('y,mm')

CD3_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD3.mat'))/125
CD4_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD4.mat'))/125
CD8_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD8.mat'))/125
CD20_Contour <- as.data.frame(readMat('./Case 5/Ctr_CD20.mat'))/125
FoxP3_Contour <- as.data.frame(readMat('./Case 5/Ctr_FoxP3.mat'))/125
# Return data order: inter_Coords,inter_Coords_large, inter_Coords_small,Band_contour

List <- infer_invasiveFront(CD3_Contour, CD4_Contour,CD8_Contour,CD20_Contour,FoxP3_Contour)

contour_Case5 <- move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = List$`1`, aes(x = List$`1`[,1], y = List$`1`[,2]),alpha = 0.3) +
  geom_path(data = List$`2`, aes(x = List$`2`[,1], y = List$`2`[,2])) +
  geom_path(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2])) +
  geom_polygon(data = List$`4`, aes(x = List$`4`[,1], y = List$`4`[,2], fill = 'Invasive front'),alpha = 0.2) +
  geom_polygon(data = List$`3`, aes(x = List$`3`[,1], y = List$`3`[,2], fill = 'Tumor'),alpha = 0.2) +
  scale_fill_manual(values  = c('#F26969','#F2F2C6')) + 
  theme_bw() +
  ylim(22,0) + xlim(0,26) +
  xlab('x,mm') +
  ylab('y,mm') +
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14)) +
  labs(fill = "Annotations") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12))
plot(contour_Case5)

