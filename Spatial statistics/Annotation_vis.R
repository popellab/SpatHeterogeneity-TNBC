# This file is to visualize annotations
library(R.matlab)
library(jpeg)
library(ggplot2)
library(ggforce)

################################
########### Case 1A #############
################################

Normal_stratified <- readRDS('./Case 0/Normal_stratified_split.rds')
Tumor_stratified <- readRDS('./Case 0/Tumor_stratified_split.rds')
Invasive_stratified <- readRDS('./Case 0/Invasive_stratified_split.rds')
Invasive_stratified$n <- as.numeric(as.character(Invasive_stratified$n))

# Invasive
Invasive_down <- Invasive_stratified[Invasive_stratified[,3] < 4,]
Reference_line <- Invasive_stratified[Invasive_stratified[,3] == 4,]
Reference_line$n <- 0
Invasive_up <- Invasive_stratified[Invasive_stratified[,3] >= 5,]
Invasive_up$n <-  8 -Invasive_up$n
# Tumor
Tumor_stratified$n <- Tumor_stratified$n + 3

# Normal
Normal_stratified$n <- Normal_stratified$n + 3

# Up
down_sec1 <- rbind(Invasive_down, Tumor_stratified)
up_sec1 <- rbind(Invasive_up, Normal_stratified)
jpeg( './Case 1/Section.jpeg',units="in", width=13.3, height=10, res=1200)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = up_sec1, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec1, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = Reference_line, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  #geom_polygon(data = Tumor_stratified, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw() +
  ylim(15,0) +
  xlim(0, 20) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text = element_text(size = 26)) 
dev.off()


################################
########### Case 1B #############
################################

Invasive_stratified <- readRDS('./Case 1/invasive_stratified.rds')
Tumor_stratified <- readRDS('./Case 1/tumor_stratified.rds')
Normal_stratified <- readRDS('./Case 1/normal_stratified.rds')

# Invasive
Invasive_down <- Invasive_stratified[Invasive_stratified[,3] < 4,]
Reference_line <- Invasive_stratified[Invasive_stratified[,3] == 4,]
Reference_line$n <- 0
Invasive_up <- Invasive_stratified[Invasive_stratified[,3] >= 5,]
Invasive_up$n <-  8 -Invasive_up$n
# Tumor
Tumor_stratified$n <- Tumor_stratified$n + 3

# Normal
Normal_stratified$n <- Normal_stratified$n + 3

# Up
down_sec <- rbind(Invasive_down, Tumor_stratified)
up_sec <- rbind(Invasive_up, Normal_stratified)
jpeg( './Case1AB_section.jpeg',units="in", width=9, height=6.675, res=300)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 22))+
  geom_polygon(data = up_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = up_sec1, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec1, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  xlab('x, mm') +
  ylab('y, mm') +
  ylim(20.824, 0) +
  xlim(0, 28.072) +
  coord_equal()
  
dev.off()#xlim(10, 30) +

################################
########### Case 2 #############
################################


Invasive_stratified <- readRDS('./Case 2/Case2_invasive_stratified(8.9).rds')
Tumor_stratified <- readRDS('./Case 2/Case2_tumor_stratified(8.31).rds')
Normal_stratified <- readRDS('./Case 2/Case2_normal_stratified.rds')

# Invasive
Invasive_down <- Invasive_stratified[Invasive_stratified[,3] < 4,]
Reference_line <- Invasive_stratified[Invasive_stratified[,3] == 4,]
Reference_line$n <- 0
Invasive_up <- Invasive_stratified[Invasive_stratified[,3] >= 5,]
Invasive_up$n <-  8 -Invasive_up$n
# Tumor
Tumor_stratified$n <- Tumor_stratified$n + 3

# Normal
Normal_stratified$n <- Normal_stratified$n + 3

# Up
down_sec <- rbind(Invasive_down, Tumor_stratified)
up_sec <- rbind(Invasive_up, Normal_stratified)
jpeg( './Case 2/Section.jpeg',units="in", width=12, height=6.7875, res=1200)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = up_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = Reference_line, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  #geom_polygon(data = Tumor_stratified, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw() +
  ylim(18.1, 0) +
  xlim(0, 32) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text = element_text(size = 26)) 
dev.off()
################################
########### Case 3 #############
################################
Invasive_stratified <- readRDS('./Case 3/invasive_stratified_regions.rds')
Tumor_stratified <- readRDS('./Case 3/tumor_stratified_regions.rds')
Normal_stratified <- readRDS('./Case 3/normal_stratified_regions.rds')

# Invasive
Invasive_down <- Invasive_stratified[Invasive_stratified[,3] < 4,]
Reference_line <- Invasive_stratified[Invasive_stratified[,3] == 4,]
Reference_line$n <- 0
Invasive_up <- Invasive_stratified[Invasive_stratified[,3] >= 5,]
Invasive_up$n <-  8 -Invasive_up$n
# Tumor
Tumor_stratified$n <- Tumor_stratified$n + 3

# Normal
Normal_stratified$n <- Normal_stratified$n + 3

# Up
down_sec <- rbind(Invasive_down, Tumor_stratified)
up_sec <- rbind(Invasive_up, Normal_stratified)
jpeg( './Case 3/Section.jpeg',units="in", width=12, height=8.015, res=1200)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = up_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = Reference_line, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  #geom_polygon(data = Tumor_stratified, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw() +
  ylim(18.76,0) +
  xlim(0,28.088) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text = element_text(size = 26)) 
dev.off()

################################
########### Case 4 #############
################################
Invasive_stratified <- readRDS('./Case 4/invasive_stratified_regions.rds')
Tumor_stratified <- readRDS('./Case 4/tumor_stratified_regions.rds')
Normal_stratified <- readRDS('./Case 4/normal_stratified_regions.rds')
names(Invasive_stratified) <- c('x', 'y', 'n')
names(Tumor_stratified) <- c('x', 'y', 'n')
# Invasive
Invasive_down <- Invasive_stratified[Invasive_stratified[,3] < 4,]
Reference_line <- Invasive_stratified[Invasive_stratified[,3] == 4,]
Reference_line$n <- 0
Invasive_up <- Invasive_stratified[Invasive_stratified[,3] >= 5,]
Invasive_up$n <-  8 -Invasive_up$n
# Tumor
Tumor_stratified$n <- Tumor_stratified$n + 3

# Normal
Normal_stratified$n <- Normal_stratified$n + 3

# Up
down_sec <- rbind(Invasive_down, Tumor_stratified)
up_sec <- rbind(Invasive_up, Normal_stratified)
jpeg( './Case 4/Section.jpeg',units="in", width=12, height=8.62, res=1200)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = up_sec, aes(x = x, y = y, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec, aes(x = x, y = y, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = Reference_line, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  #geom_polygon(data = Tumor_stratified, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw() +
  ylim(22.904,0) +
  xlim(0,31.9) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text = element_text(size = 26)) 
dev.off()
################################
########### Case 5 #############
################################
Invasive_stratified <- readRDS('./Case 5/invasive_stratified_regions.rds')
Tumor_stratified <- readRDS('./Case 5/tumor_stratified_regions.rds')
Normal_stratified <- readRDS('./Case 5/normal_stratified_regions.rds')

# Invasive
Invasive_down <- Invasive_stratified[Invasive_stratified[,3] < 4,]
Reference_line <- Invasive_stratified[Invasive_stratified[,3] == 4,]
Reference_line$n <- 0
Invasive_up <- Invasive_stratified[Invasive_stratified[,3] >= 5,]
Invasive_up$n <-  8 -Invasive_up$n
# Tumor
Tumor_stratified$n <- Tumor_stratified$n + 3

# Normal
Normal_stratified$n <- Normal_stratified$n + 3

# Up
down_sec <- rbind(Invasive_down, Tumor_stratified)
up_sec <- rbind(Invasive_up, Normal_stratified)
jpeg( './Case 5/Section.jpeg',units="in", width=12, height=10.19, res=1200)
ggplot() + #geom_path(data = poly_reorganize, aes(x = poly_reorganize[,1], y = poly_reorganize[,2], group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = up_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  geom_polygon(data = down_sec, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6, show.legend = FALSE) +
  #geom_polygon(data = Reference_line, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  #geom_polygon(data = Tumor_stratified, aes(x = V1, y = V2, group = n, fill = factor(n)), alpha = 0.6) +
  
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw() +
  ylim(21.424,0) +
  xlim(0,25.232) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text = element_text(size = 26)) 
dev.off()
