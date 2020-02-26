####################################################
### This script is used to validate registration ###
####################################################

library(R.matlab)
library(ggplot2)
library(concaveman)
library(spatstat)
library(sf)
library(rgeos)
library(raster)
library(lwgeom)

# iterate over cases
Dat <- data.frame(matrix(nrow = 0, ncol =4))
#colnames(Dat) <- c('DSC', 'Mark', 'Case')

for(case in seq(1,6,1)){
  case_name <- switch(case, './Case_0/', './Case_1/', './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  Case <- switch(case, 'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  # iterate each marker
  for(marker in seq(1, 4)){
    mark <- switch(marker, 'CD3', 'CD8', 'CD20', 'FoxP3')
    
    # read reference path
    ref_path <- data.frame(readMat(paste(case_name, 'Contour/', Case, '/reference.mat', sep = '')))/125
    ref_path <- as.data.frame(concaveman(as.matrix(ref_path[,1:2]), concavity = 2))
    sp_ref <- SpatialPolygons(list(Polygons(list(Polygon(ref_path)),1)))
    ref_Area <- gArea(sp_ref) # calculate area
  
    # read global path
    global_path <- data.frame(readMat(paste(case_name, 'Contour/',Case, '/', mark,'/', 'global.mat', sep = '')))/125
    global_path <- as.data.frame(concaveman(as.matrix(global_path[,1:2]), concavity = 2))
    
    sp_global <- SpatialPolygons(list(Polygons(list(Polygon(global_path)),1)))
    sp_global <- gBuffer(sp_global, byid=TRUE, width=0)
    
    global_Area <- gArea(sp_global) # calculate area
    
    # read local path
    local_path <- data.frame(readMat(paste(case_name, 'Contour/', Case, '/', mark,'/', 'local.mat', sep = '')))/125
    local_path <- as.data.frame(concaveman(as.matrix(local_path[,1:2]), concavity = 2))
    
    sp_local <- SpatialPolygons(list(Polygons(list(Polygon(local_path)),1)))
    
    local_Area <- gArea(sp_local) # calculate area
    
    # Calculate areas for the intersections
    intersection_lr <- gArea(intersect(sp_ref, sp_local))
    intersection_gr <- gArea(intersect(sp_ref, sp_global))
    
    # Dice Similarity Coefficient:
    DSC_local <- 2*intersection_lr/(local_Area + ref_Area)
    DSC_global <- 2*intersection_gr/(global_Area + ref_Area)
    
    # store data
    subDat <- cbind(DSC_local,DSC_global, mark, Case)
    colnames(subDat) <- c('DSC_local','DSC_global', 'Mark', 'Case')
    Dat <- rbind(Dat, subDat)
    
  }
}
# formatting
Dat$DSC_local <- as.numeric(as.character(Dat$DSC_local))
Dat$DSC_global <- as.numeric(as.character(Dat$DSC_global))

Dat$DSC_local <- Dat$DSC_local/(2 - Dat$DSC_local)

Dat$DSC_global <- Dat$DSC_global/(2 - Dat$DSC_global)


# assign type info
Dat_ROI_local <- data.frame(Dat$DSC_local)
Dat_ROI_local$criteria <- 'local'

Dat_ROI_global <- data.frame(Dat$DSC_global)
Dat_ROI_global$criteria <- 'global'
names(Dat_ROI_local) <- c('value', 'criteria')
names(Dat_ROI_global) <- c('value', 'criteria')

Dat_ROI <- rbind(Dat_ROI_local, Dat_ROI_global)

# statistical analysis

wilcox.test(Dat_ROI_local$value,Dat_ROI_global$value, paired = FALSE, alternative = 'two.sided')

#t.test(Dat_ROI_local$value,Dat_ROI_global$value,  alternative = 'greater')

