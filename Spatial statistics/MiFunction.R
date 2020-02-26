##############################
## fitting point process model
##############################

# load spatstat package
library(spatstat) 
# test CSR, return true if reject
testRejectCSR <- function(pppPattern, alpha){
  if(pppPattern$n > 3){
    # Quatrat test
    # not used
    #n = floor(sqrt(pppPattern$n + 1 / 2))
    #qTest <- quadrat.test(pppPattern, nx = n, ny = n)
    #if (qTest$p.value < alpha){
    #  return (TRUE)
    #}
    ceResult = clarkevans.test(pppPattern, alternative='clustered',nsim=100)
    if(ceResult$p.value<p_th)
      return(TRUE)
  }
  return (FALSE)
}

# fit subregion to point process model
subRegionFit <- function(points, xwindow, ywindow, alpha){
  n = nrow(points)
  intensity = n/xwindow/ywindow
  p1 = c(n, intensity)
  p2 = c(0, 0, 0, 0,0)
  if(n>0){
    pattern <- ppp(points[,1], points[,2], c(0,xwindow), c(0,ywindow))
    #print(pattern)
    # test CSR
    if(testRejectCSR(pattern,p_th)){
      p2[1] = 1
      p2[2] = 1
      tryCatch( {
        myfit = kppm(pattern, ~1, "Thomas")
        p = parameters(myfit)
        p2[3:5] = c(p$kappa, p$scale, p$mu)
      } , warning = function(w) {
      }, error = function(e) {
      }, finally = {
      })
    }
    else{
      p2[1] = 0
    }
  }
  return(c(p1, p2))
}
# process model fitting for one slide
# get subregions of the slide using a moving window
# fit point pattern in sub regions to spatial process model
# record fitted parameters
processSlide <-function(path, originfile_name,pointPath,outname, alpha, window){
  
  ################
  p_th = 5E-2
  origin_coords <- as.data.frame(read.table(originfile_name))/125 #previously '/2000'
  Point_files = list.files(path=pointPath,pattern = 'txt', recursive=FALSE, full.names = FALSE)
  count = length(Point_files)
  p_all <- matrix(nrow=0,ncol=11)
  header = c("xstart","xend", "ystart","yend","n","intensity","notCSR", "fit", "kappa", "sigma2", "mu")
  colnames(p_all) = header
  
  
  for (id in 1:count) { #i: x
    # get x,y coordinates for origins
    x <- origin_coords[id,1]
    y <- origin_coords[id,2]
    print(id)
    window = 0.4 # mm
    filename = paste(slidepath,"DigitalPathology",id-1,"_allPoint.txt", sep="")
    if(file.info(filename)$size != 0){
      mydata <- try(read.table(filename))/2000
      p = subRegionFit(mydata, window, window, alpha)
      p_all = rbind(p_all, c(x,x+window,y,y+window,p))
    }
  }
  # write results to file
  write.csv(p_all, file = paste(path, outname, sep=""),
            row.names=FALSE)
}

# fit spatial point process model and record fitted parameters
spatstatFit <- function(subdir, alpha){
  
  # path to files and load data
  path = paste(subdir, "\\", sep="")
  print(path)
  filePattern = "^\\d+_\\d+.+csv$"
  patchCoords = list.files(path=path, pattern = filePattern
                           , recursive=FALSE, include.dirs = TRUE
                           , full.names = FALSE)
  
  # 7 columns for p_all: x, y, intensity, CSR? kappa, sigma2, mu
  p_all <- matrix(nrow=0,ncol=7)
  header = c("x", "y", "intensity", "notCSR", "kappa", "sigma2", "mu")
  colnames(p_all) = header
  
  for (patch in patchCoords) { #i: x
    coord = as.integer( strsplit(patch,"[_.]")[[1]][1:2])
    i = coord[1]
    j = coord[2]
    # loop through grid
    #print(filename)
    filePath = paste(path, patch, sep="")
    mydata <- read.csv(filePath)
    if(dim(mydata[1])>0){
      mypattern <- ppp(mydata[,1], mydata[,2], c(0,window), c(0,window))
      intensity = nrow(mydata)/window^2
      p1 = c(i, j, intensity)
      p2 = c(0, 0, 0, 0)
      # test CSR
      if(testRejectCSR(mypattern, alpha=alpha)){
        p2[1] = 1
        tryCatch( {
          myfit = kppm(mypattern, ~1, "Thomas")
          p = parameters(myfit)
          p2[2:4] = c(p$kappa, p$scale, p$mu)
        } , warning = function(w) {
        }, error = function(e) {
        }, finally = {
        })
      }
      else{
        p2[4] = 0
      }
      p_all = rbind(p_all, c(p1, p2))
    }
  }
  
  # write results to file
  write.csv(p_all, file = paste(path, "fittedResult_newCSR.csv", sep=""),
            row.names=FALSE)
}

##############################
### batch shape descriptor
##############################


# packages
library(largeVis)
library(alphahull)
library(stats)

# class AshapePol:
# two member variables:
# goodAshape: TRUE/FALSE
# vert: indices of vertices
setClass(Class="AshapePol",
         representation(
           goodAshape="logical",
           vert="numeric"
         )
)

## get the polygon boundary points in right order
library(sp)
library(igraph)

getAlphaShapeInOrder <- function(ashapeX){
  returnVal=new("AshapePol", goodAshape=TRUE, vert=0)
  
  if (nrow(ashapeX$edges)==0) {
    #stop("Graph not connected")
    returnVal@goodAshape = FALSE
  }
  
  ashapeGraph = graph.edgelist(cbind(as.character(ashapeX$edges[, "ind1"]), 
                                     as.character(ashapeX$edges[, "ind2"])), directed = FALSE)
  
  
  
  #plot(ashapeGraph)
  ## check if is closed
  if (!is.connected(ashapeGraph)) {
    #stop("Graph not connected")
    returnVal@goodAshape = FALSE
  }
  if (any(degree(ashapeGraph) != 2)) {
    #stop("Graph not circular")
    returnVal@goodAshape = FALSE
  }
  if (clusters(ashapeGraph)$no > 1) {
    #stop("Graph composed of more than one circle")
    returnVal@goodAshape = FALSE
  }
  
  if(! returnVal@goodAshape){
    return(returnVal)
  }
  
  cutg = ashapeGraph - E(ashapeGraph)[1]
  # find chain end points
  ends = names(which(degree(cutg) == 1))
  path = get.shortest.paths(cutg, ends[1], ends[2])$vpath[[1]]
  # this is an index into the points
  pathX = as.numeric(V(ashapeGraph)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  
  aShapePath <- Polygon(ashapeX$x[pathX, ])  
  # polygon(X[pathX, 1], X[pathX,2], col = "gray", border = "red")
  
  ## test within polygon
  inAshape = point.in.polygon(ashapeX$x[,1],ashapeX$x[,2],ashapeX$x[pathX,1],ashapeX$x[pathX,2])
  if (any(inAshape==0)){
    #stop("point outside alpha shape")
    returnVal@goodAshape = FALSE
  }else{
    returnVal@vert = pathX
  }
  return(returnVal)
}

infer_invasiveFront <- function(subCtr1, subCtr2, subCtr3, subCtr4, subCtr5){
  
  subCtr1 <- as.matrix(subCtr1)
  subCtr2 <- as.matrix(subCtr2)
  subCtr3 <- as.matrix(subCtr3)
  subCtr4 <- as.matrix(subCtr4)
  subCtr5 <- as.matrix(subCtr5)
  # Assign names
  file_name <- c('x_Coords', 'y_Coords')
  names(subCtr1) <- file_name
  names(subCtr2) <- file_name
  names(subCtr3) <- file_name
  names(subCtr4) <- file_name
  names(subCtr5) <- file_name
  
  
  contour_all <- rbind(subCtr1,subCtr2, subCtr3, subCtr4, subCtr5)
  # Concave hull?
  concave_pts_subCtr1 <- as.data.frame(concaveman(subCtr1, concavity = 2))
  
  concave_pts_subCtr2 <- as.data.frame(concaveman(subCtr2, concavity = 2))
  
  concave_pts_subCtr3 <- as.data.frame(concaveman(subCtr3, concavity = 2))
  
  concave_pts_subCtr4 <- as.data.frame(concaveman(subCtr4, concavity = 2))
  
  concave_pts_subCtr5 <- as.data.frame(concaveman(subCtr5, concavity = 2))
  # Form Polygons
  Polygon_subCtr1 <- Polygon(concave_pts_subCtr1)
  Polygon_subCtr2 <- Polygon(concave_pts_subCtr2)
  Polygon_subCtr3 <- Polygon(concave_pts_subCtr3)
  Polygon_subCtr4 <- Polygon(concave_pts_subCtr4)
  Polygon_subCtr5 <- Polygon(concave_pts_subCtr5)
  # Turn polygon to spatialpolygon
  
  sp_subCtr1 <- SpatialPolygons(list(Polygons(list(Polygon(Polygon_subCtr1)),1)))
  
  sp_subCtr2 <- SpatialPolygons(list(Polygons(list(Polygon(Polygon_subCtr2)),1)))
  
  sp_subCtr3 <- SpatialPolygons(list(Polygons(list(Polygon(Polygon_subCtr3)),1)))
  
  sp_subCtr4 <- SpatialPolygons(list(Polygons(list(Polygon(Polygon_subCtr4)),1)))
  
  sp_subCtr5 <- SpatialPolygons(list(Polygons(list(Polygon(Polygon_subCtr5)),1)))
  
  
  # Find intersections
  sp_intsc <- gIntersection(sp_subCtr1, sp_subCtr2)
  sp_intsc <- gIntersection(sp_intsc, sp_subCtr3)
  sp_intsc <- gIntersection(sp_intsc, sp_subCtr4)
  sp_intsc <- gIntersection(sp_intsc, sp_subCtr5)
  
  # invert the data frame
  inter_Coords <- as.data.frame(sp_intsc@polygons[[1]]@Polygons[[1]]@coords)
  inter_Coords <- inter_Coords[order(nrow(inter_Coords):1),] #invert row order
  
  # Create buffer to represent invasive front
  inter_Coords_large <- gBuffer(sp_intsc, width=  0.25) 
  inter_Coords_small <- gBuffer(sp_intsc, width= -0.25) 
  inter_Coords_large <- as.data.frame(inter_Coords_large@polygons[[1]]@Polygons[[1]]@coords)
  inter_Coords_small <- as.data.frame(inter_Coords_small@polygons[[1]]@Polygons[[1]]@coords)
  Band_contour = inter_Coords_large[order(nrow(inter_Coords_large):1),] #invert row order
  Band_contour <- rbind(Band_contour, inter_Coords_small)
  List <- list('1' = inter_Coords, '2' = inter_Coords_large, '3' = inter_Coords_small, '4' =Band_contour, 'invasiveCoords' = Band_contour, 'tumor' = inter_Coords_small)
  return(List)
  #return(c(inter_Coords,inter_Coords_large, inter_Coords_small,Band_contour))
}

pointCharacterize <- function(coords_data, invasive_poly, tumor_poly){
  normal_set <-  matrix(nrow=,ncol=2)
  invasive_set <-  matrix(nrow=0,ncol=2)
  tumor_set <-  matrix(nrow=0,ncol=2)
  
  len <- nrow(coords_data)
  for (idx in seq(1,len,1)){
    x_coords <- coords_data[idx,1]
    y_coords <- coords_data[idx,2]
    Coords <- cbind(x_coords, y_coords)
    # Note: how to determine the boundary condition
    if(point.in.polygon(x_coords, y_coords, invasive_poly[,1], invasive_poly[,2]) != 0){
      invasive_set <- rbind(invasive_set,Coords)
    }
    else if(point.in.polygon(x_coords, y_coords, tumor_poly[,1], tumor_poly[,2]) == 1){
      tumor_set <- rbind(tumor_set,Coords)
    }
    else{
      normal_set <- rbind(normal_set,Coords)
    }
    
  }
  Charac_point <- list('normal' = normal_set, 'tumor' = tumor_set, 'invasive' = invasive_set)
  
  return(Charac_point)
}

infer_invasiveFrontColror <- function(subCtr1, subCtr2, subCtr3, subCtr4, subCtr5,Color_coords, threshold, row_num, band_width){
  subCtr1 <- as.matrix(subCtr1)
  subCtr2 <- as.matrix(subCtr2)
  subCtr3 <- as.matrix(subCtr3)
  subCtr4 <- as.matrix(subCtr4)
  subCtr5 <- as.matrix(subCtr5)
  Color_coords <- rbind(c(0.5, 0.5),Color_coords)
  
  # Assign names
  file_name <- c('x_Coords', 'y_Coords')
  names(subCtr1) <- file_name
  names(subCtr2) <- file_name
  names(subCtr3) <- file_name
  names(subCtr4) <- file_name
  names(subCtr5) <- file_name
  names(Color_coords) <- file_name
  
  concave_pts_CD3 <- as.data.frame(concaveman(subCtr1, concavity = 2))
  
  concave_pts_CD4 <- as.data.frame(concaveman(subCtr2, concavity = 2))
  
  concave_pts_CD8 <- as.data.frame(concaveman(subCtr3, concavity = 2))
  
  concave_pts_CD20 <- as.data.frame(concaveman(subCtr4, concavity = 2))
  
  concave_pts_FoxP3 <- as.data.frame(concaveman(subCtr5, concavity = 2))
  
  Charac_color_dist <- as.data.frame(point.in.polygon(Color_coords[,1],Color_coords[,2], concave_pts_CD3[,1], concave_pts_CD3[,2]))
  Charac_color_dist_CD4 <- as.data.frame(point.in.polygon(Color_coords[,1],Color_coords[,2], concave_pts_CD4[,1], concave_pts_CD4[,2]))
  Charac_color_dist_CD8 <- as.data.frame(point.in.polygon(Color_coords[,1],Color_coords[,2], concave_pts_CD8[,1], concave_pts_CD8[,2]))
  Charac_color_dist_CD20 <- as.data.frame(point.in.polygon(Color_coords[,1],Color_coords[,2], concave_pts_CD20[,1], concave_pts_CD20[,2]))
  Charac_color_dist_FoxP3 <- as.data.frame(point.in.polygon(Color_coords[,1],Color_coords[,2], concave_pts_FoxP3[,1], concave_pts_FoxP3[,2]))
  Charac_color_dist_int <- Charac_color_dist_CD4 + Charac_color_dist_CD8 + Charac_color_dist_CD20 + Charac_color_dist_FoxP3 + Charac_color_dist
  
  # Modify picture
  Charac_color_dist_int <- as.numeric(unlist(Charac_color_dist_int))
  
  Color_img <- matrix(Charac_color_dist_int, nrow = row_num, byrow = TRUE)
  
  #Color_img <- apply(Color_img, 2, rev)
  Color_img <- t(Color_img)
  
  
  Color_img[Color_img < threshold] <- 0
  Color_img[Color_img >= threshold] <- 5
  
  # Get index 
  indice_x <- row(Color_img)[which(!Color_img == 0)]
  indice_y <- col(Color_img)[which(!Color_img == 0)]
  indice <- as.matrix(cbind(indice_x,indice_y))
  concave_pts <- as.data.frame(concaveman(indice, concavity = 2))/125
  
  Polygon_ctr <- Polygon(concave_pts)
  sp_reference <- SpatialPolygons(list(Polygons(list(Polygon(Polygon_ctr)),1)))
  # invert the data frame
  
  ref_Coords <- as.data.frame(sp_reference@polygons[[1]]@Polygons[[1]]@coords)
  ref_Coords = ref_Coords[order(nrow(ref_Coords):1),] #invert row order
  # Create buffer to represent invasive front
  
  sp_Coords_large<- gBuffer(sp_reference, width= band_width) 
  sp_Coords_small <- gBuffer(sp_reference, width= -band_width)
  
  sp_Coords_large <- as.data.frame(sp_Coords_large@polygons[[1]]@Polygons[[1]]@coords)
  sp_Coords_small <- as.data.frame(sp_Coords_small@polygons[[1]]@Polygons[[1]]@coords)
  Band_contour = sp_Coords_large[order(nrow(sp_Coords_large):1),] #invert row order
  
  Band_contour <- rbind(Band_contour, sp_Coords_small)
  List <- list('reference' = ref_Coords, 'outter' = sp_Coords_large, 'inner' = sp_Coords_small, 'invasiveCoords' = Band_contour, 'tumor' = sp_Coords_small, 'indice' = indice)
  return(List)
}

# Region: normal? tumor? invasive?
fitResultChara <- function(path,case_idx, marker_idx, region){
  fitResults_subtype <- read.csv(path)
  # Claim path to save files
  dir <- paste(case_idx, marker_idx,'/', sep = '')
  for(flag in seq(1,4,1)){
    # get the type
    type <- switch(flag, 'kappa','sigma2','mu','intensity')
    
    data_sub <- fitResults_subtype[fitResults_subtype$notCSR != 0,]
    data_sub <- data_sub[complete.cases(data_sub), ]
    
    data_sub <- as.data.frame(data_sub[,type])
    
    # filter data
    data_sub <- data_sub[data_sub[,1] != 0,]
    
    file_name <- paste(dir, type, '_', region,'.csv', sep = '')
    write.csv(data_sub,  file = file_name)
    
  }
}
# This function characterizing clusters based on their location.
clusterCharacterize <- function(locInfo,shape_outlineFile ){
  # Get region statistics
  n_stroma <- locInfo[, 1]
  n_front <- locInfo[, 2]
  n_tumor <- locInfo[, 3]
  
  # Get the number of cluster
  Cluster_num <- shape_outlineFile[nrow(shape_outlineFile),1]
  
  shape_outlineFileDist <- data.frame(matrix(nrow = nrow(shape_outlineFile), ncol = 1))
  
  Cluster_dist <- data.frame(matrix(nrow = Cluster_num, ncol = 2))
  
  indicator <- 1
  for(id in seq(1, Cluster_num, 1)){
    Cluster_dist[id, 1] <- id
    # Extract sub clusters
    sub_Cluster <- shape_outlineFile[(shape_outlineFile[,1] == id), -1]
    step_size <- nrow(sub_Cluster)
    End_idx <- step_size + indicator - 1
    # Define end index for one assignment
    if(locInfo[id,2] != 0){
      # 1: invasive front, 2: tumor, 3: stromal
      Cluster_dist[id, 2] <- 1
      shape_outlineFileDist[c(indicator: End_idx), 1] <- 1
    }
    else if(locInfo[id,3] != 0){
      Cluster_dist[id, 2] <- 2
      shape_outlineFileDist[c(indicator : End_idx), 1] <- 2
      
    } else{
      Cluster_dist[id, 2] <- 3
      shape_outlineFileDist[c(indicator: End_idx),1] <- 3
    }
    indicator <- indicator + step_size
  }
  shape_outlineFileDist <- cbind(shape_outlineFile, shape_outlineFileDist)
  names(shape_outlineFileDist) <- c('cluster_id','x','y','location')
  Cluster_info <-  list('Cluster_dist' = Cluster_dist, 'shape_outlineFileDist' = shape_outlineFileDist)
  
  return(Cluster_info)
}


DistanceProfile <- function(Invasive_inner, Invasive_outer,Invasive_Polygon, Tumor_Pts, step_size){
  # Convet format to sp objects
  sp_Invasive <- SpatialPolygons(list(Polygons(list(Polygon(Invasive_Polygon)),1)))
  
  
  sp_TumorPts <- SpatialPoints(Tumor_Pts)
  
  # Claim data frame for data storage
  DistFile <- as.data.frame(matrix(nrow = 0, ncol = 1))
  Densi_file <- as.data.frame(matrix(nrow = 0, ncol = 3))
  names(DistFile) <- 'Distance'
  names(Densi_file) <- c('Distance','Area','Density')
  len <- nrow(Tumor_Pts)
  
  # Loop through the tumor coordinates to calculate distances
  for(ind in seq(1, len, 1)){
    Dist <- gDistance(sp_Invasive, sp_TumorPts[ind,])
    names(Dist) <- 'Distance'
    DistFile <- rbind(DistFile,Dist)
  }
  Tumor_Pts <- cbind(Tumor_Pts, DistFile)
  # Calculate section num
  n <- 1
  
  Hallmark <- TRUE
  bound_low  <- 0
  bound_high <- step_size
  while(Hallmark == TRUE){
    
    # Extract idealed points
    Tumor_sec <- as.matrix(Tumor_Pts[bound_low <= Tumor_Pts[,3] & Tumor_Pts[,3] < bound_high, -3])
    if(nrow(Tumor_sec) == 0){
      Hallmark <- FALSE
    }else{
      count <- nrow(Tumor_sec)
      #Calculate polygon area
      Tumor_sec <- as.data.frame(concaveman(Tumor_sec, concavity = 2))
      Tumor_sec <- SpatialPolygons(list(Polygons(list(Polygon(Tumor_sec)),1)))
      Area <- gArea(Tumor_sec)
      
      #Calculate density
      Density <- count/Area # Unit: cells/ mm^2
      
      sub_Density <- cbind(n*step_size, Area,Density)
      names(sub_Density) <- c('Distance','Area' ,'Density')
      Densi_file <- rbind(Densi_file, sub_Density) 
      bound_low <- bound_high
      bound_high <- bound_high + step_size
      n <- n + 1
      print(n)
      
    }
    
    
  }
  Profile <-  list('Densi_file' = Densi_file, 'Tumor_Pts' = Tumor_Pts)
  
  return(Profile)
}

bivarAnalysis <- function(from_Coords, to_Coords, InvasiveFile, TumorFile, step_size){
  
  # Declare global variable
  AUC_dist <- data.frame(matrix(nrow = 0, ncol = 2))
  names(AUC_dist) <- c('Section', 'AUC')
  # Characterize point (CD3)
  Charac_to <- pointCharacterize(to_Coords, InvasiveFile, TumorFile)
  Tumor_to <- data.frame(Charac_to$tumor)
  CD3_dist <- DistanceProfile(InvasiveFile, Tumor_to, 0.15)
  CD3_dist <- CD3_dist$Tumor_Pts
  CD3_dist$mark <- 1
  Tumor_to$mark <- 1
  
  Invasive_to <- data.frame(Charac_to$invasive)
  Invasive_to$mark <- 1
  
  Normal_to <- data.frame(Charac_to$normal)
  Normal_to$mark <- 1
  
  # Characterize point (CD8)
  Charac_from <- pointCharacterize(from_Coords, InvasiveFile, TumorFile)
  Tumor_from <- data.frame(Charac_from$tumor)
  CD8_dist <- DistanceProfile(InvasiveFile, Tumor_from, 0.15)
  CD8_dist <- CD8_dist$Tumor_Pts
  
  CD8_dist$mark <- 2
  Tumor_from$mark <- 2
  
  Invasive_from <- data.frame(Charac_from$invasive)
  Invasive_from$mark <- 2
  
  Normal_from <- data.frame(Charac_from$normal)
  Normal_from$mark <- 2
  
  # Perform bivariate analysis - Gcross function
  Tumor_comb <- rbind(Tumor_to, Tumor_from)
  Invasive_comb <- rbind(Invasive_to, Invasive_from)
  Normal_comb <- rbind(Normal_to, Normal_from)
  
  Tumor <- ppp(Tumor_comb[,1],Tumor_comb[,2], c(0,32), c(0,25),marks=factor(Tumor_comb$mark))
  Invasive <- ppp(Invasive_comb[,1],Invasive_comb[,2], c(0,32), c(0,25),marks=factor(Invasive_comb$mark))
  Normal <- ppp(Normal_comb[,1],Normal_comb[,2], c(0,32), c(0,25),marks=factor(Normal_comb$mark))
  
  Gcross_tumor <- data.frame(Gcross(Tumor, '2', '1', seq(0, 0.1,0.00002)))
  Gcross_invasive <- data.frame(Gcross(Invasive, '2', '1', seq(0, 0.1,0.00002)))
  Gcross_normal <- data.frame(Gcross(Normal, '2', '1', seq(0, 0.1,0.00002)))
  
  # AUC related to distance to boundary
  
  
  bound_low <- 0
  bound_high <- bound_low + step_size
  n <- 1
  Tumor_sec_CD8 <- CD8_dist[bound_low <= CD8_dist[,3] & CD8_dist[,3] < bound_high, ]
  Tumor_sec_CD3 <- CD3_dist[bound_low <= CD3_dist[,3] & CD3_dist[,3] < bound_high, ]
  
  while(nrow(Tumor_sec_CD3) != 0 & nrow(Tumor_sec_CD8) != 0){
    Tumor_sec_CD8 <- CD8_dist[bound_low <= CD8_dist[,3] & CD8_dist[,3] < bound_high,]
    Tumor_sec_CD3 <- CD3_dist[bound_low <= CD3_dist[,3] & CD3_dist[,3] < bound_high,]
    # get numbers to determine termination
    count_CD8 <- nrow(Tumor_sec_CD8)
    count_CD3 <- nrow(Tumor_sec_CD3)
    names(Tumor_sec_CD3) <- c('x','y','distance','mark')
    names(Tumor_sec_CD8) <- c('x','y','distance','mark')
    
    if(count_CD8 != 0 & count_CD3 != 0){
      # Gcross function for sub data
      Tumor_sub <- rbind(Tumor_sec_CD3, Tumor_sec_CD8)
      Tumor_sub <- ppp(Tumor_sub[,1],Tumor_sub[,2], c(0,32), c(0,25),marks=factor(Tumor_sub$mark))
      AUC_sub <- as.data.frame(trapz(data.frame(Gcross(Tumor_sub, '2', '1', seq(0, 0.02,0.00002)))[,1]*1000, data.frame(Gcross(Tumor_sub, '2', '1', seq(0, 0.02,0.00002)))[,5]))
      AUC_sub <- cbind(n, AUC_sub)
      names(AUC_sub) <- c('Section', 'AUC')
      
      AUC_dist <- rbind(AUC_dist, AUC_sub)
      print(AUC_sub)
      bound_low <- bound_high
      bound_high <- bound_high + step_size
      n <- n + 1
    }
  }
  AUC_tumor <- as.data.frame(trapz(data.frame(Gcross(Tumor, '2', '1', seq(0, 0.02,0.00002)))[,1]*1000, data.frame(Gcross(Tumor, '2', '1', seq(0, 0.02,0.00002)))[,2]))
  AUC_invasive <- as.data.frame(trapz(data.frame(Gcross(Invasive, '2', '1', seq(0, 0.02,0.00002)))[,1]*1000, data.frame(Gcross(Invasive, '2', '1', seq(0, 0.02,0.00002)))[,2]))
  AUC_normal <- as.data.frame(trapz(data.frame(Gcross(Normal, '2', '1', seq(0, 0.02,0.00002)))[,1]*1000, data.frame(Gcross(Normal, '2', '1', seq(0, 0.02,0.00002)))[,2]))
  
  AUC <- cbind(AUC_tumor, AUC_invasive, AUC_normal)
  names(AUC) <- c('Tumor', 'Invasive', 'Normal')
  bivarResults <-  list('Gcross_tumor' = Gcross_tumor, 'Gcross_invasive' = Gcross_invasive, 'Gcross_normal' = Gcross_normal, 'AUC' = AUC_dist)
  
}
