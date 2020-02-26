############################################################################
#### This script is used for clustering morphometric analysis (Fig. 8) #####
############################################################################

library('R.matlab')
library('ggplot2')
library('largeVis')
library('plyr')
library('reshape2')
library(magick)
library(ggplotify)
library(tidyverse)
library(gginnards)
library(grid)
library(cowplot)
library(alphahull)
library(plyr)
# for formatting grids
library(ggpubr)
source('MiFunction.R')
# Change case path here
plotStat <- data.frame(matrix(nrow = 0, ncol = 9))
names(plotStat) <- c('Density', 'a-Shape area', 'Fitted ellipse area', 'Circularity', 'Convexity', 'Eccentricity', 'location', 'marker', 'case')
segSuffix = "_Points.mat"

for(case in seq(1,6,1)){
  clusterPath <- switch(case,'./Case_0/', './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  clusterPath_name <- switch(case,'Case 1A', 'Case 1B' , 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  
  ########## requires coordinates that define regions
  if(case == 1){
    InvasiveFile <- readRDS(paste(clusterPath, 'Invasive_1a.rds', sep = ''))
    TumorFile <- readRDS(paste(clusterPath, 'Tumor_1a.rds', sep = ''))
    NormalFile <- readRDS(paste(clusterPath, 'Normal_1a.rds', sep = ''))
    }
  if(case == 2){
    InvasiveFile <- readRDS(paste(clusterPath, 'Invasive_1b.rds', sep = ''))
    
    TumorFile <- readRDS(paste(clusterPath, 'Tumor_1b.rds', sep = ''))
    
    NormalFile <- readRDS(paste(clusterPath, 'Normal_1b.rds', sep = ''))
    }
  if(case > 2){
    
    InvasiveFilePath <- paste(clusterPath, 'Invasive.rds', sep = '')
    InvasiveFile <- readRDS(InvasiveFilePath)
  
    TumorFilePath <- paste(clusterPath, 'Tumor.rds', sep = '')
    TumorFile <- readRDS(TumorFilePath)
    
    NormalFile <- readRDS(paste(clusterPath, 'Normal.rds', sep = ''))
    }
  for(flag in seq(1,5,1)){
    Current_lab <- switch(flag, 'CD3','CD4','CD8','CD20','FoxP3')
    print(paste('Processing', Current_lab))
    clusterOutPath = paste(clusterPath, Current_lab,"/", sep = '')
    filename <- paste(clusterOutPath,Current_lab,segSuffix,sep="")
    outfile = paste(clusterOutPath, Current_lab,"_shape.csv", sep="")
    outfileOutline = paste(clusterOutPath, Current_lab, "_shape_outline.csv", sep="")
  
  
  if (file.exists(filename)){
    allPoints <-  as.data.frame(readMat(filename))/125
  }else{
    print(filename)
    print('file not exist')
    next
  }
  subPoints <- allPoints

  # clustering
  minPts = 30
  K = 4
  dat <- as.matrix(subPoints[, 1:2])
  cd8vis <- largeVis(t(dat), dim=2, K = K, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)
  clusters <- hdbscan(cd8vis, K=K, minPts = minPts, verbose = FALSE, threads = 1)
  
  # for each cluster
  
  hullStats <- matrix(nrow = 0, ncol = 4) # n, perimeter, area, convexity
  ellStats <- matrix(nrow = 0, ncol = 4) # x, y, major, minor
  centStats <- matrix(nrow = 0, ncol = 3) # xc, yc, moment of inertia
  regionStats<- matrix(nrow = 0, ncol = 3) # number in tumor, stroma, front
  locInfo <- as.data.frame(matrix(nrow = 0, ncol  = 3)) # location info
  polygonInfo <- matrix(nrow = 0, ncol = 3) #clusterID, x, y
  
  alpha0 <- 0.01
  alphaStep <- 0.01
  
  maxCluster = max(as.numeric(clusters$clusters), na.rm = TRUE)
  minCluster = min(as.numeric(clusters$clusters), na.rm = TRUE)
  notNA = !is.na(clusters$clusters)
  
  for(cid_target in 1:maxCluster) {
    print(cid_target)
    cidMatch <- notNA & clusters$clusters == cid_target
    goodAlpha = FALSE
    alpha <- alpha0
    #X<-subPoints[cidMatch,]
    X <- subPoints[cidMatch,1:2]
    X <- unique(X)
    # try different alpha until shape covers all points
    while(!goodAlpha){
      print(alpha)
      alpha <- alpha + alphaStep
      ashapeX <- ashape(X, alpha = alpha)
      #plot(ashapeX, xlab = "x", ylab = "y")
      res=getAlphaShapeInOrder(ashapeX)
      goodAlpha <- res@goodAshape
    }
    
    nPt <-nrow(X)
    
    ### get alpha-perimeter
    Pconc <- ashapeX$length
    
    ### get alpha-area:
    alphaShape <- Polygon(X[res@vert,])
    Aconc <- alphaShape@area
    ### polygon outline
    polygonInfo = rbind(polygonInfo, 
                        cbind(rep(cid_target, length(alphaShape@coords[,1])),
                              alphaShape@coords))
    hullPts <- chull(X)
    hullPts <- c(hullPts, hullPts[1])
    
    polygon_set <- matrix(nrow=nrow(X[hullPts,]),ncol=3)
    
    ### get convex hull
    
    #lines(X[hullPts, ])
    convShape <- Polygon(X[hullPts, ])  
    ## plot the polygon of convex hull !!!!!!!!

    polygon_set[,1] <- X[hullPts, 1]
    polygon_set[,2] <- X[hullPts, 2]
    polygon_set[,3] <- cid_target
    polygon_set <- data.frame(polygon_set)

    ### convex hull area
    Aconv<-convShape@area
    
    hullStats = rbind(hullStats, as.matrix(t(c(nPt, Pconc, Aconc, Aconv))))
    
    ## fitting ellipse
    # we can get center and covariance matrix from the cluster points:
    ellFit <- cov.wt(X)
    # eigen vectors indicate orientations
    ellEigen <- eigen(ellFit$cov)
    # eigen values are associated with major/minor axes length (half)

    # alternative: chi square distribution
    axes <- sqrt(ellEigen$values * qchisq(.95, df=2))
    
    ellStats = rbind(ellStats, as.matrix(t(c(ellFit$center[1],ellFit$center[2],
                                             axes[1], axes[2]))))
    #centroid and moment of inersia
    centStats = rbind(centStats, as.matrix(t(c(colMeans(X), sum(sweep(X,2,colMeans(X))^2)))))
    
    regionStats = rbind(regionStats, colSums(subPoints[cidMatch,], dims=1)[3:5])
    
    # Count numbers of points associated to each cluster
    front_num <- data.frame(point.in.polygon(X[,1],X[,2],InvasiveFile[,1],InvasiveFile[,2]))
    front_num <- length(which(front_num[,1] == 1))
    tumor_num <- data.frame(point.in.polygon(X[,1],X[,2],TumorFile[,1],TumorFile[,2]))
    tumor_num <- length(which(tumor_num[,1]== 1))

    # calculate stroma number 
    stroma_num <- data.frame(point.in.polygon(X[,1],X[,2],NormalFile[,1],NormalFile[,2]))
    
    stroma_num <- length(which(stroma_num[,1]== 1))
    locInfo <- rbind(locInfo, c(stroma_num, front_num, tumor_num)) 
  }
  
  hullStats=data.frame(hullStats)
  colnames(hullStats)<-c("n", "perimConcave","areaConcave","areaConvex")
  
  shapeStats = hullStats
  shapeStats$n_tumor = locInfo[,3] #regionStats[,1]
  shapeStats$n_stroma = locInfo[,1] #regionStats[,2]
  shapeStats$n_front <- locInfo[,2] #regionStats[,3]
  shapeStats$circ = shapeStats$areaConcave*4*pi/(shapeStats$perimConcave)^2
  shapeStats$conv = shapeStats$areaConcave/shapeStats$areaConvex
  shapeStats$density = shapeStats$n/shapeStats$areaConcave
  shapeStats$x = ellStats[,1]
  shapeStats$y = ellStats[,2]
  shapeStats$major = ellStats[,3]
  shapeStats$minor = ellStats[,4]
  shapeStats$eccentricity = sqrt(1-(shapeStats$minor/shapeStats$major)^2)
  shapeStats$xc = centStats[,1]
  shapeStats$yc = centStats[,2]
  shapeStats$moi = centStats[,3]
  

  # Characterize clusters
  ptsPath <- paste(clusterPath, Current_lab, '/', Current_lab, '_Points.mat', sep = '')
  Cluster_info <- clusterCharacterize(locInfo, data.frame(polygonInfo))
  Cluster_dist <- Cluster_info$Cluster_dist
  
  # Assign location
  shapeStats$location = Cluster_dist[,2]
    
  # Data summary
  polygonDF = data.frame(polygonInfo)
  colnames(polygonDF)<-c("cluster_ID", "x","y")
  
  # 
  write.csv(shapeStats, file = outfile)
  write.csv(polygonDF, file = outfileOutline, row.names=FALSE)
  # calculate ellipse area
  areaEllipse <- pi*shapeStats$major*shapeStats$minor
  plotStat_sub <- cbind(shapeStats$density, alphaShape@area, areaEllipse, shapeStats$circ, shapeStats$conv, shapeStats$eccentricity,shapeStats$location, Current_lab, clusterPath_name)
  colnames(plotStat_sub) <- c('Density', 'a-Shape area', 'Fitted ellipse area', 'Circularity', 'Convexity', 'Eccentricity', 'location', 'marker', 'case')
  plotStat <- rbind(plotStat, plotStat_sub)
  }

}
saveRDS(plotStat, './plotStat.rds')
plotStat <- readRDS('./plotStat.rds')

# plot plotStat data
plotStat_bac <- plotStat
plotStat_bac$location <- as.numeric(as.character(plotStat_bac$location))
plotStat_bac$Density <- as.numeric(as.character(plotStat_bac$Density))
plotStat_bac$`a-Shape area` <- as.numeric(as.character(plotStat_bac$`a-Shape area`))
plotStat_bac$`Fitted ellipse area` <- as.numeric(as.character(plotStat_bac$`Fitted ellipse area`))
plotStat_bac$Circularity <- as.numeric(as.character(plotStat_bac$Circularity))
plotStat_bac$Convexity <- as.numeric(as.character(plotStat_bac$Convexity))
plotStat_bac$Eccentricity <- as.numeric(as.character(plotStat_bac$Eccentricity))


# Change x-axis name
plotStat_bac[plotStat_bac$location == 1,]$location <- 'IF'
plotStat_bac[plotStat_bac$location == 2,]$location <- 'CT'
plotStat_bac[plotStat_bac$location == 3,]$location <- 'N'

den_plotStat <- plotStat_bac[, c(1,7,8,9)]
ashape_plotStat <- plotStat_bac[, c(2,7,8,9)]
ellipse_plotStat <- plotStat_bac[, c(3,7,8,9)]
circ_plotStat <- plotStat_bac[, c(4,7,8,9)]
conv_plotStat <- plotStat_bac[, c(5,7,8,9)]
eccen_plotStat <- plotStat_bac[, c(6,7,8,9)]

# characterize the shape
mean_dat <- matrix(nrow = 0, ncol = 5)
max_dat <- matrix(nrow = 0, ncol = 5)
nodular_dat <- matrix(nrow = 0, ncol = 5)
irregular_dat <- matrix(nrow = 0, ncol = 5)
for(ind in seq(1, 4, 1)){
  criteria_name <- switch(ind, 'Density', 'Density of cluster', 'nodular count', 'Elongated count')
  label <- matrix(nrow = 0, ncol = 6)
  if(isTRUE(ind == 1)){
    for(case in seq(1, 6, 1)){
      case_name <- switch(case , 'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
      for(label in seq(1, 5, 1)){
        label_name <- switch(label , 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
        Stat_N <- (plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'N' & plotStat_bac$marker == label_name,])
        Stat_CT <- (plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'CT' & plotStat_bac$marker == label_name,])
        Stat_IF <- (plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'IF' & plotStat_bac$marker == label_name,])
        
        if(isTRUE( nrow(Stat_N) !=0) ){
          meanDen_N <- mean(Stat_N$Density)
        } 
        if(isTRUE(nrow(Stat_N) == 0)){
          meanDen_N = 0
          }
        
        if(isTRUE(nrow(Stat_CT) !=0)){
          meanDen_CT <- mean(Stat_CT$Density)
        } 
        if(isTRUE(nrow(Stat_CT) == 0)){
          meanDen_CT = 0
          }
        
        if(isTRUE(nrow(Stat_IF) !=0)){
          meanDen_IF <- mean(Stat_IF$Density)
        } 
        if(isTRUE(nrow(Stat_IF) == 0)){
          meanDen_IF = 0
          }
        
        subDat_N <- cbind(meanDen_N, criteria_name, case_name, label_name, 'N')
        subDat_CT <- cbind(meanDen_CT, criteria_name, case_name, label_name, 'CT')
        subDat_IF <- cbind(meanDen_IF, criteria_name, case_name, label_name, 'IF')
        
        caseSub_stat <- data.frame(rbind(subDat_N, subDat_CT, subDat_IF))
        colnames(caseSub_stat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        mean_dat <- data.frame(rbind(mean_dat, caseSub_stat))
        colnames(mean_dat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        
        }
    }
  }
  if(isTRUE(ind == 2)){
    for(case in seq(1, 6, 1)){
      case_name <- switch(case , 'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
      N_area <- switch(case, 27.898, 25.524, 107.2385, 4.6115, 166.5614, 113.6515)
      IF_area <- switch(case, 16.9748, 144.5, 33.2, 10.5, 56.4, 17.8)
      CT_area <- switch(case, 47.5, 25.32, 103.4, 238.3, 244.3, 90.6)
      for(label in seq(1, 5, 1)){
        label_name <- switch(label , 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
        Stat_N <- length(which(plotStat_bac$case == case_name & plotStat_bac$location == 'N' & plotStat_bac$marker == label_name))/N_area
        Stat_CT <- length(which(plotStat_bac$case == case_name & plotStat_bac$location == 'CT' & plotStat_bac$marker == label_name))/CT_area
        Stat_IF <- length(which(plotStat_bac$case == case_name & plotStat_bac$location == 'IF' & plotStat_bac$marker == label_name))/IF_area
        
        
        
        subDat_N <- cbind(Stat_N, criteria_name, case_name, label_name, 'N')
        subDat_CT <- cbind(Stat_CT, criteria_name, case_name, label_name, 'CT')
        subDat_IF <- cbind(Stat_IF, criteria_name, case_name, label_name, 'IF')
        
        caseSub_stat <- data.frame(rbind(subDat_N, subDat_CT, subDat_IF))
        colnames(caseSub_stat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        max_dat <- data.frame(rbind(max_dat, caseSub_stat))
        colnames(max_dat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        
      }
    }
  }
  if(isTRUE(ind == 3)){
    for(case in seq(1, 6, 1)){
      case_name <- switch(case , 'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
      for(label in seq(1, 5, 1)){
        label_name <- switch(label , 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
        Stat_N <- plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'N' & plotStat_bac$marker == label_name,]
        Stat_CT <- plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'CT' & plotStat_bac$marker == label_name,]
        Stat_IF <- plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'IF' & plotStat_bac$marker == label_name,]
        
        circular_N <- nrow(Stat_N[Stat_N$Eccentricity < 0.8 & Stat_N$Convexity > 0.8 & Stat_N$Circularity > 0.5,])
        circular_CT <- nrow(Stat_CT[Stat_CT$Eccentricity < 0.8 & Stat_CT$Convexity > 0.8 & Stat_CT$Circularity > 0.5,])
        circular_IF <- nrow(Stat_IF[Stat_IF$Eccentricity < 0.8 & Stat_IF$Convexity > 0.8 & Stat_IF$Circularity > 0.5,])
        
        subDat_N <- cbind(circular_N, criteria_name, case_name, label_name, 'N')
        subDat_CT <- cbind(circular_CT, criteria_name, case_name, label_name, 'CT')
        subDat_IF <- cbind(circular_IF, criteria_name, case_name, label_name, 'IF')
        
        caseSub_stat <- data.frame(rbind(subDat_N, subDat_CT, subDat_IF))
        colnames(caseSub_stat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        nodular_dat <- data.frame(rbind(nodular_dat, caseSub_stat))
        colnames(nodular_dat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        
      }
    }
  }
  if(isTRUE(ind == 4)){
    for(case in seq(1, 6, 1)){
      case_name <- switch(case , 'Case 1A', 'Case 1B', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
      for(label in seq(1, 5, 1)){
        label_name <- switch(label , 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
        Stat_N <- plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'N' & plotStat_bac$marker == label_name,]
        Stat_CT <- plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'CT' & plotStat_bac$marker == label_name,]
        Stat_IF <- plotStat_bac[plotStat_bac$case == case_name & plotStat_bac$location == 'IF' & plotStat_bac$marker == label_name,]
        
        irregular_N <- nrow(Stat_N[Stat_N$Eccentricity > 0.9 | Stat_N$Convexity < 0.3 | Stat_N$Circularity < 0.3,])
        irregular_CT <- nrow(Stat_CT[Stat_CT$Eccentricity > 0.9 | Stat_CT$Convexity < 0.3 | Stat_CT$Circularity < 0.3,])
        irregular_IF <- nrow(Stat_IF[Stat_IF$Eccentricity > 0.9 | Stat_IF$Convexity < 0.3 | Stat_IF$Circularity < 0.3,])
        
        subDat_N <- cbind(irregular_N, criteria_name, case_name, label_name, 'N')
        subDat_CT <- cbind(irregular_CT, criteria_name, case_name, label_name, 'CT')
        subDat_IF <- cbind(irregular_IF, criteria_name, case_name, label_name, 'IF')
        
        caseSub_stat <- data.frame(rbind(subDat_N, subDat_CT, subDat_IF))
        colnames(caseSub_stat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
        irregular_dat <- data.frame(rbind(irregular_dat, caseSub_stat))
        colnames(irregular_dat) <- c('value', 'Criteria', 'Case', 'Mark', 'location')
      }
    }
  }
}
mean_dat$value <- as.numeric(as.character(mean_dat$value))
max_dat$value <- as.numeric(as.character(max_dat$value))
nodular_dat$value <- as.numeric(as.character(nodular_dat$value))
irregular_dat$value <- as.numeric(as.character(irregular_dat$value))

Stat_comb <- rbind(mean_dat, max_dat, nodular_dat, irregular_dat)


# Density plot A
{
CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Density',]
CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Density',]
CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Density',]
CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Density',]
FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Density',]


c <- compare_means(value ~ location, CD3, ref.group = "N", method = 'wilcox.test', paired = FALSE)
c %<>% mutate(y_pos = c(7400, 8000), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))

p1 <- 
ggplot(CD3, aes(x = location, y = value)) + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
  ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                        manual = TRUE, textsize = 6)  +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
  scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  ylab(expression(paste("Density of cell, ", mm^{-2}))) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
plot(p1)

c <- compare_means(value ~ location, CD4, ref.group = "N", method = 'wilcox.test')
c %<>% mutate(y_pos = c(11200, 12200), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))

p2 <- 
  ggplot(CD4, aes(x = location, y = value)) + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
  ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                        manual = TRUE, textsize = 6)  +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
  scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  ylab(NULL) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
plot(p2)

c <- compare_means(value ~ location, CD8, ref.group = "N", method = 'wilcox.test')
c %<>% mutate(y_pos = c(5400, 6000), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
p3 <- 
  ggplot(CD8, aes(x = location, y = value)) + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
  ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                        manual = TRUE, textsize = 6)  +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
  scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  ylab(NULL) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

c <- compare_means(value ~ location, CD20, ref.group = "N", method = 'wilcox.test')
c %<>% mutate(y_pos = c(9000,9800), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
p4 <- 
  ggplot(CD20, aes(x = location, y = value)) + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
  ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                        manual = TRUE, textsize = 6)  +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
  scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  ylab(NULL) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
plot(p4)

c <- compare_means(value ~ location, FoxP3, ref.group = "N", method = 'wilcox.test')
c %<>% mutate(y_pos = c(7200,7800), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
p5 <- 
  ggplot(FoxP3, aes(x = location, y = value)) + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
  ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                        manual = TRUE, textsize = 6)  +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
  scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
  ylab(NULL) +
  xlab(NULL) +
  #ylab(expression(paste('Density, /', mm^2, )))+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
        strip.background = element_rect(colour="black", fill="white", size = 2),
        strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))
plot(p5)
Density_grid <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# max density plot B
{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Density of cluster',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Density of cluster',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Density of cluster',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Density of cluster',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Density of cluster',]
  
  
  c <- compare_means(value ~ location, CD3, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(5,5.4), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p1 <- 
    ggplot(CD3, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(expression(paste("Density of cluster, ", mm^{-2}))) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD4, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(6.6,7.1), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p2 <- 
    ggplot(CD4, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD8, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(5,5.4), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  p3 <- 
    ggplot(CD8, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD20, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(7.5,8.1), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  p4 <- 
    ggplot(CD20, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, FoxP3, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(4, 4.3), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))  

  p5 <- 
    ggplot(FoxP3, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  Cluseter_grid <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# max density plot B
{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'nodular count',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'nodular count',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'nodular count',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'nodular count',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'nodular count',]
  
  
  c <- compare_means(value ~ location, CD3, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(38,41), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p1 <- 
    ggplot(CD3, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab('Nodular cluster counts') +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD4, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(83, 89.5), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p2 <- 
    ggplot(CD4, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD8, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(56, 61), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  p3 <- 
    ggplot(CD8, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD20, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(30, 32.5), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  p4 <- 
    ggplot(CD20, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, FoxP3, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(15, 16.2), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p5 <- 
    ggplot(FoxP3, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab(NULL) +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  Nodular_grid <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

# max density plot B
{
  CD3 <- Stat_comb[Stat_comb$Mark == 'CD3' & Stat_comb$Criteria == 'Elongated count',]
  CD4 <- Stat_comb[Stat_comb$Mark == 'CD4' & Stat_comb$Criteria == 'Elongated count',]
  CD8 <- Stat_comb[Stat_comb$Mark == 'CD8' & Stat_comb$Criteria == 'Elongated count',]
  CD20 <- Stat_comb[Stat_comb$Mark == 'CD20' & Stat_comb$Criteria == 'Elongated count',]
  FoxP3 <- Stat_comb[Stat_comb$Mark == 'FoxP3' & Stat_comb$Criteria == 'Elongated count',]
  
  
  c <- compare_means(value ~ location, CD3, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(147,160), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p1 <- 
    ggplot(CD3, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab('Elongated cluster counts') +
    xlab('CD3') +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD4, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(310, 335), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p2 <- 
    ggplot(CD4, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('CD4') +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD8, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(200, 216), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  p3 <- 
    ggplot(CD8, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('CD8') +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, CD20, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(151, 163), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  p4 <- 
    ggplot(CD20, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('CD20') +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  c <- compare_means(value ~ location, FoxP3, ref.group = "N", method = 'wilcox.test')
  c %<>% mutate(y_pos = c(85,92), labels = ifelse(p < 0.05, sprintf("%2.2e",p),sprintf('%5.4f',p)))
  
  p5 <- 
    ggplot(FoxP3, aes(x = location, y = value)) + 
    theme_bw(base_rect_size = 1) +
    geom_boxplot(alpha = 0.3, size = 1.5, show.legend = FALSE, aes(color = location, fill = location)) +
    ggsignif::geom_signif(data = as.data.frame(c), aes(xmin=group1, xmax=group2, annotations=labels, y_position=y_pos),
                          manual = TRUE, textsize = 6)  +
    geom_jitter(width = 0.2, size = 2, show.legend = FALSE, aes(color = factor(location), fill = factor(location))) +
    scale_color_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    scale_fill_manual(values=c("#4cdee6", "#e47267", "#13ec87")) +
    ylab(NULL) +
    xlab('FoxP3') +
    #ylab(expression(paste('Density, /', mm^2, )))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 38), legend.title = element_text(size = 40), 
          strip.background = element_rect(colour="black", fill="white", size = 2),
          strip.text = element_text(margin = margin(10, 10, 10, 10), size = 40), panel.grid = element_line(size = 1.5))

  Elongated_grid <- plot_grid(p1, p2,p3, p4,p5, labels = c('', '','','',''), label_size = 30, ncol = 5)
}

jpeg('./Cluster_stat.jpeg', units="in", width=22, height=25, res=300)
plot_grid(Density_grid, Cluseter_grid, Nodular_grid, Elongated_grid, labels = c('', '','',''),  ncol = 1, align = 'hv') +
  xlab('Locations')
dev.off()






##############################################
### Count the number of clusters #############
##############################################
suffix <- '_shape.csv'

Current_mark <- switch(3, 'CD3', 'CD4', 'CD8', 'CD20', 'FoxP3')
for(case in seq(1,5,1)){
    Case <- switch(case, './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
    FilePath <- paste(Case, Current_mark, '/', Current_mark, suffix, sep = '')
    File <- read.csv(FilePath)
    print(Case)
    num_NM <- nrow(File[File$location == 3,])
    print(paste('NM:', num_NM))
    num_IM <- nrow(File[File$location == 1,])
    print(paste('IM:', num_IM))
    num_TM <- nrow(File[File$location == 2,])
    print(paste('TM:', num_TM))
}

##############################################
### Map clusters and spatial points ##########
##############################################
clusterPath <- switch(2, './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
Current_lab <- switch(3, 'CD3','CD4','CD8','CD20','FoxP3')
ptsPath <- paste(clusterPath, Current_lab, '/', Current_lab, '_Points.mat', sep = '')

clusterOutlinepath <- paste(clusterPath, Current_lab, '/', Current_lab, '_shape_outline.csv',sep='')
Pts <- as.data.frame(readMat(ptsPath))/125
Cluster_dats <- read.csv(clusterOutlinepath)

cluster_pts <- as.data.frame(clusters$clusters) 
# Characterize clusters based on their location
#Cluster_dist <- clusterCharacterize(Cluster_dats, TumorFile, InvasiveFile)

Pts <- cbind(Pts,cluster_pts)
Pts <- Pts[is.na(Pts[,3]) == FALSE, ]
names(Pts) <- c('x','y','cluster')
# Visualization for quality check

# Change to specific case here 
imgage <- readJPEG('./Case2_CD4.jpeg')
for_plot <-  ggplot() +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  31.872, -18.024, 0) +
  xlab('x,mm')+
  ylab('y,mm')
jpeg( './Case_2/Cluster_map.jpeg',units="in", width=16, height=9, res=1200)
move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_polygon(data = Cluster_dats, aes(x = Cluster_dats[,2], y = Cluster_dats[,3], group = Cluster_dats[,1], colour = as.factor(Cluster_dats[,1])), size = 1, fill = NA, show.legend = FALSE) +
  xlim(0,31.872) +
  ylim(18.024, 0) +
  geom_point(data = Pts, aes(x = Pts[,1], y = Pts[,2],colour = as.factor(Pts[,3])), size = 0.5, fill = NA,show.legend = FALSE) +
  
  #labs(title = 'Case2-CD4 clusters') +
  xlab('x, mm') +
  ylab('y, mm') +
  theme_bw()  +
  theme(axis.title = element_text(size = 35)) +
  theme(axis.text =  element_text(size = 30)) +
  theme(plot.title = element_text(size = 35)) 
dev.off()

# Shape descriptors ,Circularity
CD3_circ <- data.frame(matrix(nrow = 0, ncol = 3))
names(CD3_circ) <- c('circ', 'case', 'marker')
CD4_circ <- data.frame(matrix(nrow = 0, ncol = 3))
names(CD3_circ) <- c('circ', 'case','marker')

CD8_circ <- data.frame(matrix(nrow = 0, ncol = 3))
names(CD8_circ) <- c('circ', 'case','marker')

CD20_circ <- data.frame(matrix(nrow = 0, ncol = 3))
names(CD20_circ) <- c('circ', 'case','marker')

FoxP3_circ <- data.frame(matrix(nrow = 0, ncol = 3))
names(FoxP3_circ) <- c('circ', 'case','marker')

for(case in seq(1,5,1)){
  Case <- switch(case, './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  CD3_shape <- read.csv(paste(Case, 'CD3', '/', 'CD3_shape.csv',sep = ''))
  CD4_shape <- read.csv(paste(Case, 'CD4', '/', 'CD4_shape.csv',sep = ''))
  CD8_shape <- read.csv(paste(Case, 'CD8', '/', 'CD8_shape.csv',sep = ''))
  CD20_shape <- read.csv(paste(Case, 'CD20', '/', 'CD20_shape.csv',sep = ''))
  FoxP3_shape <- read.csv(paste(Case, 'FoxP3', '/', 'FoxP3_shape.csv',sep = ''))

  #########
  CD3_sub <- cbind(CD3_shape$circ, case, 'CD3')
  names(CD3_sub) <-  c('circ', 'case', 'marker')
  CD3_circ <- rbind(CD3_circ, CD3_sub)
  
  CD4_sub <- cbind(CD4_shape$circ, case, 'CD4')
  names(CD4_sub) <-  c('circ', 'case','marker')
  CD4_circ <- rbind(CD4_circ, CD4_sub)
  
  CD8_sub <- cbind(CD8_shape$circ, case, 'CD8')
  names(CD8_sub) <-  c('circ', 'case','marker')
  CD8_circ <- rbind(CD8_circ, CD8_sub)
  
  CD20_sub <- cbind(CD20_shape$circ, case, 'CD20')
  names(CD20_sub) <-  c('circ', 'case','marker')
  CD20_circ <- rbind(CD20_circ, CD20_sub)
  
  FoxP3_sub <- cbind(FoxP3_shape$circ, case, 'FoxP3')
  names(FoxP3_sub) <-  c('circ', 'case','marker')
  FoxP3_circ <- rbind(FoxP3_circ, FoxP3_sub)
  
}

Cluster_shape <- rbind(CD3_circ, CD4_circ,CD8_circ,CD20_circ,FoxP3_circ)
names(Cluster_shape) <- c('circ', 'case', 'marker')
Cluster_shape <- melt(Cluster_shape, value.name = "circ")
ggplot(Cluster_shape, aes(x = case, y = circ, group = case)) + 
  geom_boxplot(aes(fill = case)) +
  facet_wrap(~ marker) +
  scale_y_discrete(breaks = seq(0,1,0.2))

