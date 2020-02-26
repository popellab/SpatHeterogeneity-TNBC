######################################################################################
### This script is used to visualize subregion fitting results (Fig.S12 and Fig.7) ###
######################################################################################


library(ggplot2)
library(R.matlab)
library(jpeg)
library(grid)
library(tidyverse)
library(largeVis)
library(gginnards)
library(alphahull)
library(spatstat)
library(reshape2)
source('MiFunction.R')

# read image
image <- readJPEG('./Case_3/Cluster_local.jpeg')


for_plot <-  ggplot() +
  annotation_custom(rasterGrob(image, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  3.2, -3.2, -0) + 
  xlab('x,mm')+
  ylab('y,mm')



#Case2_CD3 <- data.frame(readMat('./Case_1/CD3_Points.mat'))/125
Case3_CD4 <- data.frame(readMat('./Case_3/CD4/CD4_Points.mat'))/125
subSet <- Case3_CD4[Case3_CD4[,1] >= 14.72 & Case3_CD4[,1] <= 17.92 & Case3_CD4[,2] <= 3.52 & Case3_CD4[,2] >= 0.32,]
subSet$Point.All.1 <- subSet$Point.All.1 - 14.72
subSet$Point.All.2 <- subSet$Point.All.2 - 0.32
# clustering
minPts = 30
K = 4
dat <- as.matrix(subSet[, 1:2])
cd8vis <- largeVis(t(dat), dim=2, K = K, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)
clusters <- hdbscan(cd8vis, K=K, minPts = minPts, verbose = FALSE, threads = 1)
cluster_pts <- as.data.frame(clusters$clusters) 
Pts <- cbind(subSet,cluster_pts)
Pts <- Pts[is.na(Pts[,3]) == FALSE, ]
names(Pts) <- c('x','y','cluster')


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
  X <- subSet[cidMatch,1:2]
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
  

}

hullStats=data.frame(hullStats)
colnames(hullStats)<-c("n", "perimConcave","areaConcave","areaConvex")

shapeStats = hullStats
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


# Data summary
polygonDF = data.frame(polygonInfo)
colnames(polygonDF)<-c("cluster_ID", "x","y")

# 
write.csv(shapeStats, file = outfile)
write.csv(polygonDF, file = outfileOutline, row.names=FALSE)
# calculate ellipse area
areaEllipse <- pi*shapeStats$major*shapeStats$minor
plotStat_sub <- cbind(shapeStats$density, alphaShape@area, areaEllipse, shapeStats$circ, shapeStats$conv, shapeStats$eccentricity,shapeStats$location, Current_lab, clusterPath_name)
colnames(plotStat_sub) <- c('Density', 'a-Shape area', 'Fitted ellipse area', 'Circularity', 'Convexity', 'Eccentricity')
plotStat <- rbind(plotStat, plotStat_sub)


############
# Fig.S12A #
############
jpeg('./Cluster_plot/Cluster_plot.jpeg', units="in", width=10, height=10, res=800)

cluster_map <- gplot(clusters, subSet[, 1:2]) +
  annotation_custom(rasterGrob(image, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  3.2, -3.2, -0) + 
  xlab('x,mm') +
  ylab('y,mm') +
  coord_fixed(ratio = 1) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  geom_hline(yintercept=0.5, color = "gray") +
  geom_hline(yintercept=1, color = "gray") +
  geom_hline(yintercept=1.5, color = "gray") +
  geom_hline(yintercept=2, color = "gray") +
  geom_hline(yintercept=2.5, color = "gray") +
  geom_hline(yintercept=3, color = "gray") +
  geom_vline(xintercept=0.5, color = "gray") +
  geom_vline(xintercept=1, color = "gray") +
  geom_vline(xintercept=1.5, color = "gray") +
  geom_vline(xintercept=2, color = "gray") +
  geom_vline(xintercept=2.5, color = "gray") +
  geom_vline(xintercept=3, color = "gray") +
  #geom_point(aes(c(23.1745, 27.12), c(6.72, 3.52)), size = 3) +
  scale_x_continuous(limits = c(0,3.2), expand = c(0, 0)) +
  scale_y_reverse(limits = c(3.2,0), expand = c(0, 0))
dev.off()

############
# Fig.S12B #
############

jpeg('./Cluster_plot/Cluster_withGlosh.jpeg', units="in", width=10, height=10, res=800)
move_layers(cluster_map,"GeomCustomAnn", position = 'bottom') +
  theme(legend.position = "none")
dev.off()

############
# Fig.S12C #
############

jpeg('./Cluster_plot/Cluster_withPoly.jpeg', units="in", width=10, height=10, res=800)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_path(data = polygonDF, aes(x, y, group = cluster_ID, colour = factor(cluster_ID)), size = 1.5) +
  geom_point(data = Pts, aes(x, y, colour = factor(cluster)), size = 0.5) +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  geom_hline(yintercept=0.5, color = "gray") +
  geom_hline(yintercept=1, color = "gray") +
  geom_hline(yintercept=1.5, color = "gray") +
  geom_hline(yintercept=2, color = "gray") +
  geom_hline(yintercept=2.5, color = "gray") +
  geom_hline(yintercept=3, color = "gray") +
  geom_vline(xintercept=0.5, color = "gray") +
  geom_vline(xintercept=1, color = "gray") +
  geom_vline(xintercept=1.5, color = "gray") +
  geom_vline(xintercept=2, color = "gray") +
  geom_vline(xintercept=2.5, color = "gray") +
  geom_vline(xintercept=3, color = "gray") +
  scale_x_continuous(limits = c(0,3.2), expand = c(0, 0)) +
  scale_y_reverse(limits = c(3.2,0), expand = c(0, 0))
  
dev.off()

############
# Fig.S12D #
############

jpeg('./Cluster_plot/Cluster_withEllipse.jpeg', units="in", width=10, height=10, res=800)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme(legend.position = "none") +
  coord_fixed(ratio = 1) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  geom_hline(yintercept=0.5, color = "gray") +
  geom_hline(yintercept=1, color = "gray") +
  geom_hline(yintercept=1.5, color = "gray") +
  geom_hline(yintercept=2, color = "gray") +
  geom_hline(yintercept=2.5, color = "gray") +
  geom_hline(yintercept=3, color = "gray") +
  geom_vline(xintercept=0.5, color = "gray") +
  geom_vline(xintercept=1, color = "gray") +
  geom_vline(xintercept=1.5, color = "gray") +
  geom_vline(xintercept=2, color = "gray") +
  geom_vline(xintercept=2.5, color = "gray") +
  geom_vline(xintercept=3, color = "gray") +
  geom_point(data = Pts, aes(x, y, colour = factor(cluster)), size = 0.5) +
  stat_ellipse(data  = Pts, aes(x, y, group = cluster, colour = factor(cluster)),size = 1.5) +
  scale_x_continuous(limits = c(0,3.2), expand = c(0, 0)) +
  scale_y_reverse(limits = c(3.2,0), expand = c(0, 0))
dev.off()

# Prepare data for individual cluster
minCluster <- min(as.numeric(as.character(Pts[,3])))
maxCluster <- max(as.numeric(as.character(Pts[,3])))
chull_all <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(chull_all) <- c('x', 'y', 'cluster')

for(cid in seq(minCluster, maxCluster)){
  cur_Cluster <- Pts[Pts$cluster == cid,]
  convex <- chull(cur_Cluster[,1], cur_Cluster[,2])
  chull_pts <- data.frame(cur_Cluster[convex,])
  chull_pts$cluster <- cid
  colnames(chull_pts) <- c('x', 'y', 'cluster')
  
  chull_all <- rbind(chull_all, chull_pts)
  
}



################################
# read one cluster #############
################################
# manipulate coords, convex hull
cluster_sub <- chull_all[chull_all$cluster == 20,]
cluster_sub$x <- cluster_sub$x + 14.72
cluster_sub$y <- cluster_sub$y + 0.32

cluster_sub$x <- cluster_sub$x - 14.96
cluster_sub$y <- cluster_sub$y - 2.16
# manipulate coords, concave hull
polygonDF_sub <- polygonDF[polygonDF$cluster_ID == 20,]
polygonDF_sub$x <- polygonDF_sub$x + 14.72
polygonDF_sub$y <- polygonDF_sub$y + 0.32

polygonDF_sub$x <- polygonDF_sub$x - 14.96
polygonDF_sub$y <- polygonDF_sub$y - 2.16

image <- readJPEG('./Case_3/cluster1.jpeg')

for_plot <-  ggplot() +
  annotation_custom(rasterGrob(image, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  0.55, -0.875, -0) + 
  xlab('x,mm')+
  ylab('y,mm')

############
# Fig.S12E #
############

jpeg('./Cluster_plot/Cluster_individual.jpeg', units="in", width=5.5, height=8.75, res=800)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme(legend.position = "none") +
  #coord_fixed(ratio = 0.875/0.55) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  
  geom_path(data = cluster_sub, aes(x, y, colour = factor(cluster)), size = 1.5, linetype = 'dashed') +
  geom_polygon(data  = polygonDF_sub, aes(x, y, colour = factor(cluster_ID), fill = factor(cluster_ID)),size = 1.5, alpha = 0.1) +
  geom_point(data  = cluster_sub, aes(x, y, colour = factor(cluster)),size = 5) +
  
  scale_x_continuous(limits = c(0,0.55), expand = c(0, 0)) +
  scale_y_reverse(limits = c(0.875,0), expand = c(0, 0))
dev.off()

# draw individual ellipse
Pts_sub <- Pts[Pts$cluster == 20,]
Pts_sub$x <- Pts_sub$x + 14.72
Pts_sub$y <- Pts_sub$y + 0.32

Pts_sub$x <- Pts_sub$x - 14.96
Pts_sub$y <- Pts_sub$y - 2.16

############
# Fig.S12G #
############

jpeg('./Cluster_plot/Cluster_indEllipse.jpeg', units="in", width=5.5, height=8.75, res=800)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme(legend.position = "none") +
  #coord_fixed(ratio = 0.875/0.55) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  stat_ellipse(data  = Pts_sub, aes(x, y, group = cluster, colour = factor(cluster)),size = 2) +
  geom_point(data  = Pts_sub, aes(x, y, colour = factor(cluster)),size = 5, alpha = 0.8) +
  
  scale_x_continuous(limits = c(0,0.55), expand = c(0, 0)) +
  scale_y_reverse(limits = c(0.875,0), expand = c(0, 0))
dev.off()

##################################
# Second individual cluster ######
##################################

image <- readJPEG('./Case_3/cluster2.jpeg')

for_plot <-  ggplot() +
  annotation_custom(rasterGrob(image, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  1.1, -0.875, -0) + 
  xlab('x,mm')+
  ylab('y,mm')


# manipulate coords, convex hull
cluster_sub <- chull_all[chull_all$cluster == 9,]
cluster_sub$x <- cluster_sub$x + 14.72
cluster_sub$y <- cluster_sub$y + 0.32

cluster_sub$x <- cluster_sub$x - 16.79
cluster_sub$y <- cluster_sub$y - 0.95
# manipulate coords, concave hull
polygonDF_sub <- polygonDF[polygonDF$cluster_ID == 9,]
polygonDF_sub$x <- polygonDF_sub$x + 14.72
polygonDF_sub$y <- polygonDF_sub$y + 0.32

polygonDF_sub$x <- polygonDF_sub$x - 16.79
polygonDF_sub$y <- polygonDF_sub$y - 0.95

############
# Fig.S12F #
############

jpeg('./Cluster_plot/Cluster_individual2.jpeg', units="in", width=11, height=8.75, res=800)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme(legend.position = "none") +
 #coord_fixed(ratio = 0.875/1.1) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  
  geom_path(data = cluster_sub, aes(x, y), colour = '#1E90FF', size = 1.5, linetype = 'dashed') +
  geom_polygon(data  = polygonDF_sub, aes(x, y),  colour = '#1E90FF', fill = '#1E90FF', size = 1.5, alpha = 0.1) +
  geom_point(data  = cluster_sub, aes(x, y ),colour = '#1E90FF', size = 5) +
  
  scale_x_continuous(limits = c(0,1.1), expand = c(0, 0)) +
  scale_y_reverse(limits = c(0.875,0), expand = c(0, 0))
dev.off()

# draw individual ellipse
Pts_sub <- Pts[Pts$cluster == 9,]
Pts_sub$x <- Pts_sub$x + 14.72
Pts_sub$y <- Pts_sub$y + 0.32

Pts_sub$x <- Pts_sub$x - 16.79
Pts_sub$y <- Pts_sub$y - 0.95

############
# Fig.S12H #
############

jpeg('./Cluster_plot/Cluster_indEllipse2.jpeg', units="in", width=11, height=8.75, res=800)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme(legend.position = "none") +
  #coord_fixed(ratio = 0.875/1.1) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  stat_ellipse(data  = Pts_sub, aes(x, y),colour = '#1E90FF', size = 2) +
  geom_point(data  = Pts_sub, aes(x, y, ),colour = '#1E90FF', size = 5, alpha = 0.8) +
  
  scale_x_continuous(limits = c(0,1.1), expand = c(0, 0)) +
  scale_y_reverse(limits = c(0.875,0), expand = c(0, 0))
dev.off()


####################################
# Figure region, fit, test (Fig.7) #
####################################

#Case2_CD3 <- data.frame(readMat('./Case_1/CD3_Points.mat'))/125

Case1_FoxP3 <- data.frame(readMat('./Case_3/CD4/CD4_Points.mat'))/125
subSet <- Case1_FoxP3[Case1_FoxP3[,1] >= 4.05 & Case1_FoxP3[,1] <= 4.2 & Case1_FoxP3[,2] <= 9.75 & Case1_FoxP3[,2] >= 9.6,]
subSet$Point.All.1 <- subSet$Point.All.1 - 4.05
subSet$Point.All.2 <- subSet$Point.All.2 - 9.6


################
# Fig. 7A or E #
################

image <- readJPEG('./CD3_cluster.jpeg')

for_plot <-  ggplot() +
  annotation_custom(rasterGrob(image, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE),
                    0,  0.15, -0.15, 0) + 
  xlab('x,mm')+
  ylab('y,mm')

################
# Fig. 7B or F #
################

jpeg('./Point_pattern_cluster.jpeg', units="in", width=12, height=12, res=800)

ggplot() +
  theme_bw(base_rect_size = 3)+
  xlab('x, mm') +
  ylab('y, mm') +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 42),
        plot.background = element_rect(size = 40), panel.grid = element_line(size = 2)) +
  ylim(0.15,0) +
  xlim(0,0.15) +
  geom_point(data  = subSet, aes(subSet[,1], subSet[,2]),size = 5)
dev.off()
  
move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme_bw()+
  theme(legend.position = "none") +
  #coord_fixed(ratio = 0.875/0.55) +
  geom_point(data  = subSet, aes(subSet[,1], subSet[,2]),size = 2)+
  xlim(0,0.15) +
  ylim(0.15,0) 
dev.off()


subSet_p <- ppp(subSet[,1], subSet[,2], window = owin(c(0, 0.15), c(0, 0.15)))

K_func <- Kest(subSet_p, correction = 'Ripley', rmax = 0.1)
colnames(K_func) <- c('r', 'Poisson process', 'K(r), Ripley/isotropic correction')
K_func <- melt(K_func, id.vars = 'r')

################
# Fig. 7C or G #
################
jpeg('./Ripley.jpeg', units="in", width=15, height=15, res=800)

ggplot(data = K_func, aes(x = r, col = variable, linetype = variable)) +
  geom_line(aes(r, value), size = 1.5) +
  theme_bw(base_rect_size = 3)+
  xlab('x, mm') +
  ylab('K(r), mm') +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 42),
        plot.background = element_rect(size = 40), panel.grid = element_line(size = 2),
        legend.text = element_text(size = 30),legend.title = element_blank(), legend.position = c(0.4, 0.85), 
        legend.background=element_blank(), legend.key.size =  unit(0.8, "in")) +
  ylim(0, 0.04) +
  xlim(0,0.1) +
 # geom_line(aes(x = r, y = value, linetype = variable), size = 1.5) +
  scale_color_manual(values = c("black",'#BB1200')) + 
  scale_linetype_manual(values = c(2, 1)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))
  
dev.off()

################################
### L(r) function ##############
################################


subSet_p <- ppp(subSet[,1], subSet[,2], window = owin(c(0, 0.15), c(0, 0.15)))
L_func <- Lest(subSet_p, correction = 'isotropic')
L_func$iso <- L_func$iso - L_func$r
L_func$theo <- L_func$theo - L_func$r

#Thomas <- kppm(subSet_p, clusters = "Thomas")
#Thomas_fit <- data.frame(Thomas$Fit$mcfit$fit)
#Thomas$Fit$mcfit$fit$fit <-sqrt(Thomas_fit$fit/pi)

#Thomas_fit$fit <- sqrt(Thomas_fit$fit/pi) - Thomas_fit$r

#L_func <- cbind(L_func, Thomas_fit$fit)

colnames(L_func) <- c('r', 'CSR, confidence', 'Observations') # don't forget add 'Model fit'

L_func <- melt(L_func, id.vars = 'r')

L_envelope <- envelope(subSet_p, fun = Lest, correction = 'isotropic')

colnames(L_envelope) <- c('r', 'Observations','CSR, confidence', 'low', 'high')
L_envelope$low <- L_envelope$low - L_envelope$r
L_envelope$high <- L_envelope$high - L_envelope$r


################
# Fig. 7D or H #
################

jpeg('./thomas.jpeg', units="in", width=15, height=15, res=800)

ggplot() +
  geom_ribbon(data = L_envelope, aes(x = r, ymin = low, ymax = high),fill = "grey70") +
  #geom_ribbon(data = Thomas_dg, aes(x = r, ymin = low, ymax = high),fill = "#4FCA53", alpha = 0.2) +
  
  geom_line(data = L_func, aes(x = r, y =value, col = variable, linetype = variable), size = 1.5) +
  theme_bw(base_rect_size = 3)+
  xlab('x, mm') +
  ylab('L(r) - r, mm') +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 42),
        plot.background = element_rect(size = 40), panel.grid = element_line(size = 2),
        legend.text = element_text(size = 30),legend.title = element_blank(), legend.position = c(0.3, 0.8), 
        legend.background=element_blank(), legend.key.size =  unit(0.8, "in")) +
  ylim(-0.005, 0.008) +
  xlim(0,0.038) +
  # geom_line(aes(x = r, y = value, linetype = variable), size = 1.5) +
  scale_color_manual(values = c("black",'#BB1200','#4FCA53')) + 
  scale_linetype_manual(values = c(2, 1, 1)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))
dev.off()

