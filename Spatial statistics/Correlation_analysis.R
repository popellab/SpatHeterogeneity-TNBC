#################################################################
##### This script performs correlation analysis (Fig. 9,10) #####
#################################################################

library(cluster)
library(Ecdat)
library(compareGroups)
library(largeVis)
library(R.matlab)
library(ggplot2)
#library(dbscan)
library(cowplot)
source('MiFunction.R')

# Workflow description
# Find find DoC score for each point, use DoC score to filter those points that are less localized

######### CD3 + CD8 pair

DoC_Score_CD8 <- read.csv('./Case_2/DoC_Score_CD8_Case2.csv', header = F)
DoC_Score_CD3 <- read.csv('./Case_2/DoC_Score_Case2_CD3.csv', header = F,)

CD8_Pts <- data.frame(readMat('./Case_2/CD8/CD8_Points.mat'))/125
CD3_Pts <- data.frame(readMat('./Case_2/CD3/CD3_Points.mat'))/125
print('#####')
print(nrow(data.frame(readMat('./Case_5/CD3/CD3_Points.mat'))/125))
print(nrow(data.frame(readMat('./Case_5/CD4/CD4_Points.mat'))/125))
print(nrow(data.frame(readMat('./Case_5/CD8/CD8_Points.mat'))/125))
print(nrow(data.frame(readMat('./Case_5/CD20/CD20_Points.mat'))/125))
print(nrow(data.frame(readMat('./Case_5/FoxP3/FoxP3_Points.mat'))/125))

CD8_Pts_move <- CD8_Pts

CD8_Pts_move$CD8.Registered.1 <- CD8_Pts_move$CD8.Registered.1 + 0.01/0.976 
write.table(CD8_Pts_move, file = './Case_2/CD8/CD8_Points_move.csv', row.names = FALSE, col.names = FALSE, sep = ',')
{# find threshold
  Thresh_CD8 <- read.csv('./Case_2/CD8/DoC_Score_Thresh_CD8.csv', header = F)
  Thresh_CD8m <- read.csv('./Case_2/CD8/DoC_Score_Thresh_CD8m.csv', header = F)
  
  
  Thresh_List <-trueDoC(CD8_Pts, CD8_Pts_move, Thresh_CD8, Thresh_CD8m, 2)
  Thresh_CD8_Score <- Thresh_List$listA
  Thresh_CD8m_Score <- Thresh_List$listB
  Thresh_CD8_order <- Thresh_CD8_Score[order(Thresh_CD8_Score$DoC),]
  
  Thresh_CD8m_Score <- Thresh_CD8m_Score[order(Thresh_CD8m_Score$DoC),]
  
  print(Thresh_CD8_order[round(0.9*nrow(Thresh_CD8_order)),3])
  print(Thresh_CD8m_Score[round(0.9*nrow(Thresh_CD8m_Score)),3])
  
  #print(Thresh_CD8_order[round(0.1*nrow(Thresh_CD8_order)),3])
  #print(Thresh_CD8m_Score[round(0.1*nrow(Thresh_CD8m_Score)),3])
}#

List <-trueDoC(CD3_Pts, CD8_Pts, DoC_Score_CD3, DoC_Score_CD8, 2)
CD3_Pts_Score <- List$listA
CD8_Pts_Score <- List$listB

CD3_Pts_Score <- CD3_Pts_Score[complete.cases(CD3_Pts_Score),]
CD8_Pts_Score <- CD8_Pts_Score[complete.cases(CD8_Pts_Score),]

Case1_DoC_CD3 <- CD3_Pts_Score
Case1_DoC_CD8 <- CD8_Pts_Score


QCoD_CD3 <- (quantile(CD3_Pts_Score$DoC, 0.75, na.rm = TRUE) - quantile(CD3_Pts_Score$DoC, 0.25, na.rm = TRUE))/ (quantile(CD3_Pts_Score$DoC, 0.75, na.rm = TRUE)+ quantile(CD3_Pts_Score$DoC, 0.25, na.rm = TRUE))
QCoD_CD8 <- (quantile(CD8_Pts_Score$DoC, 0.75) - quantile(CD8_Pts_Score$DoC, 0.25))/ (quantile(CD8_Pts_Score$DoC, 0.75) + quantile(CD8_Pts_Score$DoC, 0.25))



CD3_Pts_Score <- CD3_Pts_Score[CD3_Pts_Score$DoC > 0.84,]
CD8_Pts_Score <- CD8_Pts_Score[CD8_Pts_Score$DoC > 0.84,]

Case1_DoC_CD3_high <- CD3_Pts_Score
Case1_DoC_CD8_high <- CD8_Pts_Score


CD3_Pts_Score$label <- 'CD3'
CD8_Pts_Score$label <- 'CD8'

mexDat <- rbind(CD3_Pts_Score,CD8_Pts_Score)
mexDat <- mexDat[complete.cases(mexDat),]

print(nrow(mexDat[mexDat$label == 'CD3',]))
print(nrow(mexDat[mexDat$label == 'CD8',]))
# hdbscan
cd8vis <- largeVis(t(as.matrix(mexDat[, 1:2])),K = 4, dim=2, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)

clusters <- largeVis::hdbscan(cd8vis,K = 4, minPts = 30, verbose = FALSE, threads = 1)
mexDat$cluster <- clusters$clusters
mexDat <- mexDat[complete.cases(mexDat),]


# count number 
print(nrow(mexDat[mexDat$label == 'CD3',])*100/nrow(CD3_Pts))
print(nrow(mexDat[mexDat$label == 'CD8',])*100/nrow(CD8_Pts))


# get polygon infomation for clusters

maxCluster = max(as.numeric(mexDat$cluster), na.rm = TRUE)
minCluster = min(as.numeric(mexDat$cluster), na.rm = TRUE)
notNA = !is.na(mexDat$cluster)
polygonInfo <- matrix(nrow = 0, ncol = 3) #clusterID, x, y
alpha0 <- 0.01
alphaStep <- 0.01
# check number of clusters in each region

N_region <- readRDS('./Case_0/Normal_1a.rds')
IF_region <- readRDS('./Case_0/Invasive_1a.rds')
CT_region <- readRDS('./Case_0/Tumor_1a.rds')

X_IF <- 0
X_N <- 0
X_CT <- 0
for(cid_target in 1:maxCluster) {
  print(cid_target)
  goodAlpha = FALSE
  alpha <- alpha0
  #X<-subPoints[cidMatch,]
  if(isTRUE(nrow(mexDat[mexDat$cluster == cid_target & mexDat$label == 'CD3',1:2])*nrow(mexDat[mexDat$cluster == cid_target & mexDat$label == 'CD8',1:2]) != 0)){
    X <- mexDat[mexDat$cluster == cid_target,1:2]
    X <- unique(X)
    # try different alpha until shape covers all points
    while(!goodAlpha){
      #print(alpha)
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
    N_num <- data.frame(point.in.polygon(X[,1], X[,2], N_region[,1], N_region[,2]))
    IF_num <-data.frame(point.in.polygon(X[,1], X[,2], IF_region[,1], IF_region[,2]))
    CT_num <-data.frame(point.in.polygon(X[,1], X[,2], CT_region[,1], CT_region[,2]))
    
    N_num <- length(which(N_num[,1] != 0))
    IF_num <- length(which(IF_num[,1] != 0))
    CT_num <- length(which(CT_num[,1] != 0))
    
    
    if(isTRUE(IF_num != 0)){
      X_IF <- X_IF + 1
    }
    if(isTRUE(IF_num == 0 & CT_num != 0)){
      X_CT <- X_CT + 1
    }
    if(isTRUE(IF_num == 0 & CT_num == 0 & N_num != 0)){
      X_N <- X_N + 1
    }
    
  }
}

# check the number associated with each region


colnames(polygonInfo) <- c('cluster_id', 'x', 'y')
polygonInfo <- data.frame(polygonInfo)

###### for plot 1A and 1B together
Case0_polygon <- polygonInfo
Case0_mex <- mexDat

# 1B
Case1_polygon <- polygonInfo
Case1_mex <- mexDat

###########
# Fig. 9A #
###########

jpeg('./DoC_Case1AB_CD8.jpeg', units="in", width=14, height= 10.4, res=300)
ggplot() + geom_point(data = Case0_DoC_CD8, aes(x, y, color = DoC), alpha = 0.2, show.legend = FALSE,size = 0.5) +
  geom_point(data = Case1_DoC_CD8, aes(x, y, color = DoC), alpha = 0.2,size = 0.5) +
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 36)) +
  scale_colour_gradientn(colours = rainbow(5),limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1),
                        guide = guide_colourbar(barwidth=1, barheight=20), name = 'DoC') +
  xlim(0, 28.072) +
  ylim(20.824, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()
dev.off()

##################
# Fig. 9C (left) #
##################

jpeg('./DoC_Case1AB_mex.jpeg', units="in", width=14, height= 10.4, res=300)
ggplot() + geom_point(data = Case0_mex, aes(x, y), color = '#4FA7CA', alpha = 0.2, show.legend = FALSE,size = 0.5) +
  geom_point(data = Case1_mex, aes(x, y), color = '#4FA7CA', show.legend = FALSE, alpha = 0.2,size = 0.5) +
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 36)) +
  #scale_colour_gradientn(colours = rainbow(5),limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1),
   #                      guide = guide_colourbar(barwidth=1, barheight=20), name = 'DoC Score') +
  xlim(0, 28.072) +
  ylim(20.824, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()
dev.off()

###################
# Fig. 9C (right) #
###################
jpeg('./DoC_Cluster_Case1AB.jpeg', units="in", width=15.9, height=9, res=300)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme_bw(base_rect_size = 1) +
  geom_polygon(data = Poly1_path, aes(x = Poly1_path[,1], y = Poly1_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly2_path, aes(x = Poly2_path[,1], y = Poly2_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly3_path, aes(x = Poly3_path[,1], y = Poly3_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3, show.legend = FALSE) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  geom_polygon(data = Case0_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  geom_point(data = Case0_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  geom_point(data = Case0_mex, aes(x, y, group = cluster, color = factor(cluster)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  
  geom_polygon(data = Case1_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  geom_point(data = Case1_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  geom_point(data = Case1_mex, aes(x, y, group = cluster, color = factor(cluster)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  
  theme(text = element_text(size = 36)) +
  xlim(0, 28.072) +
  ylim(20.824, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()
dev.off()


ggplot() +   geom_polygon(data = polygonInfo, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 24)) +
  xlim(0, 31.872) +
  ylim(18.024, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()

# Histogram for CD3 CD8

p1 <- ggplot(Case1_DoC_CD3, aes(x=DoC)) + 
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 26)) +
  geom_histogram(color="black", fill="black", binwidth = 0.01) +
  xlab(NULL) +
  ylab('Frequency,CD3') +
  geom_vline(xintercept = 0.84, color = 'red', linetype = 'dashed', size =1)

p2 <- ggplot(Case1_DoC_CD8, aes(x=DoC)) + 
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 26)) +
  geom_histogram(color="black", fill="black", binwidth = 0.01) +
  xlab('DoC Score') +
  ylab('Frequency,CD8') +
  geom_vline(xintercept = 0.84, color = 'red', linetype = 'dashed', size =1)

###################
# Fig. 9B (left panel for 1A, right panel for 1B) ##
###################

jpeg('./Case1AB_CD8_histo.jpg', units="in", width=9, height=6, res=300)
   plot_grid(p1, p2, labels = c('', ''), label_size = 30, ncol = 1)
dev.off()


# Workflow description
# Find find DoC score for each point, use DoC score to filter those points that are less localized

######### CD4 + FoxP3

DoC_Score_CD4 <- read.csv('./Case_5/DoC_Score_Case5_CD4.csv', header = F)
DoC_Score_FoxP3 <- read.csv('./Case_5/DoC_Score_Case5_FoxP3.csv', header = F,)

CD4_Pts <- data.frame(readMat('./Case_5/CD4/CD4_Points.mat'))/125
FoxP3_Pts <- data.frame(readMat('./Case_5/FoxP3/FoxP3_Points.mat'))/125



print(nrow(CD4_Pts))
print(nrow(FoxP3_Pts))

List <-trueDoC(CD4_Pts, FoxP3_Pts, DoC_Score_CD4, DoC_Score_FoxP3, 2)
CD4_Pts_Score <- List$listA
FoxP3_Pts_Score <- List$listB

CD4_Pts_Score <- CD4_Pts_Score[complete.cases(CD4_Pts_Score),]
FoxP3_Pts_Score <- FoxP3_Pts_Score[complete.cases(FoxP3_Pts_Score),]



QCoD_CD4 <- (quantile(CD4_Pts_Score$DoC, 0.75, na.rm = TRUE) - quantile(CD4_Pts_Score$DoC, 0.25, na.rm = TRUE))/ (quantile(CD4_Pts_Score$DoC, 0.75, na.rm = TRUE)+ quantile(CD4_Pts_Score$DoC, 0.25, na.rm = TRUE))
QCoD_FoxP3 <- (quantile(FoxP3_Pts_Score$DoC, 0.75) - quantile(FoxP3_Pts_Score$DoC, 0.25))/ (quantile(FoxP3_Pts_Score$DoC, 0.75) + quantile(FoxP3_Pts_Score$DoC, 0.25))
print(QCoD_CD4)
print(QCoD_FoxP3)

CD4_Pts_Score <- CD4_Pts_Score[CD4_Pts_Score$DoC > 0.997,]
FoxP3_Pts_Score <- FoxP3_Pts_Score[FoxP3_Pts_Score$DoC > 0.997,]

CD4_Pts_Score$label <- 'CD4'
FoxP3_Pts_Score$label <- 'FoxP3'

mexDat <- rbind(CD4_Pts_Score,FoxP3_Pts_Score)
mexDat <- mexDat[complete.cases(mexDat),]

print(nrow(mexDat[mexDat$label == 'CD4',]))
print(nrow(mexDat[mexDat$label == 'FoxP3',]))

# hdbscan
cd8vis <- largeVis(t(as.matrix(mexDat[, 1:2])),K = 4, dim=2, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)

clusters <- largeVis::hdbscan(cd8vis,K = 4, minPts = 30, verbose = FALSE, threads = 1)
mexDat$cluster <- clusters$clusters
mexDat <- mexDat[complete.cases(mexDat),]


# count number 
print(nrow(mexDat[mexDat$label == 'CD4',])*100/nrow(CD4_Pts))
print(nrow(mexDat[mexDat$label == 'FoxP3',])*100/nrow(FoxP3_Pts))


# get polygon infomation for clusters

maxCluster = max(as.numeric(mexDat$cluster), na.rm = TRUE)
minCluster = min(as.numeric(mexDat$cluster), na.rm = TRUE)
notNA = !is.na(mexDat$cluster)
polygonInfo <- matrix(nrow = 0, ncol = 3) #clusterID, x, y
alpha0 <- 0.01
alphaStep <- 0.01
# check number of clusters in each region

N_region <- readRDS('./Case_5/Normal.rds')
IF_region <- readRDS('./Case_5/Invasive.rds')
CT_region <- readRDS('./Case_5/Tumor.rds')

X_IF <- 0
X_N <- 0
X_CT <- 0
for(cid_target in 1:maxCluster) {
  print(cid_target)
  goodAlpha = FALSE
  alpha <- alpha0
  #X<-subPoints[cidMatch,]
  if(isTRUE(nrow(mexDat[mexDat$cluster == cid_target & mexDat$label == 'CD4',1:2])*nrow(mexDat[mexDat$cluster == cid_target & mexDat$label == 'FoxP3',1:2]) != 0)){
    X <- mexDat[mexDat$cluster == cid_target,1:2]
    X <- unique(X)
    # try different alpha until shape covers all points
    while(!goodAlpha){
      #print(alpha)
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
    N_num <- data.frame(point.in.polygon(X[,1], X[,2], N_region[,1], N_region[,2]))
    IF_num <-data.frame(point.in.polygon(X[,1], X[,2], IF_region[,1], IF_region[,2]))
    CT_num <-data.frame(point.in.polygon(X[,1], X[,2], CT_region[,1], CT_region[,2]))
    
    N_num <- length(which(N_num[,1] != 0))
    IF_num <- length(which(IF_num[,1] != 0))
    CT_num <- length(which(CT_num[,1] != 0))
    
    
    if(isTRUE(IF_num != 0)){
      X_IF <- X_IF + 1
    }
    if(isTRUE(IF_num == 0 & CT_num != 0)){
      X_CT <- X_CT + 1
    }
    if(isTRUE(IF_num == 0 & CT_num == 0 & N_num != 0)){
      X_N <- X_N + 1
    }
    
  }
}

# check the number associated with each region


colnames(polygonInfo) <- c('cluster_id', 'x', 'y')
polygonInfo <- data.frame(polygonInfo)




jpeg('./DoC_Case2_FoxP3_High.jpeg', units="in", width=15.9, height=9, res=300)
ggplot() + geom_point(data = FoxP3_Pts_Score, aes(x, y, color = factor(label)), alpha = 0.2, show.legend = FALSE,size = 0.5) +
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 36)) +
  #scale_colour_gradientn(colours = rainbow(5),limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1),
  #                      guide = guide_colourbar(barwidth=1, barheight=20), name = 'DoC Score') +
  xlim(0, 31.872) +
  ylim(18.024, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()
dev.off()
# for mapping
jpeg('./DoC_Cluster.jpeg', units="in", width=15.9, height=9, res=300)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  geom_polygon(data = Poly1_path, aes(x = Poly1_path[,1], y = Poly1_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly2_path, aes(x = Poly2_path[,1], y = Poly2_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly3_path, aes(x = Poly3_path[,1], y = Poly3_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3, show.legend = FALSE) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  geom_polygon(data = polygonInfo, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  geom_point(data = polygonInfo, aes(x, y, group = cluster_id, color = factor(cluster_id)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  geom_point(data = mexDat, aes(x, y, group = cluster, color = factor(cluster)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  
  theme(text = element_text(size = 36)) +
  xlim(0, 31.872) +
  ylim(18.024, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()
dev.off()


ggplot() +   geom_polygon(data = polygonInfo, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 24)) +
  xlim(0, 31.872) +
  ylim(18.024, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()

###################
# Fig. 9B (right) ##
###################

jpeg('./Histogram_CD8.jpeg', units="in", width=6, height=2, res=300)

ggplot(CD8_Pts_Score, aes(x=DoC)) + 
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 12)) +
  geom_histogram(color="black", fill="black", binwidth = 0.01) +
  xlab(NULL) +
  ylab('Frequency,CD8') +
  geom_vline(xintercept = 0.84, color = 'red', linetype = 'dashed', size =1)
dev.off()
jpeg('./Histogram_FoxP3.jpeg', units="in", width=6, height=2, res=300)
ggplot(FoxP3_Pts_Score, aes(x=DoC)) + 
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 12)) +
  geom_histogram(color="black", fill="black", binwidth = 0.01) +
  xlab('DoC Score') +
  ylab('Frequency,FoxP3') +
  geom_vline(xintercept = 0.84, color = 'red', linetype = 'dashed', size =1)
dev.off()



########################################
######### CD8 + FoxP3 #################
########################################

DoC_Score_CD8 <- read.csv('./Case_5/DoC_Score_Case5_CD8.csv', header = F)
DoC_Score_FoxP3 <- read.csv('./Case_5/DoC_Score_Case5_FoxP3.csv', header = F,)

CD8_Pts <- data.frame(readMat('./Case_5/CD8/CD8_Points.mat'))/125
FoxP3_Pts <- data.frame(readMat('./Case_5/FoxP3/FoxP3_Points.mat'))/125

print(nrow(CD8_Pts))
print(nrow(FoxP3_Pts))

List <-trueDoC(CD8_Pts, FoxP3_Pts, DoC_Score_CD8, DoC_Score_FoxP3, 2)
CD8_Pts_Score <- List$listA
FoxP3_Pts_Score <- List$listB

CD8_Pts_Score <- CD8_Pts_Score[complete.cases(CD8_Pts_Score),]
FoxP3_Pts_Score <- FoxP3_Pts_Score[complete.cases(FoxP3_Pts_Score),]

Case1_DoC_CD8 <- CD8_Pts_Score
Case1_DoC_FoxP3 <- FoxP3_Pts_Score

QCoD_CD8 <- (quantile(CD8_Pts_Score$DoC, 0.75, na.rm = TRUE) - quantile(CD8_Pts_Score$DoC, 0.25, na.rm = TRUE))/ (quantile(CD8_Pts_Score$DoC, 0.75, na.rm = TRUE)+ quantile(CD8_Pts_Score$DoC, 0.25, na.rm = TRUE))
QCoD_FoxP3 <- (quantile(FoxP3_Pts_Score$DoC, 0.75) - quantile(FoxP3_Pts_Score$DoC, 0.25))/ (quantile(FoxP3_Pts_Score$DoC, 0.75) + quantile(FoxP3_Pts_Score$DoC, 0.25))
print(QCoD_CD8)
print(QCoD_FoxP3)

CD8_Pts_Score <- CD8_Pts_Score[CD8_Pts_Score$DoC < 0.84,]
FoxP3_Pts_Score <- FoxP3_Pts_Score[FoxP3_Pts_Score$DoC > 0.84,]

Case1_DoC_CD8_high <- CD8_Pts_Score
Case1_DoC_FoxP3_high <- FoxP3_Pts_Score

CD8_Pts_Score$label <- 'CD8'
FoxP3_Pts_Score$label <- 'FoxP3'

mexDat <- rbind(CD8_Pts_Score,FoxP3_Pts_Score)
mexDat <- mexDat[complete.cases(mexDat),]

print(nrow(mexDat[mexDat$label == 'CD8',]))
print(nrow(mexDat[mexDat$label == 'FoxP3',]))

# hdbscan
cd8vis <- largeVis(t(as.matrix(mexDat[, 1:2])),K = 4, dim=2, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)

clusters <- largeVis::hdbscan(cd8vis,K = 4, minPts = 30, verbose = FALSE, threads = 1)
mexDat$cluster <- clusters$clusters
mexDat <- mexDat[complete.cases(mexDat),]


# count number 
print(nrow(mexDat[mexDat$label == 'CD8',])/nrow(CD8_Pts))
print(nrow(mexDat[mexDat$label == 'FoxP3',])/nrow(FoxP3_Pts))


# get polygon infomation for clusters

maxCluster = max(as.numeric(mexDat$cluster), na.rm = TRUE)
minCluster = min(as.numeric(mexDat$cluster), na.rm = TRUE)
notNA = !is.na(mexDat$cluster)
polygonInfo <- matrix(nrow = 0, ncol = 3) #clusterID, x, y
alpha0 <- 0.01
alphaStep <- 0.01
# check number of clusters in each region

N_region <- readRDS('./Case_5/Normal.rds')
IF_region <- readRDS('./Case_5/Invasive.rds')
CT_region <- readRDS('./Case_5/Tumor.rds')

X_IF <- 0
X_N <- 0
X_CT <- 0
for(cid_target in 1:maxCluster) {
  print(cid_target)
  goodAlpha = FALSE
  alpha <- alpha0
  #X<-subPoints[cidMatch,]
  if(isTRUE(nrow(mexDat[mexDat$cluster == cid_target & mexDat$label == 'CD8',1:2])*nrow(mexDat[mexDat$cluster == cid_target & mexDat$label == 'FoxP3',1:2]) != 0)){
    X <- mexDat[mexDat$cluster == cid_target,1:2]
    X <- unique(X)
    # try different alpha until shape covers all points
    while(!goodAlpha){
      #print(alpha)
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
    N_num <- data.frame(point.in.polygon(X[,1], X[,2], N_region[,1], N_region[,2]))
    IF_num <-data.frame(point.in.polygon(X[,1], X[,2], IF_region[,1], IF_region[,2]))
    CT_num <-data.frame(point.in.polygon(X[,1], X[,2], CT_region[,1], CT_region[,2]))
    
    N_num <- length(which(N_num[,1] != 0))
    IF_num <- length(which(IF_num[,1] != 0))
    CT_num <- length(which(CT_num[,1] != 0))
    
    
    if(isTRUE(IF_num != 0)){
      X_IF <- X_IF + 1
    }
    if(isTRUE(IF_num == 0 & CT_num != 0)){
      X_CT <- X_CT + 1
    }
    if(isTRUE(IF_num == 0 & CT_num == 0 & N_num != 0)){
      X_N <- X_N + 1
    }
    
  }
}

 # check the number associated with each region


colnames(polygonInfo) <- c('cluster_id', 'x', 'y')
polygonInfo <- data.frame(polygonInfo)

###### for plot 1A and 1B together
Case0_polygon <- polygonInfo
Case0_mex <- mexDat

# 1B
Case1_polygon <- polygonInfo
Case1_mex <- mexDat



###########
# Fig. 10 #
###########

jpeg('./DoC_Cluster_CD8_high.jpeg', units="in", width=15.9, height=9, res=300)

move_layers(for_plot, idx = 1L, position = "bottom") + 
  theme_bw(base_rect_size = 1) +
  geom_polygon(data = Poly1_path, aes(x = Poly1_path[,1], y = Poly1_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly2_path, aes(x = Poly2_path[,1], y = Poly2_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Poly3_path, aes(x = Poly3_path[,1], y = Poly3_path[,2], fill = 'Stroma'), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = Invasive_region, aes(x = Invasive_region[,1], y = Invasive_region[,2], fill = 'Invasive front'),alpha = 0.3, show.legend = FALSE) +
  geom_polygon(data = Tumor_region, aes(x = Tumor_region[,1], y =Tumor_region[,2], fill = 'Tumor'),alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F26969','#02f527','#F2F2C6')) + 
  geom_polygon(data = Case0_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  geom_point(data = Case0_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  geom_point(data = Case0_mex, aes(x, y, group = cluster, color = factor(cluster)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  
  geom_polygon(data = Case1_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  geom_point(data = Case1_polygon, aes(x, y, group = cluster_id, color = factor(cluster_id)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  geom_point(data = Case1_mex, aes(x, y, group = cluster, color = factor(cluster)), alpha = 0.8,show.legend = FALSE, size = 0.5)+
  
  theme(text = element_text(size = 36)) +
  xlim(0, 28.072) +
  ylim(20.824, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()
dev.off()


ggplot() +   geom_polygon(data = polygonInfo, aes(x, y, group = cluster_id, color = factor(cluster_id)), fill = NA, show.legend = FALSE)+
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 24)) +
  xlim(0, 31.872) +
  ylim(18.024, 0) +
  xlab('x, mm') +
  ylab('y, mm') +
  coord_equal()

# Histogram for CD8 FoxP3
jpeg('./Histogram_CD8.jpeg', units="in", width=6, height=2, res=300)

ggplot(CD8_Pts_Score, aes(x=DoC)) + 
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 12)) +
  geom_histogram(color="black", fill="black", binwidth = 0.01) +
  xlab(NULL) +
  ylab('Frequency,CD8') +
  geom_vline(xintercept = 0.84, color = 'red', linetype = 'dashed', size =1)
dev.off()
jpeg('./Histogram_FoxP3.jpeg', units="in", width=6, height=2, res=300)
ggplot(FoxP3_Pts_Score, aes(x=DoC)) + 
  theme_bw(base_rect_size = 1) +
  theme(text = element_text(size = 12)) +
  geom_histogram(color="black", fill="black", binwidth = 0.01) +
  xlab('DoC Score') +
  ylab('Frequency,FoxP3') +
  geom_vline(xintercept = 0.84, color = 'red', linetype = 'dashed', size =1)
dev.off()

