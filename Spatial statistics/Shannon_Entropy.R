######################################################################################
## This script measures inhomogeneity using a modified version of Shannon's Entropy ##
######################################################################################

library(ggplot2)
library(R.matlab)
library(cowplot)
library(flexclust)
library(FNN)
library(LaplacesDemon)

source('MiFunction.R')

Var_suffix <- '_Points.mat'

p_all <- data.frame(matrix(nrow=0,ncol= 9))
header <- c('index',"xstart","xend", "ystart","yend",'Clara','Shannon En', 'location', 'Case')

colnames(p_all) <- header

for(case in seq(1,6,1)){
  
  case_name <- switch(case,'Case 0', 'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5')
  
  # Change case here
  sub_path <- switch(case,'./Case_0/', './Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
  # Loop to calculate statistics
    
    CD3_Pts <- as.data.frame(readMat(paste(sub_path, 'CD3/CD3_Points.mat', sep = '')))/125
    CD4_Pts <- as.data.frame(readMat(paste(sub_path, 'CD4/CD4_Points.mat', sep = '')))/125
    CD8_Pts <- as.data.frame(readMat(paste(sub_path, 'CD8/CD8_Points.mat', sep = '')))/125
    CD20_Pts <- as.data.frame(readMat(paste(sub_path, 'CD20/CD20_Points.mat', sep = '')))/125
    FoxP3_Pts <- as.data.frame(readMat(paste(sub_path, 'FoxP3/FoxP3_Points.mat', sep = '')))/125
    #######################
    if(case == 1){
      InvasiveFile <- readRDS(paste(sub_path, 'Invasive_1a.rds', sep = ''))
      
      ##########
      
      TumorFile <- readRDS(paste(sub_path, 'Tumor_1a.rds', sep = ''))
      NormalFile <- readRDS(paste(sub_path, 'Normal_1a.rds', sep = ''))
    }
    if(case == 2){
      
      InvasiveFile <- readRDS(paste(sub_path, 'Invasive_1b.rds', sep = ''))
      
      TumorFile <- readRDS(paste(sub_path, 'Tumor_1b.rds', sep = ''))
      
      NormalFile <- readRDS(paste(sub_path, 'Normal_1b.rds', sep = ''))
    }
    if(case > 2){
      
      InvasiveFilePath <- paste(sub_path, 'Invasive.rds', sep = '')
      InvasiveFile <- readRDS(InvasiveFilePath)
      
      TumorFilePath <- paste(sub_path, 'Tumor.rds', sep = '')
      TumorFile <- readRDS(TumorFilePath)
      
      NormalFile <- readRDS(paste(sub_path, 'Normal.rds', sep = ''))
      
    }
    
    window <- 0.15
    
    # Index for each window
    index <- 0
    # Create files for store fit results for current case 
    for(y_origin in seq(0, 25-window, window)){
      for(x_origin in seq(0,35-window, window)){
        window_coords <- data.frame(cbind(c(x_origin, x_origin, x_origin+window, x_origin+window), c(y_origin, y_origin + window, y_origin + window, y_origin)))
        
        # Characterize the window
        if(isTRUE(ncol(NormalFile) == 3)){
          normal_polygon_num <- NormalFile[nrow(NormalFile),3]
          normal_count <- 0
          for(num in seq(1, normal_polygon_num)){
            subNormalFile <- NormalFile[NormalFile$label == num, ]
            normal <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], subNormalFile[,1], subNormalFile[,2]))
            subnormal_count <- nrow(data.frame(normal[normal[,1] >0,]))
            normal_count <- normal_count + subnormal_count
          }
        }
        if(isTRUE(ncol(NormalFile) == 2)){
          normal_count <- 0
          normal <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], NormalFile[,1], NormalFile[,2]))
          normal_count <- nrow(data.frame(normal[normal[,1] >0,]))
        }
        
        
        invasive <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], InvasiveFile[,1], InvasiveFile[,2]))
        invasive_count <- nrow(data.frame(invasive[invasive[,1] >0,]))
        
        tumor <- data.frame(point.in.polygon(window_coords[,1], window_coords[,2], TumorFile[,1], TumorFile[,2]))
        tumor_count <- nrow(data.frame(tumor[tumor[,1] >0,]))
        
        
        if(isTRUE(invasive_count != 0)){
          loc <- 'IF'
        } 
        if(isTRUE(invasive_count == 0 & tumor_count != 0)){
          loc <- 'CT'
        }
        if(isTRUE(invasive_count == 0 & tumor_count == 0 & normal_count != 0)){
          loc <- 'N'
        }
        if(isTRUE(invasive_count == 0 & tumor_count == 0 & normal_count == 0)){
          loc <- 'BG'
        }
        # customizable area
        if(isTRUE(loc != 'BG')){
          CD3_sub <- CD3_Pts[CD3_Pts[,1] < (x_origin+window)  & CD3_Pts[,1] >= x_origin & CD3_Pts[,2] < (y_origin + window) & CD3_Pts[,2] >= y_origin,]
          CD4_sub <- CD4_Pts[CD4_Pts[,1] < (x_origin+window)  & CD4_Pts[,1] >= x_origin & CD4_Pts[,2] < (y_origin + window) & CD4_Pts[,2] >= y_origin,]
          CD8_sub <- CD8_Pts[CD8_Pts[,1] < (x_origin+window)  & CD8_Pts[,1] >= x_origin & CD8_Pts[,2] < (y_origin + window) & CD8_Pts[,2] >= y_origin,]
          CD20_sub <- CD20_Pts[CD20_Pts[,1] < (x_origin+window)  & CD20_Pts[,1] >= x_origin & CD20_Pts[,2] < (y_origin + window) & CD20_Pts[,2] >= y_origin,]
          FoxP3_sub <- FoxP3_Pts[FoxP3_Pts[,1] < (x_origin+window)  & FoxP3_Pts[,1] >= x_origin & FoxP3_Pts[,2] < (y_origin + window) & FoxP3_Pts[,2] >= y_origin,]
        
          FoxP3_sub <- FoxP3_sub[complete.cases(FoxP3_sub), ]
          CD20_sub <- CD20_sub[complete.cases(CD20_sub), ]
          CD8_sub <- CD8_sub[complete.cases(CD8_sub), ]
          CD4_sub <- CD4_sub[complete.cases(CD4_sub), ]
          CD3_sub <- CD3_sub[complete.cases(CD3_sub), ]
          
          Total <- nrow(FoxP3_sub) + nrow(CD3_sub) + nrow(CD4_sub) + nrow(CD8_sub) + nrow(CD20_sub)
          
          CD3_int <- 0
          CD3_ext <- 0
          p <- 0
          if(isTRUE(nrow(CD3_sub) != 0)){
            p <- nrow(CD3_sub)/Total
            CD3_int <- mean(as.matrix(dist(CD3_sub[,-1])))
            
            if(isTRUE(nrow(CD4_sub)!=0)){
              CD3_ext <- CD3_ext + mean(dist2(CD3_sub, CD4_sub))
            }
            
            if(isTRUE(nrow(CD8_sub)!=0)){
              CD3_ext <- CD3_ext + mean(dist2(CD3_sub, CD8_sub))
            }
            
            if(isTRUE(nrow(CD20_sub)!=0)){
              CD3_ext <- CD3_ext + mean(dist2(CD3_sub, CD20_sub))
            }
            
            if(isTRUE(nrow(FoxP3_sub)!=0)){
              CD3_ext <- CD3_ext + mean(dist2(CD3_sub, FoxP3_sub))
            }
          }
          CD3_dat <- cbind(p, CD3_int, CD3_ext/5)
          colnames(CD3_dat) <- c('p', 'int', 'ext')
          ##### CD4
          CD4_int <- 0
          CD4_ext <- 0
          p <- 0
          if(isTRUE(nrow(CD4_sub) != 0)){
            p <- nrow(CD4_sub)/Total
            CD4_int <- mean(as.matrix(dist(CD4_sub[,-1])))
            
            if(isTRUE(nrow(CD3_sub)!=0)){
              CD4_ext <- CD4_ext + mean(dist2(CD4_sub, CD3_sub))
            }
            
            if(isTRUE(nrow(CD8_sub)!=0)){
              CD4_ext <- CD4_ext + mean(dist2(CD4_sub, CD8_sub))
            }
            
            if(isTRUE(nrow(CD20_sub)!=0)){
              CD4_ext <- CD4_ext + mean(dist2(CD4_sub, CD20_sub))
            }
            
            if(isTRUE(nrow(FoxP3_sub)!=0)){
              CD4_ext <- CD4_ext + mean(dist2(CD4_sub, FoxP3_sub))
            }
          }
          CD4_dat <- cbind(p, CD4_int, CD4_ext/5)
          colnames(CD4_dat) <- c('p', 'int', 'ext')
          ##### CD8
          CD8_int <- 0
          CD8_ext <- 0
          p <- 0
          if(isTRUE(nrow(CD8_sub) != 0)){
            p <- nrow(CD8_sub)/Total
            CD8_int <- mean(as.matrix(dist(CD8_sub[,-1])))
            
            if(isTRUE(nrow(CD3_sub)!=0)){
              CD8_ext <- CD8_ext + mean(dist2(CD8_sub, CD3_sub))
            }
            
            if(isTRUE(nrow(CD4_sub)!=0)){
              CD8_ext <- CD8_ext + mean(dist2(CD8_sub, CD4_sub))
            }
            
            if(isTRUE(nrow(CD20_sub)!=0)){
              CD8_ext <- CD8_ext + mean(dist2(CD8_sub, CD20_sub))
            }
            
            if(isTRUE(nrow(FoxP3_sub)!=0)){
              CD8_ext <- CD8_ext + mean(dist2(CD8_sub, FoxP3_sub))
            }
          }
          CD8_dat <- cbind(p, CD8_int, CD8_ext/5)
          colnames(CD8_dat) <- c('p', 'int', 'ext')
          ##### CD20
          CD20_int <- 0
          CD20_ext <- 0
          p <- 0
          if(isTRUE(nrow(CD20_sub) != 0)){
            p <- nrow(CD20_sub)/Total
            CD20_int <- mean(as.matrix(dist(CD20_sub[,-1])))
            
            if(isTRUE(nrow(CD3_sub)!=0)){
              CD20_ext <- CD20_ext + mean(dist2(CD20_sub, CD3_sub))
            }
            
            if(isTRUE(nrow(CD4_sub)!=0)){
              CD20_ext <- CD20_ext + mean(dist2(CD20_sub, CD4_sub))
            }
            
            if(isTRUE(nrow(CD8_sub)!=0)){
              CD20_ext <- CD20_ext + mean(dist2(CD20_sub, CD8_sub))
            }
            
            if(isTRUE(nrow(FoxP3_sub)!=0)){
              CD20_ext <- CD20_ext + mean(dist2(CD20_sub, FoxP3_sub))
            }
          }
          CD20_dat <- cbind(p, CD20_int, CD20_ext/5)
          colnames(CD20_dat) <- c('p', 'int', 'ext')
          ##### Fpx{3}
          FoxP3_int <- 0
          FoxP3_ext <- 0
          p <- 0
          if(isTRUE(nrow(FoxP3_sub) != 0)){
            p <- nrow(FoxP3_sub)/Total
            FoxP3_int <- mean(as.matrix(dist(FoxP3_sub[,-1])))
            
            if(isTRUE(nrow(CD3_sub)!=0)){
              FoxP3_ext <- FoxP3_ext + mean(dist2(FoxP3_sub, CD3_sub))
            }
            
            if(isTRUE(nrow(CD4_sub)!=0)){
              FoxP3_ext <- FoxP3_ext + mean(dist2(FoxP3_sub, CD4_sub))
            }
            
            if(isTRUE(nrow(CD8_sub)!=0)){
              FoxP3_ext <- FoxP3_ext + mean(dist2(FoxP3_sub, CD8_sub))
            }
            
            if(isTRUE(nrow(CD20_sub)!=0)){
              FoxP3_ext <- FoxP3_ext + mean(dist2(FoxP3_sub, CD20_sub))
            }
          }
          FoxP3_dat <- cbind(p, FoxP3_int, FoxP3_ext/5)
          colnames(FoxP3_dat) <- c('p', 'int', 'ext')
          
          
          Dat_collect <- rbind(CD3_dat, CD4_dat, CD8_dat, CD20_dat, FoxP3_dat)

          Claramunt <- 0
          ShannonH <- 0
            # combine row data
            
            
            for(dat in seq(1, 5)){
              p <- Dat_collect[dat,1]
              d_int <- Dat_collect[dat, 2]
              d_ext <- Dat_collect[dat, 3]
              d_final <- d_int/d_ext
                
              if(isTRUE(d_int*d_ext == 0)){
                  d_final <- 0
              }
              if(isTRUE(p != 0)){
              Claramunt <- -p*d_final*log2(p) + Claramunt
              ShannonH <- -p*log2(p) + ShannonH
              }
          }
            index <- index + 1

            p_sub <- data.frame(cbind(index,x_origin,x_origin+window,y_origin,y_origin+window,Claramunt, ShannonH, loc, case_name ))
            names(p_sub) <- header
            p_all <- rbind(p_all, p_sub)
            print(index)
          } 
          
        }
      }
    # write results to file
    print('finish')
    
}

#saveRDS(p_all, 'Shannon_Clara_Entropy.rds')

#p_all <- readRDS('Shannon_Clara_Entropy.rds')

p_all$Clara <- as.numeric(as.character(p_all$Clara))
Shannon_N_Case0 <- p_all[p_all$location == 'N' & p_all$Case == 'Case 0',]
Shannon_IF_Case0 <- p_all[p_all$location == 'IF' & p_all$Case == 'Case 0',]
Shannon_CT_Case0 <- p_all[p_all$location == 'CT' & p_all$Case == 'Case 0',]

Shannon_N_Case1 <- p_all[p_all$location == 'N' & p_all$Case == 'Case 1',]
Shannon_IF_Case1 <- p_all[p_all$location == 'IF' & p_all$Case == 'Case 1',]
Shannon_CT_Case1 <- p_all[p_all$location == 'CT' & p_all$Case == 'Case 1',]

Shannon_N_Case2 <- p_all[p_all$location == 'N' & p_all$Case == 'Case 2',]
Shannon_IF_Case2 <- p_all[p_all$location == 'IF' & p_all$Case == 'Case 2',]
Shannon_CT_Case2 <- p_all[p_all$location == 'CT' & p_all$Case == 'Case 2',]

Shannon_N_Case3 <- p_all[p_all$location == 'N' & p_all$Case == 'Case 3',]
Shannon_IF_Case3 <- p_all[p_all$location == 'IF' & p_all$Case == 'Case 3',]
Shannon_CT_Case3 <- p_all[p_all$location == 'CT' & p_all$Case == 'Case 3',]

Shannon_N_Case4 <- p_all[p_all$location == 'N' & p_all$Case == 'Case 4',]
Shannon_IF_Case4 <- p_all[p_all$location == 'IF' & p_all$Case == 'Case 4',]
Shannon_CT_Case4 <- p_all[p_all$location == 'CT' & p_all$Case == 'Case 4',]

Shannon_N_Case5 <- p_all[p_all$location == 'N' & p_all$Case == 'Case 5',]
Shannon_IF_Case5 <- p_all[p_all$location == 'IF' & p_all$Case == 'Case 5',]
Shannon_CT_Case5 <- p_all[p_all$location == 'CT' & p_all$Case == 'Case 5',]


mean(Shannon_IF$`Shannon En`)
{
p1 <- ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_density(data = Shannon_N_Case0, aes(x= Clara,color= 'N', fill= 'N'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_IF_Case0, aes(x= Clara,color= 'IF', fill= 'IF'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_CT_Case0, aes(x= Clara,color= 'CT', fill= 'CT'), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F2F2C6','#F26969','#02f527')) +
  scale_color_manual(values = c('#F2F2C6','#F26969','#02f527')) + 
  theme(text = element_text(size = 14))+
  xlab(expression(paste("Case 1A, ", H[sc]))) +
  ylab('Density')

p2 <- ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_density(data = Shannon_N_Case1, aes(x= Clara,color= 'N', fill= 'N'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_IF_Case1, aes(x= Clara,color= 'IF', fill= 'IF'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_CT_Case1, aes(x= Clara,color= 'CT', fill= 'CT'), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F2F2C6','#F26969','#02f527')) +
  scale_color_manual(values = c('#F2F2C6','#F26969','#02f527')) + 
  theme(text = element_text(size = 14))+
  xlab(expression(paste("Case 1B, ", H[sc]))) +
  ylab('Density')

p3 <- ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_density(data = Shannon_N_Case2, aes(x= Clara,color= 'N', fill= 'N'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_IF_Case2, aes(x= Clara,color= 'IF', fill= 'IF'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_CT_Case2, aes(x= Clara,color= 'CT', fill= 'CT'), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F2F2C6','#F26969','#02f527')) +
  scale_color_manual(values = c('#F2F2C6','#F26969','#02f527')) + 
  theme(text = element_text(size = 14))+
  xlab(expression(paste("Case 2, ", H[sc]))) +
  ylab('Density')

p4 <- ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_density(data = Shannon_N_Case3, aes(x= Clara,color= 'N', fill= 'N'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_IF_Case3, aes(x= Clara,color= 'IF', fill= 'IF'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_CT_Case3, aes(x= Clara,color= 'CT', fill= 'CT'), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F2F2C6','#F26969','#02f527')) +
  scale_color_manual(values = c('#F2F2C6','#F26969','#02f527')) + 
  theme(text = element_text(size = 14))+
  xlab(expression(paste("Case 3, ", H[sc]))) +
  ylab('Density')

p5 <- ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_density(data = Shannon_N_Case4, aes(x= Clara,color= 'N', fill= 'N'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_IF_Case4, aes(x= Clara,color= 'IF', fill= 'IF'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_CT_Case4, aes(x= Clara,color= 'CT', fill= 'CT'), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F2F2C6','#F26969','#02f527')) +
  scale_color_manual(values = c('#F2F2C6','#F26969','#02f527')) + 
  theme(text = element_text(size = 14))+
  xlab(expression(paste("Case 4, ", H[sc]))) +
  ylab('Density')

p6 <- 
  ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_density(data = Shannon_N_Case5, aes(x= Clara,y=..scaled..,color= 'N', fill= 'N'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_IF_Case5, aes(x= Clara,y=..scaled..,color= 'IF', fill= 'IF'), alpha = 0.5, show.legend = FALSE) +
  geom_density(data = Shannon_CT_Case5, aes(x= Clara,y=..scaled..,color= 'CT', fill= 'CT'), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values  = c('#F2F2C6','#F26969','#02f527')) +
  scale_color_manual(values = c('#F2F2C6','#F26969','#02f527')) + 
  theme(text = element_text(size = 14))+
  xlab(expression(paste("Case 5, ", H[sc]))) +
  ylab('Density')

#####################
# Plotting (Fig. 6) #
#####################

jpeg('./Shannon_Entropy.jpeg', units="in", width=10, height=6, res=300)

  plot_grid(p1, p2,p3, p4,p5,p6, labels = c('', '','','',''), ncol = 3)
dev.off()
}

# Calculate numbers on plot
pdf_of_data1 <- density(Shannon_N_Case5$Clara)
pdf_of_data2 <- density(Shannon_CT_Case5$Clara)

