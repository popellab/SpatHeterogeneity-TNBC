#######################################################################
##### This script is used for analyzing intra-tumoral heterogeenity####
#######################################################################

# annotation.R script should be run first, to obtain polygon information

library(spatstat)
library(R.matlab)
library(ggplot2)
library(scales) # to access break formatting functions
library(dplyr)
library(kulife)
library(concaveman)
library(rgeos)
source('MiFunction.R')

p_th = 5E-2
# Define global variable
Var_suffix <- '_Points.mat'
outname = "fittedResult_"
fractal_stat_total <- data.frame(matrix(nrow = 0, ncol = 3))
names(fractal_stat_total) <- c('window_size', 'num', 'QCoD')

fractal_path <- 'Fractal_Analysis'
# Change case here
sub_path <- switch(6,'./Case_0/' ,'./Case_1/' , './Case_2/', './Case_3/', './Case_4/', './Case_5/')
# Loop to calculate statistics

# For case 2
Current_lab <- switch(5,'CD3','CD4','CD8','CD20','FoxP3')
Coords_path <- paste(sub_path, Current_lab, '/',Current_lab, Var_suffix, sep = '')
marker_Coords <- as.data.frame(readMat(Coords_path))/125

# Run before change the current case

header = c("xstart","xend", "ystart","yend","n","intensity","notCSR", "fit", "kappa", "sigma2", "mu")
# Scanning the whole space
  
    
for(window in seq(0.1, 5, 0.1)){
      x_origin <- 0
      y_origin <- 0
      p_all <- matrix(nrow=0,ncol=11)
      colnames(p_all) = header
      # load coordinates for current case
      # Create files for store fit results for current case 
      for(y_origin in seq(0,25 - window, window/2)){
        for(x_origin in seq(0,30 - window, window/2)){
          SubRegion <- marker_Coords[marker_Coords[,1] <= (x_origin+window)  & marker_Coords[,1] >= x_origin & marker_Coords[,2] <= (y_origin + window) & marker_Coords[,2] >= y_origin,]
          SubRegion <- SubRegion[complete.cases(SubRegion), ]
          if(isTRUE(row(SubRegion) == 0)){
            p = c(0,0,0,0,'NA','NA','NA')
            p_all = rbind(p_all, c(x_origin,x_origin+window,y_origin,y_origin+window,p))
          } 
          if(isTRUE(nrow(SubRegion) != 0 )){
            # fit model to data
            n = nrow(SubRegion)
            intensity = n/(window * window)
            p = subRegionFit(SubRegion, x_origin,y_origin, window, alpha)
            
            p_all = rbind(p_all, c(x_origin,x_origin+window,y_origin,y_origin+window,n, intensity, p))
            
          }
        }
      }
      p_all <- data.frame(p_all)
      p_all <- p_all[p_all$notCSR == 1 & p_all$sigma2 != 0, ]
      
      # write results to file
      count_sub <- nrow(p_all)
      sigma <- as.numeric(p_all$sigma2)
       
      cluster_size <- as.data.frame(sigma*2*pi*log(20))
      Q1 <- quantile(cluster_size[,1], 0.25)
      Q3 <- quantile(cluster_size[,1], 0.75)
      QCoD <- (Q3 - Q1)/ (Q3 + Q1)
      stat <- cbind(window, count_sub, QCoD)
      names(stat) <- c('window_size', 'num','QCoD')
      
      fractal_stat_total <- rbind(fractal_stat_total, stat)
      #write.csv(p_all, file = paste(fractal_path, '/','fitResult_', window, '.csv', sep = ''), row.names=FALSE)
      
      print('finish')
}
write.csv(fractal_stat_total, file = paste(sub_path, Current_lab,'_fractalAnalysis', '.csv', sep = ''), row.names=FALSE)



sub_path <- switch(1, 'Case_0', 'Case_1' , 'Case_2', 'Case_3', 'Case_4', 'Case_5')
Current_lab <- switch(1,'CD3','CD4','CD8','CD20','FoxP3')


fractal_stat_total <- read.csv(paste('./',sub_path,'/', Current_lab,'_fractalAnalysis', '.csv', sep = ''))


##### Plotting ####


jpeg(paste('./',sub_path, '/',Current_lab, '/', 'fractal_number.jpeg', sep =''), units="in", width=10, height=8, res=1200)
# insert ggplot code

# Plot the number of aggrevated window vs window size
ggplot() + geom_line(data = fractal_stat_total, aes(x = fractal_stat_total[,1], y = fractal_stat_total[,2])) +
  theme_bw(base_rect_size = 1) +
  geom_point(data = fractal_stat_total, aes(x = fractal_stat_total[,1], y = fractal_stat_total[,2])) +
  theme(text = element_text(size = 32)) + 
  #ggtitle(paste(sub_path, '-', Current_lab, sep = '')) +
  xlab('Window, mm') +
  ylab('N(window, aggregated)') +
  scale_x_log10(breaks = c(seq(0.1, 0.3, 0.1), seq(0.4, 0.8, 0.2), 1:5), 
                labels = c(seq(0.1, 0.3, 0.1), seq(0.4, 0.8, 0.2), 1:5)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 
dev.off()

# Plot the number of QCoD (Quartile Coefficient of Dispersion) vs window size
sub_path <- switch(6, 'Case_0', 'Case_1' , 'Case_2', 'Case_3', 'Case_4', 'Case_5')
Current_lab <- switch(5,'CD3','CD4','CD8','CD20','FoxP3')
fractal_stat_total <- read.csv(paste('./',sub_path,'/', Current_lab,'_fractalAnalysis', '.csv', sep = ''))

jpeg(paste('./',sub_path, '/',Current_lab, '/', 'fractal_QCoD.jpeg', sep =''), units="in", width=10, height=8, res=1200)
ggplot() + geom_line(data = fractal_stat_total, aes(x = fractal_stat_total[,1], y = fractal_stat_total[,3])) +
  geom_point(data = fractal_stat_total, aes(x = fractal_stat_total[,1], y = fractal_stat_total[,3])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 35)) +
  theme(axis.title.y = element_text(size = 35)) +
  theme(axis.text = element_text(size = 25)) +
  xlab('Window, mm') +
  ylab('QCoD') +
  scale_x_log10(breaks = c(seq(0.1, 0.3, 0.1), seq(0.4, 0.8, 0.2), 1:5), 
                labels = c(seq(0.1, 0.3, 0.1), seq(0.4, 0.8, 0.2), 1:5)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(plot.title = element_text(size = 30))
dev.off()