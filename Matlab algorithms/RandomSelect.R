################################################################################################
# This script is used to select random subregions to evaluate the performance of the algorithm #
# This script generate a list of points, each point is the origin of the selected subregion #
################################################################################################

library(imager)
library(lhs)
library(R.matlab)
library(sampling)
# set seed for reproducity
set.seed(1000)

allPoints <- data.frame(readMat('./Case_1/CD3_Points.mat')) # read all points as the input

allPoints$labels <- seq(1, nrow(allPoints))

sample <- sample(allPoints$labels, 100, replace = FALSE, prob = NULL)

sample <- allPoints[sample,]*16

sample <- sample[,c(1,2)]

write.table(sample, file = './Sample_Coordinates/Case 1/sample_origin.txt',sep = '\t', row.names = FALSE, col.names = FALSE)

#sample <- sample/2000
#ggplot() + geom_point(data = sample ,aes(sample[,1], sample[,2]))

