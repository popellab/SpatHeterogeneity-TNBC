'''This script calcualte the DoC score for CD3+ and CD4+ points'''

import numpy as np
from scipy.spatial import cKDTree
from scipy.io import loadmat
from scipy.stats import spearmanr
import math
import pandas as pd

'''Read points coordinate'''
CD3_Points = np.array(loadmat('./Points_Bivariate/Case 1/CD3_Points.mat')['x'])/125
CD4_Points = np.array(loadmat('./Points_Bivariate/Case 1/CD4_Points.mat')['x'])/125

'''Construct KDTree (implemented in C language)'''
CD3_Tree = cKDTree(CD3_Points)
CD4_Tree = cKDTree(CD4_Points)

CD3_toCD3 = []

CD4_toCD3 = []


rho_collect = []

for pts in range(0, len(CD3_Points)):
    Point = CD3_Points[pts, ]
    for r in np.arange(0.01, 2, 0.05):
        CD3_toCD3_den = len(CD3_Tree.query_ball_point(Point, r = r))/(math.pi*r*r) # how many CD4 around CD3 -> density
        #print(CD4_toCD3_den)
        CD4_toCD3_den = len(CD4_Tree.query_ball_point(Point, r = r))/(math.pi*r*r) # how many CD3 around CD3 -> density
        # append density value from each radius level
        CD3_toCD3.append(CD3_toCD3_den)
        CD4_toCD3.append(CD4_toCD3_den)

    # calculate spearman coefficient rho for this point
    rho, pval = spearmanr(CD3_toCD3, CD4_toCD3)
    print(pts)
    rho_collect.append(rho)
    CD3_toCD3 = []
    CD4_toCD3 = []

np.savetxt("DoC_Score_Case1_CD3.csv", rho_collect, delimiter=",")

CD3_toCD4 = []

CD4_toCD4 = []

rho_collect = []

for pts in range(0, len(CD4_Points)):
    Point = CD4_Points[pts, ]
    for r in np.arange(0.01, 2, 0.05):
        CD3_toCD4_den = len(CD3_Tree.query_ball_point(Point, r = r))/(math.pi*r*r) # how many CD4 around CD3 -> density
        #print(CD4_toCD3_den)
        CD4_toCD4_den = len(CD4_Tree.query_ball_point(Point, r = r))/(math.pi*r*r) # how many CD3 around CD3 -> density
        # append density value from each radius level
        CD3_toCD4.append(CD3_toCD4_den)
        CD4_toCD4.append(CD4_toCD4_den)

    # calculate spearman coefficient rho for this point
    rho, pval = spearmanr(CD3_toCD4, CD4_toCD4)
    print(pts)
    rho_collect.append(rho)
    CD3_toCD4 = []
    CD4_toCD4 = []

np.savetxt("DoC_Score_Case1_CD4.csv", rho_collect, delimiter=",")