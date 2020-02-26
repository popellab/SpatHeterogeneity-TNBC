'''This script is used to generate coordinates for each pixel in a downsampled WSI image'''
'''Output will be used to construct distance-density profile for immune markers'''
import numpy
import pandas as pd

# the product of the width and length of the downsampled image
Pixel_coords  = numpy.zeros(shape=(10496889,2))

print(Pixel_coords[1])
idx = 0
col_count = 0
for col in numpy.arange(0.5, 2847.5, 1):
    col_count = col_count + 1
    row_count = 0
    for row in numpy.arange(0.5, 3687.5, 1):
        #print(col)
        row_count = row_count + 1
        Pixel_coords[idx][0] = row
        Pixel_coords[idx][1] = col
        #Pixel_coords[idx][2] = col_count
        #Pixel_coords[idx][3] = row_count
        #print(Pixel_coords[idx])
        idx = idx + 1
        #print(idx)

#print(Pixel_coords[0])
#print(Pixel_coords[10496888])
df = pd.DataFrame(Pixel_coords)

# save coordinate to file
df.to_csv('Pixel_coords_Case4_FoxP3.csv', index=False, header=False)
