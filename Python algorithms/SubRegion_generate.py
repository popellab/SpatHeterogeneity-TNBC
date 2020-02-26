import openslide
import numpy
import time
import cv2

start = time.time()

'''read WSI images'''
slide = openslide.open_slide("/Users/mihaoyang/Desktop/Pathology_Slides/Case 1/350713_case1_CD3.svs")
[column0, row0] = slide.level_dimensions[0]

'''Define parameters for filtering'''
k = 0
list = [0]
d = 0
color_distance = []
R_accumulate = []
G_accumulate = []
B_accumulate = []
background_counting = 0

'''Claim a matrix to store RGB information for one pixel'''
RGB_combine_matrix = []

'''Define reference color to calculate Euclidean Distance'''
vec_ref_color_bg = numpy.array([243, 243, 243])
vec_ref_color_bl = numpy.array([5, 6, 10])

'''Loop through the image and record the origin of each segmented tile'''

with open('sub_origin_CD3.txt', 'w') as output_file:

    output_file.write('x_coordinates' + '\t' + 'y_coordinates\n')
    output_file.write('---------------------------------------\n')

for i in range(0,row0, 399):
    for j in range(0, column0, 399):
        if j > column0 and i > row0:
            break
        dark_counting = 0
        background_counting = 0
        less = 0
        high = 0
        w1 = min(row0 - i, 400)
        w2 = min(column0 - j, 400)
        list[k] = numpy.array(slide.read_region((j, i), 0, (w2, w1)))
        start = time.time()
        for m in range(w1):
            for n in range(w2):
                R = list[k][m][n][0]
                G = list[k][m][n][1]
                B = list[k][m][n][2]
                RGB_array = [R, G, B]
                R_accumulate.append(R)
                B_accumulate.append(B)

                RGB_combine_matrix.append(RGB_array)

        vec1 = numpy.array(RGB_combine_matrix)
        Dist_comb = numpy.linalg.norm(vec1 - vec_ref_color_bg, axis=1)
        Dist_comb_dk = numpy.linalg.norm(vec1 - vec_ref_color_bl, axis=1)

        '''Manually tuned threshold to rule out background or overstained subregions'''
        
        Background_counting = numpy.count_nonzero(Dist_comb < 10)
        dark_counting = numpy.count_nonzero(Dist_comb_dk < 10)
        if Background_counting/(w1*w2) < 0.25 and dark_counting < 100 and numpy.sum(B_accumulate) != 0:
            filename = '//Users//mihaoyang//Desktop//Segmented_Images//Case1//CD3//' + 'DigitalPathology' + str(d) + '.png'
            d = d + 1
            cv2.imwrite(filename, list[k])
            file = open('sub_origin_CD3.txt', 'a')
            file.write('{}\t{}'.format(j,  i) + '\n')
        end = time.time()
        print("Loop finished, cost: %f second" % (end-start))
        print(i, j)

        k = k + 1
        list.append('0')
        color_distance = []
        R_accumulate = []
        G_accumulate = []
        B_accumulate = []
        RGB_combine_matrix = []