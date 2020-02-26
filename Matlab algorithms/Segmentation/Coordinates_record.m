% This function is used to save the coordinate of the centroid of each
% selected cells.

function Coordinates_record(selected_img, num)
bw_c = regionprops(selected_img,'centroid');
centroids_coordinates = cat(1, bw_c.Centroid);

% Write coordinates to txt file
FileName = sprintf('./Sources/Case 4/CD3/Coords/DigitalPathology%d_allPoint.txt',num);
%FileName = sprintf('./CD3_Case2_Points/DigitalPathology%d_allPoint.txt',num);

dlmwrite(FileName,centroids_coordinates,'delimiter','\t');
