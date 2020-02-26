%----------------------------Function file----------------------------%
% This file is used to finalized the cell list for coordinate extraction
% for the current subregion.
% Only the cells that passed the morphometric filters (i.e. the circularity and area are both below the thresholds)
%will be recorded as
% the final algorithmic detected cells. 


function selected_region = regionSelection(BinaryImg, IndexList)

%-------------------Get the Pixels Id in each region----------------------%
[row, col] = size(BinaryImg);
PixelIdList = regionprops(BinaryImg, 'PixelIdxList');

%-------------------Get the size of PixelIdList---------------------------%

Length_Pixel_Id = length(PixelIdList);

%-------------------Get the Pixels of ruled out regions-------------------%

Index_length = length(IndexList);
for i = 1:Index_length
    RegionId = IndexList(i);
    Region_Pixels_Collective = PixelIdList(RegionId);

    Region_Pixels_Collective = struct2cell(Region_Pixels_Collective);
    Region_Pixels_Collective = cell2mat(Region_Pixels_Collective);
    %Num2handle = length(Region_Pixels_Collective);

    
    %----------Loop through the regions to be handled----------%
    Region_idx = Region_Pixels_Collective;
    %------------Get the number of pixles in this region-----------%
    Pixel_num = length(Region_idx);
    %-----------------For each pixel in this region----------------%
    for j = 1:Pixel_num
        Pixel = Region_idx(j);
        %-----------The integer part as row, remainder part as column-----% 
        coordinate_x = mod(Pixel,row);
        coordinate_y = fix(Pixel/row);
        if coordinate_x > 0 && coordinate_y > 0
            BinaryImg(coordinate_x,coordinate_y+1) = 0;
        end
        if coordinate_y == 0 
            coordinate_y = 1;
            BinaryImg(coordinate_x,coordinate_y) = 0; 
        end            
    end
    
end
%------------------------Eliminate noise----------------------------------%
selected_region = BinaryImg;
se=strel('disk',1);
selected_region=bwareaopen(selected_region,20);
selected_region=imclose(selected_region,se);

%-----------------------Redo watershed algorithm--------------------------%
D = -bwdist(~selected_region);
    
mask=imextendedmin(D,4);    
D2=imimposemin(D,mask);
Ld=watershed(D);
Water_splited=selected_region;
Water_splited(Ld==0)=0;
selected_region=Water_splited;