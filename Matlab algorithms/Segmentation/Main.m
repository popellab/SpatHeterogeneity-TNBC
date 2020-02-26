%----------------------------Image Information----------------------------%
% With treatment, breast cancer.
% In this study, the resolution of WSIs is ~0.497 microns per pixel.
% 1 Pixel = 0.497 Microns, 1 Pixel^2  = 0.247 Microns^2

%%%% Note: Don't forget to change file path in this script and:
%%%% - Coordinate_record.m
%%%% before switching the slides

format compact, close all, clear all, clc;
% alternative set of standard values (HDAB from Fiji)

% Batch file processing
PathRoot = './Sources/Case 4/CD3/subregions'; % File path where subregions are stored

% Counting numbers of tiles
fileNum=length(dir(PathRoot))-2;
Cell_tot = zeros(fileNum, 1);
CMYK_dist = zeros(1,0);

for num = 0:fileNum-1
    num; % print current loop index
    imageFileName = sprintf('DigitalPathology%d.png',num); 

    %------------Claim file paths------------%
    Current_path = pwd;
    
    Segment_image_path = './Sources/Case 4/CD3/subregions'; %'./Sources/Case5/CD20'
    cd(Segment_image_path)
    img=imread(imageFileName);
    [row, col,dim] = size(img);
    cd(Current_path)

    img_double = im2double(img)*255;
    channel_1 = img_double(:,:,1);
    channel_2 = img_double(:,:,2);
    channel_3 = img_double(:,:,3);
    Color_distance = zeros(row,col);
    split_1 = img(:,:,1);
    split_2 = img(:,:,2);
    split_3 = img(:,:,3);
    
% Introduce CMYK model, filter background images
    img_cmyk = rgb2cmyk(img);
    img_yellow = img_cmyk(:,:,3);
    CMYK_dist = [CMYK_dist;mean(mean(img_yellow))];
    % Set manually tuned thresholds to rule out background and overstained
    % tiles
    if  mean(mean(channel_3)) >= 120  && mean(mean(img_yellow)) < 35 && mean(mean(channel_3)) < 220 % 2 for yellow channel
        for i=1:row
            for j = 1:col
                % Set reference color vector
                Color_distance(i, j) = sqrt((channel_1(i,j)- .75)^2+ (channel_2(i,j)-47.93)^2 + (channel_3(i,j)-91.5543)^2);
                % Set manually tuned thresholds to allow color variability
                if Color_distance(i,j) < 155
                    Color_distance(i,j) = 1;
                else
                    Color_distance(i,j) = 0;
                end
            end
        end
    end
    Dist_max = max(Color_distance);
    Dist_logical = logical(Color_distance);
    %figure;
    %imshow(Dist_logical)
    
    % Morphological operations
    se=strel('disk',1);
    bw3=bwareaopen(Dist_logical,20);
    bw4=imclose(bw3,se);
    
    %figure
    %imshow(bw4);
    
    D = -bwdist(~bw4);

    mask=imextendedmin(D,4);
    D2=imimposemin(D,mask);
    Ld=watershed(D);
    Water_splited=bw4;
    Water_splited(Ld==0)=0;
    show_img=Water_splited;

    %-----------------------Calculate Areas-----------------------------------%
    
    Area_aggregate= regionprops(show_img, 'Area');
    Area_cell = struct2cell(Area_aggregate);
    Area = cell2mat(Area_cell);
    cell_number_initial = length(Area);
    
    %---------------------Calculate Major/Minor AxisLength--------------------%
    Major_axis_struct = regionprops(show_img,'MajorAxisLength');
    Major_axis_cell = struct2cell(Major_axis_struct);
    Major_axis = cell2mat(Major_axis_cell);
    
    Minor_axis_struct = regionprops(show_img,'MinorAxisLength');
    Minor_axis_cell = struct2cell(Minor_axis_struct);
    Minor_axis = cell2mat(Minor_axis_cell);
    Ovality = Major_axis - Minor_axis;
    Ovality_mean = mean(Ovality);

    %---------------------Calculate Perimeters--------------------------------%
    
    Perimeter_aggregate = regionprops(show_img, 'Perimeter');
    Peri_cell = struct2cell(Perimeter_aggregate);
    Perimeter = cell2mat(Peri_cell);
    
    %-------------------Find the microns area of cells------------------------%
    Area_microns = Area.*0.247; % convert the unit from pixel^2 -> um^2
    Area_mean = mean(Area_microns);
    Area_std = std(Area_microns);
    
    %-----------------Find the microns perimeters of cells--------------------%
    
    Perimeter_microns = Perimeter.*0.497;
    Peri_mean = mean(Perimeter_microns);
    Peri_std = std(Perimeter_microns);
    
    %------------------Calculate the circularity for each cell----------------%
    Circularity = zeros(cell_number_initial,1);
    for i=1:cell_number_initial
        Circularity(i) = (Perimeter_microns(i).^2)./(4.*pi.*Area_microns(i));
    end
    cell_number_correct = 0;
    IndexList = [];
    for i=1:cell_number_initial
        if Circularity(i) > 0.7 && Area_microns(i) > 10 && Area_microns(i) < 85
            cell_number_correct =  cell_number_correct + 1;
        else    
            IndexList = [IndexList;i];
        end
    end
    % Function: regionSelect is used to finalize the cell list from
    % candidate pool, based on above morphometrics
    selected_img = regionSelection(show_img, IndexList);
    s = regionprops(selected_img,'centroid');
    
    centroids = cat(1,s.Centroid);
%     if isempty(centroids) ~= 1
%         fig = figure;
%         imshow(img)
%         hold on
%         plot(centroids(:,1),centroids(:,2),'color','[0.156, 0.996, 0.535]','linestyle','none','marker','+','MarkerSize',10)
%         hold off
%         set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
%         set(gcf,'Units','pixels','Position',[100 100 100 100]);  %# Modify figure size
%         f=getframe(gcf);
%         imwrite(f.cdata,MarkFileName)
%     end
    coordinates_accumulate = zeros(0,2);
    Coordinates_record(selected_img,num);
    Circularity_mean = mean(Circularity);
    
    perim=bwperim(selected_img);
    r=img(:,:,1);
    g=img(:,:,2);
    b=img(:,:,3);
    r(perim)=255;
    g(perim)=0;
    b(perim)=0;
    img(:,:,1)=r;
    img(:,:,2)=g;
    img(:,:,3)=b;


    if cell_number_correct ~= 0
        filepath = pwd;
        cd('./Sources/Case 4/CD3/Processed subregions') %'./Sources/Case5/CD20'
        imwrite(img,imageFileName)
        cd(filepath)
        fprintf('The number of cells are: %f \n',cell_number_correct);
        fprintf('The Circularity_mean is: %f, \n The Ovality_mean is: %f, \n The Area_mean is: %f, \n, The perimeter is: %d \n', Circularity_mean, Ovality_mean, Area_mean, Peri_mean);
        
        Cell_tot(num+1) = cell_number_correct;
    end
end
Cell_total = sum(Cell_tot);