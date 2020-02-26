%%
%-------This script reconstructs whole tumor cell coords from subregions coords-------%

% Read origin of subregions from txt file

Coordinate_data = dlmread('./Sources/Case 4/CD3/sub_origin_CD3.txt','\t');
% Combine points
Full_point = Coordinate_data;
Point_All = zeros(0,0);

%-------Reconstruct point patterns from subregions--------%
Current_path = pwd;
% Read point pattern file from subreagions.
PathRoot = './Sources/Case 4/CD3/Coords/';
cd('./Sources/Case 4/CD3/Coords')
% Counting file numbers
fileNum=length(dir(pwd))-2;
%for num = 0:fileNum - 1
num = 0;
 while(num <= fileNum - 1)
    Filename = sprintf('DigitalPathology%d_allPoint.txt',num);
    Test_empty = dir(Filename);
    if Test_empty.bytes ~= 0 && exist(Filename,'file') == 2 
        Pattern_data = dlmread(Filename);
        % Pattern file in .mat format for calculation
        id = num + 1;
        origin_x = Full_point(id,1);
        origin_y = Full_point(id,2);
        
        % Convert subregion coordinates into whole slide coordinates
        Pattern_data(:,1) = (Pattern_data(:,1) + origin_x);
        Pattern_data(:,2) = (Pattern_data(:,2) + origin_y);
        % Store converted coordinates into file
        Point_All = cat(1, Point_All, Pattern_data);
        %id = num + 1;
        num = num + 1;
    else
        sprintf('This file is empty: %d', num)
        num = num + 1;
    end
    
        
end
%%
cd(Current_path)
% Point_corrected = transformPointsForward(CD3_CD4_Reg.Transformation, Point_All);
% Point_corrected(:,1) = Point_corrected(:,1)*3984/4044;
% Point_corrected(:,2) = Point_corrected(:,2)*2253/2192;
% save Point_corrected_CD4.mat Point_corrected
Point_All = Point_All/16;       
%save Point_All_CD4.mat Point_All
%%
plot(Point_All(:,1), Point_All(:,2))