%-------This script registers the pathologist's annotation based on local registration transformation matrix-------%
%------- Registered annotation is used to evaluate the local registration algorithm performance ----------%

%% Load FoxP3 points
load('./Contours/Case1_FoxP3_Contour.mat');
load('./Case_1_regInfo/FoxP3_SecInfo/globalReg.mat')
FoxP3_Pts = x*125;


% Assign coordinates to sections
[len, wid] = size(FoxP3_Pts);
Char = zeros(len,1);
[length, dim] = size(FoxP3_Pts);
FoxP3_Registered = zeros(length, dim);
for i = 1:length
    for j=1:5
            [FoxP3_Registered(i,1),FoxP3_Registered(i,2)] = transformPointsForward(movingReg.Transformation,FoxP3_Pts(i,1),FoxP3_Pts(i,2));
            %FoxP3_Pts(i,1) = FoxP3_Pts(i,1) - Var.SpatialRefObj.XWorldLimits(1);
            %FoxP3_Pts(i,2) = FoxP3_Pts(i,2) - Var.SpatialRefObj.YWorldLimits(1);
            Char(i,1) = j;
    end 
end
scatter(FoxP3_Registered(:,1),FoxP3_Registered(:,2))

%% Load CD8 points
load('./Contours/Case1_CD8_Contour.mat');
load('./Case_1_regInfo/CD8_SecInfo/globalReg.mat');

CD8_Pts = x*125;
[length, dim] = size(CD8_Pts);
CD8_Registered = zeros(length, dim);

for i = 1:length
    [CD8_Registered(i,1),CD8_Registered(i,2)] = transformPointsForward(movingReg.Transformation,CD8_Pts(i,1),CD8_Pts(i,2));
end
scatter(CD8_Registered(:,1),CD8_Registered(:,2))

%% Load CD20 points
load('./Contours/Case1_CD20_Contour.mat');
load('./Case_1_regInfo/CD20_SecInfo/globalReg.mat');
CD20_Pts = x*125;

% Assign coordinates to sections

[length, dim] = size(CD20_Pts);
CD20_Registered = zeros(length, dim);

string = 'Case1_CD20_Sec';
CD20_Pts(:,1) = 3390 - CD20_Pts(:,1);
CD20_Pts(:,2) = 2440 - CD20_Pts(:,2);
%plot(CD20_Registered(:,1),CD20_Registered(:,2))
for i = 1:length 
        [CD20_Registered(i,1),CD20_Registered(i,2)] = transformPointsForward(movingReg.Transformation,CD20_Pts(i,1),CD20_Pts(i,2));
end
scatter(CD20_Registered(:,1),CD20_Registered(:,2))

 %% CD3 registration
load('./Contours/Case1_CD3_Contour.mat');
load('./Case_1_regInfo/CD3_SecInfo/globalReg.mat');

CD3_Pts = x*125;

% Assign coordinates to sections
[length, dim] = size(CD3_Pts);
CD3_Registered = zeros(length, dim);
for i = 1:length
    for j=1:4
            [CD3_Registered(i,1),CD3_Registered(i,2)] = transformPointsForward(movingReg.Transformation,CD3_Pts(i,1),CD3_Pts(i,2));
            %FoxP3_Pts(i,1) = FoxP3_Pts(i,1) - Var.SpatialRefObj.XWorldLimits(1);
            %FoxP3_Pts(i,2) = FoxP3_Pts(i,2) - Var.SpatialRefObj.YWorldLimits(1);
    end 
end
scatter(CD3_Registered(:,1),CD3_Registered(:,2))
