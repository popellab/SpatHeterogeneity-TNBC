%-------This script is used to register cell coordinate based on transformation matrices obtained from local registraiton step-------%

%% 
%Load FoxP3 points
load('./Case_1_regInfo/FoxP3_SecInfo/Point_All_FoxP3.mat');

%Load FoxP3 points
load('./Case_1_regInfo/FoxP3_SecInfo/FoxP3_Transmatrix.mat')

% Define regions (should in accord with the definition in local registration)

FoxP3_Pts = Point_All;
Reg1 = [1,1670;1,850];
Reg2 = [1,460;850,2313];
Reg3 = [460,1670;850,2313];
Reg4 = [1670,3390;1,1100];
Reg5 = [1670,3390;1100,2313];

Section = struct('Reg1',Reg1, 'Reg2', Reg2,'Reg3',Reg3, 'Reg4',Reg4,'Reg5',Reg5);
Section = struct2cell(Section);
% Assign coordinates to sections
[len, wid] = size(FoxP3_Pts);
Char = zeros(len,1);
[length, dim] = size(FoxP3_Pts);
FoxP3_Registered = zeros(length, dim);
string = 'Case1_FoxP3_Sec';
for i = 1:length
    for j=1:5
        array = sprintf('%s%d',string,j);
        Var = eval(array);
        if FoxP3_Pts(i,1) < Section{j,1}(1,2) && FoxP3_Pts(i,1) >= Section{j,1}(1,1) && FoxP3_Pts(i,2) < Section{j,1}(2,2) && FoxP3_Pts(i,2) >= Section{j,1}(2,1)
            [FoxP3_Registered(i,1),FoxP3_Registered(i,2)] = transformPointsForward(Var.Transformation,FoxP3_Pts(i,1),FoxP3_Pts(i,2));
            %FoxP3_Pts(i,1) = FoxP3_Pts(i,1) - Var.SpatialRefObj.XWorldLimits(1);
            %FoxP3_Pts(i,2) = FoxP3_Pts(i,2) - Var.SpatialRefObj.YWorldLimits(1);
            Char(i,1) = j;
        end
    end 
end
%% Load CD8 points
load('./Case_1_regInfo/CD8_SecInfo/Point_All_CD8.mat');
load('./Case_1_regInfo/CD8_SecInfo/CD8_Transmatrix.mat');

CD8_Pts = Point_All;
[length, dim] = size(CD8_Pts);
CD8_Registered = zeros(length, dim);

for i = 1:length
    [CD8_Registered(i,1),CD8_Registered(i,2)] = transformPointsForward(Case1_CD8_Sec.Transformation,CD8_Pts(i,1),CD8_Pts(i,2));
end
%% Load CD20 points
load('./Case_1_regInfo/CD20_SecInfo/Point_All_CD20.mat');
load('./Case_1_regInfo/CD20_SecInfo/CD20_Transmatrix.mat');
CD20_Pts = Point_All;

Reg1 = [1,3390;1,2440];


Section = struct('Reg1',Reg1);
Section = struct2cell(Section);
% Assign coordinates to sections

[length, dim] = size(CD20_Pts);
CD20_Registered = zeros(length, dim);

string = 'Case1_CD20_Sec';
CD20_Pts(:,1) = 3390 - CD20_Pts(:,1);
CD20_Pts(:,2) = 2440 - CD20_Pts(:,2);
%plot(CD20_Registered(:,1),CD20_Registered(:,2))
for i = 1:length 
        [CD20_Registered(i,1),CD20_Registered(i,2)] = transformPointsForward(Case1_CD20_Sec.Transformation,CD20_Pts(i,1),CD20_Pts(i,2));
end
 %% CD3 registration
load('./Case_1_regInfo/CD3_SecInfo/Point_All_CD3.mat');
load('./Case_1_regInfo/CD3_SecInfo/CD3_Transmatrix.mat');

CD3_Pts = Point_All;
Reg1 = [1,1072;1,1041];
Reg2 = [1,1072;1041,2271];
Reg3 = [1072,3627;1,1041];
Reg4 = [1072,3627;1041,2271];


Section = struct('Reg1',Reg1, 'Reg2', Reg2,'Reg3',Reg3, 'Reg4',Reg4);
Section = struct2cell(Section);
% Assign coordinates to sections
Char = zeros(284158,1);
[length, dim] = size(CD3_Pts);
CD3_Registered = zeros(length, dim);
string = 'Case1_CD3_Sec';
for i = 1:length
    for j=1:4
        array = sprintf('%s%d',string,j);
        Var = eval(array);
        if CD3_Pts(i,1) < Section{j,1}(1,2) && CD3_Pts(i,1) >= Section{j,1}(1,1) && CD3_Pts(i,2) < Section{j,1}(2,2) && CD3_Pts(i,2) >= Section{j,1}(2,1)
            [CD3_Registered(i,1),CD3_Registered(i,2)] = transformPointsForward(Var.Transformation,CD3_Pts(i,1),CD3_Pts(i,2));
            %FoxP3_Pts(i,1) = FoxP3_Pts(i,1) - Var.SpatialRefObj.XWorldLimits(1);
            %FoxP3_Pts(i,2) = FoxP3_Pts(i,2) - Var.SpatialRefObj.YWorldLimits(1);
            Char(i,1) = j;
        end
    end 
end