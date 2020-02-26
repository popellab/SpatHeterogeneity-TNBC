% This script is for local registration -- register the image part by part
%%
Case1_CD4 = imread('Case1_CD4.png');
Case1_CD3 = imread('Case1_CD3.png');

CD3_gray = rgb2gray(Case1_CD3);
CD4_gray = rgb2gray(Case1_CD4);
%%
[m,n] = size(CD3_gray);
[p,q] = size(CD4_gray);
CD3_reg = 255*ones(m,n);
CD4_reg = 255*ones(p,q);

% The ranges of 'i' depends on how you divide the WSI into subregions.
for i = 1400:2271   
    for j = 2200:3627            
        CD3_reg(i,j) = CD3_gray(i,j);
    end
end
CD3_reg = uint8(CD3_reg);
for a = 1300:2603
    for b = 2000:3509
        CD4_reg(a,b) = CD4_gray(a,b);
    end
end
CD4_reg = uint8(CD4_reg);
figure;
imshow(CD3_gray)
figure;
imshow(CD4_reg)

imwrite(CD3_reg,'CD3_reg.jpg')
imwrite(CD4_reg,'CD4_reg.jpg')
%%
imshow(Foxp3_8th.RegisteredImage)
%%
Case1_CD3 = imread('Case_1/Case1_CD3.png');
CD3_gray = rgb2gray(Case1_CD3);

Reg = imwarp(CD3_gray,movingReg1.Transformation,'OutputView',movingReg1.SpatialRefObj);
for i = 1:2603
    for j = 1:3509
        if Reg(i,j) == 0
            Reg(i,j) = 243;
        end
    end
end
figure;
imshow(Reg)
imwrite(Reg,'CD3_Case1_reg.jpg')

figure;

imshow(Int_FoxP3)
hold on;
plot(FoxP3_Pts(:,1),FoxP3_Pts(:,2),'o')
%%
Mov_image = movingReg.RegisteredImage;
[row, col] = size(Mov_image);
    for ncol = 1: col
        for nrow = 1: row
            if Mov_image(nrow, ncol) == 0
                Mov_image(nrow, ncol) = 243;
            end
        end
    end
imshow(Mov_image)
imwrite(Mov_image,'global_reg.jpg')
