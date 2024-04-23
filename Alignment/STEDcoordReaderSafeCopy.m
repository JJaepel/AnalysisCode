clear all

bufferInUm = 4;
pixelperUm = 10;

%Step 1: read in data from excel file
coordDataFile = 'C:\Users\jaepelj\Desktop\STED_Coordinates_Test.xlsx';
sheets = sheetnames(coordDataFile);
currentSheet = '2022-9-9';
[~, xls_txt, xls_all]=xlsread(coordDataFile, currentSheet);

coordInfo = struct;
fileNameCol = find(contains(xls_txt(1,:),'File'),1);
xPosCol = find(contains(xls_txt(1,:),'xPos'),1);
yPosCol = find(contains(xls_txt(1,:),'yPos'),1);
AngleCol = find(contains(xls_txt(1,:),'Angle'),1);
SizeXCol = find(contains(xls_txt(1,:),'SizeX'),1);
SizeYCol = find(contains(xls_txt(1,:),'SizeY'),1);
PixelXCol = find(contains(xls_txt(1,:),'PixelSizeX'),1);
PixelYCol = find(contains(xls_txt(1,:),'PixelSizeY'),1); 

k = 1;
for i = 2:size(xls_txt,1)
    coordInfo(k).file = xls_all(i,fileNameCol);
    coordInfo(k).xPos = cell2mat(xls_all(i,xPosCol));
    coordInfo(k).yPos = cell2mat(xls_all(i,yPosCol));
    coordInfo(k).Angle = cell2mat(xls_all(i,AngleCol));
    coordInfo(k).SizeX = cell2mat(xls_all(i,SizeXCol));
    coordInfo(k).SizeY = cell2mat(xls_all(i,SizeYCol));
    coordInfo(k).PixelX = cell2mat(xls_all(i,PixelXCol));
    coordInfo(k).PixelY = cell2mat(xls_all(i,PixelYCol));
    k = k+1;
end

%Step 2: Read in tif data, scale them to fit and rotate
tifFolder = ['Y:\Data_JJ_Ferret_project\JJ_FerretProject_Organized\F2691\' char(currentSheet) '\Copies\'];
for ind= 1:size(coordInfo,2)
    coordInfo(ind).img = imread([char(tifFolder) char(coordInfo(ind).file) '.tif']);
    coordInfo(ind).rotImg = imrotate(coordInfo(ind).img, coordInfo(ind).Angle);
    %add scaled image of 1 pixel = 0.04 um = 25 pixel per um
    scale = coordInfo(ind).SizeX*coordInfo(ind).PixelX*pixelperUm/size(coordInfo(ind).img,2);
    coordInfo(ind).ImageScaled = imresize(coordInfo(ind).rotImg,scale);
    %rotate image using angle in file
    
end

%Step 3: Find the position of the left upper corner, position of image is
%given as the center - those are um positions!
for i =1:length(coordInfo)
    %calc size of each image)
    %coordInfo(i).SizeXum = coordInfo(i).PixelX*coordInfo(i).SizeX;
    %coordInfo(i).SizeYum = coordInfo(i).PixelY*coordInfo(i).SizeY;
    
    %calc size of each image)
    coordInfo(i).SizeXum = 1/pixelperUm*size(coordInfo(i).ImageScaled,2);
    coordInfo(i).SizeYum = 1/pixelperUm*size(coordInfo(i).ImageScaled,1);

    coordInfo(i).AbsPosLeftX = coordInfo(i).xPos + coordInfo(i).SizeXum/2;
    coordInfo(i).AbsPosRightX = coordInfo(i).AbsPosLeftX - coordInfo(i).SizeXum;
    coordInfo(i).AbsPosUpY = coordInfo(i).yPos - coordInfo(i).SizeYum/2;
    coordInfo(i).AbsPosDownY = coordInfo(i).AbsPosUpY  + coordInfo(i).SizeYum;

end

%Step 4: Get the 0/0 position and write everything into relative positions
%Left upper corner of the resulting image should be maxX, MinY

startPosLeft = max([coordInfo.AbsPosLeftX]);
startPosUp = min([coordInfo.AbsPosUpY]);
endPosRight = min([coordInfo.AbsPosRightX]);
endPosDown = max([coordInfo.AbsPosDownY]);

for i =1:length(coordInfo)
    coordInfo(i).RelPosLeftX = startPosLeft - coordInfo(i).AbsPosLeftX;
    coordInfo(i).RelPosRightX = startPosLeft - coordInfo(i).AbsPosRightX;
    coordInfo(i).RelPosUpY = coordInfo(i).AbsPosUpY - startPosUp;
    coordInfo(i).RelPosDownY = coordInfo(i).AbsPosDownY - startPosUp;
end

%Step 5: Get pixel coordinates
%define size of image in total - one z per image minus cell image
sizeX = ceil((abs(startPosLeft-endPosRight))*pixelperUm)+1; %size in um -> lets have 25 pixels per um
sizeY = ceil((abs(startPosUp-endPosDown))*pixelperUm)+1;

allImages = zeros(sizeY, sizeX, length(coordInfo));

for d = 1:length(coordInfo)
    %coordinates width
    coordInfo(d).pixelPosLeftX = ceil(coordInfo(d).RelPosLeftX*pixelperUm);
    if coordInfo(d).pixelPosLeftX == 0
        coordInfo(d).pixelPosLeftX = coordInfo(d).pixelPosLeftX+1;
    end
    coordInfo(d).pixelPosRightX = coordInfo(d).pixelPosLeftX+size(coordInfo(d).ImageScaled,2)-1;
    
    %coordinates height
    coordInfo(d).pixelPosUpY = ceil(coordInfo(d).RelPosUpY*pixelperUm);
    if coordInfo(d).pixelPosUpY == 0
        coordInfo(d).pixelPosUpY = coordInfo(d).pixelPosUpY+1;
    end
    coordInfo(d).pixelPosDownY = coordInfo(d).pixelPosUpY+size(coordInfo(d).ImageScaled,1)-1;
    
    allImages(coordInfo(d).pixelPosUpY :coordInfo(d).pixelPosDownY,coordInfo(d).pixelPosLeftX:coordInfo(d).pixelPosRightX,d) = coordInfo(d).ImageScaled;
end



figure
dend = [1:14];
for i=1:length(dend)
    subplot(4,4,i)
    imagesc(allImages(:,:,dend(i)))
    axis image
    axis off
end

b = mean(allImages(:,:,1:end),3);
figure
imagesc(b)
axis image