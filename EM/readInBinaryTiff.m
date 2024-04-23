
%load the spine file
path = 'Z:\EM-DATA\Juliane\DFJ-11\Binary tiffs for script testing for Juliane';
spineFile = 'binary SPINES for psd neg1-5 kept stack.tif';
file = [path filesep spineFile];
fileID = fopen(file);
a = fread(fileID, "uint8");

%convert it and reshape it
a = logical(a);

imageInfo = imfinfo(file);
width = imageInfo.Width;
height = imageInfo.Height;
imageDesc = imageInfo.ImageDescription;
parts = strsplit(imageDesc, '\n');
slices = str2double(extractAfter(parts{2},'='));

fileSize = width*height*slices;
dataStart = size(a,1)-fileSize;

spineData = reshape(a(dataStart+1:end),[width,height,slices]);
clear a

%now it is still on the side, need to turn it and flip it
spineData = rot90(spineData,-1);
spineData = fliplr(spineData);

%now also load the psd file
psdFile = 'psd neg1-5.tiff kept stack.tif';
file = [path filesep psdFile];
fileID = fopen(file);
p = fread(fileID, "uint8");

p = logical(p);
dataStart = size(p,1)-fileSize;
psdData = reshape(p(dataStart+1:end),[width,height,slices]);
clear p

psdData = rot90(psdData,-1);
psdData = fliplr(psdData);

%find the objects within the data
cc = bwconncomp(spineData, 26);
objectProps = regionprops3(cc, 'Volume', 'Centroid', 'BoundingBox', 'Image');
objectProps = table2struct(objectProps);

outputFolder = [path filesep 'SpineImages'];
mkdir(outputFolder);

%now get the coordinates for each object - check if x and y are correct
for s = 1:length(objectProps)
   yMin = ceil(objectProps(s).BoundingBox(1));
   yMax = objectProps(s).BoundingBox(4)+yMin-2;
   
   xMin = ceil(objectProps(s).BoundingBox(2));
   xMax = objectProps(s).BoundingBox(5)+xMin-2;
   
   zMin = ceil(objectProps(s).BoundingBox(3));
   zMax = objectProps(s).BoundingBox(6) + zMin-2;
   
   %make a 3D object and write it as tif
   newImage = spineData(xMin:xMax,yMin:yMax,zMin:zMax);
   
   newFileName = ['Spine' num2str(s, '%02d') '.tif']; %_x ' num2str(xMin) ':' num2str(xMax) ', y ' num2str(yMin) ':' num2str(yMax), ' z ' num2str(zMin) ':' num2str(zMax) '.tif'];
   fullFile = [outputFolder filesep newFileName];
   
   imwrite(newImage(:,:,1), fullFile);
   for zSlice = 2:objectProps(s).BoundingBox(6)-1
       imwrite(newImage(:,:,zSlice),fullFile,'WriteMode','append')
   end
   
   %now write the psdData
   newImage = psdData(xMin:xMax,yMin:yMax,zMin:zMax);
   
   newFileName = ['PSD' num2str(s, '%02d') '.tif']; %_x ' num2str(xMin) ':' num2str(xMax) ', y ' num2str(yMin) ':' num2str(yMax), ' z ' num2str(zMin) ':' num2str(zMax) '.tif'];
   fullFile = [outputFolder filesep newFileName];
   
   imwrite(newImage(:,:,1), fullFile);
   for zSlice = 2:objectProps(s).BoundingBox(6)-1
       imwrite(newImage(:,:,zSlice),fullFile,'WriteMode','append')
   end
end

%lastly, write a table with the values
spineTable = [];
for s = 1:length(objectProps)
    tempTable = table({num2str(s, '%02d')}, {num2str(ceil(objectProps(s).BoundingBox(2)))}, {num2str(ceil(objectProps(s).BoundingBox(2))+objectProps(s).BoundingBox(5))}, ...
        {num2str(ceil(objectProps(s).BoundingBox(1)))}, {num2str(ceil(objectProps(s).BoundingBox(1))+objectProps(s).BoundingBox(4))},...
        {num2str(ceil(objectProps(s).BoundingBox(3)))}, {num2str(ceil(objectProps(s).BoundingBox(3))+objectProps(s).BoundingBox(6))});
    spineTable = [spineTable; tempTable];
end

spineTable.Properties.VariableNames = ["Spine Number", "xMin", "xMax", "yMin", "yMax", "zMin", "zMax"];
writetable(spineTable,[outputFolder filesep 'SpineData.xls'])

