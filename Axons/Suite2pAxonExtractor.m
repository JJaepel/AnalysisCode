function Suite2pAxonExtractor
%% specify the experiment & parameters

name            = 'F2444_2020-07-09';
folderNumber    = 'red';
server          =  0;
level           =  0;

if server
    drive       = 'Z:\Juliane\';
else
    drive           = 'F:\';
end

baseDirectory   = [drive 'Data\2P_data\'];
filename        = 'Projection.tif';
if level
    dirName         = [baseDirectory name '\' folderNumber '\suite2p\combined\reg_tif\'];
    saveDir         = [baseDirectory name '\' folderNumber '\suite2p\combined\'];
else
    dirName         = [baseDirectory name '\' folderNumber '\suite2p\plane0\reg_tif\'];
    saveDir         = [baseDirectory name '\' folderNumber '\suite2p\plane0\'];
end
saveFile        = [saveDir filename];

%% read in the tiff files
tifStack = [];
tifFiles = dir([dirName '*.tif']);
for currentFile = 1:length(tifFiles)
    disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);
    
    % Specify stack name
    filePath = [dirName tifFiles(currentFile).name];

    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
end


%% make projection for axon ROIing
meanImg = uint16(mean(tifStack,3));
imwrite(meanImg, saveFile, 'tiff', 'writemode', 'overwrite', 'compression', 'none')

%% semi-automatic ROI selection 
cd(saveDir)
semi_manual_axon_ROIs(saveDir, filename)

%% convert RoiSet to mat file
dim = size(meanImg);
[ROIs, ~] = ROIconvert('RoiSet.zip', [dim(1) dim(2)]);

%% extract ROI positions
data = [];
for m = 1:length(ROIs)    
    data.roi(m).xPos = ROIs(m).x;
    data.roi(m).yPos = ROIs(m).y;
    data.roi(m).mask = [ROIs(m).perimeter];
    data.roi(m).name = m;
end

%% extract F trace
ROImasks = cell(1,length(ROIs));
for nr= 1:length(ROIs)
    ROImasks{nr} = zeros(size(tifStack,1),size(tifStack,2));
    ROIsize(nr) = size(ROIs(nr).body, 1 );
    for p = 1:ROIsize(nr)
        if ROIs(nr).body(p,1)>0 && ROIs(nr).body(p,2)>0
            ROImasks{nr}(ROIs(nr).body(p,2), ROIs(nr).body(p,1) ) = 1;
        end
    end
    xVals = find(sum(ROImasks{nr},1)>0);
    yVals = find(sum(ROImasks{nr},2)>0);
    xRange = min(xVals):max(xVals);
    yRange = min(yVals):max(yVals);
    M = repmat(int16(ROImasks{nr}(yRange,xRange)), [1 1 size(tifStack,3)] );
    data.roi(nr).rawF = squeeze(squeeze(sum( sum( int16(tifStack(yRange,xRange,:)) .* M, 1 ), 2 ) ./ ROIsize(nr)))';
end

%% save meanImg as template
data.template = meanImg;
clear meanImg
clear tifStack
end