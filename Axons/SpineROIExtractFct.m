function  data = SpineROIExtractFct(analysisParams)
%% loading experiment data
if analysisParams.server
    drive       = 'Z:\Juliane\';
else
    drive           = 'F:\';
end

baseDirectory   = [drive 'Data\2P_data\'];
filename = 'ROIs.mat';
if analysisParams.level
    ROIDirName         = [baseDirectory analysisParams.animal '\' analysisParams.name '\Registered\combined\Projection\'];
    tifDir             = [baseDirectory analysisParams.animal '\' analysisParams.name '\Registered\combined\'];
else
    ROIDirName         = [baseDirectory analysisParams.animal '\' analysisParams.name '\Registered\slice1\Projection\'];
    tifDir             = [baseDirectory analysisParams.animal '\' analysisParams.name '\Registered\slice1\'];
end
loadFile        = [ROIDirName filename];

%% read in the tiff files
tifStack = [];
tifFiles = dir([tifDir '*.tif']);
for currentFile = 1:length(tifFiles)
    disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);
    
    % Specify stack name
    filePath = [tifDir tifFiles(currentFile).name];

    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
end

%% loadROIs
load(loadFile);
ROIs = data.roi;

%% extract F trace
for nr= 1:length(ROIs)
    ROImasks = zeros(size(tifStack,1),size(tifStack,2));
    ROIsize = size(ROIs(nr).body, 1 );
    for p = 1:ROIsize
        if ROIs(nr).body(p,1)>0 && ROIs(nr).body(p,2)>0
            ROImasks(ROIs(nr).body(p,2), ROIs(nr).body(p,1) ) = 1;
        end
    end
    xVals = find(sum(ROImasks,1)>0);
    yVals = find(sum(ROImasks,2)>0);
    xRange = min(xVals):max(xVals);
    yRange = min(yVals):max(yVals);
    M = repmat(int16(ROImasks(yRange,xRange)), [1 1 size(tifStack,3)] );
    data.roi(nr).rawF = squeeze(squeeze(sum( sum( int16(tifStack(yRange,xRange,:)) .* M, 1 ), 2 ) ./ ROIsize))';
    zeroPixel = find(data.roi(nr).rawF == 0);
    data.roi(nr).rawF(zeroPixel) = 1;
end
clear tifStack
end