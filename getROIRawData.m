function [ce, meanImg] = getROIRawData(info)

%Extracts the ROI raw data from the tif files when the zip files for each
%ROI are in the correct folder
%
%Input:
% - info: contains all the information about the current event
%
%Outputs:
% - ce: structure containing in each row information about each individual
% ROI, for now it is anatomical information as well as the extracted ROI
% raw trace
% - meanImg: mean Image of the imaged area, used for example for overlays
% of ROIs with it to show position etc.

%% Step 1: Load ROIs from .zip file
ROIDir = [info.TwoPDir '\slice1\Projection\']; %wher is the file
ROIFile = [ROIDir 'RoiSet.zip'];
    
[sROI] = ReadImageJROI(ROIFile); %read the file

%% Step 2: Read in additional information for each file
    
numROIs = size(sROI,2); %total number of ROIs
ROISegNum = zeros(numROIs,1); %keep track of dendritic segment number
isDendrite = zeros(numROIs,1); %keep track of dendritic ROIs
currentDendrite = 1;
%is it a dendrite? what is the correct dendrite to this ROI?
for n = 1:numROIs
    ROISegNum(n) = currentDendrite;
    if strcmp(sROI{n}.strType,'PolyLine') %strType for straightening
        isDendrite(n) = 1;
        %assumes last ROI in each segment is the dendrite
        currentDendrite = currentDendrite + 1;
    end
end

%need to add some information: ideally from analysisParams
denType = info.denType; % basal, apical, oblique
expNameSplit = strsplit(info.animal, '_');
date2p = expNameSplit{1}; 
name = info.name;
zoom = info.zoom; %zoom
scale = 512/(1000/zoom);
depth = info.depth;

%% Step 3: Load the tifFiles
tifStack = [];
tifFiles = dir([info.TwoPDir '\slice1\*.tif']); %find all files
for currentFile = 1:length(tifFiles)  
    % Specify stack name
    filePath = [info.TwoPDir '\slice1\' tifFiles(currentFile).name];

    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
end
%make the meanImg
meanImg = uint16(mean(tifStack(:,:,1:1000),3));
meanImg = double(meanImg);

%% Step 4: Go through the ROIs and add all the ROI info to the different ROIs
disp('going through ROIs')
for i = 1:numROIs
    %save ROI and associated information
    ce(i).yPos =  median(sROI{i}.mnCoordinates(:,2));
    ce(i).xPos =  median(sROI{i}.mnCoordinates(:,1));
    ce(i).depth = depth;
    ce(i).mask = sROI{i}.mnCoordinates;
    ce(i).date = date2p;
    ce(i).file = name;

    ce(i).dendrite = isDendrite(i);
    ce(i).spine = ~isDendrite(i);
    ce(i).segment = ROISegNum(i);

    ce(i).denType = denType;
    ce(i).scale = scale; %pix/um

    if isDendrite(i)
        ce(i).img = meanImg;
    end

    Mask = poly2mask(ce(i).mask(:,1),ce(i).mask(:,2),size(tifStack,1),size(tifStack,2));

    % fast way: only area around ROI
    xVals = find(sum(Mask,1)>0);
    yVals = find(sum(Mask,2)>0);
    xRange = min(xVals):max(xVals);
    yRange = min(yVals):max(yVals);
    M = repmat(int16(Mask(yRange,xRange)), [1 1 size(tifStack,3)] );

    raw = squeeze(sum(sum(int16(tifStack(yRange,xRange,:)).*M,1),2)./sum(sum(M(:,:,1))));
    if size(raw,2)>1
        raw = raw(:,1);
    end
    ce(i).raw = raw;
end

%% Step 5: Save the meanImg and the extracted ROI file
save([info.saveDir 'ROIs.mat'], 'ce', '-mat') 
save([info.saveDir 'meanImg.mat'], 'ce', 'meanImg', '-mat') 