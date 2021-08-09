tifDir = 'Z:\Juliane\Data\2P_Data\F2501_2021-03-23\t00006\registered\slice1\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\F2501_2021-03-23\t00006\';
level = 0;
field = 'rawF';

%% read in the tiff files
tifStack = [];
tic
tifFiles = dir([tifDir '*.tif']);
for fileNum = 1:length(tifFiles)
    disp(['Reading Imaging Stack ' num2str(fileNum) ' Out Of ' num2str(length(tifFiles))]);
    imgStack = read_Tiffs([tifDir '/' tifFiles(fileNum).name], 1,50);
    tifStack = cat(3,tifStack,imgStack);
end
toc    

%% loadROIs
dim = size(tifStack);
ROIfile = [tifDir 'RoiSet.zip'];
[ROIs, ~] = ROIconvert(ROIfile, [dim(1) dim(2)]);

%% extract F trace
for nr= 1:length(ROIs)
    ROImasks = zeros(size(tifStack,1),size(tifStack,2));
    ROIsize = size(ROIs(nr).body, 1 );
    for p = 1
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
end
clear tifStack
%% load spike2data
metadata.StimParams=LoadStimParams(Sp2Dir);
metadata.TwoPhoton=LoadFrameTimes(Sp2Dir);

%% do quick analysis and plot
analysis = [];
if strcmp(field, 'dff')
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
end
[analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,level, field, 'pre', 1, 'post',2,'windowStart',0, 'windowStop',2);
metadata.StimParams.numCon = 1;
metadata.StimParams.numDirections = metadata.StimParams.uniqStims-1;
metadata.StimParams.directions = linspace(0,360,metadata.StimParams.numDirections+1); %this would give us a all orientations from 0 to 180
metadata.StimParams.directions = metadata.StimParams.directions(1:end-1);
for i = 1:length(data.roi)
    PlotAvgStimResponseOri(metadata, analysis, field, i)
    PlotIndTrialStimResponseOri(metadata, analysis, field', i)
end
