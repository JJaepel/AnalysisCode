analysisParams.animal = 'F2528_2021-05-28';
analysisParams.expID = 't00010';
analysisParams.savedir = 'Z:\Juliane\Data\ImageAnalysis\';
analysisParams.level = 0;

tifDir = ['Z:\Juliane\Data\2P_Data\' analysisParams.animal filesep analysisParams.expID filesep];
Sp2Dir = ['Z:\Juliane\Data\Spike2Data\' analysisParams.animal filesep analysisParams.expID];

latencyAna = 1;
level = 0;
p = 0.05;

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
end
clear tifStack

%% load spike2data
metadata.StimParams=LoadStimParamsRet(Sp2Dir);
metadata.TwoPhoton=LoadFrameTimes(Sp2Dir);

%% calculate analysis settings
%create stimulus settings for each patch
analysis = struct;
stimType = metadata.StimParams.type;
[analysisParams, metadata] = getStimParamsPatches(analysisParams, metadata);
if latencyAna
    [analysisParams, metadata, analysis] = calculateSignificantPatches(analysisParams, metadata, data, analysis);
else
    [analysisParams, metadata, data, analysis] = calcROIPropertiesPatches(analysisParams, metadata, data, analysis);
    switch stimType
        case 'Patch'
            RespMatrixW = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
            RespMatrixB = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
            for roi = 1:length(signROIs)
                RespMatrixW(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgStimResponse(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
                RespMatrixB(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgStimResponse(metadata.StimParams.numElevation * metadata.StimParams.numAzimuth+1:end-1), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
            end
            AvgRespMatrixW = mean(RespMatrixW, 3);
            AvgRespMatrixB = mean(RespMatrixB, 3);
            
            figure
            ax(1) = subplot(1,2,1);
            imagesc(AvgRespMatrixW)
            colormap(ax(1),parula)
            axis ('square')
            axis off
            title('Average Response RF, White','FontSize', 10')

            ax(2) = subplot(1,2,2);
            imagesc(AvgRespMatrixB)
            colormap(ax(2),hot)
            axis ('square')
            axis off
            title('Average Response RF, Black','FontSize', 10')
        otherwise
            RespMatrix = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
            for roi = 1:length(signROIs)
                RespMatrix(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgStimResponse(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
            end
            AvgRespMatrix = mean(RespMatrix, 3);
        
            figure
            imagesc(AvgRespMatrix)
            colormap(hot)
            axis ('square')
            axis off
            title('Average Response RF','FontSize', 10')
    end
end
    
