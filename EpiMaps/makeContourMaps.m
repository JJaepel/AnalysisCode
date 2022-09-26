function analysis = makeContourMaps(metadata, analysis, analysisParams)

analysisParams.clippingPercentile = 0.2;
analysisParams.mapsWithContours = 0;

%% 1). Load and z-score data, average across conditions
%load response Trace
stimResponseTrace = permute(analysis.(analysisParams.field).roi.stimResponseTrace, [4 5 3 2 1]);

%remove bv and outside window area
tempRespTrace = stimResponseTrace;
for stim=1:size(tempRespTrace,3)
    for trial = 1:size(tempRespTrace,4)
        for frame = 1:size(tempRespTrace,5)
            temp = tempRespTrace(:,:,stim, trial,frame);
            temp(~analysis.maskBV(:)) = NaN;
            tempRespTrace(:,:,stim,trial,frame)=temp;
        end
    end
end

%z scoring data
disp('Z scoring data')
mu = nanmean(tempRespTrace(:));
sd = nanstd(tempRespTrace(:));
zScoreRespTrace = (tempRespTrace - mu)/sd;

%average over the active trace
includedFrames = round(metadata.Imaging.rate * metadata.StimParams.isi/2)+1:round(metadata.Imaging.rate * metadata.StimParams.isi/2)+ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration);
analysis.(analysisParams.field).zScore = mean(zScoreRespTrace(:,:,:,:,includedFrames),5);

%average over trials to get zscored trial-averaged maps
analysis.(analysisParams.field).zScoreAvgMaps = mean(analysis.(analysisParams.field).zScore,4);

%% 2.) Get thresholded and contour maps

%test threshold
analysisParams.threshold = testThreshold(analysis.(analysisParams.field).zScore(:,:,1,1));

%making  contours on trial masks
analysis.(analysisParams.field).contourMaps = zeros(size(zScoreRespTrace,1), size(zScoreRespTrace,2), size(zScoreRespTrace,3), size(zScoreRespTrace,4));
analysis.(analysisParams.field).indStimTrialMaps = zeros(size(zScoreRespTrace,1), size(zScoreRespTrace,2), size(zScoreRespTrace,3), size(zScoreRespTrace,4));
for stim = 1:size(zScoreRespTrace,3)
   for trial=1:size(zScoreRespTrace,4)
       [ContourXData,ContourYData,tempContour, tempMaps] = makeContour(analysis.(analysisParams.field).zScore(:,:,stim,trial), analysisParams.threshold);
       analysis.(analysisParams.field).ContourXData{stim,trial} = ContourXData; clear ContourXData
       analysis.(analysisParams.field).ContourYData{stim,trial} = ContourYData; clear ContourYData
       analysis.(analysisParams.field).contourMaps(:,:,stim, trial) = tempContour(1:size(tempMaps,1),1:size(tempMaps,2)); clear tempContour
       analysis.(analysisParams.field).indStimTrialMaps(:,:,stim,trial) = tempMaps; clear tempMaps
   end
end

%making contours of average masks
analysis.(analysisParams.field).contourAvgMaps = zeros(size(zScoreRespTrace,1), size(zScoreRespTrace,2), size(zScoreRespTrace,3));
analysis.(analysisParams.field).AvgRespMaps = zeros(size(zScoreRespTrace,1), size(zScoreRespTrace,2), size(zScoreRespTrace,3));
for stim = 1:size(zScoreRespTrace,3)
    [ContourXData,ContourYData,tempContour, tempMaps] = makeContour(analysis.(analysisParams.field).zScoreAvgMaps(:,:,stim), analysisParams.threshold);
    analysis.(analysisParams.field).ContourAvgXData{stim} = ContourXData; clear ContourXData
    analysis.(analysisParams.field).ContourAvgYData{stim} = ContourYData; clear ContourYData
    analysis.(analysisParams.field).contourAvgMaps(:,:,stim) =  tempContour(1:size(tempMaps,1),1:size(tempMaps,2));clear tempContour
    analysis.(analysisParams.field).AvgRespMaps(:,:,stim) = tempMaps; clear tempMaps
end 

%% 3.) Calculate overlap, etc.

%calculating overlap
analysis = calcMapVariability(analysisParams,analysis);

%% 4.) Plotting

% set general parameters
analysisParams.conC = cbrewer('qual', 'Paired',metadata.StimParams.numTrials);
analysisParams.cocV1 = cbrewer('seq', 'RdPu',5);
analysisParams.cocA19 = cbrewer('seq', 'PuBuGn',5);
if metadata.StimParams.numTrials > 6
    if metadata.StimParams.numTrials > 10
        analysisParams.nRows = 5; analysisParams.nCol = 5;
        analysisParams.plotNrTrials = [1 2 3 6 7 8 11 12 13 16 17 18 21 22 23];
        analysisParams.plotNrAll= [4:5 9:10];
        analysisParams.plotNrAverage = 14; 
        analysisParams.plotNrQuant = [15 19 20]; 
    else
        analysisParams.nRows = 4; analysisParams.nCol = 5;
        analysisParams.plotNrTrials = [1 2 3 6 7 8 11 12 13 16 17 18];
        analysisParams.plotNrAll = [4:5 9:10];
        analysisParams.plotNrAverage = 14; 
        analysisParams.plotNrQuant = [15 19 20]; 
    end
else
    analysisParams.nRows = 4; analysisParams.nCol = 4;
    analysisParams.plotNrTrials = [1 2 5 6 9 10];
    analysisParams.plotNrAll = [3:4 7:8];
    analysisParams.plotNrAverage = 11;
    analysisParams.plotNrQuant = [12 15 16];
end

% if analysisParams.mapsWithContours
%     disp('Drawing countours onto z-scored maps for individual trials') 
%     for type = 1:2
%        plotSingleTrialContours(analysisParams, analysis, metadata, type)
%     end
% end

disp('Plotting variability data for individual stims')
plotAllContoursSingleTrial(analysisParams, analysis, metadata)

disp('Plotting variability data summary')
plotAllContoursSummary(analysisParams, analysis, metadata)
