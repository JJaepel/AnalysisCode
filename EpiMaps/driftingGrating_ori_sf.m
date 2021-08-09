animal = 'F2425_2020-03-05';
expt_id = 8;
sp2id = expt_id;

close all
EpiDir = 'F:\Data\Epi\';
Sp2Dir = 'F:\Data\Spike2Data\';
SaveDir = 'F:\Data\ImageAnalysis\';

windowStop=2;
windowStart=0;
pre=1;
field = 'rawF';


EpiDirectory = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id) filesep];
if sp2id > 9
    Sp2dDirectory = [Sp2Dir animal filesep 't000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't000' num2str(expt_id) filesep];
else
    Sp2dDirectory = [Sp2Dir animal filesep 't0000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't0000' num2str(expt_id) filesep];
end

if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end

sliceparams = struct;
sliceparams.expt_id = expt_id;
sliceparams.baseDirectory = EpiDirectory;

%% load metadata
metadata.StimParams=LoadStimParams(Sp2dDirectory);
metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
metadata.StimParams.path=fullfile(Sp2dDirectory);
metadata.StimParams.series=expt_id;

%% load tiffs
data.rawF = readingImagingData(EpiDirectory);

%% create stimCodes

numberOfConditions = metadata.StimParams.numberOfStims;
stimStartIndex = zeros(numberOfConditions,1,'double'); 
stimStopIndex  = zeros(numberOfConditions,1,'double');
for i=1:numberOfConditions
    stimStartIndex(i) = find(metadata.Imaging.time>=metadata.StimParams.StimOnTimes(2,i),1,'first');
    stimStopIndex(i) = find(metadata.Imaging.time>=metadata.StimParams.StimOnTimes(2,i)+metadata.StimParams.stimDuration,1,'first');
end
metadata.StimParams.stimStartIndex = stimStartIndex;
metadata.StimParams.stimStopIndex = stimStopIndex;

%% chop traces
analysis = struct;
disp('Chopping Traces')
[analysis, metadata] = ChopStimulusTraceEpi(analysis,metadata,data,field);
analysis.rawFMeanImg = data.rawFMeanImg;
analysis.ROI = data.ROI;
clear data

%% make timecourse

eval(sprintf('cmap = %s(%d);','hsv',metadata.StimParams.uniqStims-1));
LUT = cat(1,cmap, [0 0 0]);
numberOfFrames     = size(analysis.(field).roi.stimResponseTrace,1);

for condition = 1:metadata.StimParams.uniqStims
    t = ((0:(numberOfFrames-1))/metadata.Imaging.rate)-pre;
    y = nanmedian(analysis.(field).roi.stimResponseTrace(:,:,condition,analysis.ROI(:)),4); % Takes median response of the pixel values. Reduces noise relative to mean
    yMean     = squeeze(mean(y,2));
    yStdError = std(y,[],2)/sqrt(size(y,2));
    shadedErrorBar(t,yMean,yStdError,'lineProps',{'-','LineWidth',3,'Color',LUT(condition,:)},'patchSaturation',0.4); hold on;
end
xlabel('Time (s)');
ylabel('Response amplitude');
axis square;
set(gca,'Box','off');
set(gcf, 'color', 'w');

% Add stimulus box
stimStart = 0;
stimStop  = metadata.StimParams.stimDuration;
yLimits = get(gca,'YLim');
rectangle('Position',[stimStart yLimits(2) stimStop-stimStart 0.025*range(yLimits)],'FaceColor','k')
saveas(gcf, fullfile(saveDirectory, 'Timecourse.png'))
analysis.(field).roi.stimResponseTrace = permute(analysis.(field).roi.stimResponseTrace, [4 5 3 2 1]);

%% map analysis
includedFrames = [round(metadata.Imaging.rate * metadata.StimParams.isi/2)+1:round(metadata.Imaging.rate * metadata.StimParams.isi/2)+ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration)];
stimResponseTrace = mean(analysis.(field).roi.stimResponseTrace(:,:,1:(end-1),:,includedFrames),5);

trialAveragedMaps = squeeze(mean(stimResponseTrace,4));
trialAveragedMaps(isnan(trialAveragedMaps(:))) = 0;

numberOri = metadata.StimParams.numOrientations;
numberSf = (length(metadata.StimParams.uniqStimIds)-1)/metadata.StimParams.numOrientations;
sfAveragedMaps = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberSf);
directionSfMap = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberSf);
orientationSfMap = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberSf);
directionAveragedMaps = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberOri);
orientationAveragedMaps = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberOri/2);
for sf= 1:numberSf
    sfAveragedMaps(:,:,sf) = mean(trialAveragedMaps(:,:,1+numberOri*(sf-1):numberOri*sf),3);
    directionSfMap(:,:,sf)   = vectorSum(trialAveragedMaps(:,:,1+numberOri*(sf-1):numberOri*sf),2,3);
    orientationSfMap(:,:,sf)   = vectorSum(trialAveragedMaps(:,:,1+numberOri*(sf-1):numberOri*sf),2,3);
end
for dir= 1:numberOri
    directionAveragedMaps(:,:,dir) = mean(trialAveragedMaps(:,:,dir:numberOri:end),3);
end
for ori= 1:numberOri/2
    orientationAveragedMaps(:,:,ori) = (mean(trialAveragedMaps(:,:,ori+numberOri/2:numberOri:end),3)+mean(trialAveragedMaps(:,:,ori+numberOri/2:numberOri:end),3))/2;
end

analysis.(field).roi.sFMap   = vectorSum(sfAveragedMaps,1,3);
analysis.(field).roi.directionMap   = vectorSum(directionAveragedMaps, 1, 3);
analysis.(field).roi.orientationMap   = vectorSum(directionAveragedMaps, 2, 3);
analysis.(field).roi.orientationSfMap = orientationSfMap;
analysis.(field).roi.directionSfMap = directionSfMap;              

mapType = {'Orientation','Direction','SpatialFreq'};
showEpiRespAvgCon(analysis, metadata,trialAveragedMaps, orientationAveragedMaps, directionAveragedMaps, sfAveragedMaps, field, mapType,saveDirectory)

h = makeFigureFullScreen(figure);
k = makeFigureFullScreen(figure);
l = makeFigureFullScreen(figure);
set(h,'Name','Direction maps for different spatial frequencies')
set(l,'Name','Orientation magnitude maps for different spatial frequencies')
nStims = size(sfAveragedMaps, 3);
nRows  = ceil(sqrt(nStims));
nCols  = nRows;
clippingPercentile = 0.2;
clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]);
spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
for i = 1:nStims
    figure(h); subplot(nRows,nCols,i);
    dirMapSf = analysis.(field).roi.directionSfMap(:,:,i);
    imagesc(polarMapEpi(dirMapSf, [0 clipValue(2)]));
    axis image; axis off;
    caxis([0 360]);
    cb=colorbar; colormap(cb, hsv);
    title([spatFreq_cell{i}, ' cpd'])
    figure(k); subplot(nRows,nCols,i);
    oriMapSf = analysis.(field).roi.orientationSfMap(:,:,i);
    imagesc(polarMapEpi(oriMapSf, clipValue));
    axis image; axis off;
    caxis([0 180]);
    cb=colorbar; colormap(cb, hsv);
    title([spatFreq_cell{i}, ' cpd'])
    figure(l); subplot(nRows,nCols,i);
    imagesc(magMap(oriMapSf, clipValue));
    colormap('gray');
    axis image; axis off;
    caxis(clipValue);
    title([spatFreq_cell{i}, ' cpd'])
end
saveas(h, fullfile(saveDirectory, 'DirectionMapsSeparatedbySF.png'))
saveas(k, fullfile(saveDirectory, 'OrientationMapsSeparatedbySF.png'))
saveas(l, fullfile(saveDirectory, 'OrientationMagnitudeMapsSeparatedbySF.png'))

%% make HLS maps
analysis.(field).roi.sfMapHLS = stimulusMap(sfAveragedMaps, clipValue);
figure
imshow(analysis.(field).roi.sfMapHLS);
saveas(gcf, fullfile(saveDirectory, 'HLS_SF.png'))

analysis.(field).roi.orientationMapHLS = stimulusMap(orientationAveragedMaps, clipValue);
figure
imshow(analysis.(field).roi.orientationMapHLS);
saveas(gcf, fullfile(saveDirectory, 'HLS_Orientation.png'))

analysis.(field).roi.orientationMapHLS = stimulusMap(directionAveragedMaps, clipValue);
figure
imshow(analysis.(field).roi.orientationMapHLS);
saveas(gcf, fullfile(saveDirectory, 'HLS_Direction.png'))