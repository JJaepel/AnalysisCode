close all
clear all

animal = 'F2656_2022-03-25';
expt_id = 12;
sp2id = 11;

close all
EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
SaveDir = 'Z:\Juliane\Data\ImageAnalysis\';

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
data.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 
data.rawFMeanImg = mean(data.rawF,3);
data.baseImg = mean(data.rawF(:,:,1:50),3);
data.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);

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
numberTf = (length(metadata.StimParams.uniqStimIds)-1)/metadata.StimParams.numOrientations;
tfAveragedMaps = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberTf);
directionTfMap = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberTf);
orientationTfMap = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberTf);
directionAveragedMaps = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberOri);
orientationAveragedMaps = zeros(size(trialAveragedMaps,1), size(trialAveragedMaps, 2), numberOri/2);
for tf= 1:numberTf
    tfAveragedMaps(:,:,tf) = mean(trialAveragedMaps(:,:,1+numberOri*(tf-1):numberOri*tf),3);
    directionTfMap(:,:,tf)   = vectorSum(trialAveragedMaps(:,:,1+numberOri*(tf-1):numberOri*tf),2,3);
    orientationTfMap(:,:,tf)   = vectorSum(trialAveragedMaps(:,:,1+numberOri*(tf-1):numberOri*tf),2,3);
end
for dir= 1:numberOri
    directionAveragedMaps(:,:,dir) = mean(trialAveragedMaps(:,:,dir:numberOri:end),3);
end
for ori= 1:numberOri/2
    orientationAveragedMaps(:,:,ori) = (mean(trialAveragedMaps(:,:,ori+numberOri/2:numberOri:end),3)+mean(trialAveragedMaps(:,:,ori+numberOri/2:numberOri:end),3))/2;
end

analysis.(field).roi.tFMap   = vectorSum(tfAveragedMaps,1,3);
analysis.(field).roi.directionMap   = vectorSum(directionAveragedMaps, 1, 3);
analysis.(field).roi.orientationMap   = vectorSum(directionAveragedMaps, 2, 3);
analysis.(field).roi.orientationTfMap = orientationTfMap;
analysis.(field).roi.directionTfMap = directionTfMap;              

mapType = {'Orientation','Direction','TemporalFreq'};
showEpiRespAvgCon(analysis, metadata,trialAveragedMaps, orientationAveragedMaps, directionAveragedMaps, tfAveragedMaps, field, mapType,saveDirectory)

h = makeFigureFullScreen(figure);
k = makeFigureFullScreen(figure);
l = makeFigureFullScreen(figure);
set(h,'Name','Direction maps for different temporal frequencies')
set(k,'Name','Orientation maps for different temporal frequencies')
set(l,'Name','Orientation magnitude maps for different temporal frequencies')
nStims = size(tfAveragedMaps, 3);
nRows  = ceil(sqrt(nStims));
nCols  = nRows;
clippingPercentile = 0.2;
clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]);
tempFreq_cell = strsplit(metadata.StimParams.temporalFreq(2:end-1), ',');
for i = 1:nStims
    figure(h); subplot(nRows,nCols,i);
    dirMapTf = analysis.(field).roi.directionTfMap(:,:,i);
    imagesc(polarMapEpi(dirMapTf, [0 clipValue(2)]));
    axis image; axis off;
    caxis([0 360]);
    cb=colorbar; colormap(cb, hsv);
    title([tempFreq_cell{i}, ' Hz'])
    figure(k); subplot(nRows,nCols,i);
    oriMapTf = analysis.(field).roi.orientationTfMap(:,:,i);
    imagesc(polarMapEpi(oriMapTf, clipValue));
    axis image; axis off;
    caxis([0 180]);
    cb=colorbar; colormap(cb, hsv);
    title([tempFreq_cell{i}, ' cpd'])
    figure(l); subplot(nRows,nCols,i);
    imagesc(magMap(oriMapTf, clipValue));
    colormap('gray');
    axis image; axis off;
    caxis(clipValue);
    title([tempFreq_cell{i}, ' cpd']) 
end
saveas(h, fullfile(saveDirectory, 'DirectionMapsSeparatedbyTF.png'))
saveas(k, fullfile(saveDirectory, 'OrientationMapsSeparatedbyTF.png'))
saveas(l, fullfile(saveDirectory, 'OrientationMagnitudeMapsSeparatedbyTF.png'))

%% make HLS maps
analysis.(field).roi.tfMapHLS = stimulusMap(tfAveragedMaps, clipValue);
figure
imshow(analysis.(field).roi.tfMapHLS);
saveas(gcf, fullfile(saveDirectory, 'HLS_TF.png'))

analysis.(field).roi.orientationMapHLS = stimulusMap(orientationAveragedMaps, clipValue);
figure
imshow(analysis.(field).roi.orientationMapHLS);
saveas(gcf, fullfile(saveDirectory, 'HLS_Orientation.png'))

analysis.(field).roi.orientationMapHLS = stimulusMap(directionAveragedMaps, clipValue);
figure
imshow(analysis.(field).roi.orientationMapHLS);
saveas(gcf, fullfile(saveDirectory, 'HLS_Direction.png'))