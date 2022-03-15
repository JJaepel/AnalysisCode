close all
clear all

animal = 'F2636_2022-02-14';
expt_id = 12;
sp2id = 12;

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

metadata.StimParams.minAzi = metadata.StimParams.centerPoint(2)-((metadata.StimParams.numStimAzi*metadata.StimParams.stimSize(2)/2))+ metadata.StimParams.stimSize(2)/2;
metadata.StimParams.Azimuths = linspace(metadata.StimParams.minAzi, metadata.StimParams.minAzi + (metadata.StimParams.numStimAzi-1)*metadata.StimParams.stimSize(2), metadata.StimParams.numStimAzi);

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
stimResponseTrace = mean(analysis.(field).roi.stimResponseTrace(:,:,1:end-1,:,includedFrames),5);

trialAveragedMaps = squeeze(median(stimResponseTrace,4));
trialAveragedMaps(isnan(trialAveragedMaps(:))) = 0;

clippingPercentile = 0.95;
clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]); 

AzimuthHLS = stimulusMap(trialAveragedMaps, clipValue);
AzimuthMap = vectorSum(trialAveragedMaps,2,3);

%% show maps
nStims = size(trialAveragedMaps,3);
nRows  = ceil(sqrt(nStims+1));
nCols  = nRows;

h = makeFigureFullScreen(figure);
for i = 1:size(trialAveragedMaps,3)
    figure(h); subplot(nRows,nCols,i);
        imagesc(trialAveragedMaps(:,:,i));
        colorbar; 
        colormap('gray');
        title(metadata.StimParams.Azimuths(i))
        axis image; axis off;
        caxis(clipValue);
end
saveas(gcf, fullfile(saveDirectory, 'TrialAveragedMaps.png'))

i = makeFigureFullScreen(figure);
figure(i); subplot(1,2,1);
imagesc(polarMap(AzimuthMap)); %'responseImg', responseImg
axis image; axis off;
cb=colorbar; colormap(cb, hsv);
title('Polar azimuth map')

figure(i); subplot(1,2,2);
imshow(AzimuthHLS)
title('HLS azimuth map')
saveas(gcf, fullfile(saveDirectory, 'AzimuthMaps.png'))
