animal = 'F2537_2021-06-30';
expt_id = 18;
sp2id = 27;

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
metadata.StimParams.stimDuration = metadata.StimParams.stimDuration/2;
metadata.StimParams.isi = metadata.StimParams.isi/2;

%% load tiffs
data.rawF = readingImagingData(EpiDirectory);
data.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 
data.rawFMeanImg = mean(data.rawF,3);
data.baseImg = mean(data.rawF(:,:,1:50),3);
data.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);

%% create stimCodes

numberOfConditions = metadata.StimParams.numberOfStims;
metadata.StimParams.StimOnTimes(1,metadata.StimParams.StimOnTimes(1,:)==1)=2;
metadata.StimParams.StimOnTimes(1,metadata.StimParams.StimOnTimes(1,:)==0)=1;
metadata.StimParams.uniqStimIds = unique(metadata.StimParams.StimOnTimes(1,:));
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
stimResponseTrace = mean(analysis.(field).roi.stimResponseTrace(:,:,1:2,:,includedFrames),5);

trialAveragedMaps = squeeze(median(stimResponseTrace,4));
trialAveragedMaps(isnan(trialAveragedMaps(:))) = 0;

diffMap = -diff(-trialAveragedMaps,[],3);
diffMap = (diffMap- min(diffMap(:)))/ (max(diffMap(:))-min(diffMap(:))); 
rgbImg  = convertRGB(-trialAveragedMaps);

%% show maps
clippingPercentile = 0.95;
clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]); 

h = makeFigureFullScreen(figure);
for i = 1:size(trialAveragedMaps,3)
    figure(h); subplot(2,2,i);
        imagesc(trialAveragedMaps(:,:,i));
        colorbar; 
        colormap('gray');
        if i==1
            title('ON')
        else
            title('OFF')
        end
        axis image; axis off;
        caxis(clipValue);
end

% Show ON-OFF difference maps
set(h,'Name','ON-OFF')
figure(h); subplot(2,2,3);
    imagesc(diffMap);
    axis image; axis off; colorbar;
    title('Difference Map')
figure(h); subplot(2,2,4);
    imagesc(rgbImg);
    axis image; axis off; colorbar;
    title('Polarization Map')
saveas(gcf, fullfile(saveDirectory, 'Maps.png'))
