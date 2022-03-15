close all
clear all
addpath(genpath('F:\Code\FromMadineh\ScanImage'));

animal = 'F2635_2022-02-17';
expt_id =15;
sp2id = 16;

close all
EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
SaveDir = 'Z:\Juliane\Data\ImageAnalysis\';

analysisParams.pre = 0;
analysisParams.field = 'rawF';
analysisParams.downsample = 2;

%% make directories
analysisParams.EpiDirectory = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id) filesep];
if sp2id > 9
    analysisParams.Sp2dDirectory = [Sp2Dir animal filesep 't000' num2str(sp2id) filesep];
    analysisParams.saveDirectory = [SaveDir animal filesep 't000' num2str(expt_id) filesep];
else
    analysisParams.Sp2dDirectory = [Sp2Dir animal filesep 't0000' num2str(sp2id) filesep];
    analysisParams.saveDirectory = [SaveDir animal filesep 't0000' num2str(expt_id) filesep];
end

if ~exist(analysisParams.saveDirectory, 'dir')
    mkdir(analysisParams.saveDirectory);  
end

analysisParams.expt_id = expt_id;
analysisParams.baseDirectory = EpiDir;

%% load metadata
metadata.StimParams=LoadStimParams(analysisParams.Sp2dDirectory);
metadata.Imaging=LoadFrameTimes(analysisParams.Sp2dDirectory);
metadata.StimParams.path=fullfile(analysisParams.Sp2dDirectory);
metadata.StimParams.series=expt_id;

%% load tiffs
t0=tic;
data.rawF = readingImagingData(analysisParams.EpiDirectory,analysisParams.downsample);
toc(t0)

expParam.rawFMeanImg = mean(data.rawF,3);
expParam.baseImg = mean(data.rawF(:,:,1:50),3);
expParam.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);
expParam.ROI =true([size(data.rawF,1),size(data.rawF,2)]); 


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
[analysis, metadata] = ChopStimulusTraceEpi(analysis,metadata,data,analysisParams.field);
analysis.rawFMeanImg = expParam.rawFMeanImg;
analysis.ROI = expParam.ROI;
clear data

%% make timecourse

eval(sprintf('cmap = %s(%d);','hsv',metadata.StimParams.uniqStims-1));
LUT = cat(1,cmap, [0 0 0]);
numberOfFrames     = size(analysis.(analysisParams.field).roi.stimResponseTrace,1);

for condition = 1:metadata.StimParams.uniqStims
    t = ((0:(numberOfFrames-1))/metadata.Imaging.rate)-analysisParams.pre;
    y = nanmedian(analysis.(analysisParams.field).roi.stimResponseTrace(:,:,condition,analysis.ROI(:)),4); % Takes median response of the pixel values. Reduces noise relative to mean
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
saveas(gcf, fullfile(analysisParams.saveDirectory, 'Timecourse.png'))

%% map analysis
showEpiRespAvg(analysis, metadata,analysisParams.field, analysisParams.saveDirectory);