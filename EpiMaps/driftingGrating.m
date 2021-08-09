close all
clear all
addpath(genpath('F:\Code\FromMadineh\ScanImage'));

animal = 'F2538_2021-07-02';
expt_id =2;
sp2id = expt_id;

close all
EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
SaveDir = 'Z:\Juliane\Data\ImageAnalysis\';

windowStop=2;
windowStart=0;
pre=1;
field = 'rawF';
intrinsic = 0;
downsample = 1;
svd = 0;


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
t0=tic;
data.rawF = readingImagingData(EpiDirectory);
toc(t0)

%% if applicaple, downsample
if downsample
    downsampleFactor = 2;
    stackSize        = size(data.rawF);
    stackSize     = floor(stackSize/downsampleFactor);
    downsampledStack = zeros(stackSize,class(data.rawF));
    for i=1:downsampleFactor
        downsampledStack = downsampledStack+data.rawF(i:downsampleFactor:downsampleFactor*stackSize(1),i:downsampleFactor:downsampleFactor*stackSize(2),i:downsampleFactor:downsampleFactor*stackSize(3))/downsampleFactor;
    end
    data.rawF = downsampledStack;
    clear downsampledStack
end
if svd
    [data.rawF, mixedfilters, percent] = SVDsimple(data.rawF);
end
expParam.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 
expParam.rawFMeanImg = mean(data.rawF,3);
expParam.baseImg = mean(data.rawF(:,:,1:50),3);
expParam.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);

%% create stimCodes

numberOfConditions = metadata.StimParams.numberOfStims;
stimStartIndex = zeros(numberOfConditions,1,'double'); 
stimStopIndex  = zeros(numberOfConditions,1,'double');

metadata.StimParams.stimDuration = 5;
metadata.StimParams.isi = 5;
metadata.StimParams.numTrials = 9;

for i=1:numberOfConditions
    stimStartIndex(i) = find(metadata.Imaging.time>=metadata.StimParams.StimOnTimes(2,i),1,'first');
    stimStopIndex(i) = find(metadata.Imaging.time>=metadata.StimParams.StimOnTimes(2,i)+metadata.StimParams.stimDuration,1,'first');
end
metadata.StimParams.stimStartIndex = stimStartIndex;
metadata.StimParams.stimStopIndex = stimStopIndex;

if downsample
    metadata.Imaging.time = metadata.Imaging.time(1:downsampleFactor:stackSize(3));
    metadata.Imaging.rate = metadata.Imaging.rate/downsampleFactor;
    try
        metadata.StimParams.stimStartIndex = floor(metadata.StimParams.stimStartIndex/downsampleFactor);
        metadata.StimParams.stimStopIndex = floor(metadata.StimParams.stimStopIndex/downsampleFactor);
    catch
    end
end

%% if applicaple, downsample
if intrinsic
    downsampleFactor = 5;
    stackSize        = size(data.rawF);
    stackSize(3)     = floor(stackSize(3)/downsampleFactor);
    downsampledStack = zeros(stackSize,class(data.rawF));
    for i=1:downsampleFactor
        downsampledStack = downsampledStack+data.rawF(:,:,i:downsampleFactor:downsampleFactor*stackSize(3))/downsampleFactor;
    end
    data.rawF = downsampledStack;

    metadata.Imaging.time = metadata.Imaging.time(1:downsampleFactor:stackSize(3));
    metadata.Imaging.rate = metadata.Imaging.rate/downsampleFactor;
    try
        metadata.StimParams.stimStartIndex = floor(metadata.StimParams.stimStartIndex/downsampleFactor);
        metadata.StimParams.stimStopIndex = floor(metadata.StimParams.stimStopIndex/downsampleFactor);
    catch
    end
end

%% chop traces
analysis = struct;
disp('Chopping Traces')
[analysis, metadata] = ChopStimulusTraceEpi(analysis,metadata,data,field);
analysis.rawFMeanImg = expParam.rawFMeanImg;
analysis.ROI = expParam.ROI;
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

%% map analysis
showEpiRespAvg(analysis, metadata,field, saveDirectory);