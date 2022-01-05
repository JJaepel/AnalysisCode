close all; 
clear all

imageDirectory = 'Z:\Juliane\Data\Epi\'; 
metaDirectory = 'Z:\Juliane\Data\Spike2Data\';
analysisDir = 'Z:\Juliane\Data\ImageAnalysis\';
name = 'F2604_2021-12-15'; 
exptID = 6;
sp2ID = exptID;

%define parameters
spatialDownsamplingFactor = 5;
timecourse = 1;
responseSign = -1;
bandpassfilter = 0;
contours = 0;

%define location of imaging data and meta data
imagingDataDirectory = sprintf('%s%s\\tseries_%d\\',imageDirectory,name,exptID);
metaDataDirectory    = sprintf('%s%s\\t%05d\\',metaDirectory,name,sp2ID);
saveDataDirectory   = sprintf('%s%s\\t%05d\\',analysisDir,name,exptID);
saveDataDir = sprintf('%s%s\\',analysisDir,name);

%make savedir
if ~exist(saveDataDirectory, 'dir')
    mkdir(saveDataDirectory);
end

%% load data
%Load imaging data
data = readingImagingData(imagingDataDirectory,spatialDownsamplingFactor);

%Load metadata
ImagingTime = load(strcat(metaDataDirectory,'\twophotontimes.txt')); %alternatively: '\frametrigger.txt'

stim_matrix=[];
disp(strcat('Loading....',metaDataDirectory, '\stimontimes.txt'))
stimtimes=load(strcat(metaDataDirectory, '\stimontimes.txt'));
i=1;
while i <length(stimtimes)
    idx= stimtimes(i);
    stim_matrix(1, floor(i/2)+1)= stimtimes(i);
    stim_matrix(2, floor(i/2)+1)= stimtimes(i+1);
    i=i+2;
end

StimParams = struct;
StimParams.StimOnTimes = stim_matrix(:,1:end);
StimParams.uniqStimIds= unique(StimParams.StimOnTimes(1,1:end));
StimParams.uniqStims = length(StimParams.uniqStimIds);

%load stimulus parameters
files=dir(strcat(metaDataDirectory, '\*.py'));
files={files.name};
file = [metaDataDirectory filesep char(files{1})];
StimParams.stimDuration = FindStimulusParam(file, 'stimDuration'); 
StimParams.numOrientations = FindStimulusParam(file, 'numOrientations');
StimParams.isi = FindStimulusParam(file, 'isi');
StimParams.numTrials = FindStimulusParam(file, 'numTrials');
StimParams.numBlank = FindStimulusParam(file, 'nBlank');
StimParams.rate = 1/mean(diff(ImagingTime));
StimParams.spatialFreq = FindStimulusParam(file, 'spatialFreq');
StimParams.temporalFreq = FindStimulusParam(file, 'temporalFreq');

for i =1:size(StimParams.StimOnTimes,2)
    StimParams.stimStartIndex(i) = floor(find(ImagingTime> StimParams.StimOnTimes(2,i),1));
    StimParams.stimStopIndex(i)  = ceil(StimParams.stimStartIndex(i) + StimParams.stimDuration*StimParams.rate);
end
clear ImagingTime

%% chop stimulus trace

preTrialTime = 0.5 * StimParams.isi;
postTrialTime = 0.5 * StimParams.isi;
offsetPre = round(StimParams.rate * preTrialTime);
offsetPost = round(StimParams.rate * postTrialTime);
analysisPeriod = -offsetPre:(StimParams.stimDuration*StimParams.rate+offsetPost);
blPeriod = 1:offsetPre;
nFrames = length(analysisPeriod);
stimResponseTrace = zeros([size(data,1),size(data,2),nFrames, StimParams.uniqStims,StimParams.numTrials],'single');
avgResponseTrace = zeros(nFrames,StimParams.uniqStims,StimParams.numTrials, 'single');
for stimID = 1:StimParams.uniqStims
    %get stimID
    stimulus = StimParams.uniqStimIds(stimID);
    disp(strcat('Processing StimID: ', num2str(stimulus)))
    stimIndices = find(StimParams.StimOnTimes(1,:)==stimulus);
    
    %loop through each trial
    for trialNumber = 1:StimParams.numTrials
        startFrame = StimParams.stimStartIndex(stimIndices(trialNumber));
        stimResponseTrace(:,:,:,stimID,trialNumber) = single(data(:,:,startFrame+analysisPeriod));
    end
end

%clear data

%blank stimReponseTrace
disp('Pre-Stimulus Blank')
for stimID = 1:StimParams.uniqStims
    for trial = 1:StimParams.numTrials
        blmean = squeeze(nanmean(stimResponseTrace(:,:, blPeriod, stimID, trial),3));
        for frameNum = 1:nFrames
            stimResponseTrace(:,:, frameNum, stimID, trial) = (squeeze(stimResponseTrace(:,:,frameNum, stimID, trial)) - blmean)./blmean;
        end
        clear blmean
    end
end

if bandpassfilter
    ROI = ones(size(stimResponseTrace,1), size(stimResponseTrace,2));
    for stimID = 1:StimParams.uniqStims
        for trial = 1:StimParams.numTrials
            for frameNum=1:nFrames
                temp = stimResponseTrace(:,:, frameNum, stimID, trial);
                tempBpf=LowHighNormalize(double(squeeze(squeeze(temp))),ROI);
                stimResponseTrace(:,:, frameNum, stimID, trial)=tempBpf;
            end
        end
    end
end



if timecourse
    %show response trace
    disp('Plotting time course');
    figure;
    CM = colormap('jet');
    close;
    IX = round(linspace(1,64,StimParams.uniqStims));
    for i = 1:length(IX)
        C{i} = CM(IX(i),:);
    end
    ROI = true(size(stimResponseTrace,1), size(stimResponseTrace,2));
    responseTraces = permute(stimResponseTrace, [5,3,4,1,2]);
    for condition = 1:StimParams.uniqStims
        t = ((0:(nFrames-1))/StimParams.rate)-preTrialTime;
        y = nanmedian(responseTraces(:,:,condition,ROI(:)),4);
        yMean = squeeze(mean(y));
        yStdError = std(y,[],1)/sqrt(size(y,1));
        plot(t, yMean, '-', 'LineWidth', 3, 'Color', C{condition}); hold on
        if condition == StimParams.uniqStims
            plot(t, yMean, '-', 'LineWidth', 3, 'Color', 'k'); hold on
        end
    end
    xlabel('Time (s)');
    ylabel('Response amplitude');
    axis square;
    set(gca,'Box','off');
    clear responseTraces

    % Add stimulus box
    stimStart = 0;
    stimStop  = StimParams.stimDuration;
    yLimits = get(gca,'YLim');
    rectangle('Position',[stimStart yLimits(2) stimStop-stimStart 0.025*range(yLimits)],'FaceColor','k')
end
saveas(gcf, fullfile(saveDataDirectory, 'Timecourse.png'))

%make trialAveragedMaps
disp('Making trial averaged maps')
stimResponseTrace = permute(stimResponseTrace, [1, 2, 4, 5, 3]);
stimResponseTrace = mean(stimResponseTrace(:,:,1:end-1,:,:),5);
trialAveragedMaps = squeeze(median(stimResponseTrace,4));
trialAveragedMaps(isnan(trialAveragedMaps(:))) = 0;
nOri =(size(trialAveragedMaps,3)/2);
oriTrialAveragedMaps = (trialAveragedMaps(:,:,1:nOri)+ trialAveragedMaps(:,:,(nOri+1):end))/2;
clear stimResponseTrace

collapsedData = permute(trialAveragedMaps, [3, 1, 2]);
clippingPercentile = 0.2;
clipValue = prctile(collapsedData(:), [clippingPercentile 100-clippingPercentile]);
clear collapsedData

%compute maps
directionMap = polarMap(vectorSum(responseSign*imgaussfilt(trialAveragedMaps), 1, 3).*ROI);
orientationMap = polarMap(vectorSum(responseSign*imgaussfilt(trialAveragedMaps), 2, 3).*ROI);

%display direction maps
nRows = ceil(sqrt(StimParams.numOrientations+1));
j = makeFigureFullScreen(figure);
figure(j)
for stims = 1:StimParams.numOrientations
    subplot(nRows,nRows,stims)
    imagesc(trialAveragedMaps(:,:,stims));
    colorbar; colormap('gray')
    title(num2str(stims), 'fontsize', 5);
    caxis(clipValue);
    axis image; axis off
end

subplot(nRows, nRows, stims+1)
imagesc(directionMap)
axis image; axis off
caxis([0 360]);
cb = colorbar; colormap(cb, hsv);
title('Direction map')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDataDirectory, 'Direction.png'))

%display orientation maps
nRows = ceil(sqrt(StimParams.numOrientations/2+1));
h = makeFigureFullScreen(figure);
figure(h)
for stims = 1:StimParams.numOrientations/2
    subplot(nRows,nRows,stims)
    imagesc(imgaussfilt(oriTrialAveragedMaps(:,:,stims)));
    colorbar; colormap('gray')
    title(num2str(stims), 'fontsize', 5);
    axis image; axis off
    caxis(clipValue);
end

subplot(nRows, nRows, stims+1)
imagesc(orientationMap)
axis image; axis off
caxis([0 180]);
cb = colorbar; colormap(cb, hsv);
title('Orientation map')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDataDirectory, 'Orientation.png'))
