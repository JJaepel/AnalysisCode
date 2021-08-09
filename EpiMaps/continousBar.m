animal = 'F2380_2019-11-07';
expt_id_Up = 9;
sp2id_Up = expt_id_Up;
expt_id_Down = 9;
sp2id_Down = expt_id_Down;

close all
EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
SaveDir = 'Z:\Juliane\Data\ImageAnalysis\';

windowStop=2;
windowStart=0;
pre=1;
field = 'rawF';
intrinsic = 0;
responseSign = -1;


EpiDirectoryUp = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id_Up) filesep];
EpiDirectoryDown = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id_Down) filesep];
if sp2id_Up > 9
    Sp2dDirectoryUp = [Sp2Dir animal filesep 't000' num2str(sp2id_Up) filesep];
    saveDirectory = [SaveDir animal filesep 't000' num2str(expt_id_Up) filesep];
else
    Sp2dDirectoryUp = [Sp2Dir animal filesep 't0000' num2str(sp2id_Up) filesep];
    saveDirectory = [SaveDir animal filesep 't0000' num2str(expt_id_Up) filesep];
end
if sp2id_Down > 9
    Sp2dDirectoryDown = [Sp2Dir animal filesep 't000' num2str(sp2id_Down) filesep];
else
    Sp2dDirectoryDown = [Sp2Dir animal filesep 't0000' num2str(sp2id_Down) filesep];
end

if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end

sliceparams = struct;
sliceparams.expt_id_Up = expt_id_Up;
sliceparams.expt_id_Down = expt_id_Down;
sliceparams.baseDirectoryUp = EpiDirectoryUp;
sliceparams.baseDirectoryDown = EpiDirectoryDown;

%% load metadata
metadata.StimParamsUp=Load_stimparams(Sp2dDirectoryUp);
metadata.ImagingUp=LoadFrameTimes(Sp2dDirectoryUp);
metadata.StimParamsUp.path=fullfile(Sp2dDirectoryUp);
metadata.StimParamsUp.series=expt_id_Up;

metadata.StimParamsDown=LoadStimParams(Sp2dDirectoryDown);
metadata.ImagingDown=LoadFrameTimes(Sp2dDirectoryDown);
metadata.StimParamsDown.path=fullfile(Sp2dDirectoryDown);
metadata.StimParamsDown.series=expt_id_Down;

%% load tiffs
dataUp = LoadImagingData(EpiDirectoryUp);
dataDown = LoadImagingData(EpiDirectoryDown);

%% create stimCodes
stimTimesOn_HorUp = metadata.StimParamsUp.StimOnTimes(2,:);
stimTimesOn_HorDown = metadata.StimParamsDown.StimOnTimes(2,:);

stimStartIndex_HorUp = zeros(length(stimTimesOn_HorUp),1);
stimStartIndex_HorDown = zeros(length(stimTimesOn_HorDown),1);

for i = 1:length(stimTimesOn_HorUp)
    stimStartIndex_HorUp(i) = find(metadata.ImagingUp.time > stimTimesOn_HorUp(i), 1);
    stimStartIndex_HorDown(i) = find(metadata.ImagingDown.time > stimTimesOn_HorDown(i), 1);
end

metadata.StimParamsUp.stimStartIndex = stimStartIndex_HorUp;
metadata.StimParamsDown.stimStartIndex = stimStartIndex_HorDown;

%% if applicaple, downsample
if intrinsic
    downsampleFactor = 5;
    stackSize        = size(dataUp.rawF);
    stackSize(3)     = floor(stackSize(3)/downsampleFactor);
    downsampledStack = zeros(stackSize,class(dataUp.rawF));
    for i=1:downsampleFactor
        downsampledStack = downsampledStack+dataUp.rawF(:,:,i:downsampleFactor:downsampleFactor*stackSize(3))/downsampleFactor;
    end
    dataUp.rawF = downsampledStack;
    
    stackSize        = size(dataDown.rawF);
    stackSize(3)     = floor(stackSize(3)/downsampleFactor);
    downsampledStack = zeros(stackSize,class(dataDown.rawF));
    for i=1:downsampleFactor
        downsampledStack = downsampledStack+dataDown.rawF(:,:,i:downsampleFactor:downsampleFactor*stackSize(3))/downsampleFactor;
    end
    dataDown.rawF = downsampledStack;
    
    metadata.ImagingUp.time = metadata.ImagingUp.time(1:downsampleFactor:stackSize(3));
    metadata.ImagingUp.rate = metadata.ImagingUp.rate/downsampleFactor;
    metadata.ImagingDown.time = metadata.ImagingDown.time(1:downsampleFactor:stackSize(3));
    metadata.ImagingDown.rate = metadata.ImagingDown.rate/downsampleFactor;
    try
        metadata.StimParamsUp.stimStartIndex = floor(metadata.StimParams.stimStartIndexUp/downsampleFactor);
        metadata.StimParamsDown.stimStartIndex = floor(metadata.StimParams.stimStartIndexDown/downsampleFactor);
    catch
    end
end

%% 
time=tic;
harmonicsToTest = 1;
numberOfTrials = length(stimStartIndex_HorUp);
framesPerCycle = mode(diff(stimStartIndex_HorUp));
freqData  = zeros(size(dataDown.rawF,1),size(dataDown.rawF,2),length(harmonicsToTest),numberOfTrials*2);
avgImgResponse    = zeros(framesPerCycle,numberOfTrials);

%temporally restrict data
for trial = 1:numberOfTrials
    selectedFrames = metadata.StimParamsDown.stimStartIndex(trial):(metadata.StimParamsDown.stimStartIndex(trial)+framesPerCycle-1);
    data = zeros(size(dataDown.rawF,1),size(dataDown.rawF,2),length(selectedFrames)); 
    data(:,:,1:length(selectedFrames)) = dataDown.rawF(:,:,selectedFrames);
    avgImgResponse(:,trial) = mean(mean(data,1),2);
    clear data
end

%create temporal mask for times with low std
tempstd = std(avgImgResponse, [],2);
threshold = mean(tempstd) - 1 * std(tempstd);
indices = find(tempstd < threshold);
tempStart = indices(1);
tempStop = indices(end);
tempMask = 1;
if tempMask
    figure
    plot(tempstd)                    
    disp(['analysis time window is between image ' num2str(tempStart-5) ' and ' num2str(tempStop+5)])
end

for currentCycle = 1:numberOfTrials*2
    data = zeros(size(dataDown.rawF,1),size(dataDown.rawF,2),framesPerCycle);
    if currentCycle < numberOfTrials+1
        selectedFrames = stimStartIndex_HorUp(currentCycle):(stimStartIndex_HorUp(currentCycle)+framesPerCycle-1);
        data(:,:,1:length(selectedFrames)) = dataDown.rawF(:,:,selectedFrames);
    else
        selectedFrames = stimStartIndex_HorDown(currentCycle-numberOfTrials):(stimStartIndex_HorDown(currentCycle-numberOfTrials)+framesPerCycle-1);
        data(:,:,1:length(selectedFrames)) = flip(dataUp.rawF(:,:,selectedFrames),3);
    end
    if tempMask
        data = double(data);
        meanimg = mean(data,3);
        for i = 1:tempStart-5
            data(:,:,i) = meanimg;
        end
        for j = tempStop+5:size(data,3)
            data(:,:,j) = meanimg;
        end
    end
    data = double(data);
    if tempMask
        meanTrialFrame = mean(data(:,:,tempStart:tempStop),3); % Used as a Cocktail Blank
    else
        meanTrialFrame = mean(data,3); % Used as a Cocktail Blank
    end
    currentRelativeHarmonic = 0;
    for currentStimFreq = harmonicsToTest
        currentRelativeHarmonic = currentRelativeHarmonic+1;
            freqData(:,:,currentRelativeHarmonic,currentCycle) = freqData(:,:,currentRelativeHarmonic,currentCycle)...
                + responseSign*sum((data-repmat(meanTrialFrame,[1 1 size(data,3)]))...
                            .* repmat(reshape(exp(2*pi*1i*currentStimFreq*([1:size(data,3)]/framesPerCycle)),[1 1 size(data,3)]),[size(data,1) size(data,2) 1]),3);
    end
    disp(['  *Finished Reading and Constructing Imaging Data for Trial ' num2str(currentCycle) ' - Time Elapsed: ' num2str(toc(time)) ' seconds']);
end

%% create Mask
Mask2 = imgaussfilt(mean(dataUp.rawF, 3), 4);
Mean = nanmean(nanmean(Mask2,2));
Std = std(std(Mask2));
threshold = Mean - 0.9 * Std;
Mask2(Mask2<threshold) = nan;
Mask2(~isnan(Mask2)) = 1;
Mask2(isnan(Mask2)) = 0;
[labeledObject,nPolygons] = bwlabel(Mask2);
MaskSize = zeros(nPolygons,1);
for p = 1:nPolygons
    MaskSize(p) = sum(labeledObject(:)==p);
end
[~,LargestObjNr] = max(MaskSize);
Mask3 = zeros(size(Mask2,1), size(Mask2,2));
Mask3(labeledObject == LargestObjNr) = 1;
figure
subplot(2, 2, 1)
imagesc(Mask3)
title('Mask')
axis off

%% make maps
magnitudeMap = abs(mean(freqData(:,:,1,:),4)); 
magnitudeMap = times(Mask3, magnitudeMap);
magnitudeMap(magnitudeMap == 0) = nan;
magnitudeMap_Scaled = rescale(magnitudeMap);

subplot(2,2,2)
imagesc(magnitudeMap)
title('Response Amplitude')
axis off

phaseMap = imgaussfilt(angle(mean(freqData(:,:,1,:),4))*180/pi, 2);
ind_phaseMaps = zeros(size(phaseMap,1), size(phaseMap,2), size(freqData,4));
ind_magMap = zeros(size(phaseMap,1), size(phaseMap,2), size(freqData,4));
for trial = 1:size(freqData,4)
    ind_phaseMaps(:,:,trial) = times(imgaussfilt(angle(freqData(:,:,1,trial))*180/pi, 2), Mask3);
    ind_magMap(:,:,trial) = times(abs(freqData(:,:,1,trial)), Mask3);
end

phaseMap = times(Mask3, phaseMap);
phaseMap(phaseMap == 0) = nan;
ind_phaseMaps(ind_phaseMaps == 0) = nan;
ind_magMap(ind_magMap == 0) = nan;

positiveInput = (phaseMap > 0);
normphaseMap = mod(phaseMap, 180);
normphaseMap((normphaseMap == 0) & positiveInput) = 180;
normalizedPhaseMap = normphaseMap - 90;

subplot(2,2,3)
imagesc(normalizedPhaseMap)
title('PhaseMap')
axis off

%%
stimbounds = [-80; -60; -40; -20; 0; 20; 40; 60; 80];

figure; CM = colormap('jet'); close;
IX = round(linspace(1,64,length(stimbounds)+2));
for i = 1:length(IX)
    C{i} = CM(IX(i),:);
end
            
retpatchMap = zeros(size(normalizedPhaseMap, 1), size(normalizedPhaseMap,2));
for retpat = 1:length(stimbounds)-1
    bigger = normalizedPhaseMap > stimbounds(retpat);
    smaller = normalizedPhaseMap < stimbounds(retpat+1);
    all = logical(times(bigger, smaller));
    retpatchMap(all) = retpat+2; 
end
retpatchMap(normalizedPhaseMap > stimbounds(end)) = length(stimbounds)+2;
retpatchMap(normalizedPhaseMap < stimbounds(1)) = 2;
retpatchMap(isnan(normalizedPhaseMap)) = 1;

subplot(2,2,4)
imagesc(retpatchMap)
axis off

nRows  = ceil(sqrt(size(freqData,4)));
nCols  = nRows;

figure
for trial = 1:size(freqData,4)
    subplot(nRows,nCols,trial);
    imagesc(ind_phaseMaps(:,:,trial));
    axis off
end

figure
for trial = 1:size(freqData,4)
    subplot(nRows,nCols,trial);
    imagesc(ind_magMap(:,:,trial));
    axis off
end

RGBmap = zeros(size(retpatchMap,1), size(retpatchMap,2),3);
for y = 1:size(retpatchMap,1)
    for x = 1: size(retpatchMap,2)
        for c = 1:3
            RGBmap(y,x,c) = C{retpatchMap(y,x)}(c) .* magnitudeMap_Scaled(y,x);
        end
    end
end

tempRGB = RGBmap;
% Add color index
[yRes,xRes,~] = size(RGBmap);
xRes = xRes - 100;
for s = 1:length(stimbounds)+2
    xRange = round(50+(((s*(xRes/(length(stimbounds)+2))) - (xRes/((length(stimbounds)+2)*2))):(s*(xRes/(length(stimbounds)+2)))));
    for c = 1:3
        RGBmap(yRes:yRes+10,xRange,c) = C{s}(c);
    end
end

figure
imagesc(RGBmap)
title(['Stimulus bounderies: ' num2str(stimbounds')])
axis off