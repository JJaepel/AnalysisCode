close all; 
clear all
imageDirectory = 'Z:\Juliane\Data\Epi\'; 
metaDirectory = 'Z:\Juliane\Data\Spike2Data\';
name = 'F2501_2021-03-23'; 
exptID = 1;
sp2ID = 1;
spatialDownsamplingFactor = 4;
analysisDir = 'Z:\Juliane\Data\ImageAnalysis\';
timecourse = 1;
responseSign = 1;
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
data = LoadImagingData(imagingDataDirectory, spatialDownsamplingFactor);

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
meanImg = mean(data, 3);
Mask = makeMask(meanImg);
%clear data
clear meanImg

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
dirMapHLS = stimulusMap(imgaussfilt(trialAveragedMaps), Mask);
oriMapHLS = stimulusMap(imgaussfilt(oriTrialAveragedMaps), Mask);
directionMap = polarMap(vectorSum(responseSign*imgaussfilt(trialAveragedMaps), 1, 3).*ROI);
orientationMap = polarMap(vectorSum(responseSign*imgaussfilt(trialAveragedMaps), 2, 3).*ROI);

%display direction maps
nRows = ceil(sqrt(StimParams.numOrientations+1));
j = makeFigureFullScreen(figure);
figure(j)
for stims = 1:StimParams.numOrientations
    subplot(nRows,nRows,stims)
    imagesc(imgaussfilt(trialAveragedMaps(:,:,stims)));
    colorbar; colormap('gray')
    title(num2str(stims), 'fontsize', 5);
    axis image; axis off
    caxis(clipValue);
end

subplot(nRows, nRows, stims+1)
imagesc(directionMap)
axis image; axis off
caxis([0 360]);
cb = colorbar; colormap(cb, hsv);
title('Direction map')

% figure
% imshow(dirMapHLS)
% title('Direction')

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

% figure
% imshow(oriMapHLS)
% title('orientation')

if contours
    RetInfo = load([saveDataDir '\RetBound.mat']);
    cMap = RetInfo.cMap_all;
    cMap_all = sum(cMap,3);
    Maskc = RetInfo.Maskc;

    cOri = zeros(size(orientationMap));
    cDir = zeros(size(directionMap));
    for c = 1:3
        temp_ori = orientationMap(:,:,c) + cMap_all;
        cOri(:,:,c) = times(temp_ori, Maskc);
        temp_dir = directionMap(:,:,c) + cMap_all;
        cDir(:,:,c) = times(temp_dir, Maskc);
    end
    figure
    imagesc(cOri)
    axis image; axis off
    caxis([0 180]);
    cb = colorbar; colormap(cb, hsv);
    title('Orientation map')
    
    figure
    imagesc(cDir)
    axis image; axis off
    caxis([0 360]);
    cb = colorbar; colormap(cb, hsv);
    title('Direction map')
    BWmap = times(sum(orientationMap,3), Maskc);
    BWmap = rescale(BWmap);
    colormap('gray')

    colorcon = zeros(size(cMap,1), size(cMap,2));
    for i=1:size(cMap,3)
        temp = cMap(:,:,i);
        colorcon(temp == 1) = i;
    end
end





currentAxis=gcf;
for i = 1:currentAxis.Number
    saveas(figure(i),[saveDataDirectory 'Online map ' num2str(i) '.tif'])
end

function data = LoadImagingData(imgPath, spatialDownsamplingFactor)
    downsample = 1/spatialDownsamplingFactor;
    % Determine number of multi-page tif images are present.    
    files= dir(strcat(imgPath,'\*.tif'));
    files={files.name}';
    numberOfFiles = size(files, 1);
    
    % Determine total number of frames present, and initialize rawF array
    framesPerFile = zeros(numberOfFiles,1);
    filesToCheck = 1:numberOfFiles;
    for n = filesToCheck 
        fileName  = char(strcat(imgPath,'\',files(n)));
        imageInfo = imfinfo(fileName);
        framesPerFile(n)=numel(imageInfo);
    end
    numberOfFrames = sum(framesPerFile);
    data = zeros([imageInfo(1).Height*downsample,imageInfo(1).Width*downsample,numberOfFrames],'uint16');
    
    % Read imaging frames into MATLAB 
    frameCounter = 0;
    for n = 1:numberOfFiles
        fileName = char(strcat(imgPath,'\',files(n)));
        data(:,:,frameCounter+[1:framesPerFile(n)])=read_Tiffs(fileName,downsample,100);
        frameCounter=frameCounter+framesPerFile(n);
    end

end
function [val]= FindStimulusParam(fname,paramKey)
    fid=fopen(fname);
    try
        % get the rows for the paramKey
        fileText= textscan(fid, '%s', 'delimiter', '\n');
        fileText= fileText{1};
        fileText= strrep(fileText, ' ', ''); % delete all whitespace
        keys = strfind(fileText, strcat(paramKey, '='));
        keys= find(~cellfun(@isempty, keys));
        line = fileText(keys(1));
        line = strsplit(char(line), '#');
        line = strsplit(char(line(1)), '=');
        val= str2num(char(line(2)));
        if isempty(val)
            val = char(line(2));
        end
    catch ME
        val= '';
    end
    fclose(fid);
end
function Mask = makeMask(meanimg)
    Mean = nanmean(nanmean(meanimg,2));
    Std = std(std(meanimg));
    threshold = Mean - 0.9 * Std;
    TempMask = meanimg;
    TempMask(TempMask<threshold) = nan;
    TempMask(~isnan(TempMask)) = 1;
    TempMask(isnan(TempMask)) = 0;
    [labeledObject,nPolygons] = bwlabel(TempMask);
    MaskSize = [];
    for p = 1:nPolygons
        MaskSize(p) = sum(labeledObject(:)==p);
    end
    [~,LargestObjNr] = max(MaskSize);
    Mask = zeros(size(TempMask,1), size(TempMask,2));
    Mask(labeledObject == LargestObjNr) = 1;
end
function fig=makeFigureFullScreen(fig,constrainDimensions)
    % Makes a full-screen figure;
    if(nargin<1),fig=gcf; end
    if(nargin<2),constrainDimensions = false;  end

    screenSize = get(0,'screensize');
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    if(constrainDimensions)
        PaperPosition = get(fig,'PaperPosition');
        PaperPosition(3) = (1.25*screenSize(3)/screenSize(4))*PaperPosition(4);
        set(fig,'PaperPosition',PaperPosition);
    end
end
function HLSmap = stimulusMap(stimulusAveragedMaps, Mask)
    % Calculate max stimulus, amplitude and tuning width (circular variance of sorted map)
    NumStims = size(stimulusAveragedMaps, 3);
    dFoFmaps = zeros(size(stimulusAveragedMaps,1), size(stimulusAveragedMaps,2), size(stimulusAveragedMaps,3));
    %BSmap = mean(stimulusAveragedMaps,3);
    for stim = 1:NumStims
        I = stimulusAveragedMaps(:,:,stim);
        %dFoFmaps(:,:,stim) = (I-BSmap);
        dFoFmaps(:,:,stim) = I;
    end
    [ResponseAmplitude,PrefStim] = max(imgaussfilt(dFoFmaps), [], 3 );
    ResponseAmplitude = times(ResponseAmplitude, Mask);

    % Create HLS map
    if NumStims > 3
        figure;
        CM = colormap('jet');
        close;
        IX = round(linspace(1,64,NumStims));
        for i = 1:length(IX)
            C{i} = CM(IX(i),:);
        end
    elseif NumStims == 3
        C{1} = [1 0 0];
        C{2} = [0 1 0];
        C{3} = [0 0 1];
    elseif NumStims == 2
        C{1} = [0 0.8 1];
        C{2} = [1 0 0.4];
    end
    ResponseAmplitude_Scaled = rescale(ResponseAmplitude);
    HLSmap = zeros( size(stimulusAveragedMaps,1), size(stimulusAveragedMaps,2), 3 );            
    for y = 1:size(stimulusAveragedMaps,1)
        for x = 1:size(stimulusAveragedMaps,2)
            for c = 1:3
                if ~isnan(PrefStim(y,x)) && ResponseAmplitude(y,x) > 0
                    HLSmap(y,x,c) = (C{PrefStim(y,x)}(c) .* ResponseAmplitude_Scaled(y,x));
                end
            end
        end
    end

    % Add color index
    [yRes,xRes,~] = size(HLSmap);
    xRes = xRes - 100;
    for s = 1:NumStims
        xRange = round(50+(((s*(xRes/NumStims)) - (xRes/(NumStims*2))):(s*(xRes/NumStims))));
        for c = 1:3
            HLSmap(yRes:yRes+10,xRange,c) = C{s}(c);
        end
    end
end
function rgbMap = polarMap(z)
    %calculate magnitude Map
    magnitudeMap = abs(z);
    excludingNaNs = ~isnan(magnitudeMap(:));
    highClipVal = mean(magnitudeMap(excludingNaNs(:)))+3*std(magnitudeMap(excludingNaNs(:)));
    lowClipVal  = mean(magnitudeMap(excludingNaNs(:)))-3*std(magnitudeMap(excludingNaNs(:)));
    magnitudeMap(magnitudeMap > highClipVal) = highClipVal;
    magnitudeMap(magnitudeMap < lowClipVal ) = lowClipVal;
    offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
    normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));
    
    %calculate Phase map
    phaseMap = angle(z)*180/pi;
    normalizedPhaseMap = wrapTo360(phaseMap)./360;
    normalizedPhaseMap(normalizedPhaseMap>1) = 1;
    normalizedPhaseMap(normalizedPhaseMap<0) = 0;
    
    HueMap = normalizedPhaseMap;
    SaturationMap = ones(size(z)); % normally ignored unless responseImg is specified
    ValueMap  = normalizedMagnitudeMap;
    rgbMap = hsv2rgb(cat(3,HueMap,SaturationMap,ValueMap));
end
