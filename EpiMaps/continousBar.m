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

metadata.StimParamsDown=Load_stimparams(Sp2dDirectoryDown);
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

%% additional functions
function [twophoton] = LoadFrameTimes(path)
    %ExtractFrameTimes Reads in Frame Triggers from frametrigger.txt
    %  	Args:
    %     path: path of the metadata file to extract
    %   Returns:
    %     twophoton (struct):   .time :trigger times (s)
    %                           .rate :frame rate (Hz)

    twoPhotonFile = 'twophotontimes.txt';
    counter = 1;

    if ~exist(path, 'dir')
        error('%s does not exist', path);
    end
    
    twophoton = struct;
    twophoton.time = [];
    tpFullFile = fullfile(path, twoPhotonFile);
    infotpFullFile = dir(tpFullFile);
    if infotpFullFile.bytes < 10
        twoPhotonFile = 'frametrigger.txt';
        tpFullFile = fullfile(path, twoPhotonFile);
    end
    while ~exist(tpFullFile, 'file') || isempty(twophoton.time)
        if exist(tpFullFile, 'file')
            disp(['Loading... ', tpFullFile])
            twophoton.time = load(tpFullFile);
            twophoton.rate = 1/median(diff(twophoton.time));

        else
            if counter > size(defaultFileNames,1)
               warning ('Could not get frame acquisition triggers')
               return;
            end
            tpFullFile = fullfile(twoPhotonPath, char(defaultFileNames{counter}));
            counter = counter + 1;
        end
    end
end
function [stim_params]=Load_stimparams(exptpath)
stim_params=struct;
    stim_matrix=[];
    disp(strcat('Loading....',exptpath, '\stimontimes.txt'))
    stimtimes=load(strcat(exptpath, '\stimontimes.txt'));
    i=1;
    while i <length(stimtimes)
        idx= stimtimes(i);
        stim_matrix(1, floor(i/2)+1)= stimtimes(i);
        stim_matrix(2, floor(i/2)+1)= stimtimes(i+1);
        i=i+2;
    end
    stim_params.StimOnTimes=stim_matrix;
    stim_params.StimOnTimes=stim_params.StimOnTimes(:, 2:end); %discard the very first stim because its garbage
    files=dir(fullfile(exptpath, 'frametrigger.txt'));
    files={files.name};
    if ~isempty(files)
        frametimes=load(fullfile(exptpath, 'frametrigger.txt'));
    else
        frametimes=[];
    end
    stim_params.frameTimes=frametimes;
    if isempty(stim_params.StimOnTimes)
        stim_params.type = 'Spontaneous';
        return;
    else
        stim_params.uniqStimIds= unique(stim_params.StimOnTimes(1,1:end)); %always ignore the first stim code
        stim_params.uniqStims= length(stim_params.uniqStimIds);
        stim_params.numberOfStims=length(stim_params.StimOnTimes);
        files=dir(strcat(exptpath, '\*.py'));
        files={files.name};

        for i=1:length(files) % find the type of stimulus
            if isempty(strfind(char(files{i}), 'serialTriggerDaqOut'))
                stim_params.type= strrep(char(files{i}), '.py', '');
                stim_params.file= strcat(exptpath,'\',char(files{i}));
            end
        end
        if isempty(stim_params.type) || strcmp(stim_params.type, 'blackScreen') % check if there wasn't a
            stim_params.type= 'Spontaneous';
        end

        %extract all the relevant field for the stimulation
        stimFields= stimField(stim_params.type);
        for i=1:length(stimFields) %
            stim_params.(stimFields{i})= FindStimulusParam(stim_params.file, char(stimFields{i}));
        end
        %correct the stim ids for driftingGrating to ignore the
        if isfield(stim_params, 'flashInterval')
            stim_params.stimDuration = stim_params.flashInterval;
        elseif isfield(stim_params, 'stimDuration')
            stim_params.stimDuration = stim_params.stimDuration;
        else 
            stim_params.stimDuration = NaN;
        end
    end

    if strcmp(stim_params.type, 'driftingGrating')
        stims = stim_params.uniqStims;
        if stim_params.doBlank
            stims = stims-1;

        end
          stim_params.theta = (0:2*pi/stims:2*pi - 2*pi/stims);
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
function [fields] = stimField(Stimtype)
    switch Stimtype
        case 'Spontaneous'
            fields ={};
        case 'driftingGrating'
            fields = {'numTrials',...
                'doBlank',...
                'changeDirectionAt',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'temporalFreq',...
                'spatialFreq',...
                'contrast',...
                'textureType',...
                'maxv',...
                'minv',...
                'dutyCycle',...
                'foregroundColor',...
                'startingPhase',...
                'centerPoint',...
                'stimSize',...
                'animalOrientation'
                };
        case 'fullScreenFlash'
            fields = {'initialDelay',...
                'lum0',...
                'lum1',...
                'flashInterval',...
                'isi',...
                'numTrials'...
                };
        case 'waveletStim_forJoe'
            fields = {'Iname',...
                'repeats',...
                'interRepeatInt',...
                'interMovieInt',...
                'stimFramesPerTrial',...
                'movies'

                };
        case 'continuousEdge_withOrientationsAndTriggers'
            fields = {'numOrientations', ...
                'animalOrientation',...
                'barColor',...
                'backGroundColor',...
                'backGroundBarColor',...
                'isRandom',...
                'edgeWidth',...
                'centerPoint',...
                'numberOfTrials',...
                'initialDelay',...
                'movementPeriod',...
                'animalName'
                };
        case 'MovingDots'
            fields = {'numTrials',...
                'doBlank',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'nDots',...
                'dotSize',...
                'contrast',...
                'fieldShape',...
                'dotLife',...
                'minv',...
                'fieldSize',...
                'speed',...
                'coherence',...
                'theta',...
                };
        case 'Patch'
            fields = {'numTrials',...
                'doBlank',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'centerPoint',...
                'stimSize',...
                'numStimElev',...
                'numStimAzi',...
                };
        case 'driftingGrating_ori_sf'
            fields = {'numTrials',...
                    'doBlank',...
                    'changeDirectionAt',...
                    'stimDuration',...
                    'isi',...
                    'isRandom',...
                    'initialDelay',...
                    'temporalFreq',...
                    'spatialFreq',...
                    'contrast',...
                    'textureType',...
                    'numOrientations',...
                    'dutyCycle',...
                    'foregroundColor',...
                    'startingPhase',...
                    'centerPoint',...
                    'stimSize',...
                    };
        otherwise
            fields ={};
    end
end
function [data] = LoadImagingData(EpiDirectory)
    filename = [EpiDirectory filesep '*.tif'];
    files = dir(filename);
    numberOfFiles = size(files, 1);
    if(numberOfFiles ==0), error('The Image path does not contain any tifs'); end
    framesPerFile = zeros(numberOfFiles,1);
    filesToCheck = [1 numberOfFiles];
    for n = filesToCheck
        fileName  = char(strcat(EpiDirectory,'\',files(n).name));
        imageInfo = imfinfo(fileName);
        framesPerFile(n)=numel(imageInfo);
    end
    framesPerFile(framesPerFile==0) = framesPerFile(1);
    numberOfFrames = sum(framesPerFile);
    data.rawF  = zeros([imageInfo(1).Height,imageInfo(1).Width,numberOfFrames],'uint16');
    data.ROI   =true( [imageInfo(1).Height,imageInfo(1).Width]); 
    
    frameCounter = 0;
    for n = 1:numberOfFiles
        fileName  = char(strcat(EpiDirectory,'\',files(n).name));
        data.rawF(:,:,frameCounter+[1:framesPerFile(n)])=read_Tiffs(fileName,1,100);
        frameCounter=frameCounter+framesPerFile(n);
    end
    data.rawFMeanImg     = mean(data.rawF,3);
end
function [analysis, metadata] = ChopStimulusTrace(analysis,metadata,data,field)
    offsetPre  = round(metadata.Imaging.rate * metadata.StimParams.isi/2);
    offsetPost = round(metadata.Imaging.rate * metadata.StimParams.isi/2);
    stimDuration = ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration);
    blPeriod   = (1:offsetPre);

    % Collect the stimulus evoked traces
    imgStack = data.(field);
    [nRows,nCols,~] = size(imgStack);
    analysisPeriod = -offsetPre:(stimDuration+offsetPost);
    nFrames        = length(analysisPeriod);
    stimResponseTrace = zeros([metadata.StimParams.uniqStims,metadata.StimParams.numTrials-1,nRows,nCols,nFrames,],'single');
    for stimID= 1:metadata.StimParams.uniqStims
        % Get stim ID
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        disp(strcat('Processing StimID: ', num2str(stimulus)))
        stimIndices = find(metadata.StimParams.StimOnTimes(1,:)==stimulus);

        % Loop through each trial for the associated stim ID
        for trialNumber= 1:metadata.StimParams.numTrials-1
            startFrame = metadata.StimParams.stimStartIndex(stimIndices(trialNumber));
            stimResponseTrace(stimID,trialNumber,:,:,:) = single(data.(field)(:,:,startFrame+analysisPeriod));
        end
    end
    
    disp('Pre-Stimulus Blank')
    for stimID = 1:metadata.StimParams.uniqStims
        for trial = 1:metadata.StimParams.numTrials-1
            blmean = squeeze(nanmean(stimResponseTrace(stimID, trial, :,:, blPeriod), 5));
            for frameNum=1:nFrames
                stimResponseTrace(stimID, trial, :,:, frameNum)=(squeeze(stimResponseTrace(stimID, trial, :,:, frameNum))-blmean)./blmean;
            end
        end
    end
    
    analysis.(field).roi.stimResponseTrace = permute(stimResponseTrace,[5 2 1 3 4]);
    %defaultOrder = {'x','y','t','c','n'}; % DO NOT CHANGE this variable!!! This is the way that the code originally structures the stimResponseTrace
                
end  
function showEpiRespAvg(analysis, metadata, field, saveDirectory)
    analysis.(field).roi.stimResponseTrace = permute(analysis.(field).roi.stimResponseTrace, [4 5 3 2 1]);
    includedFrames = [round(metadata.Imaging.rate * metadata.StimParams.isi/2)+1:round(metadata.Imaging.rate * metadata.StimParams.isi/2)+ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration)];
    stimResponseTrace = mean(analysis.(field).roi.stimResponseTrace(:,:,1:(end-1),:,includedFrames),5);
    
    trialAveragedMaps = squeeze(mean(stimResponseTrace,4));
    trialAveragedMaps(isnan(trialAveragedMaps(:))) = 0;
    
    analysis.(field).roi.orientationMap = vectorSum(trialAveragedMaps,2,3);
    analysis.(field).roi.directionMap   = vectorSum(trialAveragedMaps,1,3);

    clippingPercentile = 0.2;
    clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]); 
    
    mapType = {'Orientation','Direction'};
    for j = 1:2 
        name = mapType{j};
        switch name
            case 'Orientation'
                nOri      = (size(trialAveragedMaps,3)/2);
                mapSet    = (trialAveragedMaps(:,:,1:nOri) ...
                            +trialAveragedMaps(:,:,(nOri+1):end))/2;
                maxTheta  = 180;
                fieldName = 'orientationMap';
            case 'Direction'
                mapSet    = trialAveragedMaps;
                maxTheta  = 360;
                fieldName = 'directionMap';
        end
        nStims = size(mapSet,3);
        nRows  = ceil(sqrt(nStims+1));
        nCols  = nRows;

        % Show trial-averaged single-condition maps for orientation/direction
        h = makeFigureFullScreen(figure);
        set(h,'Name',name)
        thetas=(0:maxTheta/(nStims):maxTheta-1);
        for i = 1:nStims
            figure(h); subplot(nRows,nCols,i);
                imagesc(mapSet(:,:,i));
                colorbar; 
                colormap('gray');
                title(num2str(thetas(i)))
                axis image; axis off;
                caxis(clipValue);
        end

        % Show polar map
        figure(h); subplot(nRows,nCols,nStims+1);
            %responseImg = normalizeArray(max(trialAveragedMaps,[],3));
            %responseImg = std(trialAveragedMaps,[],3);
            imagesc(polarMap(analysis.(field).roi.(fieldName), [0 clipValue(2)])); 
            axis image; axis off;
            caxis([0 maxTheta]);
            cb=colorbar; colormap(cb, hsv);
            title(sprintf('%s map',name))
       saveas(gcf, fullfile(saveDirectory, [fieldName '.png']))
       
       % Show magnitude map
       figure
       imagesc(magMap(analysis.(field).roi.(fieldName),[0 clipValue(2)]));
       colormap('gray');
       axis image; axis off;
       caxis(clipValue);
       title(sprintf('%s magnitude map',name))
       saveas(gcf, fullfile(saveDirectory, [fieldName 'Magnitude.png']))
    end
    
    analysis.(field).roi.stimResponseTrace = permute(analysis.(field).roi.stimResponseTrace, [5 4 3 1 2]);
end
function summedVector = vectorSum(array,harmonic,dim)
    % summedVector = vectorSum(stack,harmonic,dim)
    % vector sum of input array with 2nd harmonic response 
    % along the last dimension of the array (unless 
    % otherwise specified by user input). 

    if(nargin<3), dim = length(size(array)); end
    if(nargin<2), harmonic = 2;              end

    stackSize = size(array);
    phaseArraySize = 1+0*stackSize; phaseArraySize(dim) = stackSize(dim);

    vectorArraySize = stackSize;    vectorArraySize(dim) = 1;
    phaseValues = linspace(0,360,stackSize(dim)+1); 
    summedVector = sum(array.*repmat(reshape(exp(2*pi*1i*harmonic*mod(phaseValues(1:(end-1)),360/harmonic)/360),phaseArraySize),vectorArraySize),dim);
end
function [rgbMap] = polarMap(z, clippingValue)
    
    selectedRange    = [0 1];

    %%1.) Magnitude Map
    magnitudeMap = abs(z);
    
    % Clip Magnitude Map
    typeOfClipping = length(clippingValue);
    switch typeOfClipping
        case 1 % denotes clipping by std around the mean 
            excludingNaNs = ~isnan(magnitudeMap(:));
            highClipVal = mean(magnitudeMap(excludingNaNs(:)))+clippingValue*std(magnitudeMap(excludingNaNs(:)));
            lowClipVal  = mean(magnitudeMap(excludingNaNs(:)))-clippingValue*std(magnitudeMap(excludingNaNs(:)));
        case 2 % denotes absolute threshold clipping
            highClipVal = clippingValue(2);
            lowClipVal  = clippingValue(1);
    end
    
    magnitudeMap(magnitudeMap > highClipVal) = highClipVal;
    magnitudeMap(magnitudeMap < lowClipVal ) = lowClipVal;
    
    % Normalize Magnitude Map Between 0 and 1
    switch typeOfClipping
        case 1
            offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));
        case 2
            offsetMagnitudeMap = magnitudeMap - min(min(lowClipVal));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(highClipVal));
    end
    
    %% 2.) Phase Response Map
    phaseMap = angle(z)*180/pi;
    normalizedPhaseMap = wrapTo360(phaseMap)./360;
    normalizedPhaseMap = normalizedPhaseMap - selectedRange(1);
    normalizedPhaseMap = normalizedPhaseMap ./ (selectedRange(2)-selectedRange(1));
    normalizedPhaseMap(normalizedPhaseMap>1) = 1;
    normalizedPhaseMap(normalizedPhaseMap<0) = 0;
    
    %% 3.) RGB Polar Response Map
    HueMap        = normalizedPhaseMap;
    SaturationMap = normalizedMagnitudeMap; 
    ValueMap      = ones(size(z)); % normally ignored unless responseImg is specified
    rgbMap        = hsv2rgb(cat(3,HueMap,SaturationMap,ValueMap));
end
function [magnitudeMap] = magMap(z, clippingValue)
    magnitudeMap = abs(z);
    
    % Clip Magnitude Map
    typeOfClipping = length(clippingValue);
    switch typeOfClipping
        case 1 % denotes clipping by std around the mean 
            excludingNaNs = ~isnan(magnitudeMap(:));
            highClipVal = mean(magnitudeMap(excludingNaNs(:)))+clippingValue*std(magnitudeMap(excludingNaNs(:)));
            lowClipVal  = mean(magnitudeMap(excludingNaNs(:)))-clippingValue*std(magnitudeMap(excludingNaNs(:)));
        case 2 % denotes absolute threshold clipping
            highClipVal = clippingValue(2);
            lowClipVal  = clippingValue(1);
    end
    
    magnitudeMap(magnitudeMap > highClipVal) = highClipVal;
    magnitudeMap(magnitudeMap < lowClipVal ) = lowClipVal;
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
function tifStack = read_Tiffs(filePath,imgScaling,updateFrequency,useWaitBar)
    warning off;
    if(nargin<2), imgScaling = 0.5;    end
    if(nargin<3), updateFrequency = 5; end
    if(nargin<4), useWaitBar = false;  end
    tic;

    disp(['Reading Image Stack - ' filePath]);

    % Read TIF Header and Setup Information
    InfoImage=imfinfo(filePath);
        xImage=InfoImage(1).Width;
        yImage=InfoImage(1).Height;
        NumberOfImages=length(InfoImage);
        disp(['Finished Reading Image Header - ' num2str(toc) ' seconds Elapsed']);

    % use wait bar
    if(useWaitBar)
        h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
    %     currentPosition = get(h,'Position');
    %     offset = 100;
    %     set(h,'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
    %     currentPosition = get(get(h,'Children'),'Position');
    %     set(get(get(h,'Children')),'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
    else
        updateFrequency = round((updateFrequency/100)*NumberOfImages);
    end

    % Initialize MATLAB array to contain tif stack
    scaledX = round(xImage*imgScaling);
    scaledY = round(yImage*imgScaling);
    tifStack     = zeros(scaledY,scaledX,NumberOfImages,'uint16');

    codeVersion = 'alternativeMethod'; % both methods seem pretty similar in performance, but original method is a bit faster
    switch codeVersion
        case 'originalMethod'
            % uses the tifflib function to read images fast
            FileID = tifflib('open',filePath,'r');
            rps    = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
            rps    = min(rps,yImage);
            for i=1:(NumberOfImages)
                % display read progress
                if(useWaitBar)
                    waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
                else
                    if(mod(i+1,updateFrequency)==0)
                        disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                    end
                end

                % turn off warnings
                warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');

                % Go through each strip of data.
                currentImage = zeros(yImage,xImage);
                tifflib('setDirectory',FileID,i-1);
                for r = 1:rps:yImage
                  row_inds = r:min(yImage,r+rps-1);
                  stripNum = tifflib('computeStrip',FileID,r)-1;
                  currentImage(row_inds,:) = tifflib('readEncodedStrip',FileID,stripNum);
                end

                % Rescale data
                if(imgScaling ~= 1 && imgScaling>0)
                    tifStack(:,:,i) = imresize(currentImage,[scaledY scaledX]); % Scales image size
                else
                    tifStack(:,:,i) = currentImage;
                end
            end
            tifflib('close',FileID);
            disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
            warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
        case 'alternativeMethod'
            % Setup TIF object and Read-In Basic Information
            hTif = Tiff(filePath);

            warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            for i=1:NumberOfImages
                if(useWaitBar)
                    waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
                else
                    if(mod(i+1,updateFrequency)==0)
                        disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                    end
                end

                if(imgScaling ~= 1 && imgScaling>0)
                    tifStack(:,:,i) = imresize(hTif.read(),[scaledY scaledX]); % Scales image size
                else
                    tifStack(:,:,i) = hTif.read();
                end
                if(i == NumberOfImages)
                    hTif.close();
                    warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
                    disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                else
                    hTif.nextDirectory();                
                end
            end
    end

    if(useWaitBar),close(h); warning on; end
end
function varargout=shadedErrorBar(x,y,errBar,varargin)
    narginchk(3,inf)

    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('lineProps', '-k', @(x) ischar(x) | iscell(x));
    params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
    params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);
    params.addParameter('axis', gca, @(x) isa(x, 'matlab.graphics.axis.Axes'));

    params.parse(varargin{:});

    %Extract values from the inputParser
    lineProps =  params.Results.lineProps;
    transparent =  params.Results.transparent;
    patchSaturation = params.Results.patchSaturation;
    ax= params.Results.axis;
    if ~iscell(lineProps), lineProps={lineProps}; end


    %Process y using function handles if needed to make the error bar dynamically
    if iscell(errBar) 
        fun1=errBar{1};
        fun2=errBar{2};
        errBar=fun2(y);
        y=fun1(y);
    else
        y=y(:).';
    end

    if isempty(x)
        x=1:length(y);
    else
        x=x(:).';
    end


    %Make upper and lower error bars if only one was specified
    if length(errBar)==length(errBar(:))
        errBar=repmat(errBar(:)',2,1);
    else
        s=size(errBar);
        f=find(s==2);
        if isempty(f), error('errBar has the wrong size'), end
        if f==2, errBar=errBar'; end
    end

    if length(x) ~= length(errBar)
        error('length(x) must equal length(errBar)')
    end


    %Log the hold status so we don't change
    initialHoldStatus=ishold;
    if ~initialHoldStatus, hold on,  end

    H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation, ax);

    if ~initialHoldStatus, hold off, end

    if nargout==1
        varargout{1}=H;
    end
end
function H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation, ax)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot to get the parameters of the line
    H.mainLine=plot(ax, x,y,lineProps{:});


    % Work out the color of the shaded region and associated lines.
    % Here we have the option of choosing alpha or a de-saturated
    % solid colour for the patch surface.
    mainLineColor=get(H.mainLine,'color');
    edgeColor=mainLineColor+(1-mainLineColor)*0.55;

    if transparent
        faceAlpha=patchSaturation;
        patchColor=mainLineColor;
    else
        faceAlpha=1;
        patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
    end


    %Calculate the error bars
    uE=y-errBar(1,:);
    lE=y+errBar(2,:);


    %Add the patch error bar



    %Make the patch
    yP=[lE,fliplr(uE)];
    xP=[x,fliplr(x)];

    %remove nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];


    H.patch=patch(ax, xP,yP,1,'facecolor',patchColor, ...
                  'edgecolor','none', ...
                  'facealpha',faceAlpha);


%     %Make pretty edges around the patch. 
%     H.edge(1)=plot(x,lE,'-','color',edgeColor);
%     H.edge(2)=plot(x,uE,'-','color',edgeColor);



    uistack(H.mainLine,'top') % Bring the main line to the top
end