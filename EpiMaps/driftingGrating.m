close all
clear all
addpath(genpath('F:\Code\FromMadineh\ScanImage'));

animal = 'F2424_2020-03-03';
expt_id =21;
sp2id = expt_id;

close all
EpiDir = 'F:\Data\Epi\';
Sp2Dir = 'F:\Data\Spike2Data\';
SaveDir = 'F:\Data\ImageAnalysis\';

windowStop=2;
windowStart=0;
pre=1;
field = 'rawF';
intrinsic = 1;


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
metadata.StimParams=Load_stimparams(Sp2dDirectory);
metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
metadata.StimParams.path=fullfile(Sp2dDirectory);
metadata.StimParams.series=expt_id;

%% load tiffs
t0=tic;
data.rawF = readingImagingData(EpiDirectory);
toc(t0)
expParam.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 
expParam.rawFMeanImg = mean(data.rawF,3);
expParam.baseImg = mean(data.rawF(:,:,1:50),3);
expParam.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);

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
[analysis, metadata] = ChopStimulusTrace(analysis,metadata,data,field);
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
function [tifStack] = readingImagingData(EpiDirectory)
    tifStack= [];
    tifFiles = dir([EpiDirectory filesep '*.tif']);
    numberOfFiles = size(tifFiles, 1);
    if(numberOfFiles ==0), error('The Image path does not contain any tifs'); end
    for currentFile = 1:length(tifFiles)
        disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);

        % Specify stack name
        filePath = [EpiDirectory tifFiles(currentFile).name];

        % Read images into tifStack
        tifStack = cat(3,tifStack,read_TiffsNew(filePath));
    end
end
function tifstack = read_TiffsNew(filePath)
    tHeader = tic();
    disp(['Reading Image Stack - ' filePath]);
    tiffReader = ScanImageTiffReader(filePath);
    tifstack = tiffReader.data();
    tifstack = permute(tifstack,[2 1 3]);
    tiffReader.close();
    disp(['Finished Reading Image Stack - ' num2str(toc(tHeader)) ' seconds Elapsed']);
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
       
       %show different polar map
       figure
       cmap=hsv;
       imagesc(polarMap(analysis.(field).roi.(fieldName), [0 clipValue(2)]));        colormap(cmap), cbh=colorbar; axis image ;
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
    %rgbMap        = hsv2rgb(cat(3,HueMap,SaturationMap,ValueMap));
    
    LUT = hsv;
    rgbMagnitudeMap = ind2rgb(im2uint8(normalizedMagnitudeMap)*length(LUT),gray);
    for(currentChannel = 1:3)
        rgbMagnitudeMap(:,:,currentChannel) = double(normalizedMagnitudeMap);
    end  
    rgbPhaseMap = double(ind2rgb(uint8(normalizedPhaseMap*length(LUT)),LUT)); 
    rgbMap = rgbPhaseMap.*rgbMagnitudeMap;
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
function map = bandpassFermiFilter_Revised(map,wavelengthLower,wavelengthUpper,spatialResolution,T,ROI)
    %function map = bandpassFermiFilter(map,wavelengthLower,wavelengthUpper,spatialResolution,T,ROI)
    %
    % Performs a bandpass fermi-filter at the desired wavelength cutoff and
    % spatial resolution on a complex field image (both in microns). The 
    % returned image is the filtered complex field. The complex field map is 
    % filtered by performing an individual highpass filter on the respective 
    % real and complex values that makeup a complex field.
    if(nargin < 5)
        T = 0.05;
    end
    if(nargin < 6)
        ROI = logical(ones(size(map)));
    end

    % Pad Image
    imsize = size(map);
    [map startingPosition] = padImage(map);
    ROI = logical(padImage(ROI));

    % Z-Score Image and Remove Any Pixels Outside ROI
    %map = map-mean(map(ROI));
    %map = map./std(map(ROI));
    map(ROI == 0) = 0;

    % Bandpass Filter Image
    if(wavelengthLower ~= -1) %lowpass filter only
        map = fermifilter(map,wavelengthLower,spatialResolution,T,ROI);
    end
    if(wavelengthUpper ~= -1) %highpass filter only
        map = map-fermifilter(map,wavelengthUpper,spatialResolution,T,ROI);
    end

    % Z-Score Image and Remove Any Pixels Outside ROI
    %map = map-mean(map(ROI));
    %map = map./std(map(ROI));
    map(ROI == 0) = 0;

    % Crop Image To Original Proportions
    map = removePaddingToImage(map,startingPosition,imsize);
    ROI = removePaddingToImage(ROI,startingPosition,imsize);
end
%% Miscellaneous Functions
function im = removePaddingToImage(im,startLocation,imsize)
    % function im = cropIm(im,dims)
    % takes in the dimensions specified by dims and crops the image around the 
    % image's center

    im = im(startLocation(1):(startLocation(1)+imsize(1)-1),startLocation(2):(startLocation(2)+imsize(2)-1));
end
function [paddedIm startingPosition] = padImage(im,paddingValue)
    % pads image to nearest value of 2^n and with a padding value defined by paddingValue

    if(nargin<2)
        paddingValue = 0;
    end

    paddedImageSize = round(2.^(ceil(log(max(size(im)))/log(2))));
    paddedIm = paddingValue.*ones(paddedImageSize);
    startingPosition = round((paddedImageSize-size(im))/2);
    for(cdim = 1:2) %Ensure startingPosition is within map
        if(startingPosition(cdim)<1)
            startingPosition(cdim) = 1;
        elseif(startingPosition(cdim)>size(im,cdim))
            startingPosition(cdim) = size(im,cdim);
        end
    end
    paddedIm(startingPosition(1):(startingPosition(1)+size(im,1)-1),startingPosition(2):(startingPosition(2)+size(im,2)-1)) = im;
end
function filtim=fermifilter(im,mm_bound,pix_res,T,ROI)
    %This is a function in which a spatial fermi filter is applied to an image
    %in fourier space.
    %
    %USAGE
    %function filtim=fermifilter(im,mm_boundLower,mm_boundUpper,pix_res,T,ROI)
    %
    %VARIABLE DEFINITIONS
    %im - Input image to be filtered
    %mm_bound - The wavelength in µm at which the fermi filter is defined to give a 
    %			50% response.  This is effectively the border of the low pass
    %			filter to be applied ot the image.
    %pix_res - The pixel resolution for the image. 
    %filtim - The filtered image to be output.
    %T - temperature of filter (should be 0.05)
    %ROI - roi for image to be filtered

    %From the image size, pixel resolution, and desired 50% filter response, calculate
    %the appropriate pixel cutoff for the filter to be applied to the 2D fft. 
    imsize=size(im);
    ap_size=pix_res*imsize(2);
    pix_bound=ap_size./mm_bound;

    %Create a matrix that is the same size as the input image and will serve as a map 
    %of the distance from the center of the 2D fft image.
    [Y,X]=meshgrid(1:imsize(1),1:imsize(2));
    Y=Y-(1+imsize(1)/2);
    X=X-(1+imsize(2)/2);
    dist=sqrt(X.^2+Y.^2);

    %Define the fermi filter to be used
    fermi2=1./(1+exp(-((pix_bound-abs(dist))./(T*pix_bound))));
    fermi2=fftshift(fermi2);

    %calculate the 2D fft of the input image and apply the fermi filter to that image.
    %next, take the inverse transform in order to get the filtered image back. 
    fftedImage=fft2(im,imsize(1),imsize(2));
    fftedROI=fft2(ROI,imsize(1),imsize(2));
    filtim=ifft2(fftedImage.*fermi2);
    filtROI=ifft2(fftedROI.*fermi2);
    filtim(ROI)=filtim(ROI)./filtROI(ROI);
end