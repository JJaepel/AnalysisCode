function spontaneousNewFct(analysisParams)
%close all

% change things here - nowhere else!!! those are variables to decide which
% experiment to run and whether you need to reload, what filtering and how
% much to downsample
% variables of experiment
animal = analysisParams.animal;
expt_id = analysisParams.expID;
sp2id = analysisParams.sp2ID;

% variables for running - need to be 1 the first time it is run and can be
% changed afterwards to speed it up
reload = 0;
register = 1; 
newMask = 1;
field = analysisParams.field; %decide if you want to use dff, high or filtered data
downsample = analysisParams.downsample; % if you run into memory issues, change this number
verbose = 0; %change to 1 if you want to see all event locations and traces


%% 0.) Make folder structure
% choose folders and some variables depending on the computer that is
% running it
computer = getenv('COMPUTERNAME');
switch computer
    case 'DF-LAB-WS38'
        EpiDir = 'Z:\Juliane\Data\Epi\';
        Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
        SaveDir = 'Z:\Juliane\Data\ImageAnalysis\';
        lowMemory = 0;
    case 'DF-LAB-WS22'
        EpiDir = 'Z:\Juliane\Data\Epi\';
        Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
        SaveDir = 'Z:\Hannah\Data\ImageAnalysis\';
        lowMemory = 1;
end

EpiDirectory = [EpiDir filesep animal filesep num2str(expt_id) filesep];
Sp2dDirectory = [Sp2Dir animal filesep num2str(sp2id) filesep];
saveDirectory = [SaveDir animal filesep num2str(sp2id) filesep];


if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end

sliceparams = struct;
sliceparams.expt_id = expt_id;
sliceparams.baseDirectory = EpiDirectory;

%% 1.) Load tiffs
if reload %either load the data from scratch or load it from the last time you ran it and saved it
    t0=tic;
    disp('Loading data from tif files')
    rawF = readingImagingDataMaskRestricted(EpiDirectory, downsample, 1);
    toc(t0)
    
    if register %if necessary, do registration 
        numberOfImages = size(rawF,3);
        try
            disp('Loading shift data')
            load(fullfile(EpiDirectory, 'shifts.mat'),'rowShift', 'colShift')
        catch
            referenceImg = real(bandpassFermiFilter_Revised(rawF(:,:,1),-1,600,1000/172));
            rowShift = zeros(numberOfImages,1);
            colShift = zeros(numberOfImages,1);
            tic
            disp('Registering images')
            parfor(ii = 2:numberOfImages) %takes roughly 15 min
                [~, rowShift(ii), colShift(ii)] = registerImages(referenceImg,rawF(:,:,ii),[-1 600 1000/172],true)
            end
            toc
            save(fullfile(EpiDirectory, 'shifts.mat'),'rowShift', 'colShift');
        end
        disp('Shifting images')
        for ii = 2:numberOfImages
            [rawF(:,:,ii)] = shiftImages(rawF(:,:,ii),rowShift(ii), colShift(ii), true);
        end
        view_tiff_stack(rawF)
    end
    disp('Saving raw data to mat file')
    save(fullfile(saveDirectory, 'rawData.mat'),'-v7.3', 'rawF')
else
    disp('Loading raw data from mat file')
    load(fullfile(saveDirectory, 'rawData.mat'), 'rawF');
    load(fullfile(saveDirectory, 'data.mat'), 'data');
end

%% 2.) load metadata
metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
metadata.StimParams.path=fullfile(Sp2dDirectory);
expTimeTotal = mean(diff(metadata.Imaging.time))*size(rawF,3);
frameTimes = linspace(0,expTimeTotal, size(rawF,3));

%% 3.) get ROI of window and of subregions
if newMask %takes care of making the ROIs of the two areas
    disp('Making masks for analysis')
    analysis = struct;
    analysis.rawFMeanImg = mean(rawF,3);
    analysis.baseImg = mean(rawF(:,:,1:50),3);
    analysis.gaussMeanImg = imgaussfilt(mean(rawF, 3), 4);
    analysis.ROI =true( [size(rawF,1),size(rawF,2)]); 
    analysis = makeMasks(rawF, analysis, saveDirectory);
    save(fullfile(saveDirectory, 'masks.mat'), 'analysis');
else
    disp('Loading masks for analysis')
    load(fullfile(saveDirectory, 'masks.mat'), 'analysis');
end

if reload || newMask %if you have a new mask or reloaded the data, do dff and filtering again
    disp('Normalizing stack')
    data.dff = NormalizeStack(metadata,rawF);
    if lowMemory
        switch field
            case 'high'
                disp('High-pass filtering data')
                data.filt=LowHighNormalizeData(double(data.dff), analysis.mask, 1,10);
            case 'filt'
                disp('Applying bandpass filtering')
                data.high = HighNormalizeData(double(data.dff), analysis.mask,5);
        end
    else
        disp('Applying bandpass and highpass filter')
        data.filt=LowHighNormalizeData(double(data.dff), analysis.mask, 1,10);
        data.high = HighNormalizeData(double(data.dff), analysis.mask,5);
    end
    disp('Saving filtered data')
    %save(fullfile(saveDirectory, 'data.mat'),'-v7.3', 'data');
end
clear rawF

%% 4.) Get traces of each area
%data is goign to be turned into integers, but was normalized and filtered
%beforehand -> so multiply by 100 as you will otherwise cut off events
data.(field) = data.(field)*1000;

%save all important information in the analysis structure
V1BVmask = analysis.maskBV .* analysis.maskV1;
xVals = find(sum(V1BVmask,1)>0);
yVals = find(sum(V1BVmask,2)>0);
xRange = min(xVals):max(xVals);
yRange = min(yVals):max(yVals);
MV1 = repmat(int16(V1BVmask(yRange,xRange)), [1 1 size(data.(field),3)] ); %make a mask over the whole range
analysis.(field).DataV1 = int16(data.(field)(yRange,xRange,:)).* MV1; %multiplies it with the mask
V1Trace = squeeze(squeeze(sum( sum( int16(data.(field)(yRange,xRange,:)) .* MV1, 1 ), 2 ))) ./ nnz(V1BVmask);

A19BVmask = analysis.maskBV .* analysis.maskA19;
xVals = find(sum(A19BVmask,1)>0);
yVals = find(sum(A19BVmask,2)>0);
xRange = min(xVals):max(xVals);
yRange = min(yVals):max(yVals);
MA19 = repmat(int16(A19BVmask(yRange,xRange)), [1 1 size(data.(field),3)] ); %find the mask for it
analysis.(field).DataA19 = int16(data.(field)(yRange,xRange,:)) .* MA19; %saves only the data of the region to restrict memory use
A19Trace = squeeze(squeeze(sum( sum( int16(data.(field)(yRange,xRange,:)) .* MA19, 1 ), 2 ))) ./ nnz(A19BVmask);

% normalize Traces
analysis.(field).V1Trace = (V1Trace-min(V1Trace))/(max(V1Trace)-min(V1Trace));
analysis.(field).A19Trace = (A19Trace-min(A19Trace))/(max(A19Trace)-min(A19Trace));

%% 5.a) Find events in areas
analysis.(field).eventsV1 = findEvents(analysis.(field).V1Trace,'V1',verbose, saveDirectory); %function could still be a bit improved, set to verbose to look for what is goingw wrong - like separating multiple events
analysis.(field).eventsA19 = findEvents(analysis.(field).A19Trace,'A19',verbose, saveDirectory);

%% 5.b) Find event location at onset and peak
analysis.(field).eventsV1 = findLocations(analysis.(field).eventsV1, analysis.(field).DataV1,'V1',verbose, saveDirectory);
analysis.(field).eventsA19 = findLocations(analysis.(field).eventsA19, analysis.(field).DataA19,'A19',verbose, saveDirectory);

%% 5.c) Find local traces and events details
analysis.(field).eventsV1 = findLocalEvents(analysis.(field).DataV1,analysis.(field).eventsV1,'V1',verbose,saveDirectory); 
analysis.(field).eventsA19 = findLocalEvents(analysis.(field).DataA19,analysis.(field).eventsA19,'A19',verbose,saveDirectory); 

%% 5.d) Classify events based on where global onsets are
windowTime = [0.5, 1,2]; %within what time frame should the other event onset be? -> can be played around with
for i = 1:length(windowTime)
    analysis.(field).eventsV1 = findEventSequence(analysis.(field).eventsV1,analysis.(field).eventsA19, metadata, windowTime(i));
    analysis.(field).eventsA19 = findEventSequence(analysis.(field).eventsA19,analysis.(field).eventsV1, metadata, windowTime(i));
    if windowTime(i) < 1
        columnField = ['class0' num2str(windowTime(i)*10)];
    else
        columnField = ['class' num2str(windowTime(i)*10)];
    end
    
    figure
    subplot(1,2,1)
    labels = {'no other event', 'simultaneous', 'preceding events', 'following events'};
    classSizes = [sum([analysis.(field).eventsV1.(columnField)] == 0), sum([analysis.(field).eventsV1.(columnField)] == 1),sum([analysis.(field).eventsV1.(columnField)] == 2),sum([analysis.(field).eventsV1.(columnField)] == 3)];
    classPercentage = classSizes/sum(classSizes);
    bar(classPercentage)
    set(gca,'xticklabel',labels)
    ylabel('Percentage of events')
    title(['Events in V1, window: ' num2str(windowTime(i)) ' s']);
    
    subplot(1,2,2)
    classSizes = [sum([analysis.(field).eventsA19.(columnField)] == 0), sum([analysis.(field).eventsA19.(columnField)] == 1),sum([analysis.(field).eventsA19.(columnField)] == 2),sum([analysis.(field).eventsA19.(columnField)] == 3)];
    classPercentage = classSizes/sum(classSizes);
    bar(classPercentage)
    set(gca,'xticklabel',labels)
    ylabel('Percentage of events')
    title(['Events in A19, window: ' num2str(windowTime(i)) ' s']);
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, ['Event ' columnField ' s']))
end

%% 6.a) Do cross-correlogram for all events (globally, independent of whether there is an event in the other area)
doCrossCorrelation(analysis.(field).eventsV1,analysis.(field).V1Trace,analysis.(field).A19Trace,metadata, 'onset', saveDirectory, 'V1');
doCrossCorrelation(analysis.(field).eventsV1,analysis.(field).V1Trace,analysis.(field).A19Trace,metadata, 'peak', saveDirectory, 'V1');
doCrossCorrelation(analysis.(field).eventsA19,analysis.(field).A19Trace,analysis.(field).V1Trace,metadata, 'onset', saveDirectory, 'A19');
doCrossCorrelation(analysis.(field).eventsA19,analysis.(field).A19Trace,analysis.(field).V1Trace,metadata, 'peak', saveDirectory, 'A19');

%% 6.b) Do cross-correlogram for all local events that have a corresponding event in the other area
% for i = 1:length(windowTime)
doLocalCrossCorrelation(analysis.(field).eventsV1,analysis.(field).eventsA19,windowTime,metadata, 'onset', saveDirectory, 'V1'); %WORK IN PROGRESS
% end
%% get ative Framstes, compute correlations of the imaging stack and show it
% Computes correlations of the imaging stack (spontaneous, response, signal, or noise). 
            % The image dimensions can be "N" dimensions, but must be organized such that they are:
            %    *Spontaneous: Should be (x,y,t). The active frames are automatically extracted. 
            %    *Response:    Collapsed from a multidimensional array (i.e. [x,y,nCond,nTrials,t])
            %                  to (x,y,n), where n is all images
            %    *Signal:      Takes in a (x,y,nCond,nTrials,t), averages along the fourth dimension,
            %    *Noise        Computes correlations along nCond.
% [activeFrameStack,numberOfActiveEvents,eventOnset,eventDuration] = getActiveFrames(data.rawF,expParam);
% corrTable = computeCorrelationTable(activeFrameStack,expParam.ROI);
% showCorrelationStructure(corrTable,expParam, saveDirectory)