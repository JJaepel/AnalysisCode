%Test script for online 2p viewing of orientaiton responses. This was reformated on
%9/1/2020 by AG for testing across the Code Review Group

%% set datapath
tic
baseDirectory  = 'Z:\Juliane\Data\2P_Data\';
Sp2Directory = 'Z:\Juliane\Data\Spike2Data\';
name = 'F2454_2020-09-03'
filenum_epi =38;
filenum_spk2 =38;

if(filenum_epi>99)
    dirName = [baseDirectory name '\t00' num2str(filenum_epi) '\'];
elseif filenum_epi>9
    dirName = [baseDirectory name '\t000' num2str(filenum_epi) '\'];
else 
    dirName = [baseDirectory name '\t0000' num2str(filenum_epi) '\'];
end%% Read imaging data
tifStack = [];
tifFiles = dir([dirName '*.tif']);
for(currentFile = 1:length(tifFiles))
    disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);
    % Specify stack name
    filePath = [dirName tifFiles(currentFile).name];
    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs2019(filePath,1,50));
end
%% read 2p times from spike2 file
fileForTriggers = 'twophotontimes.txt';
if(filenum_spk2>99)
    currentDirectory = [Sp2Directory name '\t00' num2str(filenum_spk2) '\'];
elseif filenum_spk2>9
    currentDirectory = [Sp2Directory name '\t000' num2str(filenum_spk2) '\'];
else 
    currentDirectory = [Sp2Directory name '\t0000' num2str(filenum_spk2) '\'];
end

% Read file with frame times (if present)
CCD_FramesTimes  = load([currentDirectory fileForTriggers]); %twophotontimes.txt
startTime        = CCD_FramesTimes(1);
offsetFrameTimes = CCD_FramesTimes-startTime;
meanCCDTime      = median(diff(offsetFrameTimes));

StimData = load([currentDirectory 'stimontimes.txt']);
stimConditionLabel     = StimData(1:2:end);
stimConditionOnsetTime = StimData(2:2:end);
stimConditionOnsetTime = stimConditionOnsetTime-startTime;

%% Stimulus Parameters
frameRate       = 1/meanCCDTime;         % Frame Rate of Camera (was fixed at 1/66E-3)
stimConditionOnsetTimeInFrames = ceil(stimConditionOnsetTime*frameRate); % use this one. not the one in the MAT file
duration=input('stim duration in seconds = ');
%% Make functional maps
stimulusDuration = floor(frameRate*duration);  % Now in frames
selectedStimulusCycle = [];
for(ii = 1:length(stimConditionOnsetTimeInFrames))
    selectedStimulusCycle = cat(2,selectedStimulusCycle,...
                                  stimConditionOnsetTimeInFrames(ii):(stimConditionOnsetTimeInFrames(ii)+stimulusDuration-1));
end

% Compute Baseline for Fluorescence Changes NOTE: this can be changed based
% on preference (i.e. here it uses minimum across tifStack)
baselineImage = double(min(tifStack,[],3));

% Estimate Functional Maps (Capable of OR/Dir Maps and ON/OFF)
currentSelectedTrial  = 1:sum(stimConditionLabel==1);
% Determine frames to assemble single-condition maps
trialLength = min(diff(stimConditionOnsetTimeInFrames(stimConditionLabel~=0)));
selectedFrames = []; 
numberOfConditions = max(stimConditionLabel);
for(currentSingleCondition = 1:numberOfConditions)
    selectedFrames{currentSingleCondition}.trialframes      = [];
    selectedFrames{currentSingleCondition}.trialStartFrames = [];
    trialLabels = find(stimConditionLabel == currentSingleCondition);
    if(min(currentSelectedTrial) > 0)
        trialLabels = trialLabels(currentSelectedTrial);
    end

    for(ii = trialLabels)
        putativeTrialFrames = stimConditionOnsetTimeInFrames(ii):(stimConditionOnsetTimeInFrames(ii)+floor(trialLength)-1); %floor(trialLength/(2-isOrientationData))
        framesPresentInSelectedStimulusCycle = max(repmat(putativeTrialFrames',[1 length(selectedStimulusCycle)])==repmat(selectedStimulusCycle,[length(putativeTrialFrames) 1]),[],2);
        putativeTrialFrames = putativeTrialFrames(framesPresentInSelectedStimulusCycle);

        selectedFrames{currentSingleCondition}.trialframes      = cat(2,selectedFrames{currentSingleCondition}.trialframes     ,putativeTrialFrames);
        selectedFrames{currentSingleCondition}.trialStartFrames = cat(2,selectedFrames{currentSingleCondition}.trialStartFrames,putativeTrialFrames(1));
    end
end

numberOfConditions  = length(selectedFrames)-mod(length(selectedFrames),2); orthogonalCondition = numberOfConditions/4; useBlankCondition   = true;

singleConditionMap = zeros(size(tifStack,1),size(tifStack,2),numberOfConditions);
for(currentSingleCondition = 1:numberOfConditions)
    currentSingleCondition
    % STEP 1 - Determine frames used for functional responses 
    currentSingleConditions = mod(currentSingleCondition+numberOfConditions/2+[0 numberOfConditions/2]-1,numberOfConditions)+1;
    referenceFrames         = cat(2,selectedFrames{currentSingleConditions(1)}.trialframes     ,selectedFrames{currentSingleConditions(2)}.trialframes     );
    referenceFrames_Start   = cat(2,selectedFrames{currentSingleConditions(1)}.trialStartFrames,selectedFrames{currentSingleConditions(2)}.trialStartFrames);

    % STEP 2 - Construct Blank-Subtracted Single-Condition Maps
    preTrialFrames       = []; 
    numberOfTrials       = length(referenceFrames_Start);

    % find frames in the specified pre-stimulus period
    for(ii = 1:numberOfTrials)
        preTrialFrames  = cat(2,preTrialFrames ,referenceFrames_Start(ii) -(round(frameRate*0.5):-1:1));
    end
    preTrialFrames(preTrialFrames>5385)=[];

    % compute mean fluorescence for the specified pre-stimulus period
    firstFrame      = mean(tifStack(:,:,preTrialFrames),3);
    
    % STEP 3 - Compute Functional Maps
    % Create Delta F / F Maps For Single Condition Maps
    currentMap = (double(mean(tifStack(:,:,referenceFrames) ,3))-firstFrame)./baselineImage; 
    currentMap = (double(mean(tifStack(:,:,referenceFrames) ,3))-baselineImage)./baselineImage; 

    % Construct Blank-Subtracted Single-Condition Maps
    singleConditionMap(:,:,currentSingleCondition) = currentMap;
end

%% construct figure
orientationMap = vectorSum(singleConditionMap,2);
fixedRange     = 0.075/2;
figure; imagesc(polarMap(orientationMap,hsv,[0 2*fixedRange]));
toc
