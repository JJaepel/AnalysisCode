function calculateSTA_Patches(analysisParams)
% analysisParams.animal = 'F2359_2019-09-06';
% analysisParams.name = 't00003';
% analysisParams.sp2ID = 't00007';
% analysisParams.expID = 't00003';
% analysisParams.level = 0;

%% 0) define folders and structures
if analysisParams.server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

base2pDirectory= [TwoPhontondir analysisParams.animal];
analysisParams.baseDirectory = base2pDirectory;

Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep 'ROIs_STA' filesep];
if ~exist(saveDirectory, 'dir')
    % make new file directory
    mkdir(saveDirectory); 
end


%% 1) Grab spike2 times

metadata.StimParams=LoadStimParamsRet(Sp2dDirectory);
metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);

%define referenceFrames = first frame after onset of stimulus
StimOnTimes = metadata.StimParams.StimOnTimes;
StimOffTimes = metadata.StimParams.StimOnTimes(2,:) + metadata.StimParams.stimDuration;
frameTimes = metadata.StimParams.frameTimes;
numStims = size(StimOnTimes,2);
referenceFrames= zeros(size(frameTimes,2), 1);

for StimNr = 1:numStims
    StartFrame = find(frameTimes > StimOnTimes(2,StimNr),1);
    StopFrame = find(frameTimes > StimOffTimes(StimNr),1);
    referenceFrames(StartFrame:StopFrame) = StimOnTimes(1,StimNr);
end

referenceFrames(referenceFrames == 0) = metadata.StimParams.uniqStims;

%% 2) read in the ROIs

%load ROIs 
ROIs = LoadRoisS2p(analysisParams);
numROIs = length(ROIs.roi);
twoPhotonFrames = length(ROIs.roi(1).spks);

%% 3) Make image of the stims presented

%create images presented
StimType = metadata.StimParams.type;
switch StimType
    case 'Retinotopy_2D'
        stimSize = strsplit(metadata.StimParams.stimSize(2:end-1), ',');
        numStimX = ceil((abs(metadata.StimParams.startPointx-metadata.StimParams.endPointx))/ (str2double(stimSize{1})/2));
        numStimY = ceil((abs(metadata.StimParams.startPointy-metadata.StimParams.endPointy))/ (str2double(stimSize{2})/2));

        numUniqStims = numStimX * numStimY;
        PosX = reshape(repmat(linspace(1, numStimX, numStimX),numStimY,1),numUniqStims,[])';
        PosY = repmat(linspace(1, numStimY, numStimY),1,numStimX);

        stimImages = zeros(numStimX, numStimY, numUniqStims+1);
        for Images = 1:numUniqStims
            stimImages(PosX(Images),PosY(Images),Images) = 1;
        end
    case 'Patch'
        numStimX = metadata.StimParams.numStimElev;
        numStimY = metadata.StimParams.numStimAzi;
        
        numUniqStims = numStimX * numStimY;
        PosX = reshape(repmat(linspace(1, numStimX, numStimX),numStimY,1),numUniqStims,[])';
        PosY = repmat(linspace(1, numStimY, numStimY),1,numStimX);
        
        stimImages = zeros(numStimX, numStimY, (2*numUniqStims)+1);
        %first do black patches
        for Images = 1:numUniqStims
            stimImages(PosX(Images),PosY(Images),Images) = 1;
        end
        
        %then do white patches
        for WhiteImages = 1:numUniqStims
            stimImages(PosX(WhiteImages),PosY(WhiteImages),numUniqStims+WhiteImages) = 1;
        end
    case 'PatchGrating'
        numStimX = metadata.StimParams.numStimElev;
        numStimY = metadata.StimParams.numStimAzi;

        numUniqStims = numStimX * numStimY;
        PosX = reshape(repmat(linspace(1, numStimX, numStimX),numStimY,1),numUniqStims,[])';
        PosY = repmat(linspace(1, numStimY, numStimY),1,numStimX);

        stimImages = zeros(numStimX, numStimY, numUniqStims+1);
        for Images = 1:numUniqStims
            stimImages(PosX(Images),PosY(Images),Images) = 1;
        end
        
end

%this is your list of which image was presented at each frame
numOffSetFrames = round(metadata.TwoPhoton.rate * analysisParams.OffSet);
OffSetFrames = zeros(1, numOffSetFrames);
OffSetFrames(OffSetFrames == 0) = metadata.StimParams.uniqStims;
stimOrder= [OffSetFrames'; referenceFrames];

%this pulls the stims in the correct order they were presented
stimsUsed=zeros(numStimX,numStimY,length(stimOrder));
for i = 1:length(stimOrder)
    stimsUsed(:,:,i)=double(stimImages(:,:,stimOrder(i)));
end


%% 4. Calculate receptive field using reverse correlation for all ROIs

for ROINr = 1:numROIs
    RFData = ROIs.roi(ROINr).spks;

    %make sure that stimsUsed and 2p data have the same amount of frames for the
    %calculation
    if length(stimOrder) == twoPhotonFrames
        stimFrames = twoPhotonFrames;
    else
        if length(stimOrder) > twoPhotonFrames
            stimOrder = stimOrder(1:twoPhotonFrames);
            stimsUsed = stimsUsed(:,:,1:twoPhotonFrames);
            stimFrames = twoPhotonFrames;
        else
            RFData = RFData(1,1:length(stimOrder));
            stimFrames = length(stimOrder);
        end
    end
    
    
    %two variants for multiplier: a) taking the response at it is or b)
    %every spk is the same value
    multiplier = RFData;
    multiplier(multiplier > 0) = 1; %this is weighting every spks the same amount = b)
    
    
    %apply your multiplier to the stims and grab the average response
    stimsUsedReshaped=reshape(stimsUsed,[numStimX*numStimY,stimFrames]);
    avgResp=stimsUsedReshaped.*multiplier;
    avgResp=reshape(avgResp,[numStimX,numStimY,stimFrames]);
    figure; 
    imagesc(mean(avgResp,3))
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, ['ROI_Nr_' num2str(ROINr) '_STA_Offset_' num2str(analysisParams.OffSet) ' s.png']))
    close gcf
end