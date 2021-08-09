analysisParams.animal = 'F2359_2019-09-06';
analysisParams.name = 't00003';
analysisParams.sp2ID = 't00007';
analysisParams.expID = 't00003';
analysisParams.server = 0;
analysisParams.level = 0;
analysisParams.OffSet = 2;

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

if analysisParams.level
    tifDirectory = [base2pDirectory '\' analysisParams.name '\Registered\combined\'];
else
    tifDirectory = [base2pDirectory  '\' analysisParams.name '\Registered\slice1\'];
end

Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep 'ROIs_RF' filesep];
if ~exist(saveDirectory, 'dir')
    % make new file directory
    mkdir(saveDirectory); 
end


%% 1) Grab spike2 times

metadata.StimParams=LoadStimParamsRet(Sp2dDirectory);
metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);

%define referenceFrames = first frame after onset of stimulus
StimOnTimes = metadata.StimParams.StimOnTimes;
frameTimes = metadata.StimParams.frameTimes;
numStims = size(StimOnTimes,2);
referenceFrames_Start = zeros(numStims, 1);

for StimNr = 1:numStims
    referenceFrames_Start(StimNr) = find(frameTimes > StimOnTimes(2,StimNr),1);
end

%% 2) Make image of the stims presented

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
end

%this is your list of which image was presented; you can change this to
%upload the file you want
stimOrder= StimOnTimes(1,:);

%this pulls the stims in the correct order they were presented
stimsUsed=zeros(numStimX,numStimY,length(referenceFrames_Start));
for i = 1:length(referenceFrames_Start)
    stimsUsed(:,:,i)=double(stimImages(:,:,stimOrder(i)));
end

%% 3) read in the tiff files & ROIs
tifStack = [];
tifFiles = dir([tifDirectory '*.tif']);
for currentFile = 1:length(tifFiles)
    disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);
    
    % Specify stack name
    filePath = [tifDirectory tifFiles(currentFile).name];

    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
end

%save size of tifStack
sz=size(tifStack);

%make sure that imaging stack is at least as long as max of reference
%frames + offset
while sz(3) < max(referenceFrames_Start)+analysisParams.OffSet
    referenceFrames_Start = referenceFrames_Start(1:end-1);
    stimsUsed = stimsUsed(:,:,1:end-1);
end

FrameData=tifStack(:,:,referenceFrames_Start+analysisParams.OffSet);
clear tifStack;

%turn your framedata into a double to perform easier computations
FrameData = double(FrameData);

%load ROIs 
ROIs = LoadRoisS2p(analysisParams);
numROIs = length(ROIs.roi);

%then apply your ROI mask
for ROINr = 1:numROIs
    ROImask = zeros(sz(1),sz(2));
    for pixel = 1:size(ROIs.roi(ROINr).mask,1)
        ROImask(ROIs.roi(ROINr).mask(pixel,1),ROIs.roi(ROINr).mask(pixel,2)) = 1;
    end

%figure; imagesc maks{1} -> Control that the mask is working

%% 4. Calculate receptive field using reverse correlation for all ROIs

%finally, grab the frame you are using to indicate response, and take the
%largest value to use as a multiplier
    RFData = FrameData.*ROImask;
    multiplier=reshape(RFData,sz(1)*sz(2),length(referenceFrames_Start));
    multiplier=max(multiplier);

    %apply your multiplier to the stims and grab the average response
    stimsUsedReshaped=reshape(stimsUsed,[numStimX*numStimY,length(referenceFrames_Start)]);
    avgResp=stimsUsedReshaped.*multiplier;
    avgResp=reshape(avgResp,[numStimX,numStimY,length(referenceFrames_Start)]);
    figure; 
    imagesc(mean(avgResp,3))
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, ['ROI_Nr_' num2str(ROINr) '_ReceptiveField_.png']))
    close gcf
end