function calcSignificantPatches(analysisParams)

analysisParams.animal = 'F2526_2021-05-26';
analysisParams.name = 't00018';
analysisParams.sp2ID = 't00018';
analysisParams.expID = 't00018';
analysisParams.server = 1;
analysisParams.level = 0;
analysisParams.cvsFile = 'SparseNoise_17_06_47.cvs';

analysisParams.smoothing = 3;
analysisParams.StimOnSetDelay = 0.1;
analysisParams.StimOvershoot = 0.2;
analysisParams.threshold =2;
analysisParams.windowLength = 5;  %how many times longer than the stimperiod
analysisParams.p = 0.05;
analysisParams.reloadData = 0;
analysisParams.dataType = 1;
analysisParams.field = 'dff';

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
    ROIDirectory = [base2pDirectory '\' analysisParams.name '\Registered\combined\'];
else
    ROIDirectory = [base2pDirectory  '\' analysisParams.name '\Registered\slice1\'];
end

Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep];
figSaveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep 'cleanRF'];
figSaveDirectory2 = [savedir analysisParams.animal filesep analysisParams.expID filesep 'RFs_lat'];
figSaveDirectory3 = [savedir analysisParams.animal filesep analysisParams.expID filesep 'Summary'];
if ~exist(saveDirectory, 'dir')
    % make new file directory
    mkdir(saveDirectory); 
end
if ~exist(figSaveDirectory, 'dir')
    % make new file directory
    mkdir(figSaveDirectory); 
else
    % delete files from older runs
end
if ~exist(figSaveDirectory2, 'dir')
    % make new file directory
    mkdir(figSaveDirectory2); 
end
if ~exist(figSaveDirectory3, 'dir')
    % make new file directory
    mkdir(figSaveDirectory3); 
end
analysis = struct;

%% 1) Grab spike2 times & ROIs

%load metadata
metadata.StimParams=LoadStimParamsRet(Sp2dDirectory);
metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);

disp('Loading data')
if analysisParams.reloadData
    % make sure that you also do analysis afterwards
    analysisParams.reanalyse = 1;
    
    %load stimulus parameters and spike2data
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParamsRet(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    
    Stimtype = metadata.StimParams.type;
    switch Stimtype
        case 'SparseNoise'
            metadata.StimParams.numTrials = metadata.StimParams.trial;
        case 'rotatingGratingPatch'
            metadata.StimParams.stimDuration=metadata.StimParams.spinPeriodInSeconds;
    end
    
    % load tif data
    switch analysisParams.dataType
        case 1
            data = LoadRoisS2p(analysisParams);
            
            %do baseline filtering and compute dff
            if strcmp(analysisParams.field, 'dff')
                [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                data = computeDff(data, 'rawF', 'baseline', 'dff');
            end
        case 2
            if analysisParams.makeROIs
                if analysisParams.level == 1
                    Suite2pAxonTifCombiner(analysisParams)
                end
                data = Suite2pAxonExtractorFct(analysisParams);
            else
                if analysisParams.level == 1

                    dataFile = ([base2pDirectory filesep analysisParams.expID filesep 'suite2p\combined\data.mat']);
                else
                    dataFile = ([base2pDirectory filesep analysisParams.expID filesep 'suite2p\plane0\data.mat']);
                end
                File = load(dataFile, 'data');
                data = File.data;
                clear File
            end
            %do baseline filtering and compute dff
            if strcmp(analysisParams.field, 'dff')
                [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                data = computeDff(data, 'rawF', 'baseline', 'dff');
            end
            
            
        case 3 
            %do spine data reloading here
            data = SpineROIExtractFct(analysisParams);
            data = computeDffSpines(data);
            structureType = analysisParams.region;
            switch structureType
                case 'cells'
                    data = FitCellsToDendrite(analysisParams,metadata, data);
                case 'dendrite'
                    data = computeDendriticSubstraction(analysisParams, metadata, data);    
            end
            analysisParams.field = 'rawRes';     
    end
    save(fullfile(saveDirectory, 'Latency.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
else
    load(fullfile(saveDirectory, 'Latency.mat'), 'data', 'metadata', 'analysis');
end


% %load ROIs
% try
%   data = LoadRoisS2p(analysisParams);
% catch
%     load([base2pDirectory '\' analysisParams.name '\suite2p\combined\data.mat']);
%     ROIs = data;
% end

NumROIs = length(data.roi);

%% 2) Calculate analysis settings

stimType = metadata.StimParams.type;
%create stimulus settings for each patch
switch stimType
    case 'Patch'
        numTrials = metadata.StimParams.numTrials;
        numPatches = metadata.StimParams.numStimElev * metadata.StimParams.numStimAzi * 2;
        numCol = 2;
        numUniqStims = metadata.StimParams.numStimElev * metadata.StimParams.numStimAzi;
        metadata.StimParams.minAzi = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.maxAzi = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) - metadata.StimParams.stimSize(1)/2;
        %metadata.StimParams.PosX = repmat(linspace(minAzi, maxAzi, metadata.StimParams.numStimAzi), 1, metadata.StimParams.numStimElev);
        numX = repmat(linspace(1, metadata.StimParams.numStimAzi, metadata.StimParams.numStimAzi), 1, metadata.StimParams.numStimElev);
        metadata.StimParams.minElev = metadata.StimParams.centerPoint(2) - ((metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2))/2) + metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.maxElev = metadata.StimParams.centerPoint(2) + ((metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        numY = reshape(repmat(linspace(1, metadata.StimParams.numStimElev, metadata.StimParams.numStimElev), metadata.StimParams.numStimAzi, 1), numUniqStims, [])';
    case 'Retinotopy_2D'
        numTrials = metadata.StimParams.numberOfTrials;
        numCol = 1;
        stimSize = strsplit(metadata.StimParams.stimSize(2:end-1), ',');
        metadata.StimParams.stimSize = [str2double(stimSize{1}), str2double(stimSize{2})];
        if analysisParams.special
            metadata.StimParams.stimSize(1) = metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.stimSize(2) = metadata.StimParams.stimSize(2)/2;
        end
        if mod((metadata.StimParams.endPointy - metadata.StimParams.startPointy), (str2double(metadata.StimParams.stimSize(2)))) == 0
            metadata.StimParams.numStimElev = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(metadata.StimParams.stimSize(2)));
            metadata.StimParams.numStimAzi = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(metadata.StimParams.stimSize(1)));
        else
            metadata.StimParams.numStimElev = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(metadata.StimParams.stimSize(2)))+1;
            metadata.StimParams.numStimAzi = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(metadata.StimParams.stimSize(1)))+1;
        end
        metadata.StimParams.isi= metadata.StimParams.ISI;
        
        numPatches = metadata.StimParams.numStimElev * metadata.StimParams.numStimAzi;
        metadata.StimParams.minAzi = metadata.StimParams.startPointx;
        metadata.StimParams.maxElev = metadata.StimParams.startPointy;
        if min(metadata.StimParams.stimSize) < 2
            metadata.StimParams.maxAzi = floor(metadata.StimParams.startPointx + (metadata.StimParams.numStimAzi) * metadata.StimParams.stimSize(1) - 0.5);
            metadata.StimParams.minElev = floor(metadata.StimParams.startPointy + (metadata.StimParams.numStimElev) * metadata.StimParams.stimSize(2) - 0.5);
        else
            metadata.StimParams.maxAzi = floor(metadata.StimParams.startPointx + (metadata.StimParams.numStimAzi) * metadata.StimParams.stimSize(1) - 1);
            metadata.StimParams.minElev = floor(metadata.StimParams.startPointy + (metadata.StimParams.numStimElev) * metadata.StimParams.stimSize(2) - 1);
        end
        numX = repmat(linspace(1, metadata.StimParams.numStimAzi, metadata.StimParams.numStimAzi), 1, metadata.StimParams.numStimElev);
        numY = reshape(repmat(linspace(1, metadata.StimParams.numStimElev, metadata.StimParams.numStimElev), metadata.StimParams.numStimAzi, 1), numPatches, [])';
    case 'rotatingGratingPatch'
        numTrials = metadata.StimParams.numTrials;
        numPatches = metadata.StimParams.numStimElev * metadata.StimParams.numStimAzi;
        numCol = 1;
        numUniqStims = metadata.StimParams.numStimElev * metadata.StimParams.numStimAzi;
        metadata.StimParams.minAzi = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.maxAzi = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) - metadata.StimParams.stimSize(1)/2;
        numX = repmat(linspace(1, metadata.StimParams.numStimAzi, metadata.StimParams.numStimAzi), 1, metadata.StimParams.numStimElev);
        metadata.StimParams.minElev = metadata.StimParams.centerPoint(2) - ((metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2))/2) + metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.maxElev = metadata.StimParams.centerPoint(2) + ((metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        numY = reshape(repmat(linspace(1, metadata.StimParams.numStimElev, metadata.StimParams.numStimElev), metadata.StimParams.numStimAzi, 1), numUniqStims, [])';
        metadata.StimParams.stimDuration = metadata.StimParams.spinPeriodInSeconds;
    case 'SparseNoise'
        numCol = 2;
        numPatches = 25*2;
        metadata.StimParams.stimSize = [5 5]; metadata.StimParams.numStimElev = 5; metadata.StimParams.numStimAzi = 5; %Might have to change this!!!
        metadata.StimParams.minAzi = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.maxAzi = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) - metadata.StimParams.stimSize(1)/2;
        numX = repmat(linspace(1, metadata.StimParams.numStimAzi, metadata.StimParams.numStimAzi), 1, metadata.StimParams.numStimElev);
        metadata.StimParams.minElev = metadata.StimParams.centerPoint(2) - ((metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2))/2) + metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.maxElev = metadata.StimParams.centerPoint(2) + ((metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        numY = reshape(repmat(linspace(1, metadata.StimParams.numStimElev, metadata.StimParams.numStimElev), metadata.StimParams.numStimAzi, 1), numPatches/2, [])';
        cvsFile = [Sp2dir analysisParams.animal filesep analysisParams.cvsFile];
        stimOrder = csvread(cvsFile);
        if size(metadata.StimParams.StimOnTimes,2) < stimOrder
            missingTriggerSites = find(diff(metadata.StimParams.StimOnTimes(2,:)) > 0.65);
            addIns = 0;
            medianTriggerTime = median(diff(metadata.StimParams.StimOnTimes(2,:)));
            for m = 1:length(missingTriggerSites)
                lastTrigger = missingTriggerSites(m)+addIns;
                timeDiff = metadata.StimParams.StimOnTimes(2,lastTrigger+1)-metadata.StimParams.StimOnTimes(2,lastTrigger);
                triggersToAdd = round(timeDiff/medianTriggerTime)-1;
                matrixToAdd = ones(2,triggersToAdd);
                for t=1:triggersToAdd
                    matrixToAdd(2,t)=metadata.StimParams.StimOnTimes(2,lastTrigger)+t*medianTriggerTime;
                end
                metadata.StimParams.StimOnTimes = [metadata.StimParams.StimOnTimes(:,1:lastTrigger) matrixToAdd metadata.StimParams.StimOnTimes(:,lastTrigger+1:end)];
                addIns = addIns + tiggersToAdd;
            end
        end
        timeStampSplit = split(cvsFile,'Noise_');
        timeStamp = timeStampSplit{2};
        StimInfoFile = [Sp2dir analysisParams.animal filesep 'StimulusInformation_' timeStamp(1:end-4) '.mat'];
        imageInformation = load(StimInfoFile);
        %imageInformation = load([Sp2dir analysisParams.animal filesep 'StimulusInformation.mat']);
        stimArray = reshape(imageInformation.allUsedStims, 4, size(imageInformation.selectedFrames,2));
        orderedStimArray = stimArray(:,stimOrder);
        allStims = reshape(orderedStimArray,1,size(orderedStimArray,1)*size(orderedStimArray,2));
        for P = 1:numPatches
            temp= ceil(find(allStims == P)./4);
            stimTimesAll{P} = temp(temp <= size(metadata.StimParams.StimOnTimes,2));
        end
        numTrials = min(cellfun('size', stimTimesAll, 2));
end
%create stimulus settings for each patch
PatchNr = 1;
for C = 1:numCol
    for P =1:numPatches/numCol
        metadata.StimParams.Stimulus(PatchNr).X = numX(P);
        metadata.StimParams.Stimulus(PatchNr).Y = numY(P);
        metadata.StimParams.Stimulus(PatchNr).Color = C;
        PatchNr = PatchNr + 1;
    end
end

if analysisParams.level
    metadata.TwoPhoton.time = metadata.TwoPhoton.time(1:5:end);
end

%write StimOnSettings per Patch
for P = 1:numPatches
    switch stimType
        case 'SparseNoise'
            stimTimesNum = stimTimesAll{P};
            stimTimesNum = stimTimesNum(randperm(length(stimTimesNum)));
            stimTimesNum = stimTimesNum(1:numTrials);
        otherwise
            stimTimesNum = find(metadata.StimParams.StimOnTimes(1,:) == P);
    end
    StimOnSetTimes = metadata.StimParams.StimOnTimes(2,stimTimesNum);
    for i = 1:length(StimOnSetTimes)
        StartFrames(i) = find(metadata.TwoPhoton.time > StimOnSetTimes(i), 1);
        StopFrames(i) = find(metadata.TwoPhoton.time > (StimOnSetTimes(i)+ metadata.StimParams.stimDuration),1);
    end
    StimOnSetFrames{P} = StartFrames;
    StimOffSetFrames{P} = StopFrames;
end
SamplingFreq = metadata.TwoPhoton.rate;
StimDur_corr=round(SamplingFreq*metadata.StimParams.stimDuration);

%calculate windows
BaselineWindow = round((-0.7 * metadata.StimParams.isi) * SamplingFreq):round(0.5 * analysisParams.StimOnSetDelay * SamplingFreq);
TestingWindow = round(analysisParams.StimOnSetDelay * SamplingFreq):round((metadata.StimParams.stimDuration + analysisParams.StimOvershoot) * SamplingFreq);
DisplayWindow = round((-0.7 * metadata.StimParams.isi) * SamplingFreq):round((metadata.StimParams.stimDuration + (analysisParams.windowLength * metadata.StimParams.isi)) * SamplingFreq);

disp(' ');
disp(['Baseline period  = ' sprintf('%5.2f',min(BaselineWindow)./SamplingFreq) ' to ' sprintf('%5.2f',max(BaselineWindow)./SamplingFreq) ...
    's (frames: ' num2str(min(BaselineWindow)) ' .. ' num2str(max(BaselineWindow)) ')']);
disp(' ');
disp(['Display period   = ' sprintf('%5.2f',min(DisplayWindow)./SamplingFreq) ' to ' sprintf('%5.2f',max(DisplayWindow)./SamplingFreq) ...
    's (frames: ' num2str(min(DisplayWindow)) ' .. ' num2str(max(DisplayWindow)) ')']);
disp(['Testing period   = ' sprintf('%5.2f',min(TestingWindow)./SamplingFreq) ' to ' sprintf('%5.2f',max(TestingWindow)./SamplingFreq) ...
    's (frames: ' num2str(min(TestingWindow)) ' .. ' num2str(max(TestingWindow)) ')']);
disp(' ');

%% 3) Loop ROIs and get PSTH's and tuning curves

for nr = 1:NumROIs
    %define F
    raw_data = data.roi(nr).rawF;
    F = raw_data;
    F = smooth(F, analysisParams.smoothing);
    clear raw_data;
    
    BTrace = []; WTrace = [];
    
    %get tuning curves
    for P = 1:numPatches      
       Y = metadata.StimParams.Stimulus(P).Y;
       X = metadata.StimParams.Stimulus(P).X;
       C = metadata.StimParams.Stimulus(P).Color;    
       for t = 1:numTrials-1
           if length(StimOnSetFrames{P}) < t
               disp(['Warning: Patch (x =' num2str(X) ', y = ' num2str(Y) ', trial ' num2str(t) ') was not recorded'])
               if C == 1
                   BlackTraces{nr, Y, X}(t,:) = NaN;
                   BlackAvrTraces(nr, Y, X, t) = NaN;
                   TC1DBlack(nr, ((Y-1)*numTrials)+t, X) = NaN;
               else
                   WhiteTraces{nr, Y, X}(t,:) = NaN;
                   WhiteAvrTraces(nr, Y, X, t) = NaN;
                   TC1DWhite(nr, ((Y-1)*numTrials)+t, X) = NaN;
               end
           else
               DISPix = StimOnSetFrames{P}(t)+DisplayWindow;
               TESTix = StimOnSetFrames{P}(t)+TestingWindow;
               if C == 1
                    BlackTraces{nr,Y,X}(t,:) =F(DISPix );
                    temp=F(TESTix);
                    BTest(nr,Y,X,t,:)=temp;
                    BlackAvrTraces(nr,Y,X,t)=nanmean(temp);
                    temp2=F(DISPix);
                    TC1DBlack(nr, ((Y-1)*numTrials)+t, X) =nanmean(temp2);
               else
                    WhiteTraces{nr,Y,X}(t,:) =F(DISPix );
                    temp=F(TESTix);
                    WTest(nr,Y,X,t,:)=temp;
                    WhiteAvrTraces(nr,Y,X,t)=nanmean(temp);
                    temp2=F(DISPix);
                    TC1DWhite(nr, ((Y-1)*numTrials)+t, X) = nanmean(temp2);
               end
           end
       end       
       if C == 1
           temp3=(squeeze(nanmean(BTest(nr,Y,X,:,:),4)))';
           BTrace=[BTrace temp3];
       else
           temp4=(squeeze(nanmean(WTest(nr,Y,X,:,:),4)))';
           WTrace=[WTrace temp4];
       end      
    end
    BTestTrace(nr,:)=BTrace;
    try
        WTestTrace(nr,:)=WTrace;
    catch
        WTestTrace(nr,:) = NaN;
    end
end

AverageB= median(BTestTrace,1);
AverageW =  median(WTestTrace,1);

save([saveDirectory 'TestTraces.mat'],'BTestTrace','WTestTrace','AverageB','AverageW');
try
    save([saveDirectory 'TuningCurves.mat'], 'WhiteTraces','BlackTraces','TC1DBlack','TC1DWhite','WhiteAvrTraces','BlackAvrTraces',...
        'SamplingFreq', 'DisplayWindow', 'TestingWindow');
catch
    save([saveDirectory 'TuningCurves.mat'], 'BlackTraces','TC1DBlack','BlackAvrTraces',...
        'SamplingFreq', 'DisplayWindow', 'TestingWindow')
end

%% 4) Get receptive fields

for ind = 1:NumROIs
    
    %it will analyze from 1x(length 1 stim period) frames before stim to 5x(length 1 stim period) after start
    count = 1;
    for xnd = 1:metadata.StimParams.numStimAzi
        for ynd = 1:metadata.StimParams.numStimElev
            switch stimType
                case 'Patch'
                    ONdF{ind}(:,ynd,xnd)=nanmean(WhiteTraces{ind,ynd,xnd}(:,:)); %average over repeats of those frames in the display window shape: frames, y,x.
                case 'SparseNoise'
                    ONdF{ind}(:,ynd,xnd)=nanmean(WhiteTraces{ind,ynd,xnd}(:,:));
            end
            OFFdF{ind}(:,ynd,xnd)=nanmean(BlackTraces{ind,ynd,xnd}(:,:));%variable Black_Traces dF/F Johannes way
            for knd = 1:length(DisplayWindow)
                switch stimType
                    case 'Patch'
                        ATracesON{ind}(knd,count,:)=WhiteTraces{ind,ynd,xnd}(:,knd); %size this way will be Reps , PatchX,PatchY,DisplayWindow
                    case 'SparseNoise'
                         ATracesON{ind}(knd,count,:)=WhiteTraces{ind,ynd,xnd}(:,knd);
                end
                ATracesOFF{ind}(knd,count,:)=BlackTraces{ind,ynd,xnd}(:,knd);
            end
            count = count + 1;
        end
    end
    
    %create data matrix for heatmap
    switch stimType
        case 'Patch'
            Group_meanspatches_ON(ind,:)= nanmean(nanmean(ATracesON{ind},3),1); %mean of activity per patch over whole display window and all reps
            analysis.roi(ind).Group_meanspatches_ON=Group_meanspatches_ON(ind,:);
        case 'SparseNoise'
            Group_meanspatches_ON(ind,:)= nanmean(nanmean(ATracesON{ind},3),1); %mean of activity per patch over whole display window and all reps
            analysis.roi(ind).Group_meanspatches_ON=Group_meanspatches_ON(ind,:);
    end
    Group_meanspatches_OFF(ind,:)= nanmean(nanmean(ATracesOFF{ind},3),1);
    analysis.roi(ind).Group_meanspatches_OFF=Group_meanspatches_OFF(ind,:);
    
    % do ANOVA and Kruskal-Wallis test to assess whether response is significant at each frame over all repeats and all positions
    for knd=1:length(DisplayWindow)
        switch stimType 
            case 'Patch'
                Significance(1,knd,ind)=anova1(squeeze(ATracesON{ind}(knd,:,:))',[],'off');
                Significance(3,knd,ind)=kruskalwallis(squeeze(ATracesON{ind}(knd,:,:))',[],'off');
            case 'SparseNoise'
                Significance(1,knd,ind)=anova1(squeeze(ATracesON{ind}(knd,:,:))',[],'off');
                Significance(3,knd,ind)=kruskalwallis(squeeze(ATracesON{ind}(knd,:,:))',[],'off');
            otherwise
                Significance(1,knd,ind) = NaN;
                Significance(3,knd,ind) = NaN;
        end
        Significance(2,knd,ind)=anova1(squeeze(ATracesOFF{ind}(knd,:,:))',[],'off');
        Significance(4,knd,ind)=kruskalwallis(squeeze(ATracesOFF{ind}(knd,:,:))',[],'off');
    end
    
    %find latency of the signal
    switch stimType
        case 'SparseNoise'
            [Lat,LatON,LatOFF]=FindLatency(StimDur_corr,Significance,DisplayWindow,ind, analysisParams.p);
        case 'rotatingGratingPatch'
            [Lat,LatON,LatOFF]=FindLatencyRotPatch(StimDur_corr,Significance,DisplayWindow,analysisParams.p);
        otherwise
            [Lat,LatON,LatOFF]=FindLatency(StimDur_corr,Significance,DisplayWindow,analysisParams.p);
    end
    
    %calculate responsive cells stats over patch positions over the response positive window
    switch stimType 
        case 'Patch'
            TracesON{ind}=median(ATracesON{ind},3);
        case 'SparseNoise'
            TracesON{ind}=median(ATracesON{ind},3);
    end
    TracesOFF{ind}=median(ATracesOFF{ind},3);
    
    if isnan(LatON)==0 && isnan(LatOFF) == 0
        try
            pBlack=anova2(TracesOFF{ind}(LatOFF:(LatOFF+StimDur_corr)-1,:),1,'off');
        catch
            pBlack=anova2(TracesOFF{ind}(LatOFF:end,:),1,'off');
        end
        try
            pWhite=anova2(TracesON{ind}(LatON:(LatON+StimDur_corr)-1,:),1,'off');
        catch
            pWhite=anova2(TracesON{ind}(LatON:end,:),1,'off');
        end
        analysis.roi(ind).pWhite= min(pWhite);
        analysis.roi(ind).pBlack =min(pBlack);
        [thzRFtempONci]= CalcRFsSigLocations(ONdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        [thzRFtempOFFci]= CalcRFsSigLocations(OFFdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        
        %to calculate fields at latency
        try
            if isempty(thzRFtempONci)==0 && nnz(thzRFtempONci(:,:,LatON))~=0 
                ON_field= thzRFtempONci(:,:,LatON)+thzRFtempONci(:,:,LatON+1)+thzRFtempONci(:,:,LatON+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanON = PlotCleanRF(ON_field,metadata,analysisParams.threshold,ind, figSaveDirectory);
                ONout1 = bwperim(cleanON);
                [Ony,Onx]=find(ONout1);
                [ONperimeter,ONx,ONy,ONsize]= define_fields(ONout1,metadata);
                analysis.roi(ind).cleanON=cleanON;
                analysis.roi(ind).ONperimeter=ONperimeter;
                analysis.roi(ind).ONx=ONx;
                analysis.roi(ind).ONy=ONy;
                analysis.roi(ind).ONsize=ONsize;
            end
        catch
            disp('Latency outside analysis window')
        end
        
        try
            if  isempty(thzRFtempOFFci)==0 && nnz(thzRFtempOFFci(:,:,LatOFF))~=0
                OFF_field= thzRFtempOFFci(:,:,LatOFF)+thzRFtempOFFci(:,:,LatOFF+1)+thzRFtempOFFci(:,:,LatOFF+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanOFF=PlotCleanRF(OFF_field,metadata,analysisParams.threshold,ind, figSaveDirectory);
                %calculate perimeter
                OFFout1 = bwperim(cleanOFF);
                [Offy,Offx]=find(OFFout1);
                [OFFperimeter,OFFx,OFFy,OFFsize]= define_fields(OFFout1,metadata);
                analysis.roi(ind).cleanOFF=cleanOFF;
                analysis.roi(ind).OFFperimeter=OFFperimeter;
                analysis.roi(ind).OFFx=OFFx;
                analysis.roi(ind).OFFy=OFFy;
                analysis.roi(ind).OFFsize=OFFsize;
            end
        catch
            disp('Latency outside analysis window')
        end
        
        if isempty(cleanON) == 0 && isempty(cleanOFF == 0)
            %calculate both fields
            if max(max(thzRFtempOFFci(:,:,LatOFF)))>max(max(thzRFtempONci(:,:,LatON)))
                ONOFF=smooth2a(imabsdiff(cleanON,cleanOFF),0,0);
            else
                ONOFF=smooth2a(cleanON-cleanOFF,0,0); %subtracts the ON map matrix from the OFF one
            end
        elseif isempty(cleanON)==0
            ONOFF=smooth2a(cleanON,0,0);
        elseif  isempty(cleanOFF)==0
            ONOFF=smooth2a(cleanOFF,0,0);
        end
        
        smONOFF=smooth2a(ONOFF,5,5);
        smthONOFF=smONOFF;
        smthONOFF(smONOFF==0)=min(smONOFF(1:end))+(range(smONOFF(1:end))/2);
        %plot both fields, perimeter and statistics
        PlotResultingRFs(smthONOFF,Onx,Ony,Offx,Offy,Significance,DisplayWindow,metadata,ind,figSaveDirectory2);
        
        analysis.roi(ind).cleanRF=smONOFF;
        analysis.roi(ind).Latencies=[Lat LatON LatOFF];
        
    elseif isnan(LatON) == 0
        try
            pWhite=anova2(TracesON{ind}(LatON:(LatON+StimDur_corr)-1,:),1,'off');
        catch
            pWhite=anova2(TracesON{ind}(LatON:end,:),1,'off');
        end
        analysis.roi(ind).pWhite= min(pWhite);
        [thzRFtempONci]=CalcRFsSigLocations(ONdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        try
            if isempty(thzRFtempONci)|| nnz(thzRFtempONci(:,:,LatON))==0
                disp(['For ROI' num2str(ind) ': Empty RF matrix at latency'])
            else
                ON_field= thzRFtempONci(:,:,LatON)+thzRFtempONci(:,:,LatON+1)+thzRFtempONci(:,:,LatON+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanON = PlotCleanRF(ON_field,metadata,analysisParams.threshold, ind, figSaveDirectory);
                ONout1 = bwperim(cleanON);
                [Ony,Onx]=find(ONout1);
                [ONperimeter,ONx,ONy,ONsize]= define_fields(ONout1,metadata);
                analysis.roi(ind).cleanON=cleanON;
                analysis.roi(ind).Latencies = LatON;
                analysis.roi(ind).ONperimeter=ONperimeter;
                analysis.roi(ind).ONx=ONx;
                analysis.roi(ind).ONy=ONy;
                analysis.roi(ind).ONsize=ONsize;
                %plot field, perimeter and statistics
                PlotSingleRFs(cleanON,Significance,DisplayWindow,metadata,ind,figSaveDirectory2,1, Onx, Ony) 
            end
        catch
            disp('Latency outside analysis window')
        end
        
    elseif isnan(LatOFF) == 0
        try
            pBlack=anova2(TracesOFF{ind}(LatOFF:(LatOFF+StimDur_corr)-1,:),1,'off');
        catch
            pBlack=anova2(TracesOFF{ind}(LatOFF:end,:),1,'off');
        end
        analysis.roi(ind).pBlack =min(pBlack);
        [thzRFtempOFFci]=CalcRFsSigLocations(OFFdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        try
            if isempty(thzRFtempOFFci)|| nnz(thzRFtempOFFci(:,:,LatOFF))==0
                disp(['For ROI' num2str(ind) ': Empty RF matrix at latency'])
            else
                OFF_field= thzRFtempOFFci(:,:,LatOFF)+thzRFtempOFFci(:,:,LatOFF+1)+thzRFtempOFFci(:,:,LatOFF+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanOFF=PlotCleanRF(OFF_field,metadata,analysisParams.threshold, ind, figSaveDirectory);
                %calculate perimeter
                OFFout1 = bwperim(cleanOFF);
                [Offy,Offx]=find(OFFout1);
                [OFFperimeter,OFFx,OFFy,OFFsize]= define_fields(OFFout1,metadata);
                analysis.roi(ind).cleanOFF=cleanOFF;
                analysis.roi(ind).Latencies = LatOFF;
                analysis.roi(ind).OFFperimeter=OFFperimeter;
                analysis.roi(ind).OFFx=OFFx;
                analysis.roi(ind).OFFy=OFFy;
                analysis.roi(ind).OFFsize=OFFsize;
                %plot field, perimeter and statistics
                PlotSingleRFs(cleanOFF,Significance,DisplayWindow,metadata,ind,figSaveDirectory2, 2, Offx, Offy)            
            end
        catch
            disp('Latency outside analysis window')
        end
    end
    
    analysis.roi(ind).Significance = Significance;
    
end
switch stimType 
    case 'Patch'
        save([savedir 'Group_meanspatches.mat'], 'Group_meanspatches_ON', 'Group_meanspatches_OFF');
        save([savedir 'Atraces.mat'], 'ATracesON', 'ATracesOFF');
    case 'SparseNoise'
        save([savedir 'Group_meanspatches.mat'], 'Group_meanspatches_ON', 'Group_meanspatches_OFF');
        save([savedir 'Atraces.mat'], 'ATracesON', 'ATracesOFF');
    otherwise
        save([savedir 'Group_meanspatches.mat'], 'Group_meanspatches_OFF');
        save([savedir 'Atraces.mat'], 'ATracesOFF');
end
%save([savedir 'Patches.mat'], 'analysis');

save([saveDirectory 'Patches.mat'],'analysis');


%% 5) Make summary data
% calculate totals
All_OFFx=calculateTotals(analysis.roi,NumROIs,'OFFx');
All_OFFy=calculateTotals(analysis.roi,NumROIs,'OFFy');
All_ONx=calculateTotals(analysis.roi,NumROIs,'ONx');
All_ONy=calculateTotals(analysis.roi,NumROIs,'ONy');
All_ONsize=calculateTotals(analysis.roi,NumROIs,'ONsize');
All_OFFsize=calculateTotals(analysis.roi,NumROIs,'OFFsize');

x1 = 0; x2 = 0; y1 = 0; y2 = 0;

% plot field centers
figure;    hold on;
line([x1 x2],[y1 y1],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
line([x1 x2],[y2 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
line([x1 x1],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
line([x2 x2],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
switch stimType 
    case 'Patch'
        try
            if true((All_ONy~=0)==(All_ONx~=0))
                plot(All_ONx(All_ONx~=0),All_ONy(All_ONy~=0),'o','Color',[0.7 0 0]);    
            end
        end
    case 'SparseNoise'
        try
            if true((All_ONy~=0)==(All_ONx~=0))
                plot(All_ONx(All_ONx~=0),All_ONy(All_ONy~=0),'o','Color',[0.7 0 0]);    
            end
        end
end
try
    if true((All_OFFy~=0)==(All_OFFx~=0))
        plot(All_OFFx(All_OFFx~=0),All_OFFy(All_OFFy~=0),'x','Color',[0 0 0.7]);
    end
end
box on;
set(gca,'YDir','normal')
set(gca,'TickDir','out');
set(gca,'TickLen',[0.01 0.01]);
xlim([metadata.StimParams.minAzi metadata.StimParams.maxAzi])
try
    ylim([metadata.StimParams.minElev metadata.StimParams.maxElev])
catch
    ylim([metadata.StimParams.minElev metadata.StimParams.maxElev])
end
title('On and OFF field centers (red=ON, blue=OFF)');
saveas(gcf, fullfile(figSaveDirectory3, 'On and OFF field centers.png'))
pause(0.1);
%close(gcf);

% plot field distribution
figure;
if isfield(analysis.roi,'ONperimeter')
    subplot(1,2,1); hold on;
    line([x1 x2],[y1 y1],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x2],[y2 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x1],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x2 x2],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    % draw outline of ONfields
    Color = [0.7 0 0];
    LineWidth = 1;
    for nr = 1:NumROIs
        if  isempty(analysis.roi(nr).ONperimeter)==0
            for p = 1:size(analysis.roi(nr).ONperimeter, 1 )-1
                line(analysis.roi(nr).ONperimeter(p:p+1, 1), analysis.roi(nr).ONperimeter(p:p+1, 2), ...
                    'Color', Color, 'LineStyle', '-', 'LineWidth', LineWidth );
            end
            line( [ analysis.roi(nr).ONperimeter(1,   1) analysis.roi(nr).ONperimeter(end, 1) ], [ analysis.roi(nr).ONperimeter(1,   2) analysis.roi(nr).ONperimeter(end, 2) ], ...
                'Color',  Color, 'LineStyle', '-', 'LineWidth', LineWidth );
        end
    end
    box on;axis equal
    set(gca,'YDir','normal')
    set(gca,'TickDir','out');
    set(gca,'TickLen',[0.01 0.01]);
    xlim([metadata.StimParams.minAzi metadata.StimParams.maxAzi])
    try
        ylim([metadata.StimParams.minElev metadata.StimParams.maxElev])
    catch
        ylim([metadata.StimParams.minElev metadata.StimParams.maxElev])
    end
    title('ON fields');
end

if isfield(analysis.roi,'OFFperimeter')
    subplot(1,2,2); hold on;
    line([x1 x2],[y1 y1],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x2],[y2 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x1],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x2 x2],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    % draw outline of OFFfields
    Color = [0 0 0.7];
    LineWidth = 1;
    for nr = 1:NumROIs
        if  isempty(analysis.roi(nr).OFFperimeter)==0
            for p = 1:size( analysis.roi(nr).OFFperimeter, 1 )-1
                line( analysis.roi(nr).OFFperimeter(p:p+1, 1), analysis.roi(nr).OFFperimeter(p:p+1, 2), ...
                    'Color', Color, 'LineStyle', '-', 'LineWidth', LineWidth );
            end
            line( [ analysis.roi(nr).OFFperimeter(1,   1) analysis.roi(nr).OFFperimeter(end, 1) ], [ analysis.roi(nr).OFFperimeter(1,   2) analysis.roi(nr).OFFperimeter(end, 2) ], ...
                'Color',  Color, 'LineStyle', '-', 'LineWidth', LineWidth );
        end
    end
    box on; axis equal
    set(gca,'YDir','normal')
    set(gca,'TickDir','out');
    set(gca,'TickLen',[0.01 0.01]);
    xlim([metadata.StimParams.minAzi metadata.StimParams.maxAzi])
    try
        ylim([metadata.StimParams.minElev metadata.StimParams.maxElev])
    catch
        ylim([metadata.StimParams.minElev metadata.StimParams.maxElev])
    end
    title('OFF fields');
end
saveas(gcf, fullfile(figSaveDirectory3, 'ON and OFF field distribution.png'))
pause(0.1);
%close(gcf);

% plot field size

figure; 
hold on;
All_OFFsize(All_OFFsize ==0) = NaN;
All_ONsize(All_ONsize ==0) = NaN;

figure
try
    try
        boxplot([All_ONsize', All_OFFsize'], 'Labels', {'ON', 'OFF'})
        h = findobj(gca,'Tag','Box');
        patch(get(h(1),'XData'),get(h(1),'YData'),[0 0 0.7]);
        patch(get(h(2),'XData'),get(h(2),'YData'),[0.7 0 0]);
        hold all
        ml = findobj(gca, 'Tag', 'Median');
        line(get(ml(1),'XData'),get(ml(1),'YData'), 'Color', [0 0 0])
        line(get(ml(2),'XData'),get(ml(2),'YData'), 'Color', [0 0 0])
    catch
        boxplot(All_OFFsize')
        h = findobj(gca,'Tag','Box');
        patch(get(h(1),'XData'),get(h(1),'YData'),[0 0 0.7]);
        hold all
        ml = findobj(gca, 'Tag', 'Median');
        line(get(ml(1),'XData'),get(ml(1),'YData'), 'Color', [0 0 0])
    end
end
ylabel('Field size (deg^2)');
box off;
saveas(gcf, fullfile(figSaveDirectory3, 'Mean field size.png'))
pause(0.1);


% Receptive field size distribution
figure; 
switch stimType
    case 'Patch'
        hist(All_ONsize);hold on;
        hist(All_OFFsize);
        h = findobj(gca, 'Type', 'patch');
        set(h(2), 'facecolor', 'w', 'edgecolor', [0.7 0 0]);
        set(h(1), 'facecolor', 'w', 'edgecolor', [0 0 0.7]);
    case 'SparseNoise'
        hist(All_ONsize);hold on;
        hist(All_OFFsize);
        h = findobj(gca, 'Type', 'patch');
        set(h(2), 'facecolor', 'w', 'edgecolor', [0.7 0 0]);
        set(h(1), 'facecolor', 'w', 'edgecolor', [0 0 0.7]);
    otherwise
        hist(All_OFFsize);
        h = findobj(gca, 'Type', 'patch');
        set(h(1), 'facecolor', 'w', 'edgecolor', [0 0 0.7]);
end
xlabel('Area (deg^2)', 'FontSize', 10')
ylabel('number of boutons', 'FontSize', 10')
%set(gcf, 'color', 'w');
box off;
saveas(gcf, fullfile(figSaveDirectory3, 'Receptive field size distribution.png'))
pause(0.1);

%close(gcf);