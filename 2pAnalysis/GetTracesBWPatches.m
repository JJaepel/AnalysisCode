analysisParams.animal = 'F2359_2019-09-06';
analysisParams.name = 't00003';
analysisParams.sp2ID = 't00007';
analysisParams.expID = 't00003';
analysisParams.server = 0;
analysisParams.level = 0;

analysisParams.smoothing = 3;
analysisParams.StimOnSetDelay = 0.1;
analysisParams.StimOvershoot = 0.2;
analysisParams.threshold = 2;
analysisParams.windowLength = 5;  %how many times longer than the stimperiod
analysisParams.p = 0.05;

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


%load ROIs 
ROIs = LoadRoisS2p(analysisParams);
NumROIs = length(ROIs.roi);
twoPhotonFrames = length(ROIs.roi(1).spks);

%% 2) Calculate analysis settings
NumTrials = metadata.StimParams.numTrials;
NumPatches = metadata.StimParams.numStimElev * metadata.StimParams.numStimAzi * 2;

%create stimulus settings for each patch
numStimX = metadata.StimParams.numStimElev;
numStimY = metadata.StimParams.numStimAzi;
        
numUniqStims = numStimX * numStimY;
PosX = reshape(repmat(linspace(1, numStimX, numStimX),numStimY,1),numUniqStims,[])';
PosY = repmat(linspace(1, numStimY, numStimY),1,numStimX);

PatchNr = 1;
for C = 1:2
    for P =1:NumPatches/2
        metadata.StimParams.Stimulus(PatchNr).X = PosX(P);
        metadata.StimParams.Stimulus(PatchNr).Y = PosY(P);
        metadata.StimParams.Stimulus(PatchNr).Color = C;
        PatchNr = PatchNr + 1;
    end
end

%write StimOnSettings per Patch
for P = 1:NumPatches
    stimTimesNum = find(metadata.StimParams.StimOnTimes(1,:) == P);
    StimOnSetTimes = metadata.StimParams.StimOnTimes(2,stimTimesNum);
    for i = 1:length(StimOnSetTimes)
        StartFrames(i) = find(metadata.StimParams.frameTimes > StimOnSetTimes(i), 1);
        StopFrames(i) = find(metadata.StimParams.frameTimes > (StimOnSetTimes(i)+ metadata.StimParams.stimDuration),1);
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
    raw_data = ROIs.roi(nr).rawF;
    F = raw_data;
    F = smooth(F, analysisParams.smoothing);
    clear raw_data;
    
    BTrace = []; WTrace = [];
    
    %get tuning curves
    for P = 1:NumPatches      
       Y = metadata.StimParams.Stimulus(P).Y;
       X = metadata.StimParams.Stimulus(P).X;
       C = metadata.StimParams.Stimulus(P).Color;    
       for t = 1:NumTrials
           if length(StimOnSetFrames{P}) < t
               disp(['Warning: Patch (x =' num2str(X) ', y = ' num2str(Y) ', trial ' num2str(t) ') was not recorded'])
               if C == 1
                   BlackTraces{nr, Y, X}(t,:) = NaN;
                   BlackAvrTraces(nr, Y, X, t) = NaN;
                   TC1DBlack(nr, ((Y-1)*NumTrials)+t, X) = NaN;
               else
                   WhiteTraces{nr, Y, X}(t,:) = NaN;
                   WhiteAvrTraces(nr, Y, X, t) = NaN;
                   TC1DWhite(nr, ((Y-1)*NumTrials)+t, X) = NaN;
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
                    TC1DBlack(nr, ((Y-1)*NumTrials)+t, X) =nanmean(temp2);
               else
                    WhiteTraces{nr,Y,X}(t,:) =F(DISPix );
                    temp=F(TESTix);
                    WTest(nr,Y,X,t,:)=temp;
                    WhiteAvrTraces(nr,Y,X,t)=nanmean(temp);
                    temp2=F(DISPix);
                    TC1DWhite(nr, ((Y-1)*NumTrials)+t, X) = nanmean(temp2);
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
    WTestTrace(nr,:)=WTrace;
end

AverageB= median(BTestTrace,1);
AverageW =  median(WTestTrace,1);

save([saveDirectory 'TestTraces.mat'],'BTestTrace','WTestTrace','AverageB','AverageW');
save([saveDirectory 'TuningCurves.mat'], 'WhiteTraces','BlackTraces','TC1DBlack','TC1DWhite','WhiteAvrTraces','BlackAvrTraces',...
    'SamplingFreq', 'DisplayWindow', 'TestingWindow');

%% 4) Get receptive fields

for ind = 1:NumROIs
    
    %it will analyze from 1x(length 1 stim period) frames before stim to 5x(length 1 stim period) after start
    count = 1;
    for xnd = 1:metadata.StimParams.numStimElev
        for ynd = 1:metadata.StimParams.numStimAzi
            ONdF{ind}(:,ynd,xnd)=nanmean(WhiteTraces{ind,ynd,xnd}(:,:)); %average over repeats of those frames in the display window shape: frames, y,x.
            OFFdF{ind}(:,ynd,xnd)=nanmean(BlackTraces{ind,ynd,xnd}(:,:));%variable Black_Traces dF/F Johannes way
            for knd = 1:length(DisplayWindow)
                ATracesON{ind}(knd,count,:)=WhiteTraces{ind,ynd,xnd}(:,knd); %size this way will be Reps , PatchX,PatchY,DisplayWindow
                ATracesOFF{ind}(knd,count,:)=BlackTraces{ind,ynd,xnd}(:,knd);
            end
            count = count + 1;
        end
    end
    
    %create data matrix for heatmap
    Group_meanspatches_ON(ind,:)= nanmean(nanmean(ATracesON{ind},3),1); %mean of activity per patch over whole display window and all reps
    Group_meanspatches_OFF(ind,:)= nanmean(nanmean(ATracesOFF{ind},3),1);
    
    analysis.roi(ind).Group_meanspatches_ON=Group_meanspatches_ON(ind,:);
    analysis.roi(ind).Group_meanspatches_OFF=Group_meanspatches_OFF(ind,:);
    
    % do ANOVA and Kruskal-Wallis test to assess whether response is significant at each frame over all repeats and all positions
    for knd=1:length(DisplayWindow)
        Significance(1,knd,ind)=anova1(squeeze(ATracesON{ind}(knd,:,:))',[],'off');
        Significance(2,knd,ind)=anova1(squeeze(ATracesOFF{ind}(knd,:,:))',[],'off');
        Significance(3,knd,ind)=kruskalwallis(squeeze(ATracesON{ind}(knd,:,:))',[],'off');
        Significance(4,knd,ind)=kruskalwallis(squeeze(ATracesOFF{ind}(knd,:,:))',[],'off');
    end
    
    %find latency of the signal
    [Lat,LatON,LatOFF]=FindLatency(StimDur_corr,Significance,DisplayWindow,ind,analysisParams.p);
    
    %calculate responsive cells stats over patch positions over the response positive window
    TracesON{ind}=median(ATracesON{ind},3);
    TracesOFF{ind}=median(ATracesOFF{ind},3);
    
    if isnan(LatON)==0 && isnan(LatOFF) == 0
        pBlack=anova2(TracesOFF{ind}(LatOFF:(LatOFF+StimDur_corr)-1,:),1,'off');
        pWhite=anova2(TracesON{ind}(LatON:(LatON+StimDur_corr)-1,:),1,'off');
        analysis.roi(ind).pWhite= min(pWhite);
        analysis.roi(ind).pBlack =min(pBlack);
        [thzRFtempONci]= CalcRFsSigLocations(ONdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        [thzRFtempOFFci]= CalcRFsSigLocations(OFFdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        
        %to calculate fields at latency
        if isempty(thzRFtempONci)==0 && nnz(thzRFtempONci(:,:,LatON))~=0 
            ON_field= thzRFtempONci(:,:,LatON)+thzRFtempONci(:,:,LatON+1)+thzRFtempONci(:,:,LatON+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
            cleanON = PlotCleanRF(ON_field,metadata,analysisParams,ind, figSaveDirectory);
            ONout1 = bwperim(cleanON);
            [Ony,Onx]=find(ONout1);
            [ONperimeter,ONx,ONy,ONsize]= define_fields(ONout1,metadata);
            analysis.roi(ind).cleanON=cleanON;
            analysis.roi(ind).ONperimeter=ONperimeter;
            analysis.roi(ind).ONx=ONx;
            analysis.roi(ind).ONy=ONy;
            analysis.roi(ind).ONsize=ONsize;
        end
        
        if  isempty(thzRFtempOFFci)==0 && nnz(thzRFtempOFFci(:,:,LatOFF))~=0
            OFF_field= thzRFtempOFFci(:,:,LatOFF)+thzRFtempOFFci(:,:,LatOFF+1)+thzRFtempOFFci(:,:,LatOFF+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
            cleanOFF=PlotCleanRF(OFF_field,metadata,analysisParams,ind, figSaveDirectory);
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
        pWhite=anova2(TracesON{ind}(LatON:(LatON+StimDur_corr)-1,:),1,'off');
        analysis.roi(ind).pWhite= min(pWhite);
        [thzRFtempONci]=CalcRFsSigLocations(ONdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        
        if isempty(thzRFtempONci)|| nnz(thzRFtempONci(:,:,LatON))==0
            disp(['For ROI' num2str(ind) ': Empty RF matrix at latency'])
        else
            ON_field= thzRFtempONci(:,:,LatON)+thzRFtempONci(:,:,LatON+1)+thzRFtempONci(:,:,LatON+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
            cleanON = PlotCleanRF(ON_field,metadata,analysisParams, ind, figSaveDirectory);
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
            PlotSingleRFs(cleanON,Significance,DisplayWindow,metadata,ind,figSaveDirectory2,1, ONx, ONy) 
        end
        
    elseif isnan(LatOFF) == 0
        pBlack=anova2(TracesOFF{ind}(LatOFF:(LatOFF+StimDur_corr)-1,:),1,'off');
        analysis.roi(ind).pBlack =min(pBlack);
        [thzRFtempOFFci]=CalcRFsSigLocations(OFFdF{ind},metadata, analysisParams.threshold, SamplingFreq, DisplayWindow);
        
        if isempty(thzRFtempOFFci)|| nnz(thzRFtempOFFci(:,:,LatOFF))==0
            disp(['For ROI' num2str(ind) ': Empty RF matrix at latency'])
        else
            OFF_field= thzRFtempOFFci(:,:,LatOFF)+thzRFtempOFFci(:,:,LatOFF+1)+thzRFtempOFFci(:,:,LatOFF+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
            cleanOFF=PlotCleanRF(OFF_field,metadata,analysisParams, ind, figSaveDirectory);
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
            PlotSingleRFs(cleanOFF,Significance,DisplayWindow,metadata,ind,figSaveDirectory2, 2, OFFx, OFFy)            
        end
    end
    
    analysis.roi(ind).Significance = Significance;
    
end

save([savedir 'Group_meanspatches.mat'], 'Group_meanspatches_ON', 'Group_meanspatches_OFF');
save([savedir 'Atraces.mat'], 'ATracesON', 'ATracesOFF');
save([savedir 'Patches.mat'], 'analysis');

%% 5) Make summary data
% calculate totals
All_OFFx=calculateTotals(analysis.roi,NumROIs,'OFFx');
All_OFFy=calculateTotals(analysis.roi,NumROIs,'OFFy');
All_ONx=calculateTotals(analysis.roi,NumROIs,'ONx');
All_ONy=calculateTotals(analysis.roi,NumROIs,'ONy');
All_ONsize=calculateTotals(analysis.roi,NumROIs,'ONsize');
All_OFFsize=calculateTotals(analysis.roi,NumROIs,'OFFsize');

x1 = 0; x2 = 0; y1 = 0; y2 = 0;
StimAreaWidth = metadata.StimParams.numStimElev * metadata.StimParams.stimSize(1);
StimAreaHeight = metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(2);

% plot field centers
figure(1);    hold on;
line([x1 x2],[y1 y1],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
line([x1 x2],[y2 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
line([x1 x1],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
line([x2 x2],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
if true((All_ONy~=0)==(All_ONx~=0))
    plot(All_ONx(All_ONx~=0),All_ONy(All_ONy~=0),'o','Color',[0.7 0 0]);    
end
if true((All_OFFy~=0)==(All_OFFx~=0))
    plot(All_OFFx(All_OFFx~=0),All_OFFy(All_OFFy~=0),'x','Color',[0 0 0.7]);
end
box on;
set(gca,'YDir','reverse')
set(gca,'TickDir','out');
set(gca,'TickLen',[0.01 0.01]);
set(gca,'YLim',[-0.5 0.5]*StimAreaHeight);
set(gca,'XLim',[-0.5 0.5]*StimAreaWidth);
title('On and OFF field centers (red=ON, blue=OFF)');
saveas(gcf, fullfile(figSaveDirectory3, 'On and OFF field centers.png'))
pause(0.1);
%close(gcf);

% plot field distribution
figure(2);
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
    set(gca,'YDir','reverse')
    set(gca,'TickDir','out');
    set(gca,'TickLen',[0.01 0.01]);
    set(gca,'XLim',[-0.5 0.5]*StimAreaWidth);
    set(gca,'YLim',[-0.5 0.5]*StimAreaHeight);
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
    set(gca,'YDir','reverse')
    set(gca,'TickDir','out');
    set(gca,'TickLen',[0.01 0.01]);
    set(gca,'XLim',[-0.5 0.5]*StimAreaWidth);
    set(gca,'YLim',[-0.5 0.5]*StimAreaHeight);
    title('OFF fields');
end
saveas(gcf, fullfile(figSaveDirectory3, 'ON and OFF field distribution.png'))
pause(0.1);
%close(gcf);

% plot field size

figure(30); 
hold on;
All_OFFsize(All_OFFsize ==0) = NaN;
All_ONsize(All_ONsize ==0) = NaN;
PlotConfInts(2,nanmean(All_OFFsize), std( All_OFFsize )/ sqrt( length(All_OFFsize)),[0 0 0],0.3,'stderr','both');
bar(2,nanmean(All_OFFsize),'FaceColor',[0 0 0.7],'EdgeColor','none');
PlotConfInts(1,nanmean(All_ONsize),std( All_ONsize )/ sqrt( length(All_ONsize)),[0 0 0],0.3,'stderr','both');
bar(1,nanmean(All_ONsize),'FaceColor',[0.7 0 0],'EdgeColor','none');
title('Mean field size (deg^2)');
%set(gca,'XLim',[0.4 2.6]);
set(gca,'TickDir','out');
%set(gca,'TickLen',[0.01 0.01]);
box off;
saveas(gcf, fullfile(figSaveDirectory3, 'Mean field size.png'))
pause(0.1);

% Receptive field size distribution
figure(4); 
hist(All_ONsize);hold on;
hist(All_OFFsize);
h = findobj(gca, 'Type', 'patch');
xlabel('Area (deg^2)', 'FontSize', 10')
ylabel('number of boutons', 'FontSize', 10')
set(h(1), 'facecolor', 'w', 'edgecolor', [0.7 0 0]);
set(h(2), 'facecolor', 'w', 'edgecolor', [0 0 0.7]);
%set(gcf, 'color', 'w');
box off;
saveas(gcf, fullfile(figSaveDirectory3, 'Receptive field size distribution.png'))
pause(0.1);

%close(gcf);