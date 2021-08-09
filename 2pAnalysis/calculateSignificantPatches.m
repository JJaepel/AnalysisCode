function [analysisParams, metadata, analysis] = calculateSignificantPatches(analysisParams, metadata, data, analysis)

analysisParams.smoothing = 3;
analysisParams.StimOnSetDelay = 0.1;
analysisParams.StimOvershoot = 0.2;
analysisParams.threshold =2;
analysisParams.windowLength = 5;  %how many times longer than the stimperiod
analysisParams.p = 0.01;

saveDirectory = [analysisParams.savedir analysisParams.animal filesep analysisParams.expID filesep];
figSaveDirectory = [analysisParams.savedir analysisParams.animal filesep analysisParams.expID filesep 'cleanRF'];
figSaveDirectory2 = [analysisParams.savedir analysisParams.animal filesep analysisParams.expID filesep 'RFs_lat'];
figSaveDirectory3 = [analysisParams.savedir analysisParams.animal filesep analysisParams.expID filesep 'Summary'];

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

%% 1.) Analysis settings

stimType = metadata.StimParams.type;

%create stimulus settings for each patch
PatchNr = 1;
for C = 1:metadata.StimParams.numCol
    for P =1:metadata.StimParams.numPatches/metadata.StimParams.numCol
        metadata.StimParams.Stimulus(PatchNr).X = metadata.StimParams.numX(P);
        metadata.StimParams.Stimulus(PatchNr).Y = metadata.StimParams.numY(P);
        metadata.StimParams.Stimulus(PatchNr).Color = C;
        PatchNr = PatchNr + 1;
    end
end

if analysisParams.level
    metadata.TwoPhoton.time = metadata.TwoPhoton.time(1:5:end);
end

%write StimOnSettings per Patch
for P = 1:metadata.StimParams.numPatches
    switch stimType
        case 'SparseNoise'
            stimTimesNum = metadata.StimParams.stimTimesAll{P};
            stimTimesNum = stimTimesNum(randperm(length(stimTimesNum)));
            stimTimesNum = stimTimesNum(1:metadata.StimParams.numTrials);
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
metadata.TwoPhoton.SamplingFreq = metadata.TwoPhoton.rate;
StimDur_corr=round(metadata.TwoPhoton.SamplingFreq*metadata.StimParams.stimDuration);

%calculate windows
analysisParams.BaselineWindow = round((-0.7 * metadata.StimParams.isi) * metadata.TwoPhoton.SamplingFreq):round(0.5 * analysisParams.StimOnSetDelay * metadata.TwoPhoton.SamplingFreq);
analysisParams.TestingWindow = round(analysisParams.StimOnSetDelay * metadata.TwoPhoton.SamplingFreq):round((metadata.StimParams.stimDuration + analysisParams.StimOvershoot) * metadata.TwoPhoton.SamplingFreq);
analysisParams.DisplayWindow = round((-0.7 * metadata.StimParams.isi) * metadata.TwoPhoton.SamplingFreq):round((metadata.StimParams.stimDuration + (analysisParams.windowLength * metadata.StimParams.isi)) * metadata.TwoPhoton.SamplingFreq);

disp(' ');
disp(['Baseline period  = ' sprintf('%5.2f',min(analysisParams.BaselineWindow)./metadata.TwoPhoton.SamplingFreq) ' to ' sprintf('%5.2f',max(analysisParams.BaselineWindow)./metadata.TwoPhoton.SamplingFreq) ...
    's (frames: ' num2str(min(analysisParams.BaselineWindow)) ' .. ' num2str(max(analysisParams.BaselineWindow)) ')']);
disp(' ');
disp(['Display period   = ' sprintf('%5.2f',min(analysisParams.DisplayWindow)./metadata.TwoPhoton.SamplingFreq) ' to ' sprintf('%5.2f',max(analysisParams.DisplayWindow)./metadata.TwoPhoton.SamplingFreq) ...
    's (frames: ' num2str(min(analysisParams.DisplayWindow)) ' .. ' num2str(max(analysisParams.DisplayWindow)) ')']);
disp(['Testing period   = ' sprintf('%5.2f',min(analysisParams.TestingWindow)./metadata.TwoPhoton.SamplingFreq) ' to ' sprintf('%5.2f',max(analysisParams.TestingWindow)./metadata.TwoPhoton.SamplingFreq) ...
    's (frames: ' num2str(min(analysisParams.TestingWindow)) ' .. ' num2str(max(analysisParams.TestingWindow)) ')']);
disp(' ');

%% 2.) Loop ROIs and get PSTH's and tuning curves

for nr = 1:length(data.roi)
    %define F
    raw_data = data.roi(nr).rawF;
    F = raw_data;
    F = smooth(F, analysisParams.smoothing);
    clear raw_data;
        
    BTrace = []; WTrace = [];
    
    %get tuning curves
    for P = 1:metadata.StimParams.numPatches      
       Y = metadata.StimParams.Stimulus(P).Y;
       X = metadata.StimParams.Stimulus(P).X;
       C = metadata.StimParams.Stimulus(P).Color;    
       for t = 1:metadata.StimParams.numTrials
           if length(StimOnSetFrames{P}) < t
               disp(['Warning: Patch (x =' num2str(X) ', y = ' num2str(Y) ', trial ' num2str(t) ') was not recorded'])
               if C == 1
                   analysis.rawF.roi(nr).BlackTraces{Y, X}(t,:) = NaN;
                   analysis.rawF.roi(nr).BlackAvrTraces(Y, X, t) = NaN;
                   analysis.rawF.roi(nr).TC1DBlack(((Y-1)*metadata.StimParams.numTrials)+t, X) = NaN;
               else
                   analysis.rawF.roi(nr).WhiteTraces{Y, X}(t,:) = NaN;
                   analysis.rawF.roi(nr).WhiteAvrTraces(Y, X, t) = NaN;
                   analysis.rawF.roi(nr).TC1DWhite(((Y-1)*metadata.StimParams.numTrials)+t, X) = NaN;
               end
           else
               DISPix = StimOnSetFrames{P}(t)+analysisParams.DisplayWindow;
               TESTix = StimOnSetFrames{P}(t)+analysisParams.TestingWindow;
               if C == 1
                    analysis.rawF.roi(nr).BlackTraces{Y,X}(t,:) =F(DISPix );
                    temp=F(TESTix);
                    analysis.rawF.roi(nr).BTest(Y,X,t,:)=temp;
                    analysis.rawF.roi(nr).BlackAvrTraces(Y,X,t)=nanmean(temp);
                    temp2=F(DISPix);
                    analysis.rawF.roi(nr).TC1DBlack(((Y-1)*metadata.StimParams.numTrials)+t, X) =nanmean(temp2);
               else
                    analysis.rawF.roi(nr).WhiteTraces{Y,X}(t,:) =F(DISPix );
                    temp=F(TESTix);
                    analysis.rawF.roi(nr).WTest(Y,X,t,:)=temp;
                    analysis.rawF.roi(nr).WhiteAvrTraces(Y,X,t)=nanmean(temp);
                    temp2=F(DISPix);
                    analysis.rawF.roi(nr).TC1DWhite(((Y-1)*metadata.StimParams.numTrials)+t, X) = nanmean(temp2);
               end
           end
       end       
       if C == 1
           temp3=(squeeze(nanmean(analysis.rawF.roi(nr).BTest(Y,X,:,:),3)))';
           BTrace=[BTrace temp3];
       else
           temp4=(squeeze(nanmean(analysis.rawF.roi(nr).WTest(Y,X,:,:),3)))';
           WTrace=[WTrace temp4];
       end      
    end
    analysis.rawF.roi(nr).BTestTrace=BTrace;
    try
        analysis.rawF.roi(nr).WTestTrace=WTrace;
    catch
        analysis.rawF.roi(nr).WTestTrace = NaN;
    end
end


save([saveDirectory 'TuningCurves.mat'],'analysisParams','metadata','analysis');

%% 3.) Get receptive fields

for ind = 1:length(data.roi)
    
    %it will analyze from 1x(length 1 stim period) frames before stim to 5x(length 1 stim period) after start
    count = 1;
    for xnd = 1:metadata.StimParams.numAzimuth
        for ynd = 1:metadata.StimParams.numElevation
            switch stimType
                case 'Patch'
                    analysis.rawF.roi(ind).ONdF(:,ynd,xnd)=nanmean(analysis.rawF.roi(ind).WhiteTraces{ynd,xnd}(:,:)); %average over repeats of those frames in the display window shape: frames, y,x.
                case 'SparseNoise'
                    analysis.rawF.roi(ind).ONdF(:,ynd,xnd)=nanmean(analysis.rawF.roi(ind).WhiteTraces{ynd,xnd}(:,:));
            end
            analysis.rawF.roi(ind).OFFdF(:,ynd,xnd)=nanmean(analysis.rawF.roi(ind).BlackTraces{ynd,xnd}(:,:));%variable Black_Traces dF/F Johannes way
            for knd = 1:length(analysisParams.DisplayWindow)
                switch stimType
                    case 'Patch'
                        analysis.rawF.roi(ind).ATracesON(knd,count,:)=analysis.rawF.roi(ind).WhiteTraces{ynd,xnd}(:,knd); %size this way will be Reps , PatchX,PatchY,DisplayWindow
                    case 'SparseNoise'
                        analysis.rawF.roi(ind).ATracesON(knd,count,:)=analysis.rawF.roi(ind).WhiteTraces{ynd,xnd}(:,knd);
                end
                analysis.rawF.roi(ind).ATracesOFF(knd,count,:)=analysis.rawF.roi(ind).BlackTraces{ynd,xnd}(:,knd);
            end
            count = count + 1;
        end
    end
    
    %create data matrix for heatmap
    switch stimType
        case 'Patch'
            analysis.rawF.roi(ind).Group_meanspatches_ON= nanmean(nanmean(analysis.rawF.roi(ind).ATracesON,3),1); %mean of activity per patch over whole display window and all reps
        case 'SparseNoise'
            analysis.rawF.roi(ind).Group_meanspatches_ON= nanmean(nanmean(analysis.rawF.roi(ind).ATracesON,3),1); %mean of activity per patch over whole display window and all reps
    end
    analysis.rawF.roi(ind).Group_meanspatches_OFF= nanmean(nanmean(analysis.rawF.roi(ind).ATracesOFF,3),1);
    
    % do ANOVA and Kruskal-Wallis test to assess whether response is significant at each frame over all repeats and all positions
    for knd=1:length(analysisParams.DisplayWindow)
        switch stimType 
            case 'Patch'
                analysis.rawF.roi(ind).Significance(1,knd)=anova1(squeeze(analysis.rawF.roi(ind).ATracesON(knd,:,:))',[],'off');
                analysis.rawF.roi(ind).Significance(3,knd)=kruskalwallis(squeeze(analysis.rawF.roi(ind).ATracesON(knd,:,:))',[],'off');
            case 'SparseNoise'
                analysis.rawF.roi(ind).Significance(1,knd)=anova1(squeeze(analysis.rawF.roi(ind).ATracesON(knd,:,:))',[],'off');
                analysis.rawF.roi(ind).Significance(3,knd)=kruskalwallis(squeeze(analysis.rawF.roi(ind).ATracesON(knd,:,:))',[],'off');
            otherwise
                analysis.rawF.roi(ind).Significance(1,knd) = NaN;
                analysis.rawF.roi(ind).Significance(3,knd) = NaN;
        end
        analysis.rawF.roi(ind).Significance(2,knd)=anova1(squeeze(analysis.rawF.roi(ind).ATracesOFF(knd,:,:))',[],'off');
        analysis.rawF.roi(ind).Significance(4,knd)=kruskalwallis(squeeze(analysis.rawF.roi(ind).ATracesOFF(knd,:,:))',[],'off');
    end
    
    %find latency of the signal
    switch stimType
        case 'SparseNoise'
            [analysis.rawF.roi(ind).Lat,analysis.rawF.roi(ind).LatON,analysis.rawF.roi(ind).LatOFF]=FindLatencyRotPatch(StimDur_corr,analysis.rawF.roi(ind).Significance,analysisParams.DisplayWindow,analysisParams.p);
        case 'rotatingGratingPatch'
            [analysis.rawF.roi(ind).Lat,analysis.rawF.roi(ind).LatON,analysis.rawF.roi(ind).LatOFF]=FindLatencyRotPatch(StimDur_corr,analysis.rawF.roi(ind).Significance,analysisParams.DisplayWindow,analysisParams.p);
        otherwise
            [analysis.rawF.roi(ind).Lat,analysis.rawF.roi(ind).LatON,analysis.rawF.roi(ind).LatOFF]=FindLatencyIntegrated(StimDur_corr,analysis.rawF.roi(ind).Significance,analysisParams.DisplayWindow,analysisParams.p);
    end
    
    %calculate responsive cells stats over patch positions over the response positive window
    switch stimType 
        case 'Patch'
            analysis.rawF.roi(ind).TracesON=median(analysis.rawF.roi(ind).ATracesON,3);
        case 'SparseNoise'
            analysis.rawF.roi(ind).TracesON=median(analysis.rawF.roi(ind).ATracesON,3);
    end
    analysis.rawF.roi(ind).TracesOFF=median(analysis.rawF.roi(ind).ATracesOFF,3);
    
    if isnan(analysis.rawF.roi(ind).LatON)==0 && isnan(analysis.rawF.roi(ind).LatOFF) == 0
        try
            pBlack=anova2(analysis.rawF.roi(ind).TracesOFF(analysis.rawF.roi(ind).LatOFF:(analysis.rawF.roi(ind).LatOFF+StimDur_corr)-1,:),1,'off');
        catch
            pBlack=anova2(analysis.rawF.roi(ind).TracesOFF(analysis.rawF.roi(ind).LatOFF:end,:),1,'off');
        end
        try
            pWhite=anova2(analysis.rawF.roi(ind).TracesON(analysis.rawF.roi(ind).LatON:(analysis.rawF.roi(ind).LatON+StimDur_corr)-1,:),1,'off');
        catch
            pWhite=anova2(analysis.rawF.roi(ind).TracesON(analysis.rawF.roi(ind).LatON:end,:),1,'off');
        end
        analysis.rawF.roi(ind).pWhite= min(pWhite);
        analysis.rawF.roi(ind).pBlack =min(pBlack);
        [analysis.rawF.roi(ind).thzRFtempONci]= CalcRFsSigLocations(analysis.rawF.roi(ind).ONdF,metadata, analysisParams.threshold, metadata.TwoPhoton.SamplingFreq, analysisParams.DisplayWindow);
        [analysis.rawF.roi(ind).thzRFtempOFFci]= CalcRFsSigLocations(analysis.rawF.roi(ind).OFFdF,metadata, analysisParams.threshold, metadata.TwoPhoton.SamplingFreq, analysisParams.DisplayWindow);
        
        %to calculate fields at latency
        try
            if isempty(analysis.rawF.roi(ind).thzRFtempONci)==0 && nnz(analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON))~=0 
                analysis.rawF.roi(ind).ON_field= analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON)+analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON+1)+analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanON = PlotCleanRF(analysis.rawF.roi(ind).ON_field,metadata,analysisParams.threshold,ind, figSaveDirectory);
                ONout1 = bwperim(cleanON);
                [Ony,Onx]=find(ONout1);
                [ONperimeter,ONx,ONy,ONsize]= define_fields(ONout1,metadata);
                analysis.rawF.roi(ind).cleanON=cleanON;
                analysis.rawF.roi(ind).ONperimeter=ONperimeter;
                analysis.rawF.roi(ind).ONx=ONx;
                analysis.rawF.roi(ind).ONy=ONy;
                analysis.rawF.roi(ind).ONsize=ONsize;
            end
        catch
            disp('Latency outside analysis window')
        end
        
        try
            if  isempty(analysis.rawF.roi(ind).thzRFtempOFFci)==0 && nnz(analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF))~=0
                analysis.rawF.roi(ind).OFF_field= analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF)+analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF+1)+analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanOFF=PlotCleanRF(analysis.rawF.roi(ind).OFF_field,metadata,analysisParams.threshold,ind, figSaveDirectory);
                %calculate perimeter
                OFFout1 = bwperim(cleanOFF);
                [Offy,Offx]=find(OFFout1);
                [OFFperimeter,OFFx,OFFy,OFFsize]= define_fields(OFFout1,metadata);
                analysis.rawF.roi(ind).cleanOFF=cleanOFF;
                analysis.rawF.roi(ind).OFFperimeter=OFFperimeter;
                analysis.rawF.roi(ind).OFFx=OFFx;
                analysis.rawF.roi(ind).OFFy=OFFy;
                analysis.rawF.roi(ind).OFFsize=OFFsize;
            end
        catch
            disp('Latency outside analysis window')
        end
        
        if isempty(cleanON) == 0 && isempty(cleanOFF == 0)
            %calculate both fields
            if max(max(analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF)))>max(max(analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON)))
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
        PlotResultingRFsIntegrated(smthONOFF,Onx,Ony,Offx,Offy,analysis.rawF.roi(ind).Significance,analysisParams.DisplayWindow,metadata,ind,figSaveDirectory2);
        
        analysis.rawF.roi(ind).cleanRF=smONOFF;
        analysis.rawF.roi(ind).Latencies=[analysis.rawF.roi(ind).Lat analysis.rawF.roi(ind).LatON analysis.rawF.roi(ind).LatOFF];
        
    elseif isnan(analysis.rawF.roi(ind).LatON) == 0
        try
            pWhite=anova2(analysis.rawF.roi(ind).TracesON(analysis.rawF.roi(ind).LatON:(analysis.rawF.roi(ind).LatON+StimDur_corr)-1,:),1,'off');
        catch
            pWhite=anova2(analysis.rawF.roi(ind).TracesON(analysis.rawF.roi(ind).LatON:end,:),1,'off');
        end
        analysis.rawF.roi(ind).pWhite= min(pWhite);
        [analysis.rawF.roi(ind).thzRFtempONci]=CalcRFsSigLocations(analysis.rawF.roi(ind).ONdF,metadata, analysisParams.threshold, metadata.TwoPhoton.SamplingFreq, analysisParams.DisplayWindow);
        try
            if isempty(analysis.rawF.roi(ind).thzRFtempONci)|| nnz(analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON))==0
                disp(['For ROI' num2str(ind) ': Empty RF matrix at latency'])
            else
                analysis.rawF.roi(ind).ON_field= analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON)+analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON+1)+analysis.rawF.roi(ind).thzRFtempONci(:,:,analysis.rawF.roi(ind).LatON+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanON = PlotCleanRF(ON_field,metadata,analysisParams.threshold, ind, figSaveDirectory);
                ONout1 = bwperim(cleanON);
                [Ony,Onx]=find(ONout1);
                [ONperimeter,ONx,ONy,ONsize]= define_fields(ONout1,metadata);
                analysis.rawF.roi(ind).cleanON=cleanON;
                analysis.rawF.roi(ind).Latencies = LatON;
                analysis.rawF.roi(ind).ONperimeter=ONperimeter;
                analysis.rawF.roi(ind).ONx=ONx;
                analysis.rawF.roi(ind).ONy=ONy;
                analysis.rawF.roi(ind).ONsize=ONsize;
                %plot field, perimeter and statistics
                PlotSingleRFs(cleanON,Significance,analysisParams.DisplayWindow,metadata,ind,figSaveDirectory2,1, Onx, Ony) 
            end
        catch
            disp('Latency outside analysis window')
        end
        
    elseif isnan(analysis.rawF.roi(ind).LatOFF) == 0
        try
            pBlack=anova2(analysis.rawF.roi(ind).TracesOFF(analysis.rawF.roi(ind).LatOFF:(analysis.rawF.roi(ind).LatOFF+StimDur_corr)-1,:),1,'off');
        catch
            pBlack=anova2(analysis.rawF.roi(ind).TracesOFF(analysis.rawF.roi(ind).LatOFF:end,:),1,'off');
        end
        analysis.rawF.roi(ind).pBlack =min(pBlack);
        [analysis.rawF.roi(ind).thzRFtempOFFci]=CalcRFsSigLocations(analysis.rawF.roi(ind).OFFdF,metadata, analysisParams.threshold, metadata.TwoPhoton.SamplingFreq, analysisParams.DisplayWindow);
        try
            if isempty(analysis.rawF.roi(ind).thzRFtempOFFci)|| nnz(analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF))==0
                disp(['For ROI' num2str(ind) ': Empty RF matrix at latency'])
            else
                analysis.rawF.roi(ind).OFF_field= analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF)+analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF+1)+analysis.rawF.roi(ind).thzRFtempOFFci(:,:,analysis.rawF.roi(ind).LatOFF+2); %this line is necessary to increase the signal since the length of the transient is bigger than one stim period
                cleanOFF=PlotCleanRF(OFF_field,metadata,analysisParams.threshold, ind, figSaveDirectory);
                %calculate perimeter
                OFFout1 = bwperim(cleanOFF);
                [Offy,Offx]=find(OFFout1);
                [OFFperimeter,OFFx,OFFy,OFFsize]= define_fields(OFFout1,metadata);
                analysis.rawF.roi(ind).cleanOFF=cleanOFF;
                analysis.rawF.roi(ind).Latencies = LatOFF;
                analysis.rawF.roi(ind).OFFperimeter=OFFperimeter;
                analysis.rawF.roi(ind).OFFx=OFFx;
                analysis.rawF.roi(ind).OFFy=OFFy;
                analysis.rawF.roi(ind).OFFsize=OFFsize;
                %plot field, perimeter and statistics
                PlotSingleRFs(cleanOFF,Significance,analysisParams.DisplayWindow,metadata,ind,figSaveDirectory2, 2, Offx, Offy)            
            end
        catch
            disp('Latency outside analysis window')
        end
    end
        
end
save([saveDirectory 'Patches.mat'],'analysis');


%% 5) Make summary data
% calculate totals
All_OFFx=calculateTotals(analysis.rawF.roi,length(data.roi),'OFFx');
All_OFFy=calculateTotals(analysis.rawF.roi,length(data.roi),'OFFy');
All_ONx=calculateTotals(analysis.rawF.roi,length(data.roi),'ONx');
All_ONy=calculateTotals(analysis.rawF.roi,length(data.roi),'ONy');
All_ONsize=calculateTotals(analysis.rawF.roi,length(data.roi),'ONsize');
All_OFFsize=calculateTotals(analysis.rawF.roi,length(data.roi),'OFFsize');

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
if isfield(analysis.rawF.roi,'ONperimeter')
    subplot(1,2,1); hold on;
    line([x1 x2],[y1 y1],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x2],[y2 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x1],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x2 x2],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    % draw outline of ONfields
    Color = [0.7 0 0];
    LineWidth = 1;
    for nr = 1:length(data.roi)
        if  isempty(analysis.rawF.roi(nr).ONperimeter)==0
            for p = 1:size(analysis.rawF.roi(nr).ONperimeter, 1 )-1
                line(analysis.rawF.roi(nr).ONperimeter(p:p+1, 1), analysis.rawF.roi(nr).ONperimeter(p:p+1, 2), ...
                    'Color', Color, 'LineStyle', '-', 'LineWidth', LineWidth );
            end
            line( [ analysis.rawF.roi(nr).ONperimeter(1,   1) analysis.rawF.roi(nr).ONperimeter(end, 1) ], [ analysis.rawF.roi(nr).ONperimeter(1,   2) analysis.rawF.roi(nr).ONperimeter(end, 2) ], ...
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

if isfield(analysis.rawF.roi,'OFFperimeter')
    subplot(1,2,2); hold on;
    line([x1 x2],[y1 y1],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x2],[y2 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x1 x1],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    line([x2 x2],[y1 y2],'Color',[0.7 0.7 0.7],'linestyle',':','LineWidth',1);
    % draw outline of OFFfields
    Color = [0 0 0.7];
    LineWidth = 1;
    for nr = 1:length(data.roi)
        if  isempty(analysis.rawF.roi(nr).OFFperimeter)==0
            for p = 1:size( analysis.rawF.roi(nr).OFFperimeter, 1 )-1
                line( analysis.rawF.roi(nr).OFFperimeter(p:p+1, 1), analysis.rawF.roi(nr).OFFperimeter(p:p+1, 2), ...
                    'Color', Color, 'LineStyle', '-', 'LineWidth', LineWidth );
            end
            line( [ analysis.rawF.roi(nr).OFFperimeter(1,   1) analysis.rawF.roi(nr).OFFperimeter(end, 1) ], [ analysis.rawF.roi(nr).OFFperimeter(1,   2) analysis.rawF.roi(nr).OFFperimeter(end, 2) ], ...
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