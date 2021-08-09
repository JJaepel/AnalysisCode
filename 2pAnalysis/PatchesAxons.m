function Patches(analysisParams)

close all
if analysisParams.server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

field = analysisParams.field;
plotROIsResps = analysisParams.plotROIs;

base2pDirectory= [TwoPhontondir analysisParams.animal];
tifDirectory = [base2pDirectory filesep analysisParams.name];
Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep];
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
% if ~exist(ROIsaveDirectory, 'dir')
%     mkdir(ROIsaveDirectory);  
% end

disp('Loading data')

%% load Data and metadata
if analysisParams.reloadData
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParamsRet(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
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
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');
    metadata.ROI = struct;
    analysis = struct;
    save(fullfile(saveDirectory, 'Patches.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
    if analysisParams.checkDFF
        figure
        subplot(3,1,1)
        plot(data.roi(1).rawF)
        axis off
        subplot(3,1,2)
        plot(data.roi(1).baseline)
        axis off
        subplot(3,1,3)
        plot(data.roi(1).dff)
        axis off
    end
else
    load(fullfile(saveDirectory, 'Patches.mat'), 'data', 'metadata', 'analysis');
end

%% chop traces
disp('Chopping Traces')
Stimtype = metadata.StimParams.type;
switch Stimtype
    case 'Retinotopy_2D'
        metadata.StimParams.isi = metadata.StimParams.ISI;
        metadata.StimParams.numTrials = metadata.StimParams.numberOfTrials;
        analysisParams.windowStop=metadata.StimParams.stimDuration+0.5;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        if analysisParams.special
            if mod((metadata.StimParams.endPointy - metadata.StimParams.startPointy), (0.5*str2double(metadata.StimParams.stimSize(2)))) == 0
                metadata.StimParams.numElevation = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(0.5*str2double(metadata.StimParams.stimSize(2))));
                metadata.StimParams.numAzimuth = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(0.5* str2double(metadata.StimParams.stimSize(2))));
            else
                metadata.StimParams.numElevation = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(0.5*str2double(metadata.StimParams.stimSize(2))))+1;
                metadata.StimParams.numAzimuth = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(0.5* str2double(metadata.StimParams.stimSize(2))))+1;
            end
            metadata.StimParams.stimPosX = linspace(metadata.StimParams.startPointx,(metadata.StimParams.numAzimuth-1)*0.5*str2double(metadata.StimParams.stimSize(2))+metadata.StimParams.startPointx,metadata.StimParams.numAzimuth); 
            metadata.StimParams.stimPosY = linspace(metadata.StimParams.startPointy,(metadata.StimParams.numElevation-1)*0.5*str2double(metadata.StimParams.stimSize(2))+metadata.StimParams.startPointy,metadata.StimParams.numElevation);
        else
            if mod((metadata.StimParams.endPointy - metadata.StimParams.startPointy), (str2double(metadata.StimParams.stimSize(2)))) == 0
                metadata.StimParams.numElevation = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(str2double(metadata.StimParams.stimSize(2))));
                metadata.StimParams.numAzimuth = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(str2double(metadata.StimParams.stimSize(2))));
            else
                metadata.StimParams.numElevation = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/str2double(metadata.StimParams.stimSize(2)))+1;
                metadata.StimParams.numAzimuth = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/str2double(metadata.StimParams.stimSize(2)))+1;
            end
            metadata.StimParams.stimPosX = linspace(metadata.StimParams.startPointx,(metadata.StimParams.numAzimuth-1)*str2double(metadata.StimParams.stimSize(2))+metadata.StimParams.startPointx,metadata.StimParams.numAzimuth); 
            metadata.StimParams.stimPosY = linspace(metadata.StimParams.startPointy,(metadata.StimParams.numElevation-1)*str2double(metadata.StimParams.stimSize(2))+metadata.StimParams.startPointy,metadata.StimParams.numElevation);
        end
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        metadata.StimParams.PatchX = repmat(metadata.StimParams.stimPosX',1,metadata.StimParams.numAzimuth);
        metadata.StimParams.PatchY = repmat(metadata.StimParams.stimPosY,metadata.StimParams.numElevation,1);
        metadata.StimParams.TwoDStim = 1;
    case 'PatchGrating'
        analysisParams.windowStop=0.5;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
        metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        minAzim = 10 - ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*metadata.StimParams.numAzimuth+minAzim,metadata.StimParams.numAzimuth);
        maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*metadata.StimParams.numElevation,maxElev,metadata.StimParams.numElevation);
        metadata.StimParams.TwoDStim = 1;
    case 'RF_localbar'
        analysisParams.windowStop=0.5;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
        metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        minAzim = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*metadata.StimParams.numAzimuth+minAzim,metadata.StimParams.numAzimuth);
        maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*metadata.StimParams.numElevation,maxElev,metadata.StimParams.numElevation);
        metadata.StimParams.TwoDStim = 0;
    case 'Patch'
        analysisParams.windowStop=1;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        minAzim = -((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*(metadata.StimParams.numAzimuth-1)+minAzim,metadata.StimParams.numAzimuth);
        maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*(metadata.StimParams.numElevation-1),maxElev,metadata.StimParams.numElevation);
        metadata.StimParams.TwoDStim = 1;
    case 'azimuthBars'
        analysisParams.windowStop=1;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.numElevation = 1;
        metadata.StimParams.numStimElev = 1;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        if metadata.StimParams.minAzim == 0
            metadata.StimParams.startAzim = 0;
        else
            metadata.StimParams.startAzim = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1))/2) +metadata.StimParams. stimSize(1)/2;
        end
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        metadata.StimParams.stimPosX = linspace(metadata.StimParams.startAzim,metadata.StimParams.stimSize(2)*(metadata.StimParams.numAzimuth-1)+metadata.StimParams.startAzim,metadata.StimParams.numAzimuth);
        metadata.StimParams.stimPosY = 1;
        metadata.StimParams.TwoDStim = 0;
    case 'RF_local'
        analysisParams.windowStop=1;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.numElevation = 1;
        metadata.StimParams.numStimElev = 1;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        if metadata.StimParams.minAzim == 0
            metadata.StimParams.startAzim = 0;
        else
            metadata.StimParams.startAzim = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(2))/2) +metadata.StimParams. stimSize(2)/2;
        end
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        metadata.StimParams.stimPosX = linspace(metadata.StimParams.startAzim,metadata.StimParams.stimSize(2)*(metadata.StimParams.numAzimuth-1)+metadata.StimParams.startAzim,metadata.StimParams.numAzimuth);
        metadata.StimParams.stimPosY = 1;
        metadata.StimParams.TwoDStim = 0;
    case 'RetWedge'
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
        analysisParams.windowStop=metadata.StimParams.stimDuration;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.numPatches = metadata.StimParams.numWedges;
        metadata.StimParams.TwoDStim = 0;
    case 'Ret_Annulus'
        metadata.StimParams.TwoDStim = 0;
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
        analysisParams.windowStop=metadata.StimParams.stimDuration;
        analysisParams.windowStart=0;
        analysisParams.pre= 0.5;
        analysisParams.post = 0.5;
        metadata.StimParams.numPatches = metadata.StimParams.numSizes;
        
end
[analysis, metadata, data] = ChopStimulusTraceRet(analysis,metadata,data,analysisParams.level, field, 'pre', analysisParams.pre, 'post',analysisParams.post,'windowStart',analysisParams.windowStart, 'windowStop',analysisParams.windowStop);

%% find maxResponses and significant responses
disp('Calculating significant responses')

for i = 1:length(analysis.(field).roi)
    pretrialTime= analysis.(field).preTrialTime;
    preTrialIndex= (1:floor(pretrialTime * metadata.TwoPhoton.rate));
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    %collect our pretrial interval
    analysis.(field).roi(i).isRespSignificant = false;
    analysis.(field).roi(i).respThreshold = [];
    baselines=analysis.(field).roi(i).stimResponseTrace(:,:,preTrialIndex);
    analysis.(field).roi(i).baselineSD = std(baselines,[],3);
    analysis.(field).roi(i).baselineMean = mean(baselines,3);
    analysis.(field).roi(i).baselineMean(analysis.(field).roi(i).baselineMean < 0) = mean(mean(analysis.(field).roi(i).baselineMean,2),1);
    analysisPeriod=(analysis.(field).windowStart:analysis.(field).windowStop);
    stimResp = analysis.(field).roi(i).stimResponseTrace(:,:,analysisPeriod);
    analysis.(field).roi(i).peaks = max(stimResp,[],3);
    analysis.(field).roi(i).zscore = ([analysis.(field).roi(i).peaks]-[analysis.(field).roi(i).baselineMean])./[analysis.(field).roi(i).baselineSD];
    analysis.(field).roi(i).crosser = sum(analysis.(field).roi(i).zscore > analysisParams.zThresh,2);
    analysis.(field).roi(i).respStim = analysis.(field).roi(i).crosser >= ((metadata.StimParams.numTrials)*analysisParams.fraction);
    if sum(analysis.(field).roi(i).respStim) > 0
        analysis.(field).roi(i).isResponseSignificant = 1;
    else 
        analysis.(field).roi(i).isResponseSignificant = 0;
    end
end

%% 
disp('Calculating ROI properties')

numberOfStims = metadata.StimParams.uniqStims;
numStims=numberOfStims-1;

for i = 1:length(data.roi)
    %find preferred Patch for all cells
    medResponse = analysis.(field).roi(i).stimResponseTrace(1:end-1, :, :);
    medResponse = mean(medResponse(:,:,stimWindow),3);
    switch Stimtype
        case 'Patch'
            medResponseA = zeros(metadata.StimParams.numPatches, size(medResponse,2)*2);
            medResponseA(:,1:size(medResponse,2)) = medResponse(1:metadata.StimParams.numPatches,:);
            medResponseA(:,size(medResponse,2)+1:size(medResponse,2)*2) = medResponse(metadata.StimParams.numPatches+1:metadata.StimParams.numPatches*2,:);
            medResponse = medResponseA;
        case 'Retinotopy_2D'
           medResponse = analysis.(field).roi(i).stimResponseTrace;
           medResponse = mean(medResponse(:,:,stimWindow),3);
    end
           
    medResponsePatch= median(medResponse,2);
    [~, prefPatch] = max(medResponsePatch);
    analysis.(field).roi(i).prefPatch = prefPatch;
    
    %find preferred Elevation for all cells for all 2D stimuli
    if metadata.StimParams.TwoDStim == 1
        switch Stimtype
            case 'Patch'
                medResponseElev = zeros(metadata.StimParams.numElevation, metadata.StimParams.numAzimuth*metadata.StimParams.numTrials*2);
                for elev = 1:metadata.StimParams.numElevation
                    medResponseElev(elev,:) = reshape(medResponse((elev-1)*metadata.StimParams.numAzimuth+1:elev*metadata.StimParams.numAzimuth,:),1,metadata.StimParams.numTrials*2*metadata.StimParams.numAzimuth);
                end
            otherwise
                medResponseElev = zeros(metadata.StimParams.numElevation, metadata.StimParams.numAzimuth*metadata.StimParams.numTrials);
                for elev = 1:metadata.StimParams.numElevation
                    medResponseElev(elev,:) = reshape(medResponse((elev-1)*metadata.StimParams.numAzimuth+1:elev*metadata.StimParams.numAzimuth,:),1,metadata.StimParams.numTrials*metadata.StimParams.numAzimuth);
                end
        end
    
        medResponseElev= median(medResponseElev,2);
        [~, prefElev] = max(medResponseElev);

        %find preferred Azimuth for all cells
        switch Stimtype
            case 'Patch'
                medResponseAzi = zeros(metadata.StimParams.numAzimuth, metadata.StimParams.numElevation*metadata.StimParams.numTrials*2);
                for azi = 1:metadata.StimParams.numAzimuth
                    medResponseAzi(azi,:) = reshape(medResponse(azi:metadata.StimParams.numAzimuth:metadata.StimParams.numPatches,:),1,metadata.StimParams.numTrials*2*metadata.StimParams.numElevation);
                end
            otherwise
                medResponseAzi = zeros(metadata.StimParams.numAzimuth, metadata.StimParams.numElevation*metadata.StimParams.numTrials);
                for azi = 1:metadata.StimParams.numAzimuth
                    medResponseAzi(azi,:) = reshape(medResponse(azi:metadata.StimParams.numAzimuth:end,:),1,metadata.StimParams.numTrials*metadata.StimParams.numElevation);
                end
        end
        medResponseAzi= median(medResponseAzi,2);
        [~, prefAzi] = max(medResponseAzi);

        %write to struct
        analysis.(field).roi(i).prefPatch = prefPatch;
        analysis.(field).roi(i).prefElev = prefElev;
        analysis.(field).roi(i).prefAzi = prefAzi;
        analysis.(field).roi(i).prefElevDeg = metadata.StimParams.stimPosY(prefElev);
        analysis.(field).roi(i).prefAziDeg = metadata.StimParams.stimPosX(prefAzi);
    end
end

coc_prop = cbrewer('qual', 'Paired', 12);

%% plot Preference on top of template
for types = 1:2
    if types == 1
        rois = linspace(1,length(analysis.dff.roi),length(analysis.dff.roi));
        allPatchPrefs = [analysis.dff.roi.prefPatch];
        if metadata.StimParams.TwoDStim == 1
            allElevPrefs = [analysis.dff.roi.prefElevDeg];
            allAziPrefs = [analysis.dff.roi.prefAziDeg];
        end
    elseif types == 2
        rois = find([analysis.dff.roi.isResponseSignificant] == 1);
        allPatchPrefs = [analysis.dff.roi(rois).prefPatch];
        if metadata.StimParams.TwoDStim == 1
            allElevPrefs = [analysis.dff.roi(rois).prefElevDeg];
            allAziPrefs = [analysis.dff.roi(rois).prefAziDeg];
        end
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    
    if metadata.StimParams.TwoDStim == 1
        subplot(2,3,1)
        PlotPrefOnTemplateRet(analysis, data, metadata,2, field,data.template, rois)

        subplot(2,3,4)
        histogram(allElevPrefs,linspace(metadata.StimParams.stimPosY(1), metadata.StimParams.stimPosY(end),metadata.StimParams.numElevation), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
        ylabel('Cells');
        xlabel(sprintf('Elevation preference (%s)',char(145)));
        axis square;
        set(gca,'Box','off');

        subplot(2,3,2)
        PlotPrefOnTemplateRet(analysis, data, metadata,3, field,data.template, rois)

        subplot(2,3,5)
        histogram(allAziPrefs,linspace(metadata.StimParams.stimPosX(1), metadata.StimParams.stimPosX(end),metadata.StimParams.numAzimuth), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
        ylabel('Cells');
        xlabel(sprintf('Azimuth preference (%s)',char(145)));
        axis square;
        set(gca,'Box','off');

        subplot(2,3,3)
        PlotPrefOnTemplateRet(analysis, data, metadata,1, field,data.template, rois)

        subplot(2,3,6)
        [NprefPatch, ~] = histcounts(allPatchPrefs, metadata.StimParams.numPatches);
        NprefPatch = reshape(NprefPatch, metadata.StimParams.numElevation, metadata.StimParams.numAzimuth);
        imagesc(NprefPatch)
        colorbar
        colormap('jet')
        axis off
        axis square
    
        ylabel('Cells');
        xlabel(sprintf('Patch preference'));
        axis square;
        set(gca,'Box','off');
    else
        subplot(2,2,[1 3])
        PlotPrefOnTemplateRet(analysis, data, metadata,1, field,data.template, rois)
        
        subplot(2,2,2)
        CM = colormap('jet');
        IX = round(linspace(1,64,metadata.StimParams.numPatches));
        for cc = 1:length(IX)
            C{cc} = CM(IX(cc),:);
        end
        switch Stimtype
            case 'RetWedge'
                stimPatches = [6 1; 5 2; 4 3];
                imagesc(stimPatches)
            case 'Ret_Annulus'
                for sizes = 1:metadata.StimParams.numSizes
                    plot(0,0,'o','MarkerSize',30*sizes, 'Color', C{sizes})
                    hold on
                end
        end
        axis off
        axis square
        title('Patch positions')
        
        subplot(2,2,4)
        [NprefPatch, ~] = histcounts(allPatchPrefs, metadata.StimParams.numPatches);
        for sizes = 1:metadata.StimParams.numPatches
            h = bar(sizes, NprefPatch(sizes));
            set(h,'FaceColor', C{sizes})
            hold all
        end
        
%         switch Stimtype
%             case 'RetWedge'
%                 NprefPatch = reshape(NprefPatch([6 5 4 1 2 3]),3,2);
%                 imagesc(NprefPatch)
%                 colorbar('Location', 'SouthOutside')
%             case 'Ret_Annulus'
%                 for sizes = 1:metadata.StimParams.numSizes
%                     h = bar(sizes, NprefPatch(sizes))
%                     set(h,'FaceColor', C{sizes})
%                     hold all
%                 end             
%         end

        axis off
        axis square
        title('Patch preference')
        
    end
    
    set(gcf, 'color', 'w');
    if types == 1
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_all_cells.png'))
    elseif types == 2
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_resp_cells.png'))
    end
end

%% plot ROIs
if analysisParams.plotROIs
    if ~exist(ROIRespsaveDirectory, 'dir')
        mkdir(ROIRespsaveDirectory);
    else
        %remove old files
        filePattern = fullfile(ROIRespsaveDirectory, '*.png'); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        for k = 1 : length(theFiles)
          baseFileName = theFiles(k).name;
          fullFileName = fullfile(ROIRespsaveDirectory, baseFileName);
          delete(fullFileName);
        end
    end
    if ~analysisParams.plotRespROIsOnly
        if ~exist(ROINonRespsaveDirectory, 'dir')
            mkdir(ROINonRespsaveDirectory);
        else
            %remove old files
            filePattern = fullfile(ROIRespsaveDirectory, '*.png'); % Change to whatever pattern you need.
            theFiles = dir(filePattern);
            for k = 1 : length(theFiles)
              baseFileName = theFiles(k).name;
              fullFileName = fullfile(ROIRespsaveDirectory, baseFileName);
              delete(fullFileName);
            end
        end
    end
    for i = 1:length(data.roi)
        if analysis.(field).roi(i).isResponseSignificant == 1
            if metadata.StimParams.TwoDStim == 1
                PlotTrialStimResponse2D(metadata, analysis, field, i)
                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                close gcf
                PlotAvgStimResponse2D(metadata, analysis, field, i)
                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                close gcf
            else
                PlotTrialStimResponse1D(metadata, analysis, field, i)
                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                close gcf
                PlotAvgStimResponse1D(metadata, analysis, field, i)
                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                close gcf
            end
        else
            if ~analysisParams.plotRespROIsOnly
                if metadata.StimParams.TwoDStim == 1
                    PlotTrialStimResponse2D(metadata, analysis, field, i)
                    saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                    close gcf
                    PlotAvgStimResponse2D(metadata, analysis, field, i)
                    saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                    close gcf
                else
                    PlotTrialStimResponse1D(metadata, analysis, field, i)
                    saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                    close gcf
                    PlotAvgStimResponse1D(metadata, analysis, field, i)
                    saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                    close gcf
                end
            end
        end
    end
end

save(fullfile(saveDirectory, 'PatchesAna.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
end
%%

    
