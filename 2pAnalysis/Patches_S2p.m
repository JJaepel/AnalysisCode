function Patches_S2p(animal, expt_id, sp2id, name, reloadData, plotROIsResps, special)

close all
% TwoPhontondir = 'Z:\Juliane\Data\2P_Data\';
% Sp2dir = 'Z:\Juliane\Data\Spike2Data\';
% savedir = 'Z:\Juliane\Data\ImageAnalysis\';

TwoPhontondir = 'F:\Data\2P_Data\';
Sp2dir = 'F:\Data\Spike2Data\';
savedir = 'F:\Data\ImageAnalysis\';

windowStop=1.75;
windowStart=0;
pre= 0.5;
post = 0.5;
field = 'dff';

z_thresh = 5;
fraction = 0.5;
shufflenum = 100;
predictor = 1;

base2pDirectory= [TwoPhontondir animal];
tifDirectory = [base2pDirectory filesep name];
Sp2dDirectory = [Sp2dir animal filesep sp2id filesep];
saveDirectory = [savedir animal filesep expt_id filesep];
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
if ~exist(ROIsaveDirectory, 'dir')
    mkdir(ROIsaveDirectory);  
end
if ~exist(ROIRespsaveDirectory, 'dir')
    mkdir(ROIRespsaveDirectory);  
end
if ~exist(ROINonRespsaveDirectory, 'dir')
    mkdir(ROINonRespsaveDirectory);  
end

disp('Loading data')

%%load Data and metadata
if reloadData
    sliceparams = struct;
    sliceparams.expt_id = expt_id;
    sliceparams.baseDirectory = base2pDirectory;
    metadata.StimParams=Load_stimparams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    metadata.StimParams.series=expt_id;
    data = Load_rois(sliceparams);
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');
    metadata.ROI = struct;
    analysis = struct;
    save(fullfile(saveDirectory, 'Patches.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
else
    load(fullfile(saveDirectory, 'Patches.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
end

%% chop traces
disp('Chopping Traces')
Stimtype = metadata.StimParams.type;
metadata.StimParams.onlyAzi = 0;
switch Stimtype
    case 'Retinotopy_2D'
        metadata.StimParams.isi = metadata.StimParams.ISI;
        metadata.StimParams.numTrials = metadata.StimParams.numberOfTrials;
        windowStop=1.75;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
        if special
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
    case 'PatchGrating'
        windowStop=0.5;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
        metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        minAzim = 10 - ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*metadata.StimParams.numAzimuth+minAzim,metadata.StimParams.numAzimuth);
        maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*metadata.StimParams.numElevation,maxElev,metadata.StimParams.numElevation);
    case 'RF_localbar'
        windowStop=0.5;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
        metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        minAzim = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*metadata.StimParams.numAzimuth+minAzim,metadata.StimParams.numAzimuth);
        maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*metadata.StimParams.numElevation,maxElev,metadata.StimParams.numElevation);
    case 'Patch'
        windowStop=1;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
        metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
        metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
        metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
        minAzim = -((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
        metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*(metadata.StimParams.numAzimuth-1)+minAzim,metadata.StimParams.numAzimuth);
        maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
        metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*(metadata.StimParams.numElevation-1),maxElev,metadata.StimParams.numElevation);
    case 'azimuthBars'
        windowStop=1;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
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
        metadata.StimParams.onlyAzi = 1;
    case 'RF_local'
        windowStop=1;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
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
        metadata.StimParams.onlyAzi = 1;
    case 'RetWedge'
        windowStop=1;
        windowStart=0;
        pre= 0.5;
        post = 0.5;
        metadata.StimParams.numPatches = metadata.StimParams.numWedges;
        metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
end
[analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, 'pre', pre, 'post',post,'windowStart',windowStart, 'windowStop',windowStop);

%% find maxResponses and significant responses
disp('Calculating significant responses')
sigVector=zeros(1,length(analysis.(field).roi));

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
    analysis.(field).roi(i).crosser = sum(analysis.(field).roi(i).zscore > z_thresh,2);
    analysis.(field).roi(i).respStim = analysis.(field).roi(i).crosser >= ((metadata.StimParams.numTrials)*fraction);
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
    end
           
    medResponsePatch= median(medResponse,2);
    [~, prefPatch] = max(medResponsePatch);
    
    %find preferred Elevation for all cells
    switch Stimtype
        case 'Patch'
            medResponseElev = zeros(metadata.StimParams.numElevation, metadata.StimParams.numAzimuth*metadata.StimParams.numTrials*2);
            for elev = 1:metadata.StimParams.numElevation
                medResponseElev(elev,:) = reshape(medResponse((elev-1)*metadata.StimParams.numAzimuth+1:elev*metadata.StimParams.numAzimuth,:),1,metadata.StimParams.numTrials*2*metadata.StimParams.numAzimuth);
            end
        case 'RetWedge'
            medResponseElev = medResponse
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
        case 'RetWedge'
            medResponseAzi = medResponse;
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

coc_prop = cbrewer('qual', 'Paired', 12);

%% plot Preference on top of template
for types = 1:2
    if types == 1
        rois = linspace(1,length(analysis.dff.roi),length(analysis.dff.roi));
        allPatchPrefs = [analysis.dff.roi.prefPatch];
        allElevPrefs = [analysis.dff.roi.prefElevDeg];
        allAziPrefs = [analysis.dff.roi.prefAziDeg];
    elseif types == 2
        rois = find([analysis.dff.roi.isResponseSignificant] == 1);
        allPatchPrefs = [analysis.dff.roi(rois).prefPatch];
        allElevPrefs = [analysis.dff.roi(rois).prefElevDeg];
        allAziPrefs = [analysis.dff.roi(rois).prefAziDeg];
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    if metadata.StimParams.onlyAzi == 0
        subplot(2,3,1)
        PlotPrefOnTemplate(analysis, data, metadata,2, field,data.template, rois)

        subplot(2,3,4)
        histogram(allElevPrefs,linspace(metadata.StimParams.stimPosY(1), metadata.StimParams.stimPosY(end),metadata.StimParams.numElevation), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
        ylabel('Cells');
        xlabel(sprintf('Elevation preference (%s)',char(145)));
        axis square;
        set(gca,'Box','off');
    end

    subplot(2,3,2)
    PlotPrefOnTemplate(analysis, data, metadata,3, field,data.template, rois)
    
    subplot(2,3,5)
    histogram(allAziPrefs,linspace(metadata.StimParams.stimPosX(1), metadata.StimParams.stimPosX(end),metadata.StimParams.numAzimuth), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel(sprintf('Azimuth preference (%s)',char(145)));
    axis square;
    set(gca,'Box','off');
    
    subplot(2,3,3)
    PlotPrefOnTemplate(analysis, data, metadata,1, field,data.template, rois)
    
    subplot(2,3,6)
    [NprefPatch, ~] = histcounts(allPatchPrefs, metadata.StimParams.numPatches);
    NprefPatch = reshape(NprefPatch, metadata.StimParams.numElevation, metadata.StimParams.numAzimuth);
    colormap('jet')
    imagesc(NprefPatch)
    colorbar
    axis off
    
    ylabel('Cells');
    xlabel(sprintf('Patch preference'));
    axis square;
    set(gca,'Box','off');
    
    set(gcf, 'color', 'w');
    if types == 1
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_all_cells.png'))
    elseif types == 2
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_resp_cells.png'))
    end
end

%% plot ROIs
if plotROIsResps
    for i = 1:length(data.roi)
        if analysis.(field).roi(i).isResponseSignificant == 1
            PlotTrialStimResponse(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        else
            PlotTrialStimResponse(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        end
    end
end

save(fullfile(saveDirectory, 'Patches_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
end
%%
function [data]= Load_rois(sliceparams)
    basedirectory = sliceparams.baseDirectory;
    expt_id = sliceparams.expt_id;
    
    Suite2pFile = [basedirectory filesep expt_id filesep 'suite2p' filesep 'plane0' filesep 'Fall.mat'];
    Suite2p = load(Suite2pFile);
    load(Suite2pFile);
    disp('loaded Suite2pFile');
    cell_selector = logical(Suite2p.iscell(:,1));
    data.roi = [];
    counter = 1;
    for i = 1:length(cell_selector)
        if cell_selector(i) == 1
            data.roi(counter).xPos = Suite2p.stat{i}.med(2);
            data.roi(counter).yPos = Suite2p.stat{i}.med(1);
            data.roi(counter).mask = [Suite2p.stat{i}.xpix; Suite2p.stat{i}.ypix]';
            data.roi(counter).name = i;
            data.roi(counter).rawF = double(Suite2p.F(i,:));
            counter = counter +1;
        end
    end
    template = Suite2p.ops.meanImg(1:512,1:512);
    data.template = template./prctile(template(:),99.9);
    clear Suite2p
    disp('loaded data from Suite2p')
end
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
function y = percentileFilt1(x,percentile,n,blksz,DIM)
    narginchk(1,5);
    if nargin < 4, blksz = []; end
    if nargin < 5, DIM = []; end

    % Check the input data type. Single precision is not supported.
    % try
    %     chkinputdatatype(x,n,blksz,DIM);
    % catch ME
    %     throwAsCaller(ME);
    % end

    % Check if the input arguments are valid
    if isempty(n)
      n = 3;
    end

    if ~isempty(DIM) && DIM > ndims(x)
        error(message('signal:medfilt1:InvalidDimensions'))
    end

    % Reshape x into the right dimension.
    if isempty(DIM)
        % Work along the first non-singleton dimension
        [x, nshifts] = shiftdim(x);
    else
        % Put DIM in the first (row) dimension (this matches the order 
        % that the built-in filter function uses)
        perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
        x = permute(x,perm);
    end

    % Verify that the block size is valid.
    siz = size(x);
    if isempty(blksz)
        blksz = siz(1); % siz(1) is the number of rows of x (default)
    else
        blksz = blksz(:);
    end

    % Initialize y with the correct dimension
    y = zeros(siz); 

    % Call medfilt1D (vector)
    for i = 1:prod(siz(2:end))
        y(:,i) = prctilefilt1d(x(:,i),n,blksz,percentile);
    end

    % Convert y to the original shape of x
    if isempty(DIM)
        y = shiftdim(y, -nshifts);
    else
        y = ipermute(y,perm);
    end
end
function y = percentile(x, k)
    x = sort(x);
    n = size(x,1);

    p = 1 + (n-1) * k / 100;

    if p == fix(p)
        y = x(p);
    else
        r1 = floor(p); r2 = r1+1;
        y = x(r1) + (x(r2)-x(r1)) * k / 100;
    end
end
function y = prctilefilt1d(x,n,blksz,percentile)
    %PRCTILEFILT1D  One dimensional median filter.
    %
    % Inputs:
    %   x     - vector
    %   n     - order of the filter
    %   blksz - block size

    nx = length(x);
    if rem(n,2)~=1    % n even
        m = n/2;
    else
        m = (n-1)/2;
    end
    X = [zeros(m,1); x; zeros(m,1)];
    y = zeros(nx,1);

    % Work in chunks to save memory
    indr = (0:n-1)';
    indc = 1:nx;
    for i=1:blksz:nx
        ind = indc(ones(1,n),i:min(i+blksz-1,nx)) + ...
              indr(:,ones(1,min(i+blksz-1,nx)-i+1));
        xx = reshape(X(ind),n,min(i+blksz-1,nx)-i+1);
        y(i:min(i+blksz-1,nx)) = prctile(xx,percentile,1);
    end
end
function [zoom, setup] = getzoom(tifDirectory)
    filename = [tifDirectory filesep '*.tif'];
    files = dir(filename);
    filepath = [tifDirectory filesep files(1).name];
    InfoImage = imfinfo(filepath);
    try 
        a = InfoImage(1).Software;
        setup = 1;
    catch 
        a = InfoImage(1).ImageDescription;
        setup = 2;
    end
    zoom = regexp(a,'(?<=scanZoomFactor = )\d+\.?\d*', 'match');
    zoom = str2num(zoom{1});
end
function [fields] = stimField(Stimtype)
    switch Stimtype
        case 'Patch'
            fields = {'numTrials',...
                'doBlank',...
                'blankpercent',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'centerPoint',...
                'stimSize',...
                'numStimElev',...
                'numStimAzi',...
                };
        case 'Retinotopy_2D'
            fields = {'numberOfTrials',...
                'doBlank',...
                'stimDuration',...
                'ISI',...
                'isRandom',...
                'initialDelay',...
                'centerPoint',...
                'startPointx',...
                'endPointx',...
                'startPointy',...
                'endPointy',...
                'stimSize',...
                };
        case 'PatchGrating'
            fields = {'numTrials',...
                'doBlank',...
                'blankpercent',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'centerPos',...
                'stimPos',...
                'centerPoint',...
                'stimSize',...
                'numStimElev',...
                'numStimAzi',...
                'numOrientations',...
                };
        case 'azimuthBars'
            fields = {'numTrials',...
                'doBlank',...
                'blankpercent',...
                'stimDuration',...
                'isi',...
                'bwtime',...
                'isRandom',...
                'initialDelay',...
                'centerPos',...
                'stimPos',...
                'centerPoint',...
                'stimSize',...
                'numStimAzi',...
                'numOrientations',...
                'spatialFreq',...
                'minAzim',...
                };
        case 'RF_local'
            fields = {'numTrials',...
                'doBlank',...
                'blankpercent',...
                'stimDuration',...
                'isi',...
                'bwtime',...
                'isRandom',...
                'initialDelay',...
                'centerPos',...
                'stimPos',...
                'centerPoint',...
                'stimSize',...
                'numStimAzi',...
                'numOrientations',...
                'spatialFreq',...
                'minAzim',...
                };
        case 'RF_localbar'
            fields = {'numTrials',...
                'doBlank',...
                'blankpercent',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'centerPos',...
                'stimPos',...
                'centerPoint',...
                'stimSize',...
                'numStimElev',...
                'numStimAzi',...
                'numOrientations',...
                'minAzim',...
                'minElev',...
                };
        case 'RetWedge'
            fields = {'numTrials',...
                'doBlank',...
                'nBlank',...
                'blankpercent',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'centerPos',...
                'stimSize',...
                'numWedges',...
                'temporalFreq',...
                'spatialFreq',...
                'numOrientations',...
                };
        otherwise
            fields ={};
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
function [metadata, data]= baselinePercentileFilter(metadata,data, fieldSource,fieldTarget,filteredCutoff, desiredPercentileRank)            
    numSeries = 1;
    fps= metadata.TwoPhoton.rate; %Hz
    if ~isfield(metadata, fieldTarget)
        metadata.baseline= struct;
    end
    if ~isfield(metadata.baseline, fieldSource)
        metadata.baseline.(fieldSource)= struct;
    end
    metadata.baseline.(fieldSource).filteredCutoff=filteredCutoff;
    metadata.baseline.(fieldSource).desiredPercentileRank= desiredPercentileRank;
    for i=1:length(data.roi)
        for j = 1:numSeries
            disp({'baselining cell #',  num2str(i)});
            if isfield(metadata.StimParams, 'seriesBreaks')
                startIdx = 1+ metadata.StimParams.seriesBreaks(j);
                stopIdx = metadata.StimParams.seriesBreaks(j+1);
                inputTrace = data.roi(i).(fieldSource)(startIdx:stopIdx)';
            else
                inputTrace = data.roi(i).(fieldSource)';
            end
            paddingLength = ceil(length(inputTrace)/1);

            paddedTrace   = [inputTrace(paddingLength:-1:1); inputTrace; inputTrace(paddingLength:-1:1)];
            filteredTrace = percentileFilt1(paddedTrace, desiredPercentileRank, round(filteredCutoff*fps));
            filteredTrace = filteredTrace(paddingLength+(1:length(inputTrace)));

            butterWorthOrder = 1;
            Wn = (1/filteredCutoff) / (fps/2);
            [b,a] = butter(butterWorthOrder, Wn, 'low');
            highpassFilteredTrace = filtfilt(b,a,[filteredTrace(paddingLength:-1:1); filteredTrace; filteredTrace(1:paddingLength)]);
            highpassFilteredTrace = highpassFilteredTrace(paddingLength+[1:length(inputTrace)]);
            if j==1
                data.roi(i).(fieldTarget)=highpassFilteredTrace';
            else
                data.roi(i).(fieldTarget)= horzcat(obj.data.roi(i).(fieldTarget), highpassFilteredTrace');
            end
        end
    end
end
function [data] = computeDff(data, sourceSignal, sourceBaseline, fieldTarget)
    for i=1:length(data.roi)
        data.roi(i).(fieldTarget)= (data.roi(i).(sourceSignal)-data.roi(i).(sourceBaseline))./abs(data.roi(i).(sourceBaseline));
    end
end
function [analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, varargin)
    if mod(length(varargin),2) ~= 0
        error('You need to specify both the field and the argument')
    elseif strcmp(metadata.StimParams.type, 'Spontaneous')
        warning('Unable to chop up Spontaneous activity into meaningful chunks')
        return
    else
        % defaults for values
        windowStart= 0;  % choose the entire stim period
        windowStop = 0.5;    % choose the entire stim period
        tfield = field; % default to field name as analysis name
        pre=0;
        post=1;
        skipTrial=0;
        baseline=false;
        for i = 1:2:length(varargin)
            switch varargin{i}
                case 'analysisName'
                    tfield = varargin{i+1};
                case 'windowStart'
                    windowStart= varargin{i+1};
                case 'windowStop'
                    windowStop=varargin{i+1};
                case 'pre'
                    pre = varargin{i+1};
                case 'post'
                    post = varargin{i+1};
                case 'skipTrial'
                    skipTrial = varargin{i+1};
                case 'baseline'
                    baseline=varargin{i+1};
                otherwise
                    error(['You have specified an invalid field: ', varargin{i}])
            end
        end
    end

    TwoPhotonRate= metadata.TwoPhoton.rate;
    analysis.(tfield)=struct;
    analysis.(tfield).preTrialTime = pre;
    analysis.(tfield).postTrialTime = post;
    offsetPre = round(metadata.TwoPhoton.rate * analysis.(tfield).preTrialTime);
    offsetPost = round(metadata.TwoPhoton.rate * analysis.(tfield).postTrialTime);
    analysis.(tfield).stimStart=offsetPre+1;
    analysis.(tfield).stimStop=round(analysis.(tfield).stimStart+metadata.TwoPhoton.rate * metadata.StimParams.stimDuration());
    analysis.(tfield).roi = struct;

    windowStartIdx= round(windowStart * metadata.TwoPhoton.rate +offsetPre+1);
    windowStopIdx = round(windowStop * TwoPhotonRate + offsetPost);
    if windowStartIdx == windowStopIdx
        windowStopIdx = analysis(tfield).stimStop;
    end
    analysis.(tfield).windowStart =windowStartIdx;
    analysis.(tfield).windowStop = windowStopIdx;

    %do a little clean up here to make stim numbers match
%     while (metadata.StimParams.StimOnTimes(2,metadata.StimParams.numberOfStims)+metadata.StimParams.stimDuration > metadata.TwoPhoton.time(end))
%         if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ==0
%             metadata.StimParams.numTrials = metadata.StimParams.numTrials -1;
%         else
%             metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
%         end
%         metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
%     end

    Stimtype = metadata.StimParams.type;
    switch Stimtype
        case 'Patch'
            numVisStims = metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev *metadata.StimParams.numTrials*2;
            numBlanks = floor(metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev * metadata.StimParams.blankpercent *2) * metadata.StimParams.numTrials;
        case 'RetWedge'
            numVisStims = metadata.StimParams.numWedges * metadata.StimParams.numTrials;
            numBlanks = metadata.StimParams.nBlank * metadata.StimParams.numTrials;
        otherwise
            numBlanks = floor(metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev * metadata.StimParams.blankpercent) * metadata.StimParams.numTrials;
            numVisStims = metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev *metadata.StimParams.numTrials;
    end
    if (numBlanks+numVisStims) > metadata.StimParams.numberOfStims
        metadata.StimParams.numTrials = metadata.StimParams.numTrials-1;
    elseif max(metadata.TwoPhoton.time) < max(metadata.StimParams.StimOnTimes(2,:))
        metadata.StimParams.numTrials = metadata.StimParams.numTrials-1;
        numBlanks = floor(metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev * metadata.StimParams.blankpercent) * metadata.StimParams.numTrials;
        numVisStims = metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev *metadata.StimParams.numTrials;
    end

    for i =1:(numBlanks+numVisStims)
        %try
            metadata.StimParams.stimStartIndex(i)= find(metadata.TwoPhoton.time > metadata.StimParams.StimOnTimes(2,i),1);
            metadata.StimParams.stimStopIndex(i)= floor(metadata.StimParams.stimStartIndex(i)+ metadata.StimParams.stimDuration* metadata.TwoPhoton.rate);
        %end
    end

    stimPeriod=(analysis.(tfield).stimStart:analysis.(tfield).stimStop);
    analysisPeriod=(analysis.(tfield).windowStart:analysis.(tfield).windowStop);

    %correct for our specified window

    if windowStartIdx > max(stimPeriod) || windowStopIdx > max(stimPeriod)
        warning('You have choosen a period to average over that is longer than the stim period')
    end

    stimStarts=metadata.StimParams.stimStartIndex;
    stimStops=metadata.StimParams.stimStopIndex;

    for stimID= 1:metadata.StimParams.uniqStims
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        StimonTimes= metadata.StimParams.StimOnTimes;
        StimonTimes=StimonTimes(:,1:metadata.StimParams.numberOfStims);
        stimIndices = find(StimonTimes(1,:)==stimulus);
        for trialNumber= 1:metadata.StimParams.numTrials
            stimIndex = stimIndices(trialNumber);
            for i=1:length(data.roi)
                selectedFramesTrace = (stimStarts(stimIndex)-offsetPre):(stimStops(stimIndex)+offsetPost);
                if max(selectedFramesTrace) <= length(data.roi(i).(field))
                    analysis.(tfield).roi(i).stimResponseTrace(stimID,trialNumber,:)= data.roi(i).(field)(selectedFramesTrace);
                else
                    analysis.(tfield).roi(i).stimResponseTrace(stimID,trialNumber,:) = zeros(1,1,length(selectedFramesTrace));
                end
            end
        end

        for i=1:length(data.roi)
            if skipTrial >0 && skipTrial <=size(analysis.(tfield).roi(i).stimResponseTrace,2)
                analysis.(tfield).roi(i).stimResponseTrace(:,skipTrial,:)=[];
            end
            traces=squeeze(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:));
            bltraces=repmat(mean(traces(:,1:(stimPeriod(1)-1)),2)', size(traces,2),1)';

            if baseline
                traces=traces- bltraces;
            end
            try
                analysis.(tfield).roi(i).stimResponseTraceRaw(stimID,:,:)=analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:);
            catch ME
                ME.message
            end
            analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:)= traces;
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = mean(traces, 1);
            n= size(analysis.(tfield).roi(i).stimResponseTrace,2);
            y=analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:);
            analysis.(tfield).roi(i).SEMResponseTrace(stimID,:) = std(y,[],2)/sqrt(n);
            analysis.(tfield).roi(i).avgStimResponse(stimID,:) = mean(analysis.(tfield).roi(i).avgResponseTrace(stimID,analysisPeriod),2);
            n= size(analysis.(tfield).roi(i).stimResponseTrace,2);
            y=analysis.(tfield).roi(i).avgStimResponse(stimID,:);
            analysis.(tfield).SEMCellResponses(i,stimID)= std(y,[],2)/sqrt(n);
            analysis.(tfield).avgCellResponses(i,stimID)= mean(analysis.(tfield).roi(i).avgStimResponse(stimID,:),2);
        end
    end

    for i = 1:length(data.roi)
        for stimID = 1:metadata.StimParams.uniqStims-1
            analysis.(tfield).roi(i).avgStimResponse(stimID) = analysis.(tfield).roi(i).avgStimResponse(stimID) -analysis.(tfield).roi(i).avgStimResponse(end) ;
            analysis.(tfield).roi(i).avgResponseTraceRaw(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:);
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:) - analysis.(tfield).roi(i).avgResponseTrace(end,:);
            analysis.(tfield).roi(i).Fb(stimID,:)= mean(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,1:analysisPeriod(1)-1),3);

        end
    end
    analysis.(tfield).preTrialTime=pre;
end        
function PlotPrefOnTemplate(analysis, data, metadata,type, field,template,rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,template,template,template)/prctile(template(:),99));
    colormap(hsv) 
    if type == 1
        numPatches = metadata.StimParams.uniqStims-1;
        LUT = jet(numPatches);
        title('Patch map')
        caxis([0 numPatches]); colorbar('Location', 'southoutside');
    elseif type == 2
        title('Elevation map')
        numElev = metadata.StimParams.numElevation;
        LUT = jet(numElev);
        caxis([metadata.StimParams.stimPosY(1) metadata.StimParams.stimPosY(end)]);colorbar('Location', 'southoutside');
    elseif type == 3
        numAzi = metadata.StimParams.numAzimuth;
        LUT = jet(numAzi);
        title('Azimuth map')
        caxis([metadata.StimParams.stimPosX(1) metadata.StimParams.stimPosX(end)]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        l = rois(i);
        xpos= data.roi(l).xPos;
        ypos= data.roi(l).yPos;
        try
           if type == 1
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(l).prefPatch,:));
            elseif type == 2
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(l).prefElev,:));
            elseif type == 3
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(l).prefAzi,:));
           end
        catch
        end
    end
end
function PlotAvgStimResponse(metadata, analysis, field, roi)
    %PlotAvgStimResponse Plots Avg + SEM response for stimuli
    %  	Args:
    %     obj (Population Imaging): class
    %     field (char): analysis field key
    %     roi (int): cell index number
    %   Returns:
    %     NONE
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    stimNum = metadata.StimParams.uniqStims-1;
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    for i=1:stimNum
        ax=subplot(1,stimNum, i);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
            bl=x(stimWindow);
            patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        err=analysis.(field).roi(roi).SEMResponseTrace(i,:);
        hold on
        patch([x fliplr(x)],[y+err fliplr(y-err)],[.5 .5 .5], 'LineStyle', 'none');
        %patch([bl fliplr(bl)],[zeros(1,length(bl)), ones(1,length(bl))],[0 0 0]);
        %plot(ax,x,y, 'k-')
        plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTraceRaw(i,:)), 'k');
        hold off
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        title(metadata.StimParams.directions(i))
        axis off
    end
    hold off
    set(gcf, 'Color', 'w')
end
function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% 
% [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging), 'qual' (qualitative)
%   - cname: name of colortable. It changes depending on ctype.
%   - ncol:  number of color in the table. It changes according to ctype and
%            cname
%   - interp_method: interpolation method (see interp1.m). Default is "cubic" )
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of 
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types and
%           names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
% ------------------------------
% 18.09.2015  Minor fixes, fixed a bug where the 'spectral' color table did not appear in the preview


    % load colorbrewer data
    load('C:\Users\jaepelj\Dropbox\Work\colorbrewer.mat')
    % initialise the colormap is there are any problems
    colormap=[];
    if (~exist('interp_method', 'var'))
        interp_method='cubic';
    end

    % If no arguments
    if (~exist('ctype', 'var') | ~exist('cname', 'var') | ~exist('ncol', 'var'))
        disp(' ')
        disp('[colormap] = cbrewer(ctype, cname, ncol [, interp_method])')
        disp(' ')
        disp('INPUT:')
        disp('  - ctype: type of color table *seq* (sequential), *div* (divergent), *qual* (qualitative)')
        disp('  - cname: name of colortable. It changes depending on ctype.')
        disp('  - ncol:  number of color in the table. It changes according to ctype and cname')
        disp('  - interp_method:  interpolation method  (see interp1.m). Default is "cubic" )')

        disp(' ')
        disp('Sequential tables:')
        z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
                 'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
        disp(z')     

        disp('Divergent tables:')
        z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
        disp(z')

        disp(' ')
        disp('Qualitative tables:')
        %getfield(colorbrewer, 'qual')
        z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
        disp(z')

        plot_brewer_cmap
        return
    end

    % Verify that the input is appropriate
    ctype_names={'div', 'seq', 'qual'};
    if (~ismember(ctype,ctype_names))
        disp('ctype must be either: *div*, *seq* or *qual*')
        colormap=[];
        return
    end

    if (~isfield(colorbrewer.(ctype),cname))
        disp(['The name of the colortable of type *' ctype '* must be one of the following:'])
        getfield(colorbrewer, ctype)
        colormap=[];
        return
    end

    if (ncol>length(colorbrewer.(ctype).(cname)))
    %     disp(' ')
    %     disp('----------------------------------------------------------------------')
    %     disp(['The maximum number of colors for table *' cname '* is ' num2str(length(colorbrewer.(ctype).(cname)))])
    %     disp(['The new colormap will be extrapolated from these ' num2str(length(colorbrewer.(ctype).(cname))) ' values'])
    %     disp('----------------------------------------------------------------------')
    %     disp(' ')
        cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
        colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
        colormap=colormap./255;
        return
    end

    if (isempty(colorbrewer.(ctype).(cname){ncol}))

        while(isempty(colorbrewer.(ctype).(cname){ncol}))
            ncol=ncol+1;
        end        
        disp(' ')
        disp('----------------------------------------------------------------------')
        disp(['The minimum number of colors for table *' cname '* is ' num2str(ncol)])
        disp('This minimum value shall be defined as ncol instead')
        disp('----------------------------------------------------------------------')
        disp(' ')
    end

    colormap=(colorbrewer.(ctype).(cname){ncol})./255;
 end
function PlotTrialStimResponse(metadata, analysis, field, roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);    
    colorlevels_tr = metadata.StimParams.numTrials;
    coc_tr = cbrewer('div', 'Spectral', colorlevels_tr);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).stimResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).stimResponseTrace(:));
    i = 1;
    for elev=1:metadata.StimParams.numElevation
        for azi = 1:metadata.StimParams.numAzimuth
            ax=subplot(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, i);
            cla(ax)
            y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
            x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
            if stimWindow(1)~=0
                    bl=x(stimWindow);
                    patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
            end
            hold all
            for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                plot(ax,x, smooth(squeeze(analysis.(field).roi(roi).stimResponseTrace(i,trial,:))), 'Color', coc_tr(trial,:));
                hold all
            end
            plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(i,:)), 'k');
            hold all
            ylim([ymin, ymax])
            xlim([min(x) max(x)])
            axis off
            if metadata.StimParams.onlyAzi == 0
                title([num2str(metadata.StimParams.stimPosY(elev)) ' deg Elev, ' num2str(metadata.StimParams.stimPosX(azi)) ' deg Azi'])
            else
                title([num2str(metadata.StimParams.stimPosX(azi)) ' deg Azi'])
            end
            i = i +1;
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end