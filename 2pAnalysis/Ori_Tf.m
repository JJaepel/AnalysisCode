close all
clear all
animal = 'F2290_2019-01-17\';
TwoPhontondir = 'E:\Data\2P_Data\';
Sp2dir = 'E:\Data\Spike2Data\';
savedir = 'E:\Data\ImageAnalysis\';
expt_id='t00017';

slicenum = 1;
totalSlices=1;
channelLoad= 1;
roiFiles='';
reloadData = 0;
windowStop=2;
windowStart=0;
pre=1;
field = 'dff';

base2pDirectory= [TwoPhontondir animal];
Sp2dDirectory = [Sp2dir animal filesep expt_id filesep];
saveDirectory = [savedir animal filesep expt_id filesep];
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
if ~exist(ROIsaveDirectory, 'dir')
    mkdir(ROIsaveDirectory);  
end

%%load Data and metadata
if reloadData
    sliceparams = struct;
    sliceparams.expt_id = expt_id;
    sliceparams.channel = channelLoad ;
    sliceparams.slicenum = slicenum;
    sliceparams.totalSlices = totalSlices;
    sliceparams.baseDirectory = base2pDirectory;
    metadata.StimParams=Load_stimparams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    metadata.StimParams.series=expt_id;
    data = Load_rois(sliceparams, roiFiles);
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');
    metadata.ROI = struct;
    analysis = struct;
    save(fullfile(saveDirectory, 's1_ori_tf.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
else
    load(fullfile(saveDirectory, 's1_ori_tf.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
end

%%chop traces
disp('Chopping Traces')
[analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, 'pre', pre, 'post',metadata.StimParams.isi,'windowStart',windowStart, 'windowStop',windowStop);

%%find maxResponses and significant responses
sigVector=zeros(1,length(analysis.(field).roi));
nStims= metadata.StimParams.numOrientations;
metadata.StimParams.theta= [0:2*pi/(nStims) : 2*pi-2*pi/(nStims)];
for i = 1:length(analysis.(field).roi)
    pretrialTime= analysis.(field).preTrialTime;
    preTrialIndex= (1:floor(pretrialTime * metadata.TwoPhoton.rate));
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    %collect our pretrial interval
    analysis.(field).roi(i).isRespSignificant = false;
    analysis.(field).roi(i).respThreshold = [];
    baselines=analysis.(field).roi(i).stimResponseTrace(:,:,preTrialIndex);
    analysis.(field).roi(i).baselineSD = std(baselines(:,:),[],2);
    analysis.(field).roi(i).baselineMean = mean(baselines(:,:),2);
    stimResp = analysis.(field).roi(i).avgResponseTrace;
    analysis.(field).roi(i).maxResp = max(stimResp, [],2);
    respThreshold = NaN * ones(size(stimResp,1));
    for stim=1:size(stimResp,1)
        respThreshold(stim)= 2*analysis.(field).roi(i).baselineSD(stim) + analysis.(field).roi(i).baselineMean(stim);
        largerThanFivePercent=max(stimResp(stim,:))> (.05+ analysis.(field).roi(i).baselineMean(stim));
        largerThanTwoSDs= (max(medfilt1(stimResp(stim,:),3)) >  respThreshold(stim));
        sigFlag= largerThanFivePercent && largerThanTwoSDs;
        isRespSignificant = analysis.(field).roi(i).isRespSignificant || sigFlag;
        analysis.(field).roi(i).isRespSignificant = isRespSignificant;
    end
    analysis.(field).roi(i).respThreshold = respThreshold;
    sigVector(i)= analysis.(field).roi(i).isRespSignificant;
end
analysis.(field).significantResponses=sigVector;

%%find preferred orientation for all cells
theta = metadata.StimParams.theta;
theta=mod(theta, pi);
lastOrientation = metadata.StimParams.numOrientations;
orientations = linspace(0,180,nStims/2+1);
orientations = orientations(1:end-1);
direction = linspace(0,360, nStims+1);
metadata.StimParams.direction = direction(1:end-1);
for i = 1:length(data.roi)
    %collapse all sf for the same orientation 
    tempFreq_cell = strsplit(metadata.StimParams.temporalFreq(2:end-1), ',');
    numTf = size(tempFreq_cell,2);
    metadata.StimParams.numTf = size(tempFreq_cell,2);
    medResponse_ori = zeros(metadata.StimParams.numOrientations,metadata.StimParams.numTrials-1,size(analysis.(field).roi(i).stimResponseTrace,3));
    endnum = metadata.StimParams.numOrientations * (metadata.StimParams.numTrials-1);
    for ori = 1:metadata.StimParams.numOrientations
        medResponse_ori(ori,:,:) = mean(analysis.(field).roi(i).stimResponseTrace(ori:metadata.StimParams.numOrientations:endnum, 1:metadata.StimParams.numTrials-1, :));
    end
    %calculate parameters
    medResponse_ori = mean(medResponse_ori(:,:,stimWindow),3);
    medResponseA(1:size(medResponse_ori,1)/2, 1:size(medResponse_ori,2)) = medResponse_ori(1:size(medResponse_ori,1)/2,:);
    medResponseA(1:size(medResponse_ori,1)/2, 1+size(medResponse_ori,2):2*size(medResponse_ori,2))= medResponse_ori(1+size(medResponse_ori,1)/2:size(medResponse_ori,1),:);
    medResponseA= median(medResponseA,2);
    medResponse_ori = medResponseA(:)';
    
    theta=theta(1:lastOrientation/2);
    theta= repmat(theta, 1, length(medResponse_ori)/length(theta));
    
    [analysis.(field).roi(i).preferredOrientation,analysis.(field).roi(i).coeffOr,analysis.(field).roi(i).rsqOr] = ComputePreferredOrientations(medResponse_ori, theta);
    analysis.(field).roi(i).preferredOrientation =mod(analysis.(field).roi(i).preferredOrientation, 180);
    analysis.(field).roi(i).rsqOr= Rsquared(medResponse_ori, vonMisesFit(analysis.(field).roi(i).coeffOr, theta*2), true);
    analysis.(field).roi(i).OSIFit = computeOSIFit(analysis.(field).roi(i).coeffOr);
    [~, ind] = min(abs(orientations-analysis.(field).roi(i).preferredOrientation));
    analysis.(field).roi(i).prefOriStim = orientations(ind);
    analysis.(field).roi(i).prefOriStimInd = ind;
    medResponse_ori = mean(analysis.(field).roi(i).avgStimResponse, 2)';
    analysis.(field).vsOrientation(i)= wrapTo2Pi(angle(vectorSum(abs(medResponse_ori(1:lastOrientation)),2)))/2 * 180/pi;
    analysis.(field).orientationMag(i) = abs(vectorSum(abs(medResponse_ori(1:lastOrientation))/max(medResponse_ori(1:lastOrientation)),2));                
end

%%find preferred tF for each cell by looking at preferred Orientation
for i = 1:length(data.roi)
    Response_tf_prefOri = zeros(numTf,metadata.StimParams.numTrials-1,size(analysis.(field).roi(i).stimResponseTrace,3));
    for tf = 1:numTf
        temp = analysis.(field).roi(i).stimResponseTrace(1+(tf-1)*ori:tf*ori, 1:metadata.StimParams.numTrials-1, :);
        Response_tf_prefOri(tf, :, :) = (temp(analysis.(field).roi(i).prefOriStimInd,:,:)+temp(analysis.(field).roi(i).prefOriStimInd+nStims/2,:,:)/2);
    end
    avgResponse_tf_prefOri = mean(Response_tf_prefOri(:,:,stimWindow),3);
    medResponse_tf_prefOri = mean(avgResponse_tf_prefOri,2)';
    semResponse_tf_prefOri = (std(avgResponse_tf_prefOri,[],2)/sqrt(size(medResponse_tf_prefOri,2)))';
    PrefTFInd = find(medResponse_tf_prefOri == max(medResponse_tf_prefOri));
    analysis.(field).roi(i).prefTfStim = PrefTFInd;
    analysis.(field).roi(i).prefTf = tempFreq_cell{PrefTFInd};
    %plot(cellfun(@str2double,temp_cell),medResponse_tf_prefOri)
end

%%load template
imgFileName = [expt_id filesep 'Registered' filesep 'slice1' filesep 'slice1_c' num2str(slicenum) '_MasterTemplate.tif'];
imgFile=fullfile(base2pDirectory,imgFileName);
template=read_Tiffs(imgFile,1);
template = double(template);
data.template = template./prctile(template(:),99.9);

%plot prefOri on top of template
h=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
PlotPrefOnTemplate(analysis, data, metadata,1, field,template)
subplot(2,2,2)

alloriprefs = [analysis.dff.roi.preferredOrientation];
histogram(alloriprefs,linspace(0,180,5));
title('Histogram')
ylabel('Cells');
xlabel(sprintf('Orientation preference (%s)',char(145)));
xlim([-22.5 (360+22.5)]/2)
axis square;
set(gca,'Box','off');

%plot prefTF on top of template
subplot(2,2,3)
PlotPrefOnTemplate(analysis, data, metadata,2, field,template)

subplot(2,2,4)
allsfpref = [analysis.dff.roi.prefTfStim];
numTf = size(tempFreq_cell,2);
tf_cat = cellfun(@str2double,tempFreq_cell);
tf_counts = histcounts(allsfpref, numTf);
bar(tf_cat, tf_counts);
title('Histogram')
ylabel('Cells');
xlabel(sprintf('Temporal frequency preference'));
axis square;
set(gca,'Box','off');
saveas(gcf, fullfile(saveDirectory, 'Overlaymaps.png'))
close all

%plot all ROIs
for i = 1:length(data.roi)
    PlotAvgStimResponse(metadata, analysis, field, i)
    saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
    close gcf
end
for i = 1:length(data.roi)
    PlotTrialStimResponse(metadata, analysis, field, i)
    saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
    close gcf
end


function [data]= Load_rois(sliceparams, roiFiles)
    %IRConfigs = ConfigurationParser.ReadImageRegistrationConfig;     
    slicenum = sliceparams.slicenum;
    baseDirectory = sliceparams.baseDirectory;
    channelLoad = sliceparams.channel;
    expt_id = sliceparams.expt_id;

    RM= ij.plugin.frame.RoiManager;
    RC= RM.getInstance();
    count= RC.getCount();
    if count > 0
        RC.setSelectedIndexes([0:count-1]);
        RC.runCommand('Delete');
    end
    roiFileName = [expt_id filesep 'Registered' filesep 'slice' num2str(slicenum) filesep 'RoiSet.zip'];
    
    if isempty(roiFiles)
        roiFiles= dir(fullfile(baseDirectory, roiFileName));
        roiFiles= {roiFiles.name};
        if ~isempty(roiFiles)
            roiFile=fullfile(baseDirectory,roiFileName);
        end
    end
    disp(['Loading ROI from: ', roiFile]);

    if RC.getCount() == 0 && isempty(roiFile)
        % check if there is a zip file in the expected directory
        fprintf('Could not find ROI file in %s', roiPath, '\n');
        fprintf('Load files for %s', expt_id, '\n');
        RC.runCommand('Open','');
    elseif RC.getCount() == 0
        RC.runCommand('Open', roiFile)
    end

    if RC.getCount() == 0
        error('There are still no ROIs defined')
    end
    
    metadata.File=struct;
    data.roi=[];
            
    for chIdx=1:length(channelLoad)
        channel=channelLoad(chIdx);
        sliceparams.channel= channel;
        metadata.File.baseDirectory = baseDirectory;
        metadata.File.series = expt_id;
        metadata.File.slicenum = slicenum;
        registeredImageDirectory = [baseDirectory filesep expt_id filesep 'Registered' filesep 'slice' num2str(slicenum) filesep];
        metadata.File.ImgData = registeredImageDirectory;
        registeredFileName = 'stack*.tif';

        disp(['Loading from ...' registeredImageDirectory])
        disp(['Looking for...' fullfile(registeredImageDirectory, registeredFileName)])
        files= dir(fullfile(registeredImageDirectory, registeredFileName));
        files=sort({files.name}'); 

        for fnum = 1:size(files,1)
           files{fnum} = fullfile(registeredImageDirectory, char(files{fnum})); 
        end

        NumberofTifs= size(files, 1);
        disp(['Found ' num2str(NumberofTifs) ' files'])
        % this will save the ROIs to a specified folder
        
        count = RC.getCount();
        for i =1:count
            tempROI= RC.getRoi(i-1);
            tempPos= tempROI.getContourCentroid();
            data.roi(i).xPos= tempPos(1);
            data.roi(i).yPos= tempPos(2);
            data.roi(i).mask= tempROI.getContainedPoints();
            data.roi(i).name= tempROI.getName;
        end

        %% Export ROI fluorescence traces
        for i=1:NumberofTifs
            if i==1
                for j=1:RC.getCount()
                    if channel ==1
                        data.roi(j).rawF=[];
                    end
                end
            end
            if i < 10
                registeredFileName = ['stack_c' num2str(slicenum) '_0' num2str(i) '.tif'];
            else
                registeredFileName = ['stack_c' num2str(slicenum) '_' num2str(i) '.tif'];
            end

            filename = fullfile(registeredImageDirectory, registeredFileName);
            if(~isempty(find(ismember(files, filename), 1)))
                data.roi=ExtractRawFMIJ(RC, data, filename, channel);
            end
        end
    end
end
function [cells]=ExtractRawFMIJ(RC, data,fullfilename, channel)
    %ExtractRawF generates an array of structs that contains the centroid of
    %ROI and the Raw fluorscence
    %   Args:
    %           RC: handler for ROI Manager
    %           cells(array): array of strucs
    %   Returns:
    %           cells(array): reassignments as below
    %               cells(i).name: ROI name
    %               cells(i).x: x location
    %               cells(i).y: y location
    %               cells(i).name: ROI name
    %               cells(i).name: ROI name
    %ij.IJ.run('Image Sequence...', strrep(sprintf('open=%s sort',filename),'\','\\'));
    cells=data.roi;
    disp(strcat('Loading...', fullfilename));
    MIJ.run('Open...', strrep(sprintf('path=%s', fullfilename), '\', '\\'));
    count = RC.getCount();
    targetField='rawF';
    if channel > 1
        targetField=strcat(targetField, '_c', num2str(channel));
    end
    for i = 1:count
        RC.select(i-1)
        cells(i).name = char(RC.getName(i-1));
        % This is the way the position was extracted by using ROI names
        %this method figures out x,y position by computing the centroid
        % Pull the fluorescence plot
        if ~isfield(cells(i), targetField)
            cells(i).(targetField)=[];
        end

            a=ij.plugin.ZAxisProfiler();
            a.run('ZPlot');
            zplot=a.getPlot();

            tmp = double(zplot.getYValues());
            MIJ.run('Close','');
            cells(i).(targetField) = [cells(i).(targetField), tmp(:,1)'];
    end
    MIJ.run('Close', '');
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
        case 'driftingGrating_ori_tf'
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
    windowStopIdx = round(windowStop * TwoPhotonRate + offsetPre);
    if windowStartIdx == windowStopIdx
        windowStopIdx = analysis(tfield).stimStop;
    end
    analysis.(tfield).windowStart =windowStartIdx;
    analysis.(tfield).windowStop = windowStopIdx;

    %do a little clean up here to make stim numbers match
    while (metadata.StimParams.StimOnTimes(2,metadata.StimParams.numberOfStims)+metadata.StimParams.stimDuration > metadata.TwoPhoton.time(end))
        if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ==0
            metadata.StimParams.numTrials = metadata.StimParams.numTrials -1;
        else
            metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
        end
        metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
    end

    if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ~=0
        metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
        metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
    end

    for i =1:metadata.StimParams.numberOfStims
        metadata.StimParams.stimStartIndex(i)= find(metadata.TwoPhoton.time > metadata.StimParams.StimOnTimes(2,i),1);
        metadata.StimParams.stimStopIndex(i)= floor(metadata.StimParams.stimStartIndex(i)+ metadata.StimParams.stimDuration* metadata.TwoPhoton.rate);
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
        for trialNumber= 1:metadata.StimParams.numTrials-1 %exception!
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
        for stimID = 1:length(stimIndices)
            analysis.(tfield).roi(i).avgStimResponse(stimID) = analysis.(tfield).roi(i).avgStimResponse(stimID) -analysis.(tfield).roi(i).avgStimResponse(end) ;
            analysis.(tfield).roi(i).avgResponseTraceRaw(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:);
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:) - analysis.(tfield).roi(i).avgResponseTrace(end,:);
            analysis.(tfield).roi(i).Fb(stimID,:)= mean(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,1:analysisPeriod(1)-1),3);

        end
    end
    analysis.(tfield).preTrialTime=pre;
end        
function [prefOri, coeffOr, rsqPOr] = ComputePreferredOrientations(medResponse, theta)
    theta = theta *2;
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    coeffOr=[NaN, NaN, NaN, NaN];
    [coeffOr, rsqPOr]=fitVonMisesFunction(medResponse, theta, coeffOr, options);
    prefOri = rad2deg(coeffOr(3)/2);
    hold off;
    plot(rad2deg(theta/2), medResponse, 'bo')
    hold on;
    plot(rad2deg([0:pi/32:pi]), vonMisesFit(coeffOr, [0:pi/16:2*pi]), 'r');
    title([num2str(prefOri),':',num2str(rsqPOr)]);
    drawnow()
end
function [bestCoeff, bestrsq]=fitVonMisesFunction(ydata, tdata, bestCoeff, options)
    stims=length(ydata);
    fun=@(x, tdata) vonMisesFit(x,tdata);

    prefdir=mod(tdata(ydata == max(ydata(:))), pi);
    if length(prefdir) >1
        prefdir=prefdir(1);
    end
    %coeff0=[abs(max(ydata(:))), abs(max(ydata(:))), pi/2, pi/stims];
    
    lb=[0,0,0,0];
    ub=[abs(1.5*max(ydata(:))),abs(max(ydata(:))*2),2*pi,2.5*pi];
    coeff0=(ub-lb)./2 + lb;
    rsq=0;
    bestrsq=0;
    
    iterations=0;
    try
    if sum(isnan(lb))==0 && sum(isnan(ub)) ==0
        while rsq < 0.7 && iterations <5
            

            [currentCoeff, ~, ~,exitflag]=lsqcurvefit(fun, coeff0, tdata, ydata, lb, ub, options);

            [rsq, ~]=Rsquared(ydata, vonMisesFit(currentCoeff, tdata), true);

             
            if rsq > bestrsq
                bestCoeff= currentCoeff;
                bestrsq=rsq;
            end
            iterations=iterations+1;
            coeff0= ub-lb .* rand(1) + lb;
        end
        %disp(['Ran ', num2str(iterations),  ' iterations']);
    else
        bestCoeff=[NaN, NaN, NaN, NaN];
    end 
    catch
        bestCoeff=[NaN,NaN,NaN,NaN];
        rsq=0;
    end
        
end
function [OSIFit] = computeOSIFit(coeffOr)
    Rpref= vonMisesFit(coeffOr, coeffOr(3));
    Rnull = vonMisesFit(coeffOr, coeffOr(3) + pi);
    OSIFit = (Rpref - Rnull) / (Rpref + Rnull);
end
function [ out ] = vonMisesFit(x, tdata)
    %VONMISESFIT Linear Sum of two vonMises functions that's constrained to be
    %pi radians apart
    %   x is a vector of the parameters of the linear sum and vonMises
    %   parameters
    %   tdata is the x range over which to compute the vonMises
    %   out = [Alpha Beta mu kappa]
    if sum(isnan(x)) > 0 || length(x) <4
        out = NaN * ones(size(tdata));
        return
    end
    A= x(1);
    B= x(2);  % a dc component
    mu1=x(3); % center
    kappa1=x(4); % width
    out = A*vonMisesFunction(kappa1, mu1, tdata ) + B;
end
function PlotPrefOnTemplate(analysis, data, metadata,type, field,template)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,template,template,template)/prctile(template(:),99));
    colormap(hsv) 
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        title('Temporal frequency preference map')
        numtf = (metadata.StimParams.uniqStims-1)/metadata.StimParams.numOrientations;
        LUT = hsv(numtf);
        caxis([0 numtf]);colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(data.roi)
        xpos= data.roi(i).xPos;
        ypos= data.roi(i).yPos;
        try
            if type == 1
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(1+floor(analysis.(field).roi(i).preferredOrientation),:));
            elseif type == 2
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(i).prefTfStim,:));
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
    colorlevels = metadata.StimParams.numTf;
    coc_tf = cbrewer('qual', 'Set1', colorlevels);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    for i=1:metadata.StimParams.numOrientations
        ax=subplot(1,metadata.StimParams.numOrientations, i);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
            bl=x(stimWindow);
            patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        for tf = 1:metadata.StimParams.numTf
            id = metadata.StimParams.numOrientations*(tf-1)+i;
            err=analysis.(field).roi(roi).SEMResponseTrace(id,:);
            hold all
            patch([x fliplr(x)],[y+err fliplr(y-err)],[.5 .5 .5], 'LineStyle', 'none');
            plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', coc_tf(tf,:));
            hold all
        end
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        title(metadata.StimParams.direction(i))
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
    colorlevels = metadata.StimParams.numTf;
    colorlevels_tr = metadata.StimParams.numTrials;
    coc_tf = cbrewer('qual', 'Set1', colorlevels);
    coc_tr = cbrewer('seq', 'Greys', colorlevels_tr);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    for tf = 1:metadata.StimParams.numTf
        for i=1:metadata.StimParams.numOrientations
            id = metadata.StimParams.numOrientations*(tf-1)+i;
            ax=subplot(metadata.StimParams.numTf,metadata.StimParams.numOrientations, id);
            cla(ax)
            y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(id,:),3);
            x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
            if stimWindow(1)~=0
                    bl=x(stimWindow);
                    patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
            end
            hold all
            for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                plot(ax,x, squeeze(analysis.(field).roi(roi).stimResponseTrace(id,trial,:)), 'Color', coc_tr(trial,:));
                hold all
            end
            plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', coc_tf(tf,:), 'LineWidth', 3);
            hold all
            ylim([ymin, ymax])
            xlim([min(x) max(x)])
            title(metadata.StimParams.direction(i))
            axis off
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end