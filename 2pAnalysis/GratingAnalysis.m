function GratingAnalysis(analysisParams)

close all;

computer = getenv('COMPUTERNAME');
switch computer
    case 'DF-LAB-WS38'
        RaidDir = 'F:\';
        ServerDir = 'Z:\Juliane\';
    case 'DF-LAB-WS40'
        RaidDir = 'C:\';
        ServerDir = 'Z:\Juliane\';
end

%% 0). Defining folders and files
if analysisParams.server
    drive           = ServerDir;
else
    drive           = RaidDir;
end


TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

base2pDirectory= [TwoPhontondir analysisParams.animal];
tifDirectory = [base2pDirectory filesep analysisParams.name];
Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];

if ~exist(saveDirectory, 'dir')
    % make new file directory
    mkdir(saveDirectory); 
end

metadata = struct;
metadata.ROI = struct;
analysis = struct;

%% 1.) load data and metadata
disp('Loading data')
if analysisParams.reloadData
    % make sure that you also do analysis afterwards
    analysisParams.reanalyse = 1;
    
    % load stimulus parameters and spike2data
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    
    % load tif data
    switch analysisParams.dataType
        case 1
            if analysisParams.manual
                data = SpineROIExtractFct(analysisParams);    
            else
                data = LoadRoisS2p(analysisParams);
            end
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
                %etadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                %data = computeDff(data, 'rawF', 'baseline', 'dff');[m
                data = computeDffAxons(metadata, data);
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
    
    
    %save file
    switch analysisParams.stimType
        case 1
            save(fullfile(saveDirectory, 'oriGrating.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
        case 2
            save(fullfile(saveDirectory, 'oriSf.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
        case 5
            save(fullfile(saveDirectory, 'oriBino.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
    end

else
    switch analysisParams.stimType
        case 1
            load(fullfile(saveDirectory, 'oriGrating.mat'), 'data', 'metadata', 'analysis');
        case 2
            try
                load(fullfile(saveDirectory, 'oriSf.mat'), 'data', 'metadata', 'analysis');
            catch
                load(fullfile(saveDirectory, 's1_ori_sf.mat'), 'data', 'metadata', 'analysis');
            end
        case 5
            save(fullfile(saveDirectory, 'oriBino.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
    end
end

%% 2.) Calculate single cell metrics
if analysisParams.reanalyse
    % chop traces
    disp('Chopping Traces')
    [analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,analysisParams.level, analysisParams.field);
    
    %clean up properties
    metadata.StimParams.numStims=metadata.StimParams.uniqStims-1; %we have -1 number of stims (=Blank)
    try
        spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
        metadata.StimParams.numSf = size(spatFreq_cell,2);
        metadata.StimParams.spatialFreq = cellfun(@str2double,spatFreq_cell);
    catch
        metadata.StimParams.numSf = length(metadata.StimParams.spatialFreq);
    end
    try
        tempFreq_cell = strsplit(metadata.StimParams.tempFreq(2:end-1), ',');
        metadata.StimParams.numTf = size(tempFreq_cell,2);
        metadata.StimParams.temporalFreq = cellfun(@str2double,tempFreq_cell);
    catch
        metadata.StimParams.numTf = length(metadata.StimParams.temporalFreq);
    end
    if analysisParams.stimType == 5
        metadata.StimParams.numDirections = metadata.StimParams.numOrientations;
        metadata.StimParams.numCon=3;
    else
        metadata.StimParams.numDirections = metadata.StimParams.numStims/(metadata.StimParams.numSf*metadata.StimParams.numTf);
        metadata.StimParams.numCon = metadata.StimParams.numSf * metadata.StimParams.numTf;
    end
    metadata.StimParams.orientations = linspace(0,180,metadata.StimParams.numDirections/2+1); %this would give us a all orientations from 0 to 180
    metadata.StimParams.orientations = metadata.StimParams.orientations(1:end-1); %180 == 0, so we need to substract that one
    metadata.StimParams.directions = linspace(0,360,metadata.StimParams.numDirections+1); %this would give us a all orientations from 0 to 180
    metadata.StimParams.directions = metadata.StimParams.directions(1:end-1); %180 == 0, so we need to substract htat one
    metadata.StimParams.numOrientations = length(metadata.StimParams.orientations); %how many oris do we have?
    
    %calculate indices
    disp('Calculating indices')
    figure(9999)
    for i = 1:length(data.roi) %for each ROI
        [analysis, metadata] = calculateGratingMetrices(analysisParams, metadata, analysis, i);
    end

    % find maxResponses and significant responses
    disp('Calculating significant responses')
    disp(['Criteria: threshold (z  = ' num2str(analysisParams.zThresh) ') needs to be crossed ' num2str(analysisParams.fraction *100) ' % of all trials'])
    
    for i = 1:length(data.roi)
       analysis = classifyRespROIs(analysisParams, analysis, metadata, i);
    end
    
    % save file
    disp('Saving analyzed data')
    save(fullfile(saveDirectory, 'AnaData.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
else
    disp('Loading analyzed data')
    load(fullfile(saveDirectory, 'AnaData.mat'), 'data', 'metadata', 'analysis');
end

%% 3.) Plot single cell data
analysisParams.coc_prop = cbrewer('qual', 'Paired', 12);

% 3.1 Plot ROI positions
figure(101)
plotROIpositions(data, saveDirectory)

% 3.2 Plot prefOri on top of template
for types = 1:3
    %type 1: all ROIs, type 2: significant responders, type 3: selective
    %ones, OSI/DSI/SFSI/TFSI > 0.2
    makeOrientationPreferencePlots(analysisParams, analysis, data, metadata,types, saveDirectory)
end

if analysisParams.stimType == 2 || analysisParams.stimType == 3 || analysisParams.stimType == 5
    %plot preferred ori/dir for all conditions
    %type 1: ori pref, type 2: dir pref
    for types = 1:2
       makeGratingPreferencePlots(analysisParams, analysis, data, metadata,types, saveDirectory)
    end 
end

% 3.3 Plot percentages of all ROIs
plotPercentagesPiesOri(analysisParams, analysis, saveDirectory)

% 3.4 Plot OSI/DSI metrices
% for all cells
for types = 1:2
    %type 1: all ROIs, type 2: significant responders
    plotGratingMetrices(analysisParams, analysis, metadata,types, saveDirectory)
end

%% 4.) Calculate population indices
PopAnalysis = struct;
try
    if analysisParams.reanalysePop
        %first get um per pixel for the experiment
        [umperpixel] = getzoom(tifDirectory); 

        % calculate deltaOri vs. distance of ROIs
        disp('Calculating distance between ROIs, deltaOri and HI')
        PopAnalysis = calcDistDeltaOri(analysisParams, data, analysis, umperpixel, analysisParams.level);
        PopAnalysis.(analysisParams.field).distRespROIs = calcDistRespROIs(analysisParams, data, analysis, umperpixel, analysisParams.level) ;
        PopAnalysis.(analysisParams.field).distROIs = calcDistROIs(analysisParams, data, analysis, umperpixel, analysisParams.level) ;   
        PopAnalysis = calcMinDistances(PopAnalysis,analysisParams);
        if metadata.StimParams.numTf * metadata.StimParams.numSf > 1 || analysisParams.stimType == 5
            PopAnalysis.(analysisParams.field).HI_SF = calculateHI_SF(PopAnalysis,analysis, metadata, analysisParams.level, data, analysisParams.field);
        end

        % calculate HeterogeneityIndex
        PopAnalysis.(analysisParams.field).ori_cells.HomeogeneityIndex = calculateHI(PopAnalysis,analysis, analysisParams.level, data, analysisParams.field);

        disp('Calculating trial pattern')
        [PopAnalysis.(analysisParams.field).cellPatternSorted, PopAnalysis.(analysisParams.field).corrCoeff, PopAnalysis.(analysisParams.field).corrTrialsMatched, PopAnalysis.(analysisParams.field).corrTrialsOrtho] = trialPatternCorrelation(analysis, metadata, analysisParams.field);

        disp('calculating distance vs. correlation')
        [PopAnalysis.(analysisParams.field).shortDistanceCorr, PopAnalysis.(analysisParams.field).longDistanceCorr] = cellCorrelationDistance(analysis, data, metadata, analysisParams.field, umperpixel);

        % try to predict stimulus from the data
        if analysisParams.predictor
            disp('Predicting stimulus from cell responses')
            [PopAnalysis.(analysisParams.field).decoderData] = decoderOri(analysis, data, analysisParams.field);
        end

        % save file
        disp('Saving analyzed data')
        save(fullfile(saveDirectory, 'PopData.mat'), 'PopAnalysis');
    else
        disp('Loading analyzed data')
        load(fullfile(saveDirectory, 'PopData.mat'), 'PopAnalysis');
    end

    %% 5.) Plot population metrices
    plotPopulationMetrices(analysisParams, PopAnalysis, saveDirectory)
end

%% 6.) Plot ROIs
if analysisParams.plotROIs
    if ~exist(ROIRespsaveDirectory, 'dir') % make new file directory        
        mkdir(ROIRespsaveDirectory); 
    else
        removeOldFiles(ROIRespsaveDirectory) %remove old files
    end
    if ~analysisParams.plotRespROIsOnly
        if ~exist(ROINonRespsaveDirectory, 'dir')
            mkdir(ROINonRespsaveDirectory);
        else
             removeOldFiles(ROINonRespsaveDirectory)
        end
    end
    for i = 1:length(data.roi)
        if analysis.(analysisParams.field).roi(i).isResponseSignificant == 1
                PlotAvgStimResponseOri(metadata, analysis, analysisParams.field, i, ROIRespsaveDirectory)
                Plot2DResponseTrialTraces(metadata, analysis, analysisParams.field, i, ROIRespsaveDirectory)
        else
            if ~analysisParams.plotRespROIsOnly
                PlotAvgStimResponseOri(metadata, analysis, analysisParams.field, i, ROINonRespsaveDirectory)
                Plot2DResponseTrialTraces(metadata, analysis, analysisParams.field, i, ROINonRespsaveDirectory)
            end
        end
    end

end