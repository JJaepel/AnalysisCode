function SpontaneousCells(analysisParams)

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
    try
        metadata.StimParams=LoadStimParams(Sp2dDirectory);
        metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
        metadata.StimParams.path=fullfile(Sp2dDirectory);
    catch
        metadata.TwoPhoton.rate = 30;
        metadata.StimParams = [];
    end
    
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
    save(fullfile(saveDirectory, 'spontaneous.mat'), 'data', 'metadata', 'analysisParams', 'analysis');

else
    load(fullfile(saveDirectory, 'spontaneous.mat'), 'data', 'metadata', 'analysis');
end

%% 2.) Calculate single cell metrics

if analysisParams.reanalyse
    for i = 1:length(data.roi)
        %get active events
        globalThreshold = 0.5;
        analysis = findActivity(analysis, data, analysisParams.field, i, globalThreshold);
        %cells active or not
        if analysis.(analysisParams.field).roi(i).numEvents < 5
            analysis.(analysisParams.field).roi(i).isResp = 0;
        else
            analysis.(analysisParams.field).roi(i).isResp = 1;
        end
        
        %normalize responses for correlation
        analysis = normResponse(analysis, data, analysisParams.field, i);
    end
    
    % save file
    disp('Saving analyzed data')
    save(fullfile(saveDirectory, 'AnaData.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
else
    disp('Loading analyzed data')
    load(fullfile(saveDirectory, 'AnaData.mat'), 'data', 'metadata', 'analysis');
end

%% 3.) Calculate population indices
PopAnalysis = struct;
if analysisParams.reanalysePop
    %first get um per pixel for the experiment
    [umperpixel] = getzoom(tifDirectory); 
    
    % calculate distance vs. correlation
    PopAnalysis.(analysisParams.field).distROIs = calcDistROIs(analysisParams, data, analysis, umperpixel, analysisParams.level) ;
    PopAnalysis.(analysisParams.field).corrROIs = calcCorrROIs(analysisParams, analysis);
    [PopAnalysis] = calcCorrVsDist(PopAnalysis, analysisParams); 

    % save file
    disp('Saving analyzed population data')
    save(fullfile(saveDirectory, 'PopData.mat'), 'PopAnalysis');
else
    disp('Loading analyzed population data')
    load(fullfile(saveDirectory, 'PopData.mat'), 'PopAnalysis');
end

%% 4.) Plot
analysisParams.coc_prop = cbrewer('qual', 'Paired', 12);
showCorrelationROIs(analysisParams, PopAnalysis, saveDirectory)
plotCorrVsDistance(analysisParams, PopAnalysis, saveDirectory)
