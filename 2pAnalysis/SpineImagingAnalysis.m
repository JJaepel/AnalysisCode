function SpineImagingAnalysis(varargin)
% function SpineImagingAnalysis(type)
% 
% This function goes through all the steps of the spine imaging analysis,
% from the registration until the individual imaging analysis, including:
% - 0. Define paramenters on what to do
% - 1. Listing all experiments
% - 2. Registration -> Registered tiffs
% - 3. Drawing ROIs -> RoiSet.zip
% - 4. Go through one experiment after the other:
%    - a. Get the info about the experiment -> expInfo.mat
%    - b. Extract the traces -> ROIs.mat
%    - c. Chopping traces accordint to simulus information -> StimInfos.mat
%    - (d. Dendritic substraction -> ROIsChopped & events -> ROIsEvents.mat if spine)
%    - e. Calculate responses, plot & evalute if spines -> ROIsAna.mat
%    - f. Plot population data: ROIs on template, prefOri/dir on template, population metrics 
%    - g. Copy to a differnt folder if necessary
%
% Inputs
% - (type): which analysis steps do you want to do?

if nargin > 0
    if ~isempty(varargin{1})
        type = varargin{1};
    else
        type = 'analysis';
    end
else
    type = 'analysis';
end
if nargin > 1
    animal = varargin{2};
else 
    animal = [];
end

%% 0.) Switch board for analysis variables based on run type
%switch the analysis to fit for different runs -> this is basic settings
%for a normal run 
analysisParams = struct;

%where does it run?
analysisParams.server =1; %load from the server (1) or the raid (0)

%what does it run?
analysisParams.select = 1;
analysisParams.allInclude = 1;

%what should it do with them?
analysisParams.reregister =0; %should you reregister the data
analysisParams.reloadData = 0; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse = 1; %should you reanalyse the data or just plot?
analysisParams.plotROIs = 0;  %should you plot traces for all resp ROIs?
switch type
    case 'first'
        %we assume that this one we will run from scratch, so run from raid
        %and do everything starting from registration
        analysisParams.server = 0;
        %analysisParams.reregister =1;
        analysisParams.reloadData = 1;
        analysisParams.plotROIs = 1; 
    case 'reload'
        %let's reload - assumes that it is already copied on the server
        analysisParams.reloadData = 1;
        analysisParams.plotROIs = 1; 
    case 'allData'
        analysisParams.select = 0; %run through all data
        analysisParams.allInclude = 0;
    case 'plotROIs'
        analysisParams.reanalyse = 0; %don't reanalyze, just plot
        analysisParams.plotROIs = 1; %plot the ROIs
    case 'allExp'
        analysisParams.select = 0; %quickly run through all experiment inlcuded in the dataset
    case 'figures'
        analysisParams.reloadData = 0;
        analysisParams.plotROIs = 0;
end
analysisParams.type = type;

%General settings, independent from stim type: Stim infos
analysisParams.prestimPeriod = 0;    
analysisParams.stimDur = 2.5; %stimduration
analysisParams.postPeriod = 0.5; 

%% 1.) List all experiments
%where do we find the experiments?
filePath =  'Z:\Juliane\Organization\Animals\';
file = 'SpinePerAnimal.xlsx';
sheetName = 'bimodal';

%get all the information about the experiments
expInfo = findExpInfo([filePath file], sheetName);

%which experiments do we want to run through?
if analysisParams.select
    if ~isempty(animal) %run through all experiments from one animal
        ind = find(contains(string(expInfo.animal), animal));
    else
        ind = find(expInfo.run); %run through selected experiments
    end
else
    if analysisParams.allInclude
        ind = find(expInfo.include); %run through all included experiments
    else
        ind = 1:1:length(expInfo.animal); %run through all experiments at once
    end
end

%% 2.) Registration:
% Check if data is registered or if you want to re-register -> run as
% 'first'
% do for all data together as it is time intensive
getSubcellularRegistration(analysisParams, expInfo, ind)

%% 3.) Make ROIs if not there already
% do for all data together as it needs manual input
drawROIsSpines(analysisParams, expInfo, ind)

%% 4.) Go through one experiment after the other
for j = ind
    %% 4.a) Get the info about the experiment
    info = makeInfoStructExp(analysisParams.server, expInfo, j);
    
    %% 4.b) Extract the traces or load them from other folder
    if analysisParams.reloadData
        % a.) Get raw responses and save everything in a structure
        [ce, info.template] = getROIRawData(info);

        % b.) get deltaF/F responses
        ce = computeDffSpinesNew(info, ce);
    else
        temp = load([info.saveDir 'ROIs.mat']);
        temp2 = load([info.saveDir 'meanImg.mat']);
        ce = temp.ce; clear temp
        info.template = temp2.meanImg; clear temp2
    end

    %% 4.c) Chop the traces according to stimulus infomation
    % - Get the stimulus information
    metadata = getMetadataAndTimes(analysisParams,info);

    % - Do the chopping
    ce = chopSpineTraces(ce, metadata);
    
    %% 4.d) Remove dendritic signal from traces if it is a spine
    if ~info.isCell
        analysisParams.field = 'rawRes';
        ce = dendriteSubstraction(ce,info);
        ce = extractCaEvents(ce, info);
    end
    
    %% 4.e) Calculate responses, plot them & evaluate spines
    ce = calcOriParams(ce,metadata,info,analysisParams.plotROIs);
    if ~info.isCell
        ce = evalSpines(ce, info);
    end
    
    %% 4.f) Plot population data
    plotPosROIs(ce, info) %plots position of  ROIs on template
    plotPrefROIs(ce,info) %plot oriPref & dirPref on template
    plotMetricesROIs(ce, info)
    close all
    
    %% 4.g) If you need to copy the data to other cells, do this here
    if info.copyToOtherCells
        for i = 2:length(info.addSaveDirs)
            copyfile(info.addSaveDirs{1}, info.addSaveDirs{i});
        end
    end
end