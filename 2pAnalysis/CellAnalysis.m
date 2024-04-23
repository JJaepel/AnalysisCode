function SpineImagingAnalysis(spines, varargin)
% function SpineImagingAnalysis(type)
% 
% This function draws a polar plot.
%
% Inputs
% - (type): which analysis steps do you want to do?


if nargin > 1
    type = varargin{1};
else
    type = 'analysis';
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
        analysisParams.reregister =1;
        analysisParams.plotROIs = 1; 
    case 'reload'
        %let's reload - assumes that it is already copied on the server
        analysisParams.reloadData = 1;
        analysisParams.plotROIs = 1; 
    case 'allData'
        analysisParams.select = 0; %run through all data
        analysisParams.allInclude = 0;
    case 'plotROIs'
        analysisParams.plotROIs = 1; %analyze and plot the ROIs
    case 'allExp'
        analysisParams.select = 0; %quickly run through all experiment inlcuded in the dataset
end

%general settings, independent from stim type: stim infos
analysisParams.prestimPeriod = 0;    
analysisParams.stimDur = 2.5; %stimduration
analysisParams.postPeriod = 0.5; 

%% 1.) List all experiments
filePath =  'Z:\Juliane\Organization\Animals\';
file = 'SpinePerAnimal.xlsx';
sheetName = 'bimodal';
expInfo = findExpInfo([filePath file], sheetName);

if analysisParams.select
    ind = find(expInfo.run); %run through selected experiments
else
    if analysisParams.allInclude
        ind = find(expInfo.include); %run through all included experiments
    else
        ind = 1:1:length(expInfo.animal); %run through all experiments at once
    end
end

%only do anlysis on spine data -> remove cell imaging
onlyCells = cellfun(@(x) find(contains(x, 'cells')),expInfo.region,'UniformOutput',false);
onlySpines = cellfun(@isempty,onlyCells); 
CellInd = setdiff(ind,find(onlySpines)); 

%% 2.) Registration:
% Check if data is registered or if you want to re-register -> set analysisParams.reregister = 1
% do for all data together as it is time intensive
getSubcellularRegistration(analysisParams, expInfo, ind)

%% 3.) Make ROIs if not there already
% do for all data together as it needs manual input
drawROIsSpines(analysisParams, expInfo, ind)

%% 4.) Get ROI traces and do dff
for j = CellInd
    %save most important info
    info = [];
    if analysisParams.server
        info.drive = 'Z:\Juliane\Data\';
    else
        info.drive = 'Z:\Juliane\Data\';
    end
    
    info.isCell = 1;
    info.animal = char(expInfo.animal{j});
    info.name = char(expInfo.name{j});
    info.sp2 = char(expInfo.sp2_id{j});
    info.TwoPDir = [info.drive '2P_data\' char(info.animal) filesep char(info.name) '\Registered\'];
    info.Sp2Dir = [info.drive 'Spike2Data\' char(info.animal) filesep char(info.sp2) filesep];
    info.include =  num2str(expInfo.include(j));
    
    if str2double(info.include)
        %run the animal parser
        
        animal = char(info.animal);
        info.saveDir = ['Z:\Juliane\InputAnalysis\' animal(1:5) 'A - 2p Imaging\'  filesep char(info.name) filesep];
        info.TwoPDirROIs = info.saveDir;
    else
        info.saveDir = [info.drive 'ImageAnalysis\' char(info.animal) filesep char(info.name) filesep];
        info.TwoPDirROIs = [info.TwoPDir filesep 'slice1' filesep 'Projection\'];
    end
    info.saveDir2 = [info.drive 'ImageAnalysis\' char(info.animal) filesep char(info.name) filesep];
    info.saveDirROI = [info.saveDir 'ROIs\'];
    %info.zoom = getZoomSpines([info.drive '2P_data\' char(info.animal) filesep char(info.name)]);
    info.zoom = 3;
    info.depth = num2str(expInfo.Depth(j));
    info.denType = char(expInfo.denType{j});
    info.Modality = char(expInfo.Modality{j});
    info.ROINr = num2str(expInfo.cellROIInd(j));

    if ~exist(info.saveDir, 'dir')% make new file directory 
        mkdir(info.saveDir); 
    end
    if ~exist(info.saveDirROI, 'dir')% make new file directory 
        mkdir(info.saveDirROI); 
    end
    save([info.saveDir 'expInfo.mat'], 'info', '-mat') 
    
    if analysisParams.reloadData
        % a.) Get raw responses and save everything in a structure
        [ce, info.template] = getROIRawData(info);

        % b.) get deltaF/F responses
        ce = computeDffSpinesNew(info, ce);
    else
        temp = load([info.TwoPDirROIs 'ROIs.mat']);
        temp2 = load([info.saveDir2 'meanImg.mat']);
        ce = temp.ce; clear temp
        info.template = temp2.meanImg; clear temp2
    end

    %% 5.) Chop the traces according to stimulus infomation
    % a.) Get the stimulus information
    metadata = getMetadataAndTimes(analysisParams,info);

    % b.) Do the chopping
    ce = chopSpineTraces(ce, metadata);
    
    %% 6.) Calculate responses, plot them & evaluate spines
    ce = calcOriParamsCells(ce,metadata,info,analysisParams.plotROIs);
    
    %% 7.) Plot population data
    plotPosSpines(ce, info) %plots position of cell ROIs on template
    plotPrefCells(ce,info)
end