%% switch board for analysis variables
analysisParams = struct;
analysisParams.select = 1;
analysisParams.server =1; %load from the server (1) or the raid (0)
analysisParams.field = 'rawRes';

%stim infos
analysisParams.prestimPeriod = 0;    
analysisParams.stimDur = 2.5; %stimduration
analysisParams.postPeriod = 0.5; 

%what should it do?
analysisParams.reregister =0; %should you reregister the data
analysisParams.reloadData = 0; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse = 1; %should you reanalyse the data or just plot?
analysisParams.plotROIs = 0;   %should you plot traces for all resp ROIs?

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% 0.) List all experiments
filePath =  'Z:\Juliane\Organization\Animals\';
file = 'SpinePerAnimal.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'bimodal');

exp_info = findExpInfo(xls_txt, xls_all);

if analysisParams.select
    ind = find(exp_info.run);           %run through selected experiments
else
    ind = 1:1:length(exp_info.animal); %run through all experiments at once
end

%only do anlysis on spine data -> remove cell imaging
onlyCells = cellfun(@(x) find(contains(x, 'cells')),exp_info.region,'UniformOutput',false);
onlySpines = cellfun(@isempty,onlyCells); 
SpineInd = intersect(ind,find(onlySpines)); 

%% 1.) Registration:
% Check if data is registered or if you want to re-register -> set analysisParams.reregister = 1
% do for all data together as it is time intensive
getSubcellularRegistration(analysisParams, exp_info, ind)

%% 2.) Make ROIs if not there already
% do for all data together as it needs manual input
drawROIsSpines(analysisParams, exp_info, ind)

%% 3.) Get ROI traces and do dff
for j = SpineInd
    %save most important info
    info = [];
    if analysisParams.server
        info.drive = 'Z:\Juliane\Data\';
    else
        info.drive = 'F:\Data\';
    end
    info.animal = char(exp_info.animal{j});
    info.name = char(exp_info.name{j});
    info.sp2 = char(exp_info.sp2_id{j});
    info.TwoPDir = [info.drive '2P_data\' char(info.animal) filesep char(info.name) '\Registered\'];
    info.Sp2Dir = [info.drive 'Spike2Data\' char(info.animal) filesep char(info.sp2) filesep];
    info.saveDir = [info.drive 'ImageAnalysis\' char(info.animal) filesep char(info.name) filesep];
    info.saveDirROI = [info.saveDir 'ROIs\'];
    info.zoom = getZoomSpines([info.drive '2P_data\' char(info.animal) filesep char(info.name)]);
    info.depth = num2str(exp_info.Depth(j));
    info.denType = char(exp_info.denType{j});
    info.Modality = char(exp_info.Modality{j});
    info.ROINr = num2str(exp_info.dendROIInd(j));

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
        temp = load([info.saveDir 'ROIs.mat']);
        temp2 = load([info.saveDir 'meanImg.mat']);
        ce = temp.ce; clear temp
        info.template = temp2.meanImg; clear temp2
    end

    %% 4.) Chop the traces according to stimulus infomation
    % a.) Get the stimulus information
    metadata = getMetadataAndTimes(analysisParams,info);

    % b.) Do the chopping
    ce = chopSpineTraces(ce, metadata);
    
    %% 5.) Remove dendrite signal from traces
    ce = dendriteSubstraction(ce,info);
    ce = extractCaEvents(ce, info);

    
    %% 6.) Calculate responses, plot them & evaluate spines
    ce = calcOriParamsSpines(ce,metadata,info,analysisParams.plotROIs);
    ce = evalSpines(ce, info);
    
    %% 7.) Plot population data
    plotPosSpines(ce, info) %plots position of spine ROIs on template
    plotPrefSpines(ce,info) %plot oriPref & dirPref on template
    plotMetricesSpines(ce, info)
    close all
end