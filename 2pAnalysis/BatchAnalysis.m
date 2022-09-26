%% switch board for analysis variable
analysisParams = struct;
%which type of stimulus should it run
analysisParams.stimType = 1; % 1 = driftingGrating, 2 = OriSf, 3 = OriTf, 4 = Patch, 5 = bino, 6 = spontaneous;
analysisParams.dataType = 3; %data type: 1 = cells, 2 = axons, 3 = spines

%what should it do?
analysisParams.reloadData =1; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse = 1; %should you reanalyse the data or just plot?
analysisParams.reanalysePop = 1;
if analysisParams.reanalyse 
    analysisParams.reanalysePop = 1;
end
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.plotROIs = 1;   %should you plot traces for all resp ROIs?
analysisParams.plotRespROIsOnly = 0; %should you also plot traces for all non-resp ROIs?
analysisParams.server =0; %load from the server (1) or the raid (0)
analysisParams.makeROIs = 1;

%analysisParameters
analysisParams.zThresh = 4;
analysisParams.fraction = 0.5;
analysisParams.predictor = 0;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff'; %rawRes for spines
analysisParams.manual = 0; %do manual ROIs instead of suite2p registration for cells
analysisParams.pre = 0;

%% list all experiments 

computer = getenv('COMPUTERNAME');
switch computer
    case 'DF-LAB-WS38'
        filePath = 'F:\Organization\Animals\';
    case 'DF-LAB-WS40'
        filePath = 'Z:\Hannah\';
end

switch analysisParams.dataType
    case 1
        file = '2pExpByStimulus.xlsx';
        if analysisParams.manual == 1
            Miji
        end
    case 2
        file = '2pExpByStimulusAxon.xlsx';
        Miji
    case 3
        file = '2pExpByStimulusSpines.xlsx';
        analysisParams.field = 'rawRes';
        Miji
end

switch analysisParams.stimType 
    case 1
        [~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
    case 2
        [~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating_ori_sf');
    case 3
        [~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating_ori_tf');
    case 4
        [~, xls_txt, xls_all]=xlsread([filePath file], 'Retinotopy'); %switch here between Patches and Retinotopy!
    case 5
        [~, xls_txt, xls_all]=xlsread([filePath file], 'binocular');
    case 6
        [~, xls_txt, xls_all]=xlsread([filePath file], 'spontaneous');
end

exp_info = findExpInfo(xls_txt, xls_all);

if analysisParams.select
    ind = find(exp_info.run);           %run through selected experiments
else
    ind = 1:1:length(exp_info.animal); %run through all experiments at once
end

%% for spine experiments, check beforehand if you have registered ROI data as that is manual input
if analysisParams.dataType == 3
    getSubcellularRegistration(analysisParams, exp_info, ind)
    getROIsSpines(analysisParams, exp_info, ind)
end

if analysisParams.dataType == 1 && analysisParams.manual == 1
    getSubcellularRegistration(analysisParams, exp_info, ind)
    getROIsCells(analysisParams, exp_info, ind)
end

for i = ind
    disp(['Currently analyzing: Ferret ' char(exp_info.animal{i}) ', Experiment ' char(exp_info.exp_id{i})])
    if exp_info.vol{i} == 1
        analysisParams.level = 1;
    else
        analysisParams.level = 0;
    end
    analysisParams.animal = char(exp_info.animal{i});
    analysisParams.expID = char(exp_info.exp_id{i});
    analysisParams.sp2ID = char(exp_info.sp2_id{i});
    analysisParams.name = char(exp_info.name{i});
    try
        analysisParams.region = char(exp_info.region{i});
    end
    try
        analysisParams.special = logical(exp_info.special{i});
    catch
        analysisParams.special = 0;
    end
    try
        analysisParams.cvsFile = char(exp_info.comments{i});
    catch
        analysisParams.cvsFile = [];
    end
%     if exp_info.flag{i} == 1
%         continue
%     end
    if analysisParams.stimType == 4
        PatchType = char(exp_info.stimulus{i});
        Patches(analysisParams)
    elseif analysisParams.stimType == 6
        SpontaneousCells(analysisParams);
    else
        GratingAnalysis(analysisParams);
    end
end