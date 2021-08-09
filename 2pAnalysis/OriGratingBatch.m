%% switch board for analysis variables
analysisParams = struct;
analysisParams.plotROIs = 0;   %should you plot traces for all ROIs?
analysisParams.reloadData = 1; %should you reload from suite2p and do baselining?
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.server = 0; %load from the server (1) or the raid (0)

analysisParams.zThresh = 4;
analysisParams.fraction = 0.5;
analysisParams.predictor = 1;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff';
analysisParams.windowStart = 0;
analysisParams.windowStop = 2;
analysisParams.pre = 1;

%% list all experiments 
if analysisParams.server
    [~, xls_txt, xls_all]=xlsread('Z:\Juliane\Organization\Animals\2pExpByStimulus.xlsx', 'driftingGrating');
else
    [~, xls_txt, xls_all]=xlsread('F:\Organization\Animals\2pExpByStimulus.xlsx', 'driftingGrating');
end

exp_info = findExpInfo(xls_txt, xls_all);


if analysisParams.select
    ind = find(exp_info.run);
else
    ind = 1:1:length(exp_info.animal);
end

%% run through all experiments at once
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
    OriGrating(analysisParams);
end
