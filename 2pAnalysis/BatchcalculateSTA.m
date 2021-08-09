%% switch board for analysis variables
analysisParams.server = 0;
analysisParams.select = 1;
analysisParams.OffSet = 0.150; % in s

%% list all experiments 
if analysisParams.server
    filePath = 'Z:\Juliane\Organization\Animals\';
else
    filePath = 'F:\Organization\Animals\';
end

file = '2pExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'Patches');

exp_info = findExpInfo(xls_txt, xls_all);

if analysisParams.select
    ind = find(exp_info.run);
else
    ind = 1:1:length(exp_info.animal); %run through all experiments at once
end

%% run through all selected experiments

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
        analysisParams.special = exp_info.special(i);
    end
    
    %try
        %calculateSTA_Patches(analysisParams)
        calcSignificantPatches(analysisParams)
        close all
%     catch
%         disp('Aborted experiment analysis')
%     end
end