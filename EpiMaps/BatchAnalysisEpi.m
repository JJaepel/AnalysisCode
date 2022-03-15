%close all
clear all
analysisParams = struct;

%% switch board for analysis variable

%analysisParameters
analysisParams.stimType = 5;  % 1 = driftingGrating, 2 = OriSf, 3 = OriTf, 4 = flash, 5 = spontaneous
analysisParams.server = 0;
analysisParams.select = 1;
analysisParams.downsample = 2;
analysisParams.intrinsic = 0;
analysisParams.pre = 1;
analysisParams.changeThreshold = 0;
analysisParams.clean = 0;

analysisParams.field = 'high';
analysisParams.register = 0;
analysisParams.bandpass = 0;
%% list all experiments

filePath = 'F:\Organization\Animals\';
file = 'EpiExpByStimulus.xlsx';

switch analysisParams.stimType 
    case 1
        [~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
    case 2
        [~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating_ori_sf');
    case 3
        [~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating_ori_tf');
    case 4
        [~, xls_txt, xls_all]=xlsread([filePath file], 'fullfieldFlash');
    case 5
        [~, xls_txt, xls_all]=xlsread([filePath file], 'spontaneous');
end

exp_info = findExpInfo(xls_txt, xls_all);

if analysisParams.select
    ind = find(exp_info.run);           %run through selected experiments
else
    ind = 1:1:length(exp_info.animal); %run through all experiments at once
end

for i = ind
    disp(['Currently analyzing: Ferret ' char(exp_info.animal{i}) ', Experiment ' char(exp_info.exp_id{i})])
    analysisParams.animal = char(exp_info.animal{i});
    analysisParams.expID = char(exp_info.exp_series{i});
    analysisParams.sp2ID = char(exp_info.sp2_id{i});
    switch analysisParams.stimType 
        case 1
            %try
                driftingGrating_Fct(analysisParams);
            %end
        case 5
            %try
                spontaneous_Fct(analysisParams);
            %end
    end
end