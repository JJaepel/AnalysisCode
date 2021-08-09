%% switch board for analysis variables
analysisParams = struct;
analysisParams.plotROIs = 0;   %should you plot traces for all ROIs?
analysisParams.plotRespROIsOnly = 1;
analysisParams.reloadData = 1; %should you reload from suite2p and do baselining?
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.server = 0; %load from the server (1) or the raid (0)
analysisParams.makeROIs = 1;
analysisParams.checkDFF = 0;
analysisParams.plotRespROIsOnly = 1;

analysisParams.zThresh = 5;
analysisParams.fraction = 0.5;
analysisParams.predictor = 1;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff';
analysisParams.windowStart = 0;
analysisParams.windowStop = 2;
analysisParams.pre = 1;

%% list all experiments
if analysisParams.server
    [~, xls_txt, xls_all]=xlsread('J:\Juliane\Organization\Animals\2pExpByStimulusAxon.xlsx', 'driftingGrating_ori_sf');
else
    [~, xls_txt, xls_all]=xlsread('F:\Organization\Animals\2pExpByStimulusAxon.xlsx', 'driftingGrating_ori_sf');
end

flagCol = find(contains(xls_txt(1,:),'Flag'),1);
animalCol = find(contains(xls_txt(1,:),'animal'),1);
expCol = find(contains(xls_txt(1,:),'expNumber'),1);
spk2Col = find(contains(xls_txt(1,:),'spk2'),1);
volCol = find(contains(xls_txt(1,:),'volume'),1);
nameCol = find(contains(xls_txt(1,:),'name'),1);
runCol = find(contains(xls_txt(1,:),'Run'),1);

exp_info = struct;

k = 1;
for i = 2:size(xls_txt,1)
    flag = xls_all(i,flagCol);
    if flag{1} == 0
        exp_info.animal{k} = xls_all(i,animalCol);
        exp = xls_all(i,expCol);
        if exp{1} > 9
            if exp{1} > 99
                exp_info.exp_id{k} = ['t00' num2str(exp{1})];
            else
                exp_info.exp_id{k} = ['t000' num2str(exp{1})];
            end
        else
            exp_info.exp_id{k} = ['t0000' num2str(exp{1})];
        end
        
        sp2 = xls_all(i,spk2Col);
        if sp2{1} > 9
            if sp2{1} > 99
                exp_info.sp2_id{k} = ['t00' num2str(sp2{1})];
            else
                exp_info.sp2_id{k} = ['t000' num2str(sp2{1})];
            end
        else
            exp_info.sp2_id{k} = ['t0000' num2str(sp2{1})];
        end
        vol = xls_all(i,volCol);
        if contains(vol{1}, 'yes')
            exp_info.vol{k} = 1;
        elseif contains(vol{1}, 'no')
            exp_info.vol{k} = 0;
        end
        exp_info.name{k} = xls_all(i,nameCol);
        exp_info.run(k) = xls_all{i,runCol};
        k = k+1;
    end
end

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
    OriSfAxon(analysisParams);
end