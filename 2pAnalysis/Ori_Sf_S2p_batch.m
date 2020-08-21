% list all experiments 
[~, xls_txt, xls_all]=xlsread('Z:\Juliane\Organization\Animals\2pExpByStimulus.xlsx', 'driftingGrating_ori_sf');

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

% switch board for variables
plotROIs = 0;   %should you plot traces for all ROIs?
reloadData = 1; %should you reload from suite2p and do baselining?
select = 1;

if select
    ind = find(exp_info.run);
else
    ind = 1:1:length(exp_info.animal);
end

%run through all experiments at once
for i = ind
    disp(['Currently analyzing: Ferret ' char(exp_info.animal{i}) ', Experiment ' char(exp_info.exp_id{i})])
    if exp_info.vol{i} == 1
        Ori_Sf_S2p_level(char(exp_info.animal{i}), char(exp_info.exp_id{i}), char(exp_info.sp2_id{i}), char(exp_info.name{i}), reloadData, plotROIs);
    else
        Ori_Sf_S2p(char(exp_info.animal{i}), char(exp_info.exp_id{i}), char(exp_info.sp2_id{i}), char(exp_info.name{i}), reloadData, plotROIs);
    end
end