%% switch board for analysis variable
analysisParams = struct;

analysisParams.stimType = 1;  % 1 = driftingGrating, 2 = OriSf, 3 = OriTf
analysisParams.server = 1;
analysisParams.select = 1;
analysisParams.downsample = 2;
analysisParams.intrinsic = 0;
if analysisParams.intrinsic
    analysisParams.downsample = 5;
end
analysisParams.pre = 0;
analysisParams.changeThreshold = 0;
analysisParams.clean = 1;

analysisParams.field = 'rawF';

cocV3 = cbrewer('seq', 'RdPu',30);
cocNaive = cocV3(13:17,:);
cocEarly = cocV3(19:24,:);
cocAdult = cocV3(26:30,:);

cocV1 = cbrewer('seq', 'PuBuGn',30);
cocNaiveV1 = cocV1(13:17,:);
cocEarlyV1 = cocV1(19:24,:);
cocAdultV1 = cocV1(26:30,:);

%% 0.) Set folders, list all experiments and reanalyze if necessary
adata_dir = 'Z:\Juliane\Data\ImageAnalysis\';
save_dir = [adata_dir filesep 'MapDevelopment' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

filePath = 'F:\Organization\Animals\';
file = 'EpiExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
exp_info = findExpInfo(xls_txt, xls_all);

allExpIncluded = find(cell2mat(exp_info.included));
allExpInd = find(exp_info.run);

NaiveAge =find(cell2mat(exp_info.EO) == 0); %EO+0
EarlyAge = find(cell2mat(exp_info.EO) == 2 | cell2mat(exp_info.EO)== 3); %EO2-3
AdultAge = find(cell2mat(exp_info.EO) > 8); %EO>8

EarlyInd = intersect(allExpIncluded, EarlyAge);
NaiveInd = intersect(allExpIncluded, NaiveAge);
AdultInd = intersect(allExpIncluded, AdultAge);


%reanalyze if necessary
for i = allExpInd
    disp(['Currently analyzing: Ferret ' char(exp_info.animal{i}) ', Experiment ' char(exp_info.exp_id{i})])
    analysisParams.animal = char(exp_info.animal{i});
    analysisParams.expID = char(exp_info.exp_series{i});
    analysisParams.sp2ID = char(exp_info.sp2_id{i});
    driftingGrating_Fct(analysisParams);
end

%% 1.) Load all data into three master files
for ilNaive = 1:length(NaiveInd)
    datapath = [adata_dir char(exp_info.animal{NaiveInd(ilNaive)}) filesep char(exp_info.exp_series{NaiveInd(ilNaive)}) filesep];
    masterNaive{ilNaive} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    masterNaive{ilNaive}.metadata.ferret = char(exp_info.animal{NaiveInd(ilNaive)});
    masterNaive{ilNaive}.metadata.expID = exp_info.exp_series{NaiveInd(ilNaive)};    
end
for ilEarly = 1:length(EarlyInd)
    datapath = [adata_dir char(exp_info.animal{EarlyInd(ilEarly)}) filesep char(exp_info.exp_series{EarlyInd(ilEarly)}) filesep];
    masterEarly{ilEarly} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    masterEarly{ilEarly}.metadata.ferret = char(exp_info.animal{EarlyInd(ilEarly)});
    masterEarly{ilEarly}.metadata.expID = exp_info.exp_series{EarlyInd(ilEarly)};
end
for ilAdult= 1:length(AdultInd)
    datapath = [adata_dir char(exp_info.animal{AdultInd(ilAdult)}) filesep char(exp_info.exp_series{AdultInd(ilAdult)}) filesep];
    masterAdult{ilAdult} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    masterAdult{ilAdult}.metadata.ferret = char(exp_info.animal{AdultInd(ilAdult)});
    masterAdult{ilAdult}.metadata.expID = exp_info.exp_series{AdultInd(ilAdult)};
end

%% 2). Data consolidation
overlapA19Naive = [];
overlapV1Naive = [];
for ferret =1:length(NaiveInd)
    overlapA19Naive = [overlapA19Naive mean(masterNaive{ferret}.analysis.overlapA19(:))];
    overlapV1Naive = [overlapV1Naive mean(masterNaive{ferret}.analysis.overlapV1(:))];
end

overlapA19Early = [];
overlapV1Early= [];
for ferret =1:length(EarlyInd)
    overlapA19Early = [overlapA19Early mean(masterNaive{ferret}.analysis.overlapA19(:))];
    overlapV1Early = [overlapV1Early mean(masterNaive{ferret}.analysis.overlapV1(:))];
end

overlapA19Adult = [];
overlapV1Adult= [];
for ferret =1:length(AdultInd)
    overlapA19Adult = [overlapA19Adult mean(masterNaive{ferret}.analysis.overlapA19(:))];
    overlapV1Adult= [overlapV1Adult mean(masterNaive{ferret}.analysis.overlapV1(:))];
end

%% 3.) Plot results
figure
subplot(1,2,1)
boxplot([overlapA19Naive(:), overlapA19Early(:), overlapA19Adult(:)], 'Labels',{'Naive','Early' 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),cocNaive(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarly(4,:),'FaceAlpha',.5)
patch(get(h(3),'XData'),get(h(3),'YData'),cocAdult(4,:),'FaceAlpha',.5)
box off
ylabel('percentage overlap')
title('A19')


subplot(1,2,2)
boxplot([overlapV1Naive(:), overlapV1Early(:), overlapV1Adult(:)], 'Labels',{'Naive','Early' 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),cocNaiveV1(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarlyV1(4,:),'FaceAlpha',.5)
patch(get(h(3),'XData'),get(h(3),'YData'),cocAdultV1(4,:),'FaceAlpha',.5)
box off
ylabel('percentage overlap')
title('V1')
set(gcf, 'color', 'w');