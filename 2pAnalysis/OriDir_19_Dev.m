%% switch board for analysis variable
analysisParams = struct;
%which type of stimulus should it run
analysisParams.dataType = 1; %data type: 1 = cells, 2 = axons, 3 = spines
analysisParams.stimType = 1;

%what should it do?
analysisParams.reloadData = 0; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse =0; %should you reanalyse the data or just plot?
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.plotROIs = 0;   %should you plot traces for all resp ROIs?
analysisParams.plotRespROIsOnly = 0; %should you also plot traces for all non-resp ROIs?
analysisParams.server = 0; %load from the server (1) or the raid (0)
analysisParams.makeROIs = 1;

%analysisParameters
analysisParams.zThresh = 4;
analysisParams.fraction = 0.5;
analysisParams.predictor = 0;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff';
analysisParams.windowStart = 0;
analysisParams.windowStop = 2;
analysisParams.pre = 1;

cocV3 = cbrewer('seq', 'RdPu',30);
cocNaive = cocV3(13:17,:);
cocEarly = cocV3(19:24,:);
cocAdult = cocV3(26:30,:);

cocV1 = cbrewer('seq', 'PuBuGn',30);
cocNaiveV1 = cocV1(13:17,:);
cocEarlyV1 = cocV1(19:24,:);
cocAdultV1 = cocV1(26:30,:);

%% 0.) Set folders, list all experiments and reanalyze if necessary

adata_dir = 'F:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep '19_Dev_Ori_Dir' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

filePath = 'F:\Organization\Animals\';
file = '2pExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
exp_info = findExpInfo(xls_txt, xls_all);

allExpInd = find(exp_info.run); %all experiments that need to be reanalyzed
allV3Ind = find(cell2mat(exp_info.region) == 3);%all V3 exp
allV1Ind = find(cell2mat(exp_info.region) == 1);%all V1 exp


NaiveAge =find(cell2mat(exp_info.EO) == 0); %EO+0
EarlyAge = find(cell2mat(exp_info.EO) == 2 | cell2mat(exp_info.EO)== 3); %EO2-3
AdultAge = find(cell2mat(exp_info.EO) > 8); %EO>8

EarlyInd = intersect(allV3Ind, EarlyAge);
NaiveInd = intersect(allV3Ind, NaiveAge);
AdultInd = intersect(allV3Ind, AdultAge);

EarlyIndV1 = intersect(allV1Ind, EarlyAge);
NaiveIndV1 = intersect(allV1Ind, NaiveAge);
AdultIndV1 = intersect(allV1Ind, AdultAge);

%reanalyze if necessary
for i = allExpInd
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
    GratingAnalysis(analysisParams);
end

%% 1.) Load all data into three master files and do data consolidation
for ilNaive = 1:length(NaiveInd)
    datapath = [adata_dir char(exp_info.animal{NaiveInd(ilNaive)}) filesep char(exp_info.exp_id{NaiveInd(ilNaive)}) filesep];
    try
        masterNaive{ilNaive} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterNaive{ilNaive} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'metadata','analysis');
    end
    masterNaive{ilNaive}.metadata.ferret = char(exp_info.animal{NaiveInd(ilNaive)});
    masterNaive{ilNaive}.metadata.expID = exp_info.exp_id{NaiveInd(ilNaive)};
    respNaive{ilNaive} = find([masterNaive{ilNaive}.analysis.dff.roi.isResponseSignificant] == 1);
    
end
for ilEarly = 1:length(EarlyInd)
    datapath = [adata_dir char(exp_info.animal{EarlyInd(ilEarly)}) filesep char(exp_info.exp_id{EarlyInd(ilEarly)}) filesep];
    try
        masterEarly{ilEarly} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterEarly{ilEarly} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterEarly{ilEarly}.metadata.ferret = char(exp_info.animal{EarlyInd(ilEarly)});
    masterEarly{ilEarly}.metadata.expID = exp_info.exp_id{EarlyInd(ilEarly)};
    respEarly{ilEarly} = find([masterEarly{ilEarly}.analysis.dff.roi.isResponseSignificant] == 1);
end
for ilAdult= 1:length(AdultInd)
    datapath = [adata_dir char(exp_info.animal{AdultInd(ilAdult)}) filesep char(exp_info.exp_id{AdultInd(ilAdult)}) filesep];
    try
        masterAdult{ilAdult} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterAdult{ilAdult} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterAdult{ilAdult}.metadata.ferret = char(exp_info.animal{AdultInd(ilAdult)});
    masterAdult{ilAdult}.metadata.expID = exp_info.exp_id{AdultInd(ilAdult)};
    respAdult{ilAdult} = find([masterAdult{ilAdult}.analysis.dff.roi.isResponseSignificant] == 1);
end

for ilNaiveV1 = 1:length(NaiveIndV1)
    datapath = [adata_dir char(exp_info.animal{NaiveIndV1(ilNaiveV1)}) filesep char(exp_info.exp_id{NaiveIndV1(ilNaiveV1)}) filesep];
    try
        masterNaiveV1{ilNaiveV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterNaiveV1{ilNaiveV1} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'metadata','analysis');
    end
    masterNaiveV1{ilNaiveV1}.metadata.ferret = char(exp_info.animal{NaiveIndV1(ilNaiveV1)});
    masterNaiveV1{ilNaiveV1}.metadata.expID = exp_info.exp_id{NaiveIndV1(ilNaiveV1)};
    respNaiveV1{ilNaiveV1} = find([masterNaiveV1{ilNaiveV1}.analysis.dff.roi.isResponseSignificant] == 1);
    
end
for ilEarlyV1 = 1:length(EarlyIndV1)
    datapath = [adata_dir char(exp_info.animal{EarlyIndV1(ilEarlyV1)}) filesep char(exp_info.exp_id{EarlyIndV1(ilEarlyV1)}) filesep];
    try
        masterEarlyV1{ilEarlyV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterEarlyV1{ilEarlyV1} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterEarlyV1{ilEarlyV1}.metadata.ferret = char(exp_info.animal{EarlyIndV1(ilEarlyV1)});
    masterEarlyV1{ilEarlyV1}.metadata.expID = exp_info.exp_id{EarlyIndV1(ilEarlyV1)};
    respEarlyV1{ilEarlyV1} = find([masterEarlyV1{ilEarlyV1}.analysis.dff.roi.isResponseSignificant] == 1);
end
for ilAdultV1= 1:length(AdultIndV1)
    datapath = [adata_dir char(exp_info.animal{AdultIndV1(ilAdultV1)}) filesep char(exp_info.exp_id{AdultIndV1(ilAdultV1)}) filesep];
    try
        masterAdultV1{ilAdultV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterAdultV1{ilAdultV1} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterAdultV1{ilAdultV1}.metadata.ferret = char(exp_info.animal{AdultIndV1(ilAdultV1)});
    masterAdultV1{ilAdultV1}.metadata.expID = exp_info.exp_id{AdultIndV1(ilAdultV1)};
    respAdultV1{ilAdultV1} = find([masterAdultV1{ilAdultV1}.analysis.dff.roi.isResponseSignificant] == 1);
end
%%
%data consolidation
OSIFitNaive = [];
OriCircVarNaive = [];
cohensDNaive = [];
DSINaive = [];
DirCircVarNaive = [];
isRespNaive = [];
for ferret =1:length(NaiveInd)
    OSIFitNaive = [OSIFitNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).OSIFit];
    OriCircVarNaive = [OriCircVarNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).OriCircVar];
    cohensDNaive = [cohensDNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).cohensD];
    DSINaive = [DSINaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).DSI];
    DirCircVarNaive = [DirCircVarNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).DirCircVar];
    isRespNaive = [isRespNaive masterNaive{ferret}.analysis.dff.roi.isResponseSignificant];
end

OSIFitEarly = [];
OriCircVarEarly = [];
cohensDEarly = [];
DSIEarly= [];
DirCircVarEarly = [];
isRespEarly = [];
for ferret =1:length(EarlyInd)
    OSIFitEarly = [OSIFitEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).OSIFit];
    OriCircVarEarly = [OriCircVarEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).OriCircVar];
    cohensDEarly = [cohensDEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).cohensD];
    DSIEarly = [DSIEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).DSI];
    DirCircVarEarly = [DirCircVarEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).DirCircVar];
    isRespEarly = [isRespEarly masterEarly{ferret}.analysis.dff.roi.isResponseSignificant];
end

OSIFitAdult = [];
OriCircVarAdult  = [];
cohensDAdult  = [];
DSIAdult = [];
DirCircVarAdult  = [];
isRespAdult  = [];
for ferret =1:length(AdultInd)
    OSIFitAdult  = [OSIFitAdult  masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).OSIFit];
    OriCircVarAdult  = [OriCircVarAdult  masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).OriCircVar];
    cohensDAdult  = [cohensDAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).cohensD];
    DSIAdult  = [DSIAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).DSI];
    DirCircVarAdult  = [DirCircVarAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).DirCircVar];
    isRespAdult = [isRespAdult masterAdult{ferret}.analysis.dff.roi.isResponseSignificant];
end

OSIFitNaiveV1 = [];
OriCircVarNaiveV1 = [];
cohensDNaiveV1 = [];
DSINaiveV1 = [];
DirCircVarNaiveV1 = [];
isRespNaiveV1 = [];
for ferret =1:length(NaiveIndV1)
    OSIFitNaiveV1 = [OSIFitNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).OSIFit];
    OriCircVarNaiveV1 = [OriCircVarNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).OriCircVar];
    cohensDNaiveV1 = [cohensDNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).cohensD];
    DSINaiveV1 = [DSINaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).DSI];
    DirCircVarNaiveV1 = [DirCircVarNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).DirCircVar];
    isRespNaiveV1 = [isRespNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi.isResponseSignificant];
end

OSIFitEarlyV1 = [];
OriCircVarEarlyV1 = [];
cohensDEarlyV1 = [];
DSIEarlyV1= [];
DirCircVarEarlyV1 = [];
isRespEarlyV1 = [];
for ferret =1:length(EarlyIndV1)
    OSIFitEarlyV1 = [OSIFitEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).OSIFit];
    OriCircVarEarlyV1 = [OriCircVarEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).OriCircVar];
    cohensDEarlyV1 = [cohensDEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).cohensD];
    DSIEarlyV1 = [DSIEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).DSI];
    DirCircVarEarlyV1 = [DirCircVarEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).DirCircVar];
    isRespEarlyV1 = [isRespEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi.isResponseSignificant];
end

OSIFitAdultV1 = [];
OriCircVarAdultV1  = [];
cohensDAdultV1  = [];
DSIAdultV1 = [];
DirCircVarAdultV1  = [];
isRespAdultV1  = [];
for ferret =1:length(AdultIndV1)
    OSIFitAdultV1  = [OSIFitAdultV1  masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).OSIFit];
    OriCircVarAdultV1  = [OriCircVarAdultV1  masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).OriCircVar];
    cohensDAdultV1  = [cohensDAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).cohensD];
    DSIAdultV1  = [DSIAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).DSI];
    DirCircVarAdultV1  = [DirCircVarAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).DirCircVar];
    isRespAdultV1 = [isRespAdultV1 masterAdultV1{ferret}.analysis.dff.roi.isResponseSignificant];
end

%% 2.) Calculate population indices
%% 3.) Plot results
% a)responsiveness and selectivity portions
figure
subplot(3, 3, 1)
allNaive = length(isRespNaive);
nonRespNaive = length(find([isRespNaive] == 0)) ./allNaive;
percRespNaive = length(find([isRespNaive] == 1)) ./allNaive;
h = pie([nonRespNaive percRespNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaive(1,:));
set(hp(2), 'FaceColor', cocNaive(4,:));
title('Responsive Naive')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 4)
allEarly = length(isRespEarly);
nonRespEarly = length(find([isRespEarly] == 0)) ./allEarly;
percRespEarly = length(find([isRespEarly] == 1)) ./allEarly;
h = pie([nonRespEarly percRespEarly]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocEarly(1,:));
set(hp(2), 'FaceColor', cocEarly(4,:));
title('Responsive Early')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 7)
allAdult = length(isRespAdult);
nonRespAdult = length(find([isRespAdult] == 0)) ./allAdult;
percRespAdult = length(find([isRespAdult] == 1)) ./allAdult;
h = pie([nonRespAdult percRespAdult]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAdult(1,:));
set(hp(2), 'FaceColor', cocAdult(4,:));
title('Responsive Adult')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 2)
oriNaive = length(find([OSIFitNaive] > 0.2)) ./ length(OSIFitNaive);
% ori_V1 = length(find([OriCircVar_V1] > 0.2)) ./ length(OriCircVar_V1);
nonOriNaive = 1- oriNaive;
h = pie([nonOriNaive oriNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaive(1,:));
set(hp(2), 'FaceColor', cocNaive(4,:));
title('Naive')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 5)
oriEarly = length(find([OSIFitEarly] > 0.2)) ./length(OSIFitEarly);
%ori_V3 = length(find([OriCircVar_V3] > 0.2)) ./length(OriCircVar_V3);
nonOriEarly = 1- oriEarly;
h = pie([nonOriEarly oriEarly]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocEarly(1,:));
set(hp(2), 'FaceColor', cocEarly(4,:));
title('Early')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 8)
oriAdult = length(find([OSIFitAdult] > 0.2)) ./length(OSIFitAdult);
%ori_V3 = length(find([OriCircVar_V3] > 0.2)) ./length(OriCircVar_V3);
nonOriAdult = 1- oriAdult;
h = pie([nonOriAdult oriAdult]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAdult(1,:));
set(hp(2), 'FaceColor', cocAdult(4,:));
title('Adult')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 3)
dirNaive = length(find([DSINaive] > 0.2)) ./length(DSINaive);
nonDirNaive = 1- dirNaive;
h = pie([nonDirNaive dirNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaive(1,:));
set(hp(2), 'FaceColor', cocNaive(4,:));
title('Direction-selective Naive')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 6)
dirEarly = length(find([DSIEarly] > 0.2)) ./length(DSIEarly);
nonDirEarly = 1- dirEarly;
h = pie([nonDirEarly dirEarly]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocEarly(1,:));
set(hp(2), 'FaceColor', cocEarly(4,:));
title('Early')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 9)
dirAdult= length(find([DSIAdult] > 0.2)) ./length(DSIAdult);
nonDirAdult= 1- dirAdult;
h = pie([nonDirAdult dirAdult]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAdult(1,:));
set(hp(2), 'FaceColor', cocAdult(4,:));
title('Adult')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'resp_cells.png'))

%for V1
figure
subplot(3, 3, 1)
allNaiveV1 = length(isRespNaiveV1);
nonRespNaiveV1 = length(find([isRespNaiveV1] == 0)) ./allNaiveV1;
percRespNaiveV1 = length(find([isRespNaiveV1] == 1)) ./allNaiveV1;
h = pie([nonRespNaiveV1 percRespNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaiveV1(1,:));
set(hp(2), 'FaceColor', cocNaiveV1(4,:));
title('Responsive Naive V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 4)
allEarlyV1 = length(isRespEarlyV1);
nonRespEarlyV1 = length(find([isRespEarlyV1] == 0)) ./allEarlyV1;
percRespEarlyV1 = length(find([isRespEarlyV1] == 1)) ./allEarlyV1;
h = pie([nonRespEarlyV1 percRespEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocEarlyV1(1,:));
set(hp(2), 'FaceColor', cocEarlyV1(4,:));
title('Responsive Early V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 7)
allAdultV1 = length(isRespAdultV1);
nonRespAdultV1 = length(find([isRespAdultV1] == 0)) ./allAdultV1;
percRespAdultV1 = length(find([isRespAdultV1] == 1)) ./allAdultV1;
h = pie([nonRespAdultV1 percRespAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAdultV1(1,:));
set(hp(2), 'FaceColor', cocAdultV1(4,:));
title('Responsive Adult V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 2)
oriNaiveV1 = length(find([OSIFitNaiveV1] > 0.2)) ./ length(OSIFitNaiveV1);
nonOriNaiveV1 = 1- oriNaiveV1;
h = pie([nonOriNaiveV1 oriNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaiveV1(1,:));
set(hp(2), 'FaceColor', cocNaiveV1(4,:));
title('Naive V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 5)
oriEarlyV1 = length(find([OSIFitEarlyV1] > 0.2)) ./length(OSIFitEarlyV1);
nonOriEarlyV1 = 1- oriEarlyV1;
h = pie([nonOriEarlyV1 oriEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocEarlyV1(1,:));
set(hp(2), 'FaceColor', cocEarlyV1(4,:));
title('Early V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 8)
oriAdultV1 = length(find([OSIFitAdultV1] > 0.2)) ./length(OSIFitAdultV1);
nonOriAdultV1 = 1- oriAdultV1;
h = pie([nonOriAdultV1 oriAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAdultV1(1,:));
set(hp(2), 'FaceColor', cocAdultV1(4,:));
title('Adult V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 3)
dirNaiveV1 = length(find([DSINaiveV1] > 0.2)) ./length(DSINaiveV1);
nonDirNaiveV1 = 1- dirNaiveV1;
h = pie([nonDirNaiveV1 dirNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaiveV1(1,:));
set(hp(2), 'FaceColor', cocNaiveV1(4,:));
title('Direction-selective Naive V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 6)
dirEarlyV1 = length(find([DSIEarlyV1] > 0.2)) ./length(DSIEarlyV1);
nonDirEarlyV1 = 1- dirEarlyV1;
h = pie([nonDirEarlyV1 dirEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocEarlyV1(1,:));
set(hp(2), 'FaceColor', cocEarlyV1(4,:));
title('Early V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 9)
dirAdultV1= length(find([DSIAdultV1] > 0.2)) ./length(DSIAdultV1);
nonDirAdultV1= 1- dirAdultV1;
h = pie([nonDirAdultV1 dirAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAdultV1(1,:));
set(hp(2), 'FaceColor', cocAdultV1(4,:));
title('Adult V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'resp_cellsV1.png'))

%% b) OSIFit and DSI distribution
figure
subplot(2,3,1)
plot(nanmedian(OSIFitNaive),1.05 * max(histcounts(OSIFitNaive)),'v','MarkerSize', 8','MarkerEdgeColor',cocNaive(5,:),'MarkerFaceColor',cocNaive(4,:)); hold on
text(nanmedian(OSIFitNaive),1.2 * max(histcounts(OSIFitNaive)),num2str(round(100*nanmedian(OSIFitNaive))/100),'HorizontalAlignment','center', 'Color', cocNaive(5,:), 'FontSize', 12')
histogram(OSIFitNaive, 20, 'FaceColor', cocNaive(4,:), 'EdgeColor', cocNaive(5,:));
ylabel('Cells');
xlabel('OSI');
title('Naive')
set(gca,'Box','off');

subplot(2,3,2)
plot(nanmedian(OSIFitEarly),1.05 * max(histcounts(OSIFitEarly)),'v','MarkerSize', 8','MarkerEdgeColor',cocEarly(5,:),'MarkerFaceColor',cocEarly(4,:)); hold on
text(nanmedian(OSIFitEarly),1.2 * max(histcounts(OSIFitEarly)),num2str(round(100*nanmedian(OSIFitEarly))/100),'HorizontalAlignment','center', 'Color', cocEarly(5,:), 'FontSize', 12')
histogram(OSIFitEarly, 20, 'FaceColor', cocEarly(4,:), 'EdgeColor', cocEarly(5,:));
ylabel('Cells');
xlabel('OSI');
title('Early')
set(gca,'Box','off');

subplot(2,3,3)
plot(nanmedian(OSIFitAdult),1.05 * max(histcounts(OSIFitAdult)),'v','MarkerSize', 8','MarkerEdgeColor',cocAdult(5,:),'MarkerFaceColor',cocAdult(4,:)); hold on
text(nanmedian(OSIFitAdult),1.2 * max(histcounts(OSIFitAdult)),num2str(round(100*nanmedian(OSIFitAdult))/100),'HorizontalAlignment','center', 'Color', cocAdult(5,:), 'FontSize', 12')
histogram(OSIFitAdult, 20, 'FaceColor', cocAdult(4,:), 'EdgeColor', cocAdult(5,:));
ylabel('Cells');
xlabel('OSI');
title('Adult')
set(gca,'Box','off');

subplot(2,3,4)
plot(nanmedian(DSINaive),1.05 * max(histcounts(DSINaive)),'v','MarkerSize', 8','MarkerEdgeColor',cocNaive(5,:),'MarkerFaceColor',cocNaive(4,:)); hold on
text(nanmedian(DSINaive),1.2 * max(histcounts(DSINaive)),num2str(round(100*nanmedian(DSINaive))/100),'HorizontalAlignment','center', 'Color', cocNaive(5,:), 'FontSize', 12')
histogram(DSINaive, 20, 'FaceColor', cocNaive(4,:), 'EdgeColor', cocNaive(5,:));
ylabel('Cells');
xlabel('DSI');
title('Naive')
set(gca,'Box','off');


subplot(2,3,5)
plot(nanmedian(DSIEarly),1.05 * max(histcounts(DSIEarly)),'v','MarkerSize', 8','MarkerEdgeColor',cocEarly(5,:),'MarkerFaceColor',cocEarly(4,:)); hold on
text(nanmedian(DSIEarly),1.2 * max(histcounts(DSIEarly)),num2str(round(100*nanmedian(DSIEarly))/100),'HorizontalAlignment','center', 'Color', cocEarly(5,:), 'FontSize', 12')
histogram(DSIEarly, 20, 'FaceColor', cocEarly(4,:), 'EdgeColor', cocEarly(5,:));
ylabel('Cells');
xlabel('DSI');
xlim([0 1])
title('Early')
set(gca,'Box','off');

subplot(2,3,6)
plot(nanmedian(DSIAdult),1.05 * max(histcounts(DSIAdult)),'v','MarkerSize', 8','MarkerEdgeColor',cocAdult(5,:),'MarkerFaceColor',cocAdult(4,:)); hold on
text(nanmedian(DSIAdult),1.2 * max(histcounts(DSIAdult)),num2str(round(100*nanmedian(DSIAdult))/100),'HorizontalAlignment','center', 'Color', cocAdult(5,:), 'FontSize', 12')
histogram(DSIAdult, 20, 'FaceColor', cocAdult(4,:), 'EdgeColor', cocAdult(5,:));
ylabel('Cells');
xlabel('DSI');
xlim([0 1])
title('Adult')
set(gca,'Box','off');

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_DSI_distribution.png'))

%for V1
figure
subplot(2,3,1)
plot(nanmedian(OSIFitNaiveV1),1.05 * max(histcounts(OSIFitNaiveV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocNaiveV1(5,:),'MarkerFaceColor',cocNaiveV1(4,:)); hold on
text(nanmedian(OSIFitNaiveV1),1.2 * max(histcounts(OSIFitNaiveV1)),num2str(round(100*nanmedian(OSIFitNaiveV1))/100),'HorizontalAlignment','center', 'Color', cocNaiveV1(4,:), 'FontSize', 12')
histogram(OSIFitNaiveV1, 20, 'FaceColor', cocNaiveV1(4,:), 'EdgeColor', cocNaiveV1(5,:));
ylabel('Cells');
xlabel('OSI');
title('Naive')
set(gca,'Box','off');

subplot(2,3,2)
plot(nanmedian(OSIFitEarlyV1),1.05 * max(histcounts(OSIFitEarlyV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocEarlyV1(5,:),'MarkerFaceColor',cocEarlyV1(4,:)); hold on
text(nanmedian(OSIFitEarlyV1),1.2 * max(histcounts(OSIFitEarlyV1)),num2str(round(100*nanmedian(OSIFitEarlyV1))/100),'HorizontalAlignment','center', 'Color', cocEarlyV1(5,:), 'FontSize', 12')
histogram(OSIFitEarlyV1, 20, 'FaceColor', cocEarlyV1(4,:), 'EdgeColor', cocEarlyV1(5,:));
ylabel('Cells');
xlabel('OSI');
title('Early')
set(gca,'Box','off');

subplot(2,3,3)
plot(nanmedian(OSIFitAdultV1),1.05 * max(histcounts(OSIFitAdultV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocAdultV1(5,:),'MarkerFaceColor',cocAdultV1(4,:)); hold on
text(nanmedian(OSIFitAdultV1),1.2 * max(histcounts(OSIFitAdultV1)),num2str(round(100*nanmedian(OSIFitAdultV1))/100),'HorizontalAlignment','center', 'Color', cocAdultV1(5,:), 'FontSize', 12')
histogram(OSIFitAdultV1, 20, 'FaceColor', cocAdultV1(4,:), 'EdgeColor', cocAdultV1(5,:));
ylabel('Cells');
xlabel('OSI');
title('Adult')
set(gca,'Box','off');

subplot(2,3,4)
plot(nanmedian(DSINaiveV1),1.05 * max(histcounts(DSINaiveV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocNaiveV1(5,:),'MarkerFaceColor',cocNaiveV1(4,:)); hold on
text(nanmedian(DSINaiveV1),1.2 * max(histcounts(DSINaiveV1)),num2str(round(100*nanmedian(DSINaiveV1))/100),'HorizontalAlignment','center', 'Color', cocNaiveV1(5,:), 'FontSize', 12')
histogram(DSINaiveV1, 20, 'FaceColor', cocNaiveV1(4,:), 'EdgeColor', cocNaiveV1(5,:));
ylabel('Cells');
xlabel('DSI');
title('Naive')
set(gca,'Box','off');


subplot(2,3,5)
plot(nanmedian(DSIEarlyV1),1.05 * max(histcounts(DSIEarlyV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocEarlyV1(5,:),'MarkerFaceColor',cocEarlyV1(4,:)); hold on
text(nanmedian(DSIEarlyV1),1.2 * max(histcounts(DSIEarlyV1)),num2str(round(100*nanmedian(DSIEarlyV1))/100),'HorizontalAlignment','center', 'Color', cocEarlyV1(5,:), 'FontSize', 12')
histogram(DSIEarlyV1, 20, 'FaceColor', cocEarlyV1(4,:), 'EdgeColor', cocEarlyV1(5,:));
ylabel('Cells');
xlabel('DSI');
xlim([0 1])
title('Early')
set(gca,'Box','off');

subplot(2,3,6)
plot(nanmedian(DSIAdultV1),1.05 * max(histcounts(DSIAdultV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocAdultV1(5,:),'MarkerFaceColor',cocAdultV1(4,:)); hold on
text(nanmedian(DSIAdultV1),1.2 * max(histcounts(DSIAdultV1)),num2str(round(100*nanmedian(DSIAdultV1))/100),'HorizontalAlignment','center', 'Color', cocAdultV1(5,:), 'FontSize', 12')
histogram(DSIAdultV1, 20, 'FaceColor', cocAdultV1(4,:), 'EdgeColor', cocAdultV1(5,:));
ylabel('Cells');
xlabel('DSI');
xlim([0 1])
title('Adult')
set(gca,'Box','off');

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_DSI_distributionV1.png'))

%% c) circVar & dirCircVar
figure
subplot(1,6,1)
distributionPlot(OSIFitNaive','color', cocNaive(4,:)); hold on
boxplot(OSIFitNaive,'Label', {'Naive'})
% distributionPlot(OriCircVar_V1','color', coc_prop(1,:)); hold on
% boxplot(OriCircVar_V1,'Label', {'V1'})
ylim([0 1])
xlabel('OSI')
set(gca,'Box','off');

subplot(1,6,2)
% distributionPlot(OriCircVar_V3','color', coc_prop(2,:)); hold on
% boxplot(OriCircVar_V3, 'Label', {'V3'})
distributionPlot(OSIFitEarly','color', cocEarly(4,:)); hold on
boxplot(OSIFitEarly, 'Label', {'Early'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
% distributionPlot(OriCircVar_V3','color', coc_prop(2,:)); hold on
% boxplot(OriCircVar_V3, 'Label', {'V3'})
distributionPlot(OSIFitAdult','color', cocAdult(4,:)); hold on
boxplot(OSIFitAdult, 'Label', {'Adult'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
% distributionPlot(DirCircVar_V1','color', coc_prop(3,:)); hold on
% boxplot(DirCircVar_V1, 'Label', {'V1'})
distributionPlot(DSINaive','color', cocNaive(4,:)); hold on
boxplot(DSINaive, 'Label', {'Naive'})
ylim([0 1])
xlabel('DSI')
set(gca,'Box','off');

subplot(1, 6, 5)
% distributionPlot(DirCircVar_V3','color', coc_prop(4,:)); hold on
% boxplot(DirCircVar_V3, 'Label', {'V3'})
distributionPlot(DSIEarly','color', cocEarly(4,:)); hold on
boxplot(DSIEarly, 'Label', {'Early'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1, 6, 6)
% distributionPlot(DirCircVar_V3','color', coc_prop(4,:)); hold on
% boxplot(DirCircVar_V3, 'Label', {'V3'})
distributionPlot(DSIAdult','color', cocAdult(4,:)); hold on
boxplot(DSIAdult, 'Label', {'Adult'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSIDSIViolin.png'))

figure
subplot(1,6,1)
distributionPlot(OSIFitNaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(OSIFitNaiveV1,'Label', {'Naive V1'})
ylim([0 1])
xlabel('OSI')
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(OSIFitEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(OSIFitEarlyV1, 'Label', {'Early V1'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(OSIFitAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(OSIFitAdultV1, 'Label', {'Adult V1'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(DSINaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(DSINaiveV1, 'Label', {'Naive V1'})
ylim([0 1])
xlabel('DSI')
set(gca,'Box','off');

subplot(1, 6, 5)
distributionPlot(DSIEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(DSIEarlyV1, 'Label', {'Early V1'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1, 6, 6)
distributionPlot(DSIAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(DSIAdultV1, 'Label', {'Adult V1'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSIDSIViolinV1.png'))

%[pOriCircVar, ~, stats_pOriCircVar] = ranksum(OriCircVarNaive, OriCircVarEarly);
%[pDirCircVar, ~, stats_pDirCircVar] = ranksum(DirCircVarNaive, DirCircVarEarly);
%[pOSIFit, ~, stats_pOSIFit] = ranksum(OSIFitNaive, OSIFitEarly);
%[pDSI, ~, stats_pDSI] = ranksum(DSINaive, DSIEarly);

%% d) cohensD
figure
hdl(1) = cdfplot(cohensDNaive); hold all
hdl(2) = cdfplot(cohensDEarly); hold all
hdl(3) = cdfplot(cohensDAdult); hold all
set(hdl(1), 'Color', cocNaive(4,:), 'LineWidth', 4);
set(hdl(2), 'Color', cocEarly(4,:), 'LineWidth', 4);
set(hdl(3), 'Color', cocAdult(4,:), 'LineWidth', 4);
grid off;
xlim([0 5]);
legend('Naive', 'Early', 'Adult', 'Location', 'SouthEast'); legend('boxoff')
ylabel('Cumulative probability', 'FontSize', 12)
xlabel('Cohens D', 'FontSize', 12);

h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', 'w');
set(gcf, 'color', 'w');
set(gca, 'box', 'off')

saveas(gcf, fullfile(save_dir, 'cohensD.png'))

figure
hdl(1) = cdfplot(cohensDNaiveV1); hold all
hdl(2) = cdfplot(cohensDEarlyV1); hold all
hdl(3) = cdfplot(cohensDAdultV1); hold all
set(hdl(1), 'Color', cocNaiveV1(4,:), 'LineWidth', 4);
set(hdl(2), 'Color', cocEarlyV1(4,:), 'LineWidth', 4);
set(hdl(3), 'Color', cocAdultV1(4,:), 'LineWidth', 4);
grid off;
xlim([0 5]);
legend('Naive V1', 'Early V1', 'Adult V1', 'Location', 'SouthEast'); legend('boxoff')
ylabel('Cumulative probability', 'FontSize', 12)
xlabel('Cohens D', 'FontSize', 12);

h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', 'w');
set(gcf, 'color', 'w');
set(gca, 'box', 'off')

saveas(gcf, fullfile(save_dir, 'cohensDV1.png'))