%% switch board for analysis variable
analysisParams = struct;
%which type of stimulus should it run
analysisParams.dataType = 1; %data type: 1 = cells, 2 = axons, 3 = spines
analysisParams.stimType = 2; % 1 = driftingGrating, 2 = OriSf, 3 = OriTf, 4 = Patch, 5 = bino

% what should it do?
analysisParams.reloadData = 1; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse = 1; %should you reanalyse the data or just plot?
analysisParams.reanalysePop = 1;
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.plotROIs = 0;   %should you plot traces for all resp ROIs?
analysisParams.plotRespROIsOnly = 0; %should you also plot traces for all non-resp ROIs?
analysisParams.server = 0; %load from the server (1) or the raid (0)

% analysisParameters
analysisParams.zThresh = 4;
analysisParams.fraction = 0.5;
analysisParams.predictor = 0;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff';
analysisParams.windowStart = 0;
analysisParams.windowStop = 2;
analysisParams.pre = 1;
analysisParams.manual = 0;

% color schemes
cocV3 = cbrewer('seq', 'RdPu',30);
cocNaive = cocV3(13:17,:);
cocEarly = cocV3(19:24,:);
cocAdult = cocV3(26:30,:);
grey = [0.5 0.5 0.5];

cocV1 = cbrewer('seq', 'PuBuGn',30);
cocNaiveV1 = cocV1(13:17,:);
cocEarlyV1 = cocV1(19:24,:);
cocAdultV1 = cocV1(26:30,:);

%% 0.) Set folders, list all experiments and reanalyze if necessary

adata_dir = 'F:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep '19_Dev_Ori_SF' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

filePath = 'F:\Organization\Animals\';
file = '2pExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating_ori_sf');
exp_info = findExpInfo(xls_txt, xls_all);

% sort experiments by age and region
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

% remove bad experiments
Flag = find(exp_info.flag == 1);
EarlyInd = setdiff(EarlyInd, Flag);
NaiveInd = setdiff(NaiveInd, Flag);
AdultInd = setdiff(AdultInd, Flag);
EarlyIndV1 = setdiff(EarlyIndV1, Flag);
NaiveIndV1 = setdiff(NaiveIndV1, Flag);
AdultIndV1 = setdiff(AdultIndV1, Flag);

%reanalyze everything if necessary
if ~analysisParams.select  
    allExpInd =[EarlyInd NaiveInd AdultInd EarlyIndV1 NaiveIndV1 AdultIndV1];
end

%reanalyze only selected ones
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

%% 1.) Load all data into three master files and do data consolidation for population information
shortDistanceCorrNaive= [];
longDistanceCorrNaive= [];
corrTrialsMatchedNaive= [];
corrTrialsOrthoNaive= [];

for ilNaive = 1:length(NaiveInd)
    datapath = [adata_dir char(exp_info.animal{NaiveInd(ilNaive)}) filesep char(exp_info.exp_id{NaiveInd(ilNaive)}) filesep];
    try
        masterNaive{ilNaive} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterNaive{ilNaive} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'metadata','analysis');
    end
    masterNaive{ilNaive}.metadata.ferret = char(exp_info.animal{NaiveInd(ilNaive)});
    masterNaive{ilNaive}.metadata.expID = exp_info.exp_id{NaiveInd(ilNaive)};
    temp = load(fullfile(datapath, 'PopData.mat'), 'PopAnalysis');
    masterNaive{ilNaive}.PopAnalysis = temp.PopAnalysis;
    clear temp;
    respNaive{ilNaive} = find([masterNaive{ilNaive}.analysis.dff.roi.isResponseSignificant] == 1);
%     shortDistanceCorrNaive = [shortDistanceCorrNaive  masterNaive{ilNaive}.analysis.dff.shortDistanceCorr'];
%     longDistanceCorrNaive = [longDistanceCorrNaive  masterNaive{ilNaive}.analysis.dff.longDistanceCorr'];
%     corrTrialsMatchedNaive = [corrTrialsMatchedNaive  masterNaive{ilNaive}.analysis.dff.corrTrialsMatched(:)'];
%     corrTrialsOrthoNaive = [corrTrialsOrthoNaive  masterNaive{ilNaive}.analysis.dff.corrTrialsOrtho(:)'];    
end

% shortDistanceCorrEarly= [];longDistanceCorrEarly= [];
% corrTrialsMatchedEarly= [];corrTrialsOrthoEarly= [];
for ilEarly = 1:length(EarlyInd)
    datapath = [adata_dir char(exp_info.animal{EarlyInd(ilEarly)}) filesep char(exp_info.exp_id{EarlyInd(ilEarly)}) filesep];
    try
        masterEarly{ilEarly} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterEarly{ilEarly} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterEarly{ilEarly}.metadata.ferret = char(exp_info.animal{EarlyInd(ilEarly)});
    masterEarly{ilEarly}.metadata.expID = exp_info.exp_id{EarlyInd(ilEarly)};
    temp = load(fullfile(datapath, 'PopData.mat'), 'PopAnalysis');
    masterEarly{ilEarly}.PopAnalysis = temp.PopAnalysis;
    clear temp;
    respEarly{ilEarly} = find([masterEarly{ilEarly}.analysis.dff.roi.isResponseSignificant] == 1);
%     shortDistanceCorrEarly = [shortDistanceCorrEarly  masterEarly{ilEarly}.analysis.dff.shortDistanceCorr'];
%     longDistanceCorrEarly = [longDistanceCorrEarly  masterEarly{ilEarly}.analysis.dff.longDistanceCorr'];
%     corrTrialsMatchedEarly = [corrTrialsMatchedEarly  masterEarly{ilEarly}.analysis.dff.corrTrialsMatched(:)'];
%     corrTrialsOrthoEarly = [corrTrialsOrthoEarly masterEarly{ilEarly}.analysis.dff.corrTrialsOrtho(:)'];  
end

shortDistanceCorrAdult= []; longDistanceCorrAdult= [];
corrTrialsMatchedAdult= []; corrTrialsOrthoAdult= [];
for ilAdult= 1:length(AdultInd)
    datapath = [adata_dir char(exp_info.animal{AdultInd(ilAdult)}) filesep char(exp_info.exp_id{AdultInd(ilAdult)}) filesep];
    try
        masterAdult{ilAdult} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterAdult{ilAdult} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterAdult{ilAdult}.metadata.ferret = char(exp_info.animal{AdultInd(ilAdult)});
    masterAdult{ilAdult}.metadata.expID = exp_info.exp_id{AdultInd(ilAdult)};
    temp = load(fullfile(datapath, 'PopData.mat'), 'PopAnalysis');
    masterAdult{ilAdult}.PopAnalysis = temp.PopAnalysis;
    clear temp;
    respAdult{ilAdult} = find([masterAdult{ilAdult}.analysis.dff.roi.isResponseSignificant] == 1);
%     shortDistanceCorrAdult = [shortDistanceCorrAdult masterAdult{ilAdult}.analysis.dff.shortDistanceCorr'];
%     longDistanceCorrAdult= [longDistanceCorrAdult masterAdult{ilAdult}.analysis.dff.longDistanceCorr'];
%     corrTrialsMatchedAdult = [corrTrialsMatchedAdult  masterAdult{ilAdult}.analysis.dff.corrTrialsMatched(:)'];
%     corrTrialsOrthoAdult = [corrTrialsOrthoAdult masterAdult{ilAdult}.analysis.dff.corrTrialsOrtho(:)'];  
end

shortDistanceCorrNaiveV1= []; longDistanceCorrNaiveV1= [];
corrTrialsMatchedNaiveV1= []; corrTrialsOrthoNaiveV1= [];
for ilNaiveV1 = 1:length(NaiveIndV1)
    datapath = [adata_dir char(exp_info.animal{NaiveIndV1(ilNaiveV1)}) filesep char(exp_info.exp_id{NaiveIndV1(ilNaiveV1)}) filesep];
    try
        masterNaiveV1{ilNaiveV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterNaiveV1{ilNaiveV1} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'metadata','analysis');
    end
    masterNaiveV1{ilNaiveV1}.metadata.ferret = char(exp_info.animal{NaiveIndV1(ilNaiveV1)});
    masterNaiveV1{ilNaiveV1}.metadata.expID = exp_info.exp_id{NaiveIndV1(ilNaiveV1)};
    temp = load(fullfile(datapath, 'PopData.mat'), 'PopAnalysis');
    masterNaiveV1{ilNaiveV1}.PopAnalysis = temp.PopAnalysis;
    clear temp;
    respNaiveV1{ilNaiveV1} = find([masterNaiveV1{ilNaiveV1}.analysis.dff.roi.isResponseSignificant] == 1);
%     shortDistanceCorrNaiveV1 = [shortDistanceCorrNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.shortDistanceCorr'];
%     longDistanceCorrNaiveV1 = [longDistanceCorrNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.longDistanceCorr'];
%     corrTrialsMatchedNaiveV1 = [corrTrialsMatchedNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.corrTrialsMatched(:)'];
%     corrTrialsOrthoNaiveV1 = [corrTrialsOrthoNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.corrTrialsOrtho(:)'];   
end

% shortDistanceCorrEarlyV1= []; longDistanceCorrEarlyV1= [];
% corrTrialsMatchedEarlyV1= []; corrTrialsOrthoEarlyV1= [];
for ilEarlyV1 = 1:length(EarlyIndV1)
    datapath = [adata_dir char(exp_info.animal{EarlyIndV1(ilEarlyV1)}) filesep char(exp_info.exp_id{EarlyIndV1(ilEarlyV1)}) filesep];
    try
        masterEarlyV1{ilEarlyV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterEarlyV1{ilEarlyV1} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterEarlyV1{ilEarlyV1}.metadata.ferret = char(exp_info.animal{EarlyIndV1(ilEarlyV1)});
    masterEarlyV1{ilEarlyV1}.metadata.expID = exp_info.exp_id{EarlyIndV1(ilEarlyV1)};
    temp = load(fullfile(datapath, 'PopData.mat'), 'PopAnalysis');
    masterEarlyV1{ilEarlyV1}.PopAnalysis = temp.PopAnalysis;
    clear temp;
    respEarlyV1{ilEarlyV1} = find([masterEarlyV1{ilEarlyV1}.analysis.dff.roi.isResponseSignificant] == 1);
%     shortDistanceCorrEarlyV1 = [shortDistanceCorrEarlyV1  masterEarlyV1{ilEarlyV1}.analysis.dff.shortDistanceCorr'];
%     longDistanceCorrEarlyV1 = [longDistanceCorrEarlyV1  masterEarlyV1{ilEarlyV1}.analysis.dff.longDistanceCorr'];
%     corrTrialsMatchedEarlyV1 = [corrTrialsMatchedEarlyV1  masterEarlyV1{ilEarlyV1}.analysis.dff.corrTrialsMatched(:)'];
%     corrTrialsOrthoEarlyV1 = [corrTrialsOrthoEarlyV1 masterEarlyV1{ilEarlyV1}.analysis.dff.corrTrialsOrtho(:)']; 
end

% shortDistanceCorrAdultV1= []; longDistanceCorrAdultV1= [];
% corrTrialsMatchedAdultV1= []; corrTrialsOrthoAdultV1= [];
for ilAdultV1= 1:length(AdultIndV1)
    datapath = [adata_dir char(exp_info.animal{AdultIndV1(ilAdultV1)}) filesep char(exp_info.exp_id{AdultIndV1(ilAdultV1)}) filesep];
    try
        masterAdultV1{ilAdultV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterAdultV1{ilAdultV1} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterAdultV1{ilAdultV1}.metadata.ferret = char(exp_info.animal{AdultIndV1(ilAdultV1)});
    masterAdultV1{ilAdultV1}.metadata.expID = exp_info.exp_id{AdultIndV1(ilAdultV1)};
    temp = load(fullfile(datapath, 'PopData.mat'), 'PopAnalysis');
    masterAdultV1{ilAdultV1}.PopAnalysis = temp.PopAnalysis;
    clear temp;
    respAdultV1{ilAdultV1} = find([masterAdultV1{ilAdultV1}.analysis.dff.roi.isResponseSignificant] == 1);
%     shortDistanceCorrAdultV1 = [shortDistanceCorrAdultV1 masterAdultV1{ilAdultV1}.analysis.dff.shortDistanceCorr'];
%     longDistanceCorrAdultV1= [longDistanceCorrAdultV1 masterAdultV1{ilAdultV1}.analysis.dff.longDistanceCorr'];
%     corrTrialsMatchedAdultV1 = [corrTrialsMatchedAdultV1  masterAdultV1{ilAdultV1}.analysis.dff.corrTrialsMatched(:)'];
%     corrTrialsOrthoAdultV1 = [corrTrialsOrthoAdultV1 masterAdultV1{ilAdultV1}.analysis.dff.corrTrialsOrtho(:)'];  
end
%% Data consolidation for single cells
OSIFitNaive = []; OriCircVarNaive = [];DSINaive = [];isRespNaive = [];
prefSfNaive = []; sfSINaive = []; SFVarNaive = [];
HI_SfNaive = double.empty(0,6); HI_Naive = double.empty(0,6); oriPairsNaive = [];

for ferret =1:length(NaiveInd)
    OSIFitNaive = [OSIFitNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).OSIFit];
    OriCircVarNaive = [OriCircVarNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).OriCircVar];
    DSINaive = [DSINaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).DSI];
    isRespNaive = [isRespNaive masterNaive{ferret}.analysis.dff.roi.isResponseSignificant];
    prefSfNaive = [prefSfNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).prefSf];
    sfSINaive = [sfSINaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).SFSI];
    SFVarNaive = [SFVarNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).SFVar]; 
    HI_Naive = [HI_Naive; masterNaive{ferret}.PopAnalysis.dff.ori_cells.HomeogeneityIndex];
    HI_SfNaive = [HI_SfNaive; masterNaive{ferret}.PopAnalysis.dff.HI_SF];
    oriPairsNaive = [oriPairsNaive masterNaive{ferret}.PopAnalysis.dff.ori_pair];
end

OSIFitEarly = []; OriCircVarEarly = [];DSIEarly= []; isRespEarly = [];
prefSfEarly = []; sfSIEarly= []; SFVarEarly = [];
HI_SfEarly = double.empty(0,6); HI_Early = double.empty(0,6);oriPairsEarly = [];

for ferret =1:length(EarlyInd)
    OSIFitEarly = [OSIFitEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).OSIFit];
    OriCircVarEarly = [OriCircVarEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).OriCircVar];
    DSIEarly = [DSIEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).DSI];
    isRespEarly = [isRespEarly masterEarly{ferret}.analysis.dff.roi.isResponseSignificant];
    prefSfEarly = [prefSfEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).prefSf];
    sfSIEarly = [sfSIEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).SFSI];
    SFVarEarly = [SFVarEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).SFVar];
    HI_Early = [HI_Early; masterEarly{ferret}.PopAnalysis.dff.ori_cells.HomeogeneityIndex];
    HI_SfEarly = [HI_SfEarly; masterEarly{ferret}.PopAnalysis.dff.HI_SF];
    oriPairsEarly = [oriPairsEarly masterEarly{ferret}.PopAnalysis.dff.ori_pair];
end

OSIFitAdult = []; OriCircVarAdult  = []; DSIAdult = []; isRespAdult  = [];
prefSfAdult= []; sfSIAdult= []; SFVarAdult = [];
HI_SfAdult = double.empty(0,6); HI_Adult = double.empty(0,6);oriPairsAdult = [];

for ferret =1:length(AdultInd)
    OSIFitAdult  = [OSIFitAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).OSIFit];
    OriCircVarAdult  = [OriCircVarAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).OriCircVar];
    DSIAdult  = [DSIAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).DSI];
    isRespAdult = [isRespAdult masterAdult{ferret}.analysis.dff.roi.isResponseSignificant];
    prefSfAdult = [prefSfAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).prefSf];
    sfSIAdult = [sfSIAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).SFSI];
    SFVarAdult = [SFVarAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).SFVar];
    HI_Adult = [HI_Adult; masterAdult{ferret}.PopAnalysis.dff.ori_cells.HomeogeneityIndex];
    HI_SfAdult = [HI_SfAdult; masterAdult{ferret}.PopAnalysis.dff.HI_SF];
    oriPairsAdult = [oriPairsAdult masterAdult{ferret}.PopAnalysis.dff.ori_pair];
end

OSIFitNaiveV1 = []; OriCircVarNaiveV1 = []; DSINaiveV1 = []; isRespNaiveV1 = [];
prefSfNaiveV1 = []; sfSINaiveV1 = []; SFVarNaiveV1 = [];
HI_SfNaiveV1 = double.empty(0,6); HI_NaiveV1 = double.empty(0,6); oriPairsNaiveV1 = [];

for ferret =1:length(NaiveIndV1)
    OSIFitNaiveV1 = [OSIFitNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).OSIFit];
    OriCircVarNaiveV1 = [OriCircVarNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).OriCircVar];
    DSINaiveV1 = [DSINaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).DSI];
    isRespNaiveV1 = [isRespNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi.isResponseSignificant];
    prefSfNaiveV1 = [prefSfNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).prefSf];
    sfSINaiveV1 = [sfSINaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).SFSI];
    SFVarNaiveV1 = [SFVarNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).SFVar];
    HI_NaiveV1 = [HI_NaiveV1; masterNaiveV1{ferret}.PopAnalysis.dff.ori_cells.HomeogeneityIndex];
    HI_SfNaiveV1 = [HI_SfNaiveV1; masterNaiveV1{ferret}.PopAnalysis.dff.HI_SF];
    oriPairsNaiveV1 = [oriPairsNaiveV1 masterNaiveV1{ferret}.PopAnalysis.dff.ori_pair];
end

OSIFitEarlyV1 = []; OriCircVarEarlyV1 = []; DSIEarlyV1= []; isRespEarlyV1 = [];
prefSfEarlyV1 = []; sfSIEarlyV1= []; SFVarEarlyV1 = [];
HI_SfEarlyV1 = double.empty(0,6); HI_EarlyV1 = double.empty(0,6); oriPairsEarlyV1 = [];

for ferret =1:length(EarlyIndV1)
    OSIFitEarlyV1 = [OSIFitEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).OSIFit];
    OriCircVarEarlyV1 = [OriCircVarEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).OriCircVar];
    DSIEarlyV1 = [DSIEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).DSI];
    isRespEarlyV1 = [isRespEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi.isResponseSignificant];
    prefSfEarlyV1 = [prefSfEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).prefSf];
    sfSIEarlyV1 = [sfSIEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).SFSI];
    SFVarEarlyV1 = [SFVarEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).SFVar];
    HI_EarlyV1 = [HI_EarlyV1; masterEarlyV1{ferret}.PopAnalysis.dff.ori_cells.HomeogeneityIndex];
    HI_SfEarlyV1 = [HI_SfEarlyV1; masterEarlyV1{ferret}.PopAnalysis.dff.HI_SF];
    oriPairsEarlyV1 = [oriPairsEarlyV1 masterEarlyV1{ferret}.PopAnalysis.dff.ori_pair];
end

OSIFitAdultV1 = []; OriCircVarAdultV1  = []; DSIAdultV1 = []; isRespAdultV1  = [];
prefSfAdultV1= []; sfSIAdultV1= []; SFVarAdultV1 = [];
HI_SfAdultV1 = double.empty(0,6); HI_AdultV1 = double.empty(0,6); oriPairsAdultV1 = [];

for ferret =1:length(AdultIndV1)
    OSIFitAdultV1  = [OSIFitAdultV1  masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).OSIFit];
    OriCircVarAdultV1  = [OriCircVarAdultV1  masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).OriCircVar];
    DSIAdultV1  = [DSIAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).DSI];
    isRespAdultV1 = [isRespAdultV1 masterAdultV1{ferret}.analysis.dff.roi.isResponseSignificant];
    prefSfAdultV1 = [prefSfAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).prefSf];
    sfSIAdultV1 = [sfSIAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).SFSI];
    SFVarAdultV1 = [SFVarAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).SFVar];
    HI_AdultV1 = [HI_AdultV1; masterAdultV1{ferret}.PopAnalysis.dff.ori_cells.HomeogeneityIndex];
    HI_SfAdultV1 = [HI_SfAdultV1; masterAdultV1{ferret}.PopAnalysis.dff.HI_SF];
    oriPairsAdultV1 = [oriPairsAdultV1 masterAdultV1{ferret}.PopAnalysis.dff.ori_pair];
end


%% 2) Plot results single cells
% a)responsiveness and selectivity portions
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 3, 1)
allNaive = length(isRespNaive);
nonRespNaive = length(find([isRespNaive] == 0)) ./allNaive;
percRespNaive = length(find([isRespNaive] == 1)) ./allNaive;
h = pie([nonRespNaive percRespNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocNaive(4,:));
title('Responsive Naive')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 2)
allEarly = length(isRespEarly);
nonRespEarly = length(find([isRespEarly] == 0)) ./allEarly;
percRespEarly = length(find([isRespEarly] == 1)) ./allEarly;
h = pie([nonRespEarly percRespEarly]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocEarly(4,:));
title('Responsive Early')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 3)
allAdult = length(isRespAdult);
nonRespAdult = length(find([isRespAdult] == 0)) ./allAdult;
percRespAdult = length(find([isRespAdult] == 1)) ./allAdult;
h = pie([nonRespAdult percRespAdult]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocAdult(4,:));
title('Responsive Adult')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 4)
allNaiveV1 = length(isRespNaiveV1);
nonRespNaiveV1 = length(find([isRespNaiveV1] == 0)) ./allNaiveV1;
percRespNaiveV1 = length(find([isRespNaiveV1] == 1)) ./allNaiveV1;
h = pie([nonRespNaiveV1 percRespNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocNaiveV1(4,:));
title('Responsive Naive V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 5)
allEarlyV1 = length(isRespEarlyV1);
nonRespEarlyV1 = length(find([isRespEarlyV1] == 0)) ./allEarlyV1;
percRespEarlyV1 = length(find([isRespEarlyV1] == 1)) ./allEarlyV1;
h = pie([nonRespEarlyV1 percRespEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocEarlyV1(4,:));
title('Responsive Early V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 6)
allAdultV1 = length(isRespAdultV1);
nonRespAdultV1 = length(find([isRespAdultV1] == 0)) ./allAdultV1;
percRespAdultV1 = length(find([isRespAdultV1] == 1)) ./allAdultV1;
h = pie([nonRespAdultV1 percRespAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocAdultV1(4,:));
title('Responsive Adult V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')


set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'resp_cells.png'))

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


%% c) prefSf
figure
subplot(2,3,1)
plot(nanmedian(prefSfNaive),1.05 * max(histcounts(prefSfNaive)),'v','MarkerSize', 8','MarkerEdgeColor',cocNaive(5,:),'MarkerFaceColor',cocNaive(4,:)); hold on
text(nanmedian(prefSfNaive),1.2 * max(histcounts(prefSfNaive)),num2str(round(100*nanmedian(prefSfNaive))/100),'HorizontalAlignment','center', 'Color', cocNaive(4,:), 'FontSize', 12')
histogram(prefSfNaive, 20, 'FaceColor', cocNaive(4,:), 'EdgeColor', cocNaive(5,:));
ylabel('Cells');
title('Naive')
xlim([0 0.2])
set(gca,'Box','off');

subplot(2,3,2)
plot(nanmedian(prefSfEarly),1.05 * max(histcounts(prefSfEarly)),'v','MarkerSize', 8','MarkerEdgeColor',cocEarly(5,:),'MarkerFaceColor',cocEarly(4,:)); hold on
text(nanmedian(prefSfEarly),1.2 * max(histcounts(prefSfEarly)),num2str(round(100*nanmedian(prefSfEarly))/100),'HorizontalAlignment','center', 'Color', cocEarly(5,:), 'FontSize', 12')
histogram(prefSfEarly, 20, 'FaceColor', cocEarly(4,:), 'EdgeColor', cocEarly(5,:));
ylabel('Cells');
title('Early')
xlim([0 0.2])
set(gca,'Box','off');

subplot(2,3,3)
plot(nanmedian(prefSfAdult),1.05 * max(histcounts(prefSfAdult)),'v','MarkerSize', 8','MarkerEdgeColor',cocAdult(5,:),'MarkerFaceColor',cocAdult(4,:)); hold on
text(nanmedian(prefSfAdult),1.2 * max(histcounts(prefSfAdult)),num2str(round(100*nanmedian(prefSfAdult))/100),'HorizontalAlignment','center', 'Color', cocAdult(5,:), 'FontSize', 12')
histogram(prefSfAdult, 20, 'FaceColor', cocAdult(4,:), 'EdgeColor', cocAdult(5,:));
ylabel('Cells');
title('Adult')
xlim([0 0.2])
set(gca,'Box','off');

subplot(2,3,4)
plot(nanmedian(prefSfNaiveV1),1.05 * max(histcounts(prefSfNaiveV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocNaiveV1(5,:),'MarkerFaceColor',cocNaiveV1(4,:)); hold on
text(nanmedian(prefSfNaiveV1),1.2 * max(histcounts(prefSfNaiveV1)),num2str(round(100*nanmedian(prefSfNaiveV1))/100),'HorizontalAlignment','center', 'Color', cocNaiveV1(5,:), 'FontSize', 12')
histogram(prefSfNaiveV1, 20, 'FaceColor', cocNaiveV1(4,:), 'EdgeColor', cocNaiveV1(5,:));
ylabel('Cells');
xlabel('prefSF');
xlim([0 0.2])
title('Naive')
set(gca,'Box','off');


subplot(2,3,5)
plot(nanmedian(prefSfEarlyV1),1.05 * max(histcounts(prefSfEarlyV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocEarlyV1(5,:),'MarkerFaceColor',cocEarlyV1(4,:)); hold on
text(nanmedian(prefSfEarlyV1),1.2 * max(histcounts(prefSfEarlyV1)),num2str(round(100*nanmedian(prefSfEarlyV1))/100),'HorizontalAlignment','center', 'Color', cocEarlyV1(5,:), 'FontSize', 12')
histogram(prefSfEarlyV1, 20, 'FaceColor', cocEarlyV1(4,:), 'EdgeColor', cocEarlyV1(5,:));
ylabel('Cells');
xlabel('pref SF');
xlim([0 0.2])
title('Early')
set(gca,'Box','off');

subplot(2,3,6)
plot(nanmedian(prefSfAdultV1),1.05 * max(histcounts(prefSfAdultV1)),'v','MarkerSize', 8','MarkerEdgeColor',cocAdultV1(5,:),'MarkerFaceColor',cocAdultV1(4,:)); hold on
text(nanmedian(prefSfAdultV1),1.2 * max(histcounts(prefSfAdultV1)),num2str(round(100*nanmedian(prefSfAdultV1))/100),'HorizontalAlignment','center', 'Color', cocAdultV1(5,:), 'FontSize', 12')
histogram(prefSfAdultV1, 20, 'FaceColor', cocAdultV1(4,:), 'EdgeColor', cocAdultV1(5,:));
ylabel('Cells');
xlabel('prefSF');
xlim([0 0.2])
title('Adult')
set(gca,'Box','off');

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'prefSf.png'))

%pie charts of adult pref
perLowSfNaive = length(find(prefSfNaive <= 0.011)) /length(prefSfNaive)+0.000001;
perHighSfNaive = length(find(prefSfNaive > 0.031))/length(prefSfNaive)+0.000001;
middleSfNaive = 1-perLowSfNaive-perHighSfNaive+0.000001;
perLowSfNaiveV1 = length(find(prefSfNaiveV1 <= 0.011)) /length(prefSfNaiveV1)+0.000001;
perHighSfNaiveV1 = length(find(prefSfNaiveV1 > 0.031))/length(prefSfNaiveV1)+0.000001;
middleSfNaiveV1 = 1-perLowSfNaiveV1-perHighSfNaiveV1+0.000001;
perLowSfEarly = length(find(prefSfEarly <= 0.011)) /length(prefSfEarly)+0.000001;
perHighSfEarly = length(find(prefSfEarly > 0.031))/length(prefSfEarly)+0.000001;
middleSfEarly = 1-perLowSfEarly-perHighSfEarly+0.000001;
perLowSfEarlyV1 = length(find(prefSfEarlyV1 <= 0.011)) /length(prefSfEarlyV1)+0.000001;
perHighSfEarlyV1 = length(find(prefSfEarlyV1 > 0.031))/length(prefSfEarlyV1)+0.000001;
middleSfEarlyV1 = 1-perLowSfEarlyV1-perHighSfEarlyV1+0.000001;
perLowSfAdult = length(find(prefSfAdult <= 0.011)) /length(prefSfAdult)+0.000001;
perHighSfAdult = length(find(prefSfAdult > 0.031))/length(prefSfAdult)+0.000001;
middleSfAdult = 1-perLowSfAdult-perHighSfAdult+0.000001;
perLowSfAdultV1 = length(find(prefSfAdultV1 <= 0.011)) /length(prefSfAdultV1)+0.000001;
perHighSfAdultV1 = length(find(prefSfAdultV1 > 0.031))/length(prefSfAdultV1)+0.000001;
middleSfAdultV1 = 1-perLowSfAdultV1-perHighSfAdultV1+0.000001;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1)
h = pie([perLowSfNaive middleSfNaive perHighSfNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaive(2,:));
set(hp(2), 'FaceColor', cocEarly(3,:));
set(hp(3), 'FaceColor', cocAdult(4,:));
title('Naive')
legend({'<= 0.01 cpd','0.01-0.03 cpd', '>0.03 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2,3,2)
h = pie([perLowSfEarly middleSfEarly perHighSfEarly]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaive(2,:));
set(hp(2), 'FaceColor', cocEarly(3,:));
set(hp(3), 'FaceColor', cocAdult(4,:));
title('Early')
legend({'<= 0.01 cpd','0.01-0.03 cpd', '>0.03 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2,3,3)
h = pie([perLowSfAdult middleSfAdult perHighSfAdult]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaive(2,:));
set(hp(2), 'FaceColor', cocEarly(3,:));
set(hp(3), 'FaceColor', cocAdult(4,:));
title('Adult')
legend({'<= 0.01 cpd','0.01-0.03 cpd', '>0.03 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2,3,4)
h = pie([perLowSfNaiveV1 middleSfNaiveV1 perHighSfNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaiveV1(2,:));
set(hp(2), 'FaceColor', cocEarlyV1(3,:));
set(hp(3), 'FaceColor', cocAdultV1(4,:));
title('Naive V1')
legend({'<= 0.01 cpd','0.01-0.03 cpd', '>0.03 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2,3,5)
h = pie([perLowSfEarlyV1 middleSfEarlyV1 perHighSfEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaiveV1(2,:));
set(hp(2), 'FaceColor', cocEarlyV1(3,:));
set(hp(3), 'FaceColor', cocAdultV1(4,:));
title('Early V1')
legend({'<= 0.01 cpd','0.01-0.03 cpd', '>0.03 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2,3,6)
h = pie([perLowSfAdultV1 middleSfAdultV1 perHighSfAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocNaiveV1(2,:));
set(hp(2), 'FaceColor', cocEarlyV1(3,:));
set(hp(3), 'FaceColor', cocAdultV1(4,:));
title('Adult V1')
legend({'<= 0.01 cpd','0.01-0.03 cpd', '>0.03 cpd'}, 'Location', 'southoutside')
legend('boxoff')
%% d) sfSI
figure
subplot(1,6,1)
distributionPlot(sfSINaive','color', cocNaive(4,:)); hold on
boxplot(sfSINaive,'Label', {'Naive'})
set(gca,'Box','off');
ylim([0 5])

subplot(1,6,2)
distributionPlot(sfSIEarly','color', cocEarly(4,:)); hold on
boxplot(sfSIEarly, 'Label', {'Early'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

subplot(1,6,3)
distributionPlot(sfSIAdult','color', cocAdult(4,:)); hold on
boxplot(sfSIAdult, 'Label', {'Adult'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

subplot(1,6,4)
distributionPlot(sfSINaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(sfSINaiveV1, 'Label', {'Naive V1'})
set(gca,'Box','off');
ylim([0 5])

subplot(1, 6, 5)
distributionPlot(sfSIEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(sfSIEarlyV1, 'Label', {'Early V1'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

subplot(1, 6, 6)
distributionPlot(sfSIAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(sfSIAdultV1, 'Label', {'Adult V1'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

set(gcf, 'color', 'w');
title('SF Selectivity Index')
saveas(gcf, fullfile(save_dir, 'SFSI.png'))

%% e) SFVar
figure
subplot(1,6,1)
distributionPlot(SFVarNaive','color', cocNaive(4,:)); hold on
boxplot(SFVarNaive,'Label', {'Naive'})
set(gca,'Box','off');
ylim([0 5])

subplot(1,6,2)
distributionPlot(SFVarEarly','color', cocEarly(4,:)); hold on
boxplot(SFVarEarly, 'Label', {'Early'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

subplot(1,6,3)
distributionPlot(SFVarAdult','color', cocAdult(4,:)); hold on
boxplot(SFVarAdult, 'Label', {'Adult'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

subplot(1,6,4)
distributionPlot(SFVarNaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(SFVarNaiveV1, 'Label', {'Naive V1'})
set(gca,'Box','off');
ylim([0 5])

subplot(1, 6, 5)
distributionPlot(SFVarEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(SFVarEarlyV1, 'Label', {'Early V1'})
set(gca,'box','off','ycolor','w')
ylim([0 5])

subplot(1, 6, 6)
distributionPlot(SFVarAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(SFVarAdultV1, 'Label', {'Adult V1'})
set(gca,'box','off','ycolor','w')
title('SF Variability Index')
ylim([0 5])

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'SFVar.png'))

%% 3.) Plot results pop analysis
% a) homogeneity index
figure
subplot(1,2,1)
radius = linspace(0,250,6);
errorbar(radius(2:end), nanmean(HI_Naive(:, 1:end-1)), nanstd(HI_Naive(:, 1:end-1))/sqrt(length(HI_Naive)),'o-','Color', cocNaive(4,:), 'MarkerFaceColor', cocNaive(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_Early(:, 1:end-1)), nanstd(HI_Early(:, 1:end-1))/sqrt(length(HI_Early)),'o-', 'Color', cocEarly(4,:), 'MarkerFaceColor', cocEarly(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_Adult(:, 1:end-1)), nanstd(HI_Adult(:, 1:end-1))/sqrt(length(HI_Adult)),'o-', 'Color', cocAdult(4,:), 'MarkerFaceColor', cocAdult(4,:))
xlabel('Maximal radius in \mum')
ylabel('Homeogeneity Index')
xlim([0 250])
ylim([0 1])
legend('Naive', 'Early', 'Adult', 'Location', 'SouthEast')
legend('boxoff');
set(gca,'Box','off');

subplot(1,2,2)
errorbar(radius(2:end), nanmean(HI_NaiveV1(:, 1:end-1)), nanstd(HI_NaiveV1(:, 1:end-1))/sqrt(length(HI_NaiveV1)),'o-','Color', cocNaiveV1(4,:), 'MarkerFaceColor', cocNaiveV1(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_EarlyV1(:, 1:end-1)), nanstd(HI_EarlyV1(:, 1:end-1))/sqrt(length(HI_EarlyV1)),'o-', 'Color', cocEarlyV1(4,:), 'MarkerFaceColor', cocEarlyV1(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_AdultV1(:, 1:end-1)), nanstd(HI_AdultV1(:, 1:end-1))/sqrt(length(HI_AdultV1)),'o-', 'Color', cocAdultV1(4,:), 'MarkerFaceColor', cocAdultV1(4,:))
xlabel('Maximal radius in \mum')
ylabel('Homeogeneity Index')
xlim([0 250])
ylim([0 1])
legend('Naive V1', 'Early V1', 'Adult V1', 'Location', 'SouthEast')
legend('boxoff');

set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'HomeogeneityIndex.png'))

% b) homogeneity index for sf
figure
subplot(1,2,1)
radius = linspace(0,250,6);
errorbar(radius(2:end), nanmean(HI_SfNaive(:, 1:end-1)), nanstd(HI_SfNaive(:, 1:end-1))/sqrt(length(HI_SfNaive)),'o-','Color', cocNaive(4,:), 'MarkerFaceColor', cocNaive(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_SfEarly(:, 1:end-1)), nanstd(HI_SfEarly(:, 1:end-1))/sqrt(length(HI_SfEarly)),'o-', 'Color', cocEarly(4,:), 'MarkerFaceColor', cocEarly(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_SfAdult(:, 1:end-1)), nanstd(HI_SfAdult(:, 1:end-1))/sqrt(length(HI_SfAdult)),'o-', 'Color', cocAdult(4,:), 'MarkerFaceColor', cocAdult(4,:))
xlabel('Maximal radius in \mum')
ylabel('Proportion cells with same Sf')
xlim([0 250])
ylim([0 1])
legend('Naive', 'Early', 'Adult', 'Location', 'SouthEast')
legend('boxoff');
set(gca,'Box','off');

subplot(1,2,2)
errorbar(radius(2:end), nanmean(HI_SfNaiveV1(:, 1:end-1)), nanstd(HI_SfNaiveV1(:, 1:end-1))/sqrt(length(HI_SfNaiveV1)),'o-','Color', cocNaiveV1(4,:), 'MarkerFaceColor', cocNaiveV1(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_SfEarlyV1(:, 1:end-1)), nanstd(HI_SfEarlyV1(:, 1:end-1))/sqrt(length(HI_SfEarlyV1)),'o-', 'Color', cocEarlyV1(4,:), 'MarkerFaceColor', cocEarlyV1(4,:))
hold on
errorbar(radius(2:end), nanmean(HI_SfAdultV1(:, 1:end-1)), nanstd(HI_SfAdultV1(:, 1:end-1))/sqrt(length(HI_SfAdultV1)),'o-', 'Color', cocAdultV1(4,:), 'MarkerFaceColor', cocAdultV1(4,:))
xlabel('Maximal radius in \mum')
ylabel('Proportion cells with same Sf')
xlim([0 250])
ylim([0 1])
legend('Naive V1', 'Early V1', 'Adult V1', 'Location', 'SouthEast')
legend('boxoff');

set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'ProportioncellswithsameSf.png'))

% c) delta ori vs. distance
maxDis = floor(max([oriPairsNaive.distance oriPairsEarly.distance oriPairsAdult.distance]));
edges = linspace(0, maxDis, 15); edge = edges(1:end-1)+25;

[n, edges, binDist] = histcounts([oriPairsNaive.distance],edges);
meanDeltaOriNaive = zeros(length(edge),1);
semDeltaOriNaive = zeros(length(edge),1);
for bin = 1:length(edge)
    meanDeltaOriNaive(bin) = nanmean([oriPairsNaive(binDist == bin).deltaOri]);
    semDeltaOriNaive(bin) = nanstd([oriPairsNaive(binDist == bin).deltaOri])/sqrt(n(bin));
end
[n, edges, binDist] = histcounts([oriPairsEarly.distance],edges);
meanDeltaOriEarly = zeros(length(edge),1);
semDeltaOriEarly = zeros(length(edge),1);
for bin = 1:length(edge)
    meanDeltaOriEarly(bin) = nanmean([oriPairsEarly(binDist == bin).deltaOri]);
    semDeltaOriEarly(bin) = nanstd([oriPairsEarly(binDist == bin).deltaOri])/sqrt(n(bin));
end
[n, edges, binDist] = histcounts([oriPairsAdult.distance],edges);
meanDeltaOriAdult = zeros(length(edge),1);
semDeltaOriAdult = zeros(length(edge),1);
for bin = 1:length(edge)
    meanDeltaOriAdult(bin) = nanmean([oriPairsAdult(binDist == bin).deltaOri]);
    semDeltaOriAdult(bin) = nanstd([oriPairsAdult(binDist == bin).deltaOri])/sqrt(n(bin));
end

maxDis = floor(max([oriPairsNaiveV1.distance oriPairsEarlyV1.distance oriPairsAdultV1.distance]));
edgesV1 = linspace(0, maxDis, 15); edgeV1 = edgesV1(1:end-1)+25;

[n, edgesV1, binDist] = histcounts([oriPairsNaiveV1.distance],edgesV1);
meanDeltaOriNaiveV1 = zeros(length(edgeV1),1);
semDeltaOriNaiveV1 = zeros(length(edgeV1),1);
for bin = 1:length(edgeV1)
    meanDeltaOriNaiveV1(bin) = nanmean([oriPairsNaiveV1(binDist == bin).deltaOri]);
    semDeltaOriNaiveV1(bin) = nanstd([oriPairsNaiveV1(binDist == bin).deltaOri])/sqrt(n(bin));
end
[n, edgesV1, binDist] = histcounts([oriPairsEarlyV1.distance],edgesV1);
meanDeltaOriEarlyV1 = zeros(length(edgeV1),1);
semDeltaOriEarlyV1 = zeros(length(edgeV1),1);
for bin = 1:length(edgeV1)
    meanDeltaOriEarlyV1(bin) = nanmean([oriPairsEarlyV1(binDist == bin).deltaOri]);
    semDeltaOriEarlyV1(bin) = nanstd([oriPairsEarlyV1(binDist == bin).deltaOri])/sqrt(n(bin));
end
[n, edgesV1, binDist] = histcounts([oriPairsAdultV1.distance],edgesV1);
meanDeltaOriAdultV1 = zeros(length(edgeV1),1);
semDeltaOriAdultV1 = zeros(length(edgeV1),1);
for bin = 1:length(edgeV1)
    meanDeltaOriAdultV1(bin) = nanmean([oriPairsAdultV1(binDist == bin).deltaOri]);
    semDeltaOriAdultV1(bin) = nanstd([oriPairsAdultV1(binDist == bin).deltaOri])/sqrt(n(bin));
end

figure
subplot(1,2,1)
errorbar(edges(2:end),meanDeltaOriNaive, semDeltaOriNaive, 'o-', 'Color', cocNaive(4,:), 'MarkerFaceColor', cocNaive(4,:))
hold all
errorbar(edges(2:end),meanDeltaOriEarly, semDeltaOriEarly, 'o-', 'Color', cocEarly(4,:), 'MarkerFaceColor', cocEarly(4,:))
hold all
errorbar(edges(2:end),meanDeltaOriAdult, semDeltaOriAdult, 'o-', 'Color', cocAdult(4,:), 'MarkerFaceColor', cocAdult(4,:))
xlabel('Distance in \mum')
ylabel('\DeltaOrientation preference (\circ)')
ylim([0 90])
xlim([0 edges(end)])
legend('Naive', 'Early', 'Adult', 'Location', 'SouthEast')
legend('boxoff');
set(gca,'Box','off');

subplot(1,2,2)
errorbar(edgesV1(2:end),meanDeltaOriNaiveV1, semDeltaOriNaiveV1, 'o-', 'Color', cocNaiveV1(4,:), 'MarkerFaceColor', cocNaiveV1(4,:))
hold all
errorbar(edgesV1(2:end),meanDeltaOriEarlyV1, semDeltaOriEarlyV1, 'o-', 'Color', cocEarlyV1(4,:), 'MarkerFaceColor', cocEarlyV1(4,:))
hold all
errorbar(edgesV1(2:end),meanDeltaOriAdultV1, semDeltaOriAdultV1, 'o-', 'Color', cocAdultV1(4,:), 'MarkerFaceColor', cocAdultV1(4,:))
xlabel('Distance in \mum')
ylabel('\DeltaOrientation preference (\circ)')
ylim([0 90])
xlim([0 edgesV1(end)])
legend('Naive V1', 'Early V1', 'Adult V1', 'Location', 'SouthEast')
legend('boxoff');
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSIvsdistance.png'))

