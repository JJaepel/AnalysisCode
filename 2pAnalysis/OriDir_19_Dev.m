%% switch board for analysis variable
analysisParams = struct;
%which type of stimulus should it run
analysisParams.dataType = 1; %data type: 1 = cells, 2 = axons, 3 = spines
analysisParams.stimType = 1;

% what should it do?
analysisParams.reloadData = 0; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse =1; %should you reanalyse the data or just plot?
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.plotROIs = 0;   %should you plot traces for all resp ROIs?
analysisParams.plotRespROIsOnly = 0; %should you also plot traces for all non-resp ROIs?
analysisParams.server = 0; %load from the server (1) or the raid (0)
analysisParams.makeROIs = 1;

% analysisParameters
analysisParams.zThresh = 4;
analysisParams.fraction = 0.5;
analysisParams.predictor = 0;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff';
analysisParams.windowStart = 0;
analysisParams.windowStop = 2;
analysisParams.pre = 1;

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
save_dir = [adata_dir filesep '19_Dev_Ori_Dir' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

filePath = 'F:\Organization\Animals\';
file = '2pExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
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
Flag = find(cell2mat(exp_info.flag) == 1);
EarlyInd = setdiff(EarlyInd, Flag);
NaiveInd = setdiff(NaiveInd, Flag);
AdultInd = setdiff(AdultInd, Flag);
EarlyIndV1 = setdiff(EarlyIndV1, Flag);
NaiveIndV1 = setdiff(NaiveIndV1, Flag);
AdultIndV1 = setdiff(AdultIndV1, Flag);

%reanalyze everything if necessary
if ~analysisParams.select  
    allExpInd = [EarlyInd EarlyIndV1 NaiveInd NaiveIndV1 AdultInd AdultIndV1];
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
    respNaive{ilNaive} = find([masterNaive{ilNaive}.analysis.dff.roi.isResponseSignificant] == 1);
    shortDistanceCorrNaive = [shortDistanceCorrNaive  masterNaive{ilNaive}.analysis.dff.shortDistanceCorr'];
    longDistanceCorrNaive = [longDistanceCorrNaive  masterNaive{ilNaive}.analysis.dff.longDistanceCorr'];
    corrTrialsMatchedNaive = [corrTrialsMatchedNaive  masterNaive{ilNaive}.analysis.dff.corrTrialsMatched(:)'];
    corrTrialsOrthoNaive = [corrTrialsOrthoNaive  masterNaive{ilNaive}.analysis.dff.corrTrialsOrtho(:)'];    
end

shortDistanceCorrEarly= [];
longDistanceCorrEarly= [];
corrTrialsMatchedEarly= [];
corrTrialsOrthoEarly= [];
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
    shortDistanceCorrEarly = [shortDistanceCorrEarly  masterEarly{ilEarly}.analysis.dff.shortDistanceCorr'];
    longDistanceCorrEarly = [longDistanceCorrEarly  masterEarly{ilEarly}.analysis.dff.longDistanceCorr'];
    corrTrialsMatchedEarly = [corrTrialsMatchedEarly  masterEarly{ilEarly}.analysis.dff.corrTrialsMatched(:)'];
    corrTrialsOrthoEarly = [corrTrialsOrthoEarly masterEarly{ilEarly}.analysis.dff.corrTrialsOrtho(:)'];  
end

shortDistanceCorrAdult= [];
longDistanceCorrAdult= [];
corrTrialsMatchedAdult= [];
corrTrialsOrthoAdult= [];
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
    shortDistanceCorrAdult = [shortDistanceCorrAdult masterAdult{ilAdult}.analysis.dff.shortDistanceCorr'];
    longDistanceCorrAdult= [longDistanceCorrAdult masterAdult{ilAdult}.analysis.dff.longDistanceCorr'];
    corrTrialsMatchedAdult = [corrTrialsMatchedAdult  masterAdult{ilAdult}.analysis.dff.corrTrialsMatched(:)'];
    corrTrialsOrthoAdult = [corrTrialsOrthoAdult masterAdult{ilAdult}.analysis.dff.corrTrialsOrtho(:)'];  
end

shortDistanceCorrNaiveV1= [];
longDistanceCorrNaiveV1= [];
corrTrialsMatchedNaiveV1= [];
corrTrialsOrthoNaiveV1= [];
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
    shortDistanceCorrNaiveV1 = [shortDistanceCorrNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.shortDistanceCorr'];
    longDistanceCorrNaiveV1 = [longDistanceCorrNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.longDistanceCorr'];
    corrTrialsMatchedNaiveV1 = [corrTrialsMatchedNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.corrTrialsMatched(:)'];
    corrTrialsOrthoNaiveV1 = [corrTrialsOrthoNaiveV1  masterNaiveV1{ilNaiveV1}.analysis.dff.corrTrialsOrtho(:)'];   
end

shortDistanceCorrEarlyV1= [];
longDistanceCorrEarlyV1= [];
corrTrialsMatchedEarlyV1= [];
corrTrialsOrthoEarlyV1= [];
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
    shortDistanceCorrEarlyV1 = [shortDistanceCorrEarlyV1  masterEarlyV1{ilEarlyV1}.analysis.dff.shortDistanceCorr'];
    longDistanceCorrEarlyV1 = [longDistanceCorrEarlyV1  masterEarlyV1{ilEarlyV1}.analysis.dff.longDistanceCorr'];
    corrTrialsMatchedEarlyV1 = [corrTrialsMatchedEarlyV1  masterEarlyV1{ilEarlyV1}.analysis.dff.corrTrialsMatched(:)'];
    corrTrialsOrthoEarlyV1 = [corrTrialsOrthoEarlyV1 masterEarlyV1{ilEarlyV1}.analysis.dff.corrTrialsOrtho(:)']; 
end

shortDistanceCorrAdultV1= [];
longDistanceCorrAdultV1= [];
corrTrialsMatchedAdultV1= [];
corrTrialsOrthoAdultV1= [];
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
    shortDistanceCorrAdultV1 = [shortDistanceCorrAdultV1 masterAdultV1{ilAdultV1}.analysis.dff.shortDistanceCorr'];
    longDistanceCorrAdultV1= [longDistanceCorrAdultV1 masterAdultV1{ilAdultV1}.analysis.dff.longDistanceCorr'];
    corrTrialsMatchedAdultV1 = [corrTrialsMatchedAdultV1  masterAdultV1{ilAdultV1}.analysis.dff.corrTrialsMatched(:)'];
    corrTrialsOrthoAdultV1 = [corrTrialsOrthoAdultV1 masterAdultV1{ilAdultV1}.analysis.dff.corrTrialsOrtho(:)'];  
end
%% Data consolidation for single cells
OSIFitNaive = [];
OriCircVarNaive = [];
cohensDNaive = [];
DSINaive = [];
DirCircVarNaive = [];
isRespNaive = [];
BandthWidthNaive = [];
VINaive = [];
FanoFactorNaive = [];
for ferret =1:length(NaiveInd)
    OSIFitNaive = [OSIFitNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).OSIFit];
    OriCircVarNaive = [OriCircVarNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).OriCircVar];
    cohensDNaive = [cohensDNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).cohensD];
    DSINaive = [DSINaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).DSI];
    DirCircVarNaive = [DirCircVarNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).DirCircVar];
    isRespNaive = [isRespNaive masterNaive{ferret}.analysis.dff.roi.isResponseSignificant];
    BandthWidthNaive = [BandthWidthNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).Bandwidth];
    VINaive = [VINaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).VI];
    FanoFactorNaive = [FanoFactorNaive masterNaive{ferret}.analysis.dff.roi(respNaive{ferret}).fanoFactor]; 
end

OSIFitEarly = [];
OriCircVarEarly = [];
cohensDEarly = [];
DSIEarly= [];
DirCircVarEarly = [];
isRespEarly = [];
BandthWidthEarly = [];
VIEarly = [];
FanoFactorEarly = [];
for ferret =1:length(EarlyInd)
    OSIFitEarly = [OSIFitEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).OSIFit];
    OriCircVarEarly = [OriCircVarEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).OriCircVar];
    cohensDEarly = [cohensDEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).cohensD];
    DSIEarly = [DSIEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).DSI];
    DirCircVarEarly = [DirCircVarEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).DirCircVar];
    isRespEarly = [isRespEarly masterEarly{ferret}.analysis.dff.roi.isResponseSignificant];
    BandthWidthEarly = [BandthWidthEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).Bandwidth];
    VIEarly = [VIEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).VI];
    FanoFactorEarly = [FanoFactorEarly masterEarly{ferret}.analysis.dff.roi(respEarly{ferret}).fanoFactor]; 
end

OSIFitAdult = [];
OriCircVarAdult  = [];
cohensDAdult  = [];
DSIAdult = [];
DirCircVarAdult  = [];
isRespAdult  = [];
BandthWidthAdult= [];
VIAdult = [];
FanoFactorAdult = [];
for ferret =1:length(AdultInd)
    OSIFitAdult  = [OSIFitAdult  masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).OSIFit];
    OriCircVarAdult  = [OriCircVarAdult  masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).OriCircVar];
    cohensDAdult  = [cohensDAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).cohensD];
    DSIAdult  = [DSIAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).DSI];
    DirCircVarAdult  = [DirCircVarAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).DirCircVar];
    isRespAdult = [isRespAdult masterAdult{ferret}.analysis.dff.roi.isResponseSignificant];
    BandthWidthAdult = [BandthWidthAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).Bandwidth];
    VIAdult = [VIAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).VI];
    FanoFactorAdult = [FanoFactorAdult masterAdult{ferret}.analysis.dff.roi(respAdult{ferret}).fanoFactor]; 
end

OSIFitNaiveV1 = [];
OriCircVarNaiveV1 = [];
cohensDNaiveV1 = [];
DSINaiveV1 = [];
DirCircVarNaiveV1 = [];
isRespNaiveV1 = [];
BandthWidthNaiveV1 = [];
VINaiveV1 = [];
FanoFactorNaiveV1 = [];
for ferret =1:length(NaiveIndV1)
    OSIFitNaiveV1 = [OSIFitNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).OSIFit];
    OriCircVarNaiveV1 = [OriCircVarNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).OriCircVar];
    cohensDNaiveV1 = [cohensDNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).cohensD];
    DSINaiveV1 = [DSINaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).DSI];
    DirCircVarNaiveV1 = [DirCircVarNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).DirCircVar];
    isRespNaiveV1 = [isRespNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi.isResponseSignificant];
    BandthWidthNaiveV1 = [BandthWidthNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).Bandwidth];
    VINaiveV1 = [VINaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).VI];
    FanoFactorNaiveV1 = [FanoFactorNaiveV1 masterNaiveV1{ferret}.analysis.dff.roi(respNaiveV1{ferret}).fanoFactor]; 
end

OSIFitEarlyV1 = [];
OriCircVarEarlyV1 = [];
cohensDEarlyV1 = [];
DSIEarlyV1= [];
DirCircVarEarlyV1 = [];
isRespEarlyV1 = [];
BandthWidthEarlyV1 = [];
VIEarlyV1 = [];
FanoFactorEarlyV1 = [];
for ferret =1:length(EarlyIndV1)
    OSIFitEarlyV1 = [OSIFitEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).OSIFit];
    OriCircVarEarlyV1 = [OriCircVarEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).OriCircVar];
    cohensDEarlyV1 = [cohensDEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).cohensD];
    DSIEarlyV1 = [DSIEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).DSI];
    DirCircVarEarlyV1 = [DirCircVarEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).DirCircVar];
    isRespEarlyV1 = [isRespEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi.isResponseSignificant];
    BandthWidthEarlyV1 = [BandthWidthEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).Bandwidth];
    VIEarlyV1 = [VIEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).VI];
    FanoFactorEarlyV1 = [FanoFactorEarlyV1 masterEarlyV1{ferret}.analysis.dff.roi(respEarlyV1{ferret}).fanoFactor];
end

OSIFitAdultV1 = [];
OriCircVarAdultV1  = [];
cohensDAdultV1  = [];
DSIAdultV1 = [];
DirCircVarAdultV1  = [];
isRespAdultV1  = [];
BandthWidthAdultV1= [];
VIAdultV1 = [];
FanoFactorAdultV1 = [];
for ferret =1:length(AdultIndV1)
    OSIFitAdultV1  = [OSIFitAdultV1  masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).OSIFit];
    OriCircVarAdultV1  = [OriCircVarAdultV1  masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).OriCircVar];
    cohensDAdultV1  = [cohensDAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).cohensD];
    DSIAdultV1  = [DSIAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).DSI];
    DirCircVarAdultV1  = [DirCircVarAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).DirCircVar];
    isRespAdultV1 = [isRespAdultV1 masterAdultV1{ferret}.analysis.dff.roi.isResponseSignificant];
    BandthWidthAdultV1 = [BandthWidthAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).Bandwidth];
    VIAdultV1 = [VIAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).VI];
    FanoFactorAdultV1 = [FanoFactorAdultV1 masterAdultV1{ferret}.analysis.dff.roi(respAdultV1{ferret}).fanoFactor]; 
end

%% 2.) Plot resuls populations
% a) short vs. long range correlation of cells
figure
subplot(1,6,1)
distributionPlot(shortDistanceCorrNaive','color', cocNaive(4,:)); hold on
boxplot(shortDistanceCorrNaive,'Label', {'Naive'})
xlabel('< 100 um')
ylim([-1 1])
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(longDistanceCorrNaive','color', cocNaive(3,:)); hold on
boxplot(longDistanceCorrNaive, 'Label', {'Naive'})
xlabel('300-400 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(shortDistanceCorrEarly','color', cocEarly(4,:)); hold on
boxplot(shortDistanceCorrEarly,'Label', {'Early'})
xlabel('< 100 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(longDistanceCorrEarly','color', cocEarly(3,:)); hold on
boxplot(longDistanceCorrEarly, 'Label', {'Early'})
xlabel('300-400 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,5)
distributionPlot(shortDistanceCorrAdult','color', cocAdult(4,:)); hold on
boxplot(shortDistanceCorrAdult,'Label', {'Adult'})
xlabel('< 100 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,6)
distributionPlot(longDistanceCorrAdult','color', cocAdult(3,:)); hold on
boxplot(longDistanceCorrAdult, 'Label', {'Adult'})
xlabel('300-400 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'corrCellsDistance.png'))

%for V1
figure
subplot(1,6,1)
distributionPlot(shortDistanceCorrNaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(shortDistanceCorrNaiveV1,'Label', {'Naive'})
xlabel('< 100 um')
ylim([-1 1])
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(longDistanceCorrNaiveV1','color', cocNaiveV1(3,:)); hold on
boxplot(longDistanceCorrNaiveV1, 'Label', {'Naive'})
xlabel('300-400 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(shortDistanceCorrEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(shortDistanceCorrEarlyV1,'Label', {'Early'})
xlabel('< 100 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(longDistanceCorrEarlyV1','color', cocEarlyV1(3,:)); hold on
boxplot(longDistanceCorrEarlyV1, 'Label', {'Early'})
xlabel('300-400 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,5)
distributionPlot(shortDistanceCorrAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(shortDistanceCorrAdultV1,'Label', {'Adult'})
xlabel('< 100 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,6)
distributionPlot(longDistanceCorrAdultV1','color', cocAdultV1(3,:)); hold on
boxplot(longDistanceCorrAdultV1, 'Label', {'Adult'})
xlabel('300-400 um')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'corrCellsDistanceV1.png'))

% b( Template matching ortho vs. matched
figure
subplot(1,6,1)
distributionPlot(corrTrialsMatchedNaive','color', cocNaive(4,:)); hold on
boxplot(corrTrialsMatchedNaive,'Label', {'Naive'})
xlabel('Matched')
ylim([-1 1])
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(longDistanceCorrNaive','color', cocNaive(3,:)); hold on
boxplot(longDistanceCorrNaive, 'Label', {'Naive'})
xlabel('Ortho')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(corrTrialsMatchedEarly','color', cocEarly(4,:)); hold on
boxplot(corrTrialsMatchedEarly,'Label', {'Early'})
xlabel('Matched')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(longDistanceCorrEarly','color', cocEarly(3,:)); hold on
boxplot(longDistanceCorrEarly, 'Label', {'Early'})
xlabel('Ortho')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,5)
distributionPlot(corrTrialsMatchedAdult','color', cocAdult(4,:)); hold on
boxplot(corrTrialsMatchedAdult,'Label', {'Adult'})
xlabel('Matched')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

subplot(1,6,6)
distributionPlot(longDistanceCorrAdult','color', cocAdult(3,:)); hold on
boxplot(longDistanceCorrAdult, 'Label', {'Adult'})
xlabel('Ortho')
ylim([-1 1])
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'TrialPatternCorr.png'))

%for V1
figure
subplot(1,6,1)
distributionPlot(corrTrialsMatchedNaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(corrTrialsMatchedNaiveV1,'Label', {'Naive'})
xlabel('Matched')
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(corrTrialsOrthoNaiveV1','color', cocNaiveV1(3,:)); hold on
boxplot(corrTrialsOrthoNaiveV1, 'Label', {'Naive'})
xlabel('Ortho')
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(corrTrialsMatchedEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(corrTrialsMatchedEarlyV1,'Label', {'Early'})
xlabel('Matched')
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(corrTrialsOrthoEarlyV1','color', cocEarlyV1(3,:)); hold on
boxplot(corrTrialsOrthoEarlyV1, 'Label', {'Early'})
xlabel('Ortho')
set(gca,'box','off','ycolor','w')

subplot(1,6,5)
distributionPlot(corrTrialsMatchedAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(corrTrialsMatchedAdultV1,'Label', {'Adult'})
xlabel('Matched')
set(gca,'box','off','ycolor','w')

subplot(1,6,6)
distributionPlot(corrTrialsOrthoAdultV1','color', cocAdultV1(3,:)); hold on
boxplot(corrTrialsOrthoAdultV1, 'Label', {'Adult'})
xlabel('Ortho')
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'TrialPatternCorrV1.png'))

%% 3.) Plot results single cells
% a)responsiveness and selectivity portions
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3, 3, 1)
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

subplot(3, 3, 4)
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

subplot(3, 3, 7)
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

subplot(3, 3, 2)
oriNaive = length(find([OSIFitNaive] > 0.2)) ./ length(OSIFitNaive);
% ori_V1 = length(find([OriCircVar_V1] > 0.2)) ./ length(OriCircVar_V1);
nonOriNaive = 1- oriNaive;
h = pie([nonOriNaive oriNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
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
set(hp(1), 'FaceColor', grey);
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
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocAdult(4,:));
title('Adult')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 3)
dirNaive = length(find([DSINaive] > 0.2)) ./length(DSINaive);
nonDirNaive = 1- dirNaive;
h = pie([nonDirNaive dirNaive]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocNaive(4,:));
title('Direction-selective Naive')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 6)
dirEarly = length(find([DSIEarly] > 0.2)) ./length(DSIEarly);
nonDirEarly = 1- dirEarly;
h = pie([nonDirEarly dirEarly]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocEarly(4,:));
title('Early')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 9)
dirAdult= length(find([DSIAdult] > 0.2)) ./length(DSIAdult);
nonDirAdult= 1- dirAdult;
h = pie([nonDirAdult dirAdult]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocAdult(4,:));
title('Adult')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'resp_cells.png'))

%for V1
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3, 3, 1)
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

subplot(3, 3, 4)
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

subplot(3, 3, 7)
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

subplot(3, 3, 2)
oriNaiveV1 = length(find([OSIFitNaiveV1] > 0.2)) ./ length(OSIFitNaiveV1);
nonOriNaiveV1 = 1- oriNaiveV1;
h = pie([nonOriNaiveV1 oriNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocNaiveV1(4,:));
title('Naive V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 5)
oriEarlyV1 = length(find([OSIFitEarlyV1] > 0.2)) ./length(OSIFitEarlyV1);
nonOriEarlyV1 = 1- oriEarlyV1;
h = pie([nonOriEarlyV1 oriEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocEarlyV1(4,:));
title('Early V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 8)
oriAdultV1 = length(find([OSIFitAdultV1] > 0.2)) ./length(OSIFitAdultV1);
nonOriAdultV1 = 1- oriAdultV1;
h = pie([nonOriAdultV1 oriAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocAdultV1(4,:));
title('Adult V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 3)
dirNaiveV1 = length(find([DSINaiveV1] > 0.2)) ./length(DSINaiveV1);
nonDirNaiveV1 = 1- dirNaiveV1;
h = pie([nonDirNaiveV1 dirNaiveV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocNaiveV1(4,:));
title('Direction-selective Naive V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 6)
dirEarlyV1 = length(find([DSIEarlyV1] > 0.2)) ./length(DSIEarlyV1);
nonDirEarlyV1 = 1- dirEarlyV1;
h = pie([nonDirEarlyV1 dirEarlyV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
set(hp(2), 'FaceColor', cocEarlyV1(4,:));
title('Early V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 9)
dirAdultV1= length(find([DSIAdultV1] > 0.2)) ./length(DSIAdultV1);
nonDirAdultV1= 1- dirAdultV1;
h = pie([nonDirAdultV1 dirAdultV1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', grey);
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

%% e) VI
figure
hdl(1) = cdfplot(VINaive); hold all
hdl(2) = cdfplot(VIEarly); hold all
hdl(3) = cdfplot(VIAdult); hold all
set(hdl(1), 'Color', cocNaive(4,:), 'LineWidth', 4);
set(hdl(2), 'Color', cocEarly(4,:), 'LineWidth', 4);
set(hdl(3), 'Color', cocAdult(4,:), 'LineWidth', 4);
grid off;
legend('Naive', 'Early', 'Adult', 'Location', 'SouthEast'); legend('boxoff')
ylabel('Cumulative probability', 'FontSize', 12)
xlabel('Variability index', 'FontSize', 12);

h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', 'w');
set(gcf, 'color', 'w');
set(gca, 'box', 'off')

saveas(gcf, fullfile(save_dir, 'VI.png'))

figure
hdl(1) = cdfplot(VINaiveV1); hold all
hdl(2) = cdfplot(VIEarlyV1); hold all
hdl(3) = cdfplot(VIAdultV1); hold all
set(hdl(1), 'Color', cocNaiveV1(4,:), 'LineWidth', 4);
set(hdl(2), 'Color', cocEarlyV1(4,:), 'LineWidth', 4);
set(hdl(3), 'Color', cocAdultV1(4,:), 'LineWidth', 4);
grid off;
legend('Naive V1', 'Early V1', 'Adult V1', 'Location', 'SouthEast'); legend('boxoff')
ylabel('Cumulative probability', 'FontSize', 12)
xlabel('VI', 'FontSize', 12);

h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', 'w');
set(gcf, 'color', 'w');
set(gca, 'box', 'off')

saveas(gcf, fullfile(save_dir, 'VIV1.png'))

%% f) Fanofactor and Bandthwidth
figure
subplot(1,6,1)
distributionPlot(FanoFactorNaive','color', cocNaive(4,:)); hold on
boxplot(FanoFactorNaive,'Label', {'Naive'})
xlabel('FanoFactor')
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(FanoFactorEarly','color', cocEarly(4,:)); hold on
boxplot(FanoFactorEarly, 'Label', {'Early'})
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(FanoFactorAdult','color', cocAdult(4,:)); hold on
boxplot(FanoFactorAdult, 'Label', {'Adult'})
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(BandthWidthNaive','color', cocNaive(4,:)); hold on
boxplot(BandthWidthNaive, 'Label', {'Naive'})
xlabel('BandthWidth')
set(gca,'Box','off');

subplot(1, 6, 5)
distributionPlot(BandthWidthEarly','color', cocEarly(4,:)); hold on
boxplot(BandthWidthEarly, 'Label', {'Early'})
set(gca,'box','off','ycolor','w')

subplot(1, 6, 6)
distributionPlot(BandthWidthAdult','color', cocAdult(4,:)); hold on
boxplot(BandthWidthAdult, 'Label', {'Adult'})
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'FanoFactorBandthWidth.png'))

figure
subplot(1,6,1)
distributionPlot(FanoFactorNaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(FanoFactorNaiveV1,'Label', {'Naive V1'})
xlabel('FanoFactor')
set(gca,'Box','off');

subplot(1,6,2)
distributionPlot(FanoFactorEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(FanoFactorEarlyV1, 'Label', {'Early V1'})
set(gca,'box','off','ycolor','w')

subplot(1,6,3)
distributionPlot(FanoFactorAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(FanoFactorAdultV1, 'Label', {'Adult V1'})
set(gca,'box','off','ycolor','w')

subplot(1,6,4)
distributionPlot(BandthWidthNaiveV1','color', cocNaiveV1(4,:)); hold on
boxplot(BandthWidthNaiveV1, 'Label', {'Naive V1'})
xlabel('BandthWidth')
set(gca,'Box','off');

subplot(1, 6, 5)
distributionPlot(BandthWidthEarlyV1','color', cocEarlyV1(4,:)); hold on
boxplot(DSIEarlyV1, 'Label', {'Early V1'})
set(gca,'box','off','ycolor','w')

subplot(1, 6, 6)
distributionPlot(BandthWidthAdultV1','color', cocAdultV1(4,:)); hold on
boxplot(BandthWidthAdultV1, 'Label', {'Adult V1'})
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'FanoFactorBandthWidthV1.png'))