function MultiModalAnalysisNew(animal, varargin)

%Load the confocal to STED and confocal to 2p data to combine everthing
%together, plotting the functional characteristics of the inputs, analyzing
%them as a whole and plotting this analysis

%Inputs:
% - animal: which animal are we analysing?
% - cell: if there are multiple cells per animal, which one are we looking
% at?

%STEPS:
% - 0.) Define and if necessary, add the save dir
% - 1.) Load the confocal to STED and confocal to 2p data 
% - 2.) Get the confocal input IDs
% - 3.) Specify inputs and functional groups
% - 4.) Analysis
%   - a) Get synaptic aggregates
% - 5.) Plot
% - 6.) Save

%Outputs:
% - Folder 04_FunctionalInputAnalysis with Plots (see details at the
% beginning of Step 5)
% - 04_FunctionalInputReconstruction.mat with variables:
%   - Spines: contains all the information for all the spines, including
%   whether it is an input, their functional data, anatomical data etc.

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023

if nargin > 1
    cellName = varargin{1};
else
    cellName = [];
end

coc_prop = cbrewer('div', 'RdGy', 15);

%% 0.) %define and if necessary, add the confocal and save dir 
baseDir = 'Z:\Juliane\InputAnalysis\';
cellDir = [baseDir char(animal) filesep char(cellName)];

saveDirAna = [cellDir '\E - Analysis\'];
saveDir = ([saveDirAna filesep '04_FunctionalInputAnalysis']);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%% 1.) Load the confocal to STED and confocal to 2p data 

%conf to STED
load([saveDirAna filesep '03_InputCellReconstruction.mat']);
%conf to 2p
functInfo = load([saveDirAna filesep '02_FunctionalCellReconstruction.mat']);

%% 2.) Get the confocal input IDs
Soma = functInfo.confData.Soma;
Spines = functInfo.confData.ROIs;
Dendrites = functInfo.confData.allDendrites;
Selectors = functInfo.Selectors;
Spines = arrayfun(@(x) setfield(x, 'Input', zeros(1,1)), Spines); %adds in a column with all zeros
pwMeasures = functInfo.pwMeasures;

%get the inputs
for i = 1:length(inputROINr)
    row = inputROINr(i);
    % check the info in the input Data
    Branch = InputData(row).Branch;
    ROINrONBranch = InputData(row).ROINrOnBranch;
    
    % now find all the ones that match in the Spine data
    sameBranch = find([Spines.Branch] == Branch);
    sameROINr = find([Spines.ROINrOnBranch] == ROINrONBranch);
    spineRow = intersect(sameBranch, sameROINr);
    
    %set that input to 1
    Spines(spineRow).Input = 1;
end

% add the info to the pW Measurements
pwMeasures = arrayfun(@(x) setfield(x, 'inputPair', []), pwMeasures);
for r = 1:length(pwMeasures)
    %need to find the pairs & add them to 
    ROIsOnDend = find([Spines.Dendrite] == pwMeasures(r).Dendrite);
    spineAName = find(cellfun(@(x) strcmp(x, pwMeasures(r).spineA),{Spines.Nr}));
    spineBName = find(cellfun(@(x) strcmp(x, pwMeasures(r).spineB),{Spines.Nr}));
    
    %now combine them
    pwMeasures(r).inputPair = Spines(intersect(ROIsOnDend,spineAName)).Input + Spines(intersect(ROIsOnDend,spineBName)).Input;
end

%add the STED Tracing to the confocal data
%which spines were done with STED?
Spines = arrayfun(@(x) setfield(x, 'STEDdone', zeros(1,1)), Spines); %adds in a column with all zeros
for s = 1:length(Spines)
    Spines(s).STEDdone = InputData(s).STEDdone;
end

%what parts of the dendrite were done with STED and what is the outcome of the analysis?
Dendrites = arrayfun(@(x) setfield(x, 'STEDDone', []), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'STEDTrace', []), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'STEDPerc', []), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'STEDLength', []), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'numInputs', []), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'InputDensity', []), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'percOfAllSpines', []), Dendrites);
for d = 1:length(Dendrites)
    Dendrites(d).STEDDone = STEDTracing(d).STEDdone;
    Dendrites(d).STEDTrace = STEDTracing(d).STEDTrace;
    Dendrites(d).STEDPerc = STEDTracing(d).STEDPerc;
    Dendrites(d).STEDLength = STEDTracing(d).STEDLength;
    Dendrites(d).numInputs = STEDTracing(d).numInputs;
    Dendrites(d).InputDensity = STEDTracing(d).InputDensity;
    Dendrites(d).percOfAllSpines = STEDTracing(d).percOfAllSpines;
end

%add ROIs to the dendrite
Dendrites = arrayfun(@(x) setfield(x, 'ROIs', []), Dendrites);
for d = 1:length(Dendrites)
    ROIsOnDendrite = find([Spines.Dendrite] == d);
    if ~isempty(ROIsOnDendrite)
        Dendrites(d).ROIs = Spines(ROIsOnDendrite);
    end
end

clear InputData
clear functInfo 

%% 3.) Specify inputs and functional groups

%for Spines
%load old selectors
%logical
TwoPROIs = Selectors.TwoPROIs; %all invivo spines
respROIs = Selectors.respROIs; %all resp spines
oriSelectROIs = Selectors.oriSelectROIs; %all resp & ori select spines
dirSelectROIs = Selectors.dirSelectROIs; %all resp & dir select spines

%ROI NRs
TwoPROIsNR = Selectors.TwoPROIsNR; %all invivo spines
respROINrs = Selectors.respROINrs; %all resp spines
oriSelectROIsNr = Selectors.oriSelectROIsNr; %all resp & ori select spines
dirSelectROIsNr = Selectors.dirSelectROIsNr; %all resp & dir select spines

%add selectors for inputs and no inputs
%logical
TwoPInputLog = logical([Spines.Input] == 1) & TwoPROIs; %all twoP characterized inputs
respInputLog = logical([Spines.Input] == 1) & logical([Spines.good] == 1);%all resp inputs
respNoInputLog =  logical([Spines.Input] == 0) & logical([Spines.good] == 1); %all resp non-inputs
oriSelectInputLog = logical([Spines.Input] == 1) & oriSelectROIs; %all ori-select inputs
oriSelectNoInputLog =  logical([Spines.Input] == 0) & oriSelectROIs; %all ori-select non-inputs
dirSelectInputLog = logical([Spines.Input] == 1) & dirSelectROIs; %all dir-select inputs
dirSelectNoInputLog =  logical([Spines.Input] == 0) & dirSelectROIs; %all dir-select non-inputs

%ROI NRs
TwoPInputsNr = intersect(find([Spines.Input] == 1), TwoPROIsNR); %all twoP characterized inputs
respInputsNr = intersect(find([Spines.Input] == 1), respROINrs); %all resp inputs
respNoInputsNr = intersect(find([Spines.Input] == 0), respROINrs); %all resp non-inputs
oriSelectInputNr = intersect(find([Spines.Input] == 1), oriSelectROIsNr); %all ori-sel inputs
oriSelectNoInputNr = intersect(find([Spines.Input] == 0), oriSelectROIsNr); %all ori-sel inputs
dirSelectInputNr = intersect(find([Spines.Input] == 1), dirSelectROIsNr); %all dir-sel inputs
dirSelectNoInputNr = intersect(find([Spines.Input] == 0), dirSelectROIsNr); %all dir-sel inputs

%for branches
% get the input branches
inputBranchesNrs = unique([Spines(inputROINr).Dendrite]); %find all the ones that have at least one input
%now convert it to the rows in case of dendrNr is not row of dendrites
inputBranches = [];
for inB = 1:length(inputBranchesNrs)
    inputBranches(inB) = find([Dendrites.dendNr]  == inputBranchesNrs(inB));
end


%% 4.) Analysis
% a) Get synaptic aggregates
%get angles
angles = linspace(0,180,9); angles = angles(1:end-1);
anglesAll = linspace(0,360,17); anglesAll = anglesAll(1:end-1);

%get soma response
SomaResp = (Soma.meanResp(1:8)+Soma.meanResp(9:16))/2; %fold in ori space
SomaRespNorm = (SomaResp)/(max(SomaResp));%normalize
SomaDirNorm = Soma.meanResp/(max(Soma.meanResp));

%in orientation space
[PrefAngleFunct, curveSumFunc, eventSumFunc] = calcSumSpineResponses(Spines, oriSelectROIsNr);
if ~isempty(oriSelectInputNr)
    [PrefAngleInput, curveSumInput, eventSumInput] = calcSumSpineResponses(Spines, oriSelectInputNr);
else
    curveSumInput = []; 
    eventSumInput = [];
end
[PrefAngleNoInput, curveSumNoInput, eventSumNoInput] = calcSumSpineResponses(Spines, oriSelectNoInputNr);

%in direction space
[~, curveSumDirFunc, eventSumDirFunc] = calcSumSpineDirResponses(Spines, oriSelectROIsNr);
if ~isempty(oriSelectInputNr)
    [~, curveSumDirInput, eventSumDirInput] = calcSumSpineDirResponses(Spines, oriSelectInputNr);
else
    curveSumDirInput = []; 
    eventSumDirInput = [];
end
[~, curveSumDirNoInput, eventSumDirNoInput] = calcSumSpineDirResponses(Spines, oriSelectNoInputNr);

%shift to prefAngle
[~, prefInd] = max(SomaRespNorm);
anglesAligned = linspace(-90,90,9);
SomaRespNormShifted = shiftToPrefAngleSoma(prefInd, SomaRespNorm);
curveSumFuncShift = shiftToPrefAngleSoma(prefInd, curveSumFunc');
if ~isempty(curveSumInput)
    curveSumInputShift = shiftToPrefAngleSoma(prefInd, curveSumInput');
else
    curveSumInputShift = [];
end
curveSumNoInputShift = shiftToPrefAngleSoma(prefInd, curveSumNoInput');
eventSumFuncShift = shiftToPrefAngleSoma(prefInd, eventSumFunc');
if ~isempty(eventSumInput)
    eventSumInputShift = shiftToPrefAngleSoma(prefInd, eventSumInput');
else
    eventSumInputShift = [];
end
eventSumNoInputShift = shiftToPrefAngleSoma(prefInd, eventSumNoInput');

%% 5.) Plot
%Figure Overview
%10s - ROIs on top of tracing
%Fig 10: Which functional ROIs were inputs?
%Fig 11: Which responsive ROIs were inputs?

%20s - What are the properties of the inputs on top of tracing?
%Fig 20: What is the ori pref of the inputs?
%Fig 21: What is the dir pref of the inputs?
%Fig 22: What is the OSI of the inputs?
%Fig 23: What is the delta Ori of the inputs?
%Fig 24: What is the DSI of the inputs?
%Fig 25: What is the bandwidth of the spines?

%30s - What is the distribution of the inputs?
%Fig 30: Distribution of prefOri as histogram
%Fig 31: Distribution of prefOri as polarPlotOri
%Fig 32/320: Boxplot plot & distribtuion plot of OSI
%Fig 33: Distribution of prefDir as histogram
%Fig 34: Distribution of prefDir as polarPlot
%Fig 35/350: Boxplot plot & distribution plot of DSI
%Fig 36/360: Boxplot plot & distribution plot of DSIvect
%Fig 37: Boxplot plot & distribution plot of bandwidth

%40s - Synaptic Aggregats
%Fig 40: Synaptic aggregate vs. cell in orientation space based on curve
%Fig 41: Synaptic aggregate vs. cell in orientation space based on event
%Fig 42: Synaptic aggregate vs. cell in direction space based on curve
%Fig 43: Synaptic aggregate vs. cell in direction space based on event
 
%50s & 60s: Local environment ori(50s) and dir (60s)
%Fig 50/60: local (dir) Dispersion
%Fig 51/61: HI (dir)
%Fig 52/62: local (dir) dispersion vs. HI (dir)
%Fig 53/63: local deltaOri/Dir
%Fig 54/64: local delteOriSoma/deltaDirSOma
%Fig 55/65: nearest SimilarPrefOri/dir
%Fig 56/66: local OSI/DSI

%70s: Properties of dendritic segments
%Fig 70: Branch circular dispersion ori/dir
%Fig 71: Branch selectivity
%Fig 72: delta Ori/delta Dir
%Fig 73: Branch circular dispersion vs. deltaOri
%Fig 74: Branch circular dispersion dir vs. deltaDir

%80s: Pairwise distances
%Fig 80: Pref Ori
%Fig 81: Pref Dir
%Fig 82: OSI
%Fig 83: DSI
%Fig 84: DSIVect

%% Start Plotting
%--------------------------------------------------------------------------
%10s: Where are the ROIs
%Which ROIs do we have? 
%Fig 10: Which functional ROIs were inputs?
plotIDsOnCellTracing(Dendrites, Spines, TwoPInputLog, 10, 'blue', 0)
saveas(gcf, fullfile(saveDir, '10_CellTracing_FunctionalInputs.png'))
%Fig 11: Which responsive ROIs were inputs?
plotIDsOnCellTracing(Dendrites, Spines, respInputLog, 11, 'blue', 0)
saveas(gcf, fullfile(saveDir, '11_CellTracing_ResponsiveInputs.png'))
%Fig 12: Include on STED tracing: which functional inputs and which one of
%those were responsive and ori/selective
plotFunctInputsOnCellTracing(Dendrites, Spines, TwoPInputLog,12)
saveas(gcf, fullfile(saveDir, '12_CellTracing_TypeOfMatchedInputs.png'))


%--------------------------------------------------------------------------
%20s - What are the properties of the inputs on top of tracing?
%What is the ori pref, delta ori, OSI

%Fig 20: What is the ori pref of the inputs?
plotFieldsOnCellTracing(Dendrites, Spines, respInputLog, Soma, 'prefOri', 1, hsv(180), 20, 0)
saveas(gcf, fullfile(saveDir, '20_CellTracing_FunctionalInputs.png'))
plotCompareFieldsOnSTEDTracing(Dendrites, Spines, oriSelectROIs, oriSelectInputLog, Soma, 'prefOri', 1, hsv(180), 200, 0)
saveas(gcf, fullfile(saveDir, '200_CellTracing_FunctionalInputs.png'))

%Fig 21: What is the dir pref of the inputs?
plotFieldsOnCellTracing(Dendrites, Spines, respInputLog, Soma, 'prefDir', 1, hsv(360), 21, 0)
saveas(gcf, fullfile(saveDir, '21_CellTracing_prefDir.png'))
plotCompareFieldsOnSTEDTracing(Dendrites, Spines, dirSelectROIs, dirSelectInputLog, Soma, 'prefOri', 1, hsv(180), 210, 0)
saveas(gcf, fullfile(saveDir, '210_CellTracing_FunctionalInputs.png'))

%Fig 22: What is the OSI of the inputs?
plotFieldsOnCellTracing(Dendrites, Spines, respInputLog, Soma,  'OSI', 0, inferno(100), 22, 0)
saveas(gcf, fullfile(saveDir, '22_CellTracing_OSI.png'))
plotCompareFieldsOnSTEDTracing(Dendrites, Spines, oriSelectROIs, oriSelectInputLog, Soma, 'prefOri', 1, hsv(180), 220, 0)
saveas(gcf, fullfile(saveDir, '220_CellTracing_FunctionalInputs.png'))

%Fig 23: What is the delta Ori of the inputs?
plotFieldsOnCellTracing(Dendrites, Spines, respInputLog, Soma,  'deltaOri', 0, inferno(9000), 23, 0)
saveas(gcf, fullfile(saveDir, '23_CellTracing_deltaOri.png'))
plotCompareFieldsOnSTEDTracing(Dendrites, Spines, oriSelectROIs, oriSelectInputLog, Soma, 'prefOri', 1, hsv(180), 230, 0)
saveas(gcf, fullfile(saveDir, '230_CellTracing_FunctionalInputs.png'))

%Fig 24: What is the DSI of the inputs?
plotFieldsOnCellTracing(Dendrites, Spines, respInputLog, Soma,  'DSI', 0, viridis(100), 24, 0)
saveas(gcf, fullfile(saveDir, '24_CellTracing_DSI.png'))
plotCompareFieldsOnSTEDTracing(Dendrites, Spines, dirSelectROIs, dirSelectInputLog, Soma, 'prefDir', 1, hsv(360), 240, 0)
saveas(gcf, fullfile(saveDir, '240_CellTracing_FunctionalInputs.png'))

%Fig 25: What is the bandwidth of the spines?
plotFieldsOnCellTracing(Dendrites, Spines, respInputLog, Soma,  'Bandwidth', 0, jet(9000), 25, 0)
saveas(gcf, fullfile(saveDir, '25_CellTracing_Bandwidth.png'))
plotCompareFieldsOnSTEDTracing(Dendrites, Spines, oriSelectROIs, oriSelectInputLog, Soma, 'Bandwidth', 1, hsv(90), 250, 0)
saveas(gcf, fullfile(saveDir, '250_CellTracing_FunctionalInputs.png'))



%--------------------------------------------------------------------------
%30s - What is the distribution of the inputs?
%what is the distribution of ori prefs? 

%Fig 30: Distribution of prefOri as histogram
[counts, ~] = histcounts([Spines(oriSelectROIsNr).prefOri], [0 linspace(22.5/2,180-22.5/2,8) 180]); %do counts that have the oris as center of the bins
countsAdjusted = [counts(1)+counts(end) counts(2:end-1) counts(1)+counts(end)];
countsNormalized = countsAdjusted/max(countsAdjusted);

[countsInputs, ~] = histcounts([Spines(oriSelectInputNr).prefOri], [0 linspace(22.5/2,180-22.5/2,8) 180]); %do counts that have the oris as center of the bins
countsAdjustedInputs = [countsInputs(1)+countsInputs(end) countsInputs(2:end-1) countsInputs(1)+countsInputs(end)];
countsNormalizedInputs = countsAdjustedInputs/max(countsAdjustedInputs);

figure(30)
plot(angles, SomaRespNorm, 'Color', 'black', 'LineWidth', 5);
hold on
plot(angles, countsNormalized(1:end-1), 'Color', coc_prop(13,:), 'LineWidth', 2);
hold on
plot(angles, countsNormalizedInputs(1:end-1), 'Color', coc_prop(3,:), 'LineWidth', 2);
set(gcf, 'color', 'w');
xlabel('Pref ori in deg')
ylabel('Probability')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '30_PrefOriDistribution.png'))

%Fig 31: Distribution of prefOri as polarPlotOri
figure(31)
polarPlotOri([angles 180], [SomaRespNorm SomaRespNorm(1)], [], 'black')
hold on
polarPlotOri([angles 180], countsNormalized, [], coc_prop(13,:), '-',0)
hold on
if sum(~isnan(countsNormalizedInputs)) > 0
    polarPlotOri([angles 180], countsNormalizedInputs, [], coc_prop(3,:), '-',0)
end
saveas(gcf, fullfile(saveDir, '31_PrefOriAll_PolarPlot.png'))

%Fig 32: Distribution plot of OSI
if ~isempty(respInputsNr)
    figure(32)
    subplot(1,2,1)
    distributionPlot([Spines(respROINrs).OSI]', 'color', coc_prop(13,:))
    ylim([0 1])
    ylabel('Orientation selectivity index')
    xticklabels('all spines')
    subplot(1,2,2)
    distributionPlot([Spines(respInputsNr).OSI]', 'color', coc_prop(3,:))
    ylim([0 1])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'32_OSIDistribution.png'))

    %Fig 320: Boxplot plot of OSI
    figure(320)
    subplot(1,2,1)
    boxplot([Spines(respROINrs).OSI], 'color',coc_prop(15,:))
    ylim([0 1])
    hold on
    scatter(repmat(1:1,size([Spines(respROINrs).OSI],1),1), [Spines(respROINrs).OSI], 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Orientation selectivity index')
    xticklabels('all spines')
    box off
    subplot(1,2,2)
    boxplot([Spines(respInputsNr).OSI], 'color',coc_prop(3,:))
    hold on
    scatter(repmat(1:1,size([Spines(respInputsNr).OSI],1),1), [Spines(respInputsNr).OSI], 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylim([0 1])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'320_OSIDataPoints.png'))

    %Fig 33: Distribution of prefDir as histogram
    %test with both respROIs and dirSelectROIs!
    [countsDir, ~] = histcounts([Spines(dirSelectROIsNr).prefDir], [0 linspace(22.5/2,360-22.5/2,16) 360]); %do counts that have the oris as center of the bins
    countsAdjustedDir = [countsDir(1)+countsDir(end) countsDir(2:end-1) countsDir(1)+countsDir(end)];
    countsNormalizedDir = countsAdjustedDir/max(countsAdjustedDir);

    [countsDirInput, ~] = histcounts([Spines(dirSelectInputNr).prefDir], [0 linspace(22.5/2,360-22.5/2,16) 360]); %do counts that have the oris as center of the bins
    countsAdjustedDirInput = [countsDirInput(1)+countsDirInput(end) countsDirInput(2:end-1) countsDirInput(1)+countsDirInput(end)];
    countsNormalizedDirInput = countsAdjustedDirInput/sum(countsAdjustedDirInput);

    countsNormalizedDirRand = NaN(17,1000);

    for r = 1:1000
        newOrder = randperm(length(dirSelectROIsNr));
        newRandInputs = newOrder(1:length(dirSelectInputNr));
        [countsDirRand, ~] = histcounts([Spines(newRandInputs).prefDir], [0 linspace(22.5/2,360-22.5/2,16) 360]); %do counts that have the oris as center of the bins
        countsAdjustedDirRand = [countsDir(1)+countsDir(end) countsDir(2:end-1) countsDir(1)+countsDir(end)];
        countsNormalizedDirRand(:,r) = countsAdjustedDirRand/max(countsAdjustedDirRand);
    end
    avgRandPrefDir = mean(countsNormalizedDirRand,2)/(max(mean(countsNormalizedDirRand,2)));


    figure(33)
    plot(anglesAll, SomaDirNorm, 'Color', 'black', 'LineWidth', 5);
    hold on
    plot(anglesAll, countsNormalizedDir(1:end-1), 'Color', coc_prop(13,:), 'LineWidth', 2);
    hold on
    plot(anglesAll, countsNormalizedDirInput(1:end-1), 'Color', coc_prop(3,:), 'LineWidth', 2);
    hold on
    plot(anglesAll, avgRandPrefDir(1:end-1), 'Color', 'blue', 'LineWidth', 2);
    set(gcf, 'color', 'w');
    xlabel('Pref dir in deg')
    ylabel('Probability')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '33_PrefDirAll.png'))

    %Fig 34: Distribution of prefDir as polarPlot
    figure(34)
    polarPlot([anglesAll 360], [SomaDirNorm SomaDirNorm(1)], [], 'black')
    hold on
    polarPlot([anglesAll 360], countsNormalizedDir, [], coc_prop(13,:), '-',0)
    hold on
    polarPlot([anglesAll 360], countsNormalizedDirInput, [], coc_prop(3,:), '-',0)
    hold on
    %polarPlot([anglesAll 360], avgRandPrefDir, [], 'blue')
    saveas(gcf, fullfile(saveDir, '34_PrefDirAll_PolarPlot.png'))

    %Fig 35: Distribution plot of DSI
    figure(35)
    subplot(1,2,1)
    distributionPlot([Spines(respROINrs).DSI]', 'color', coc_prop(13,:))
    ylim([0 1])
    ylabel('Direction selectivity index')
    xticklabels('all spines')
    subplot(1,2,2)
    distributionPlot([Spines(respInputsNr).DSI]', 'color', coc_prop(3,:))
    ylim([0 1])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'35_DSIDistribution.png'))

    %Fig 350: Boxplot plot of DSI
    figure(350)
    subplot(1,2,1)
    boxplot([Spines(respROINrs).DSI], 'color',coc_prop(15,:))
    ylim([0 1])
    hold on
    scatter(repmat(1:1,size([Spines(respROINrs).DSI],1),1), [Spines(respROINrs).DSI], 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Direction selectivity index')
    xticklabels('all spines')
    box off
    subplot(1,2,2)
    boxplot([Spines(respInputsNr).DSI], 'color',coc_prop(3,:))
    hold on
    scatter(repmat(1:1,size([Spines(respInputsNr).DSI],1),1), [Spines(respInputsNr).DSI], 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylim([0 1])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'350_DSIDataPoints.png'))

    %Fig 36: Distribution plot of DSIvect
    figure(36)
    subplot(1,2,1)
    distributionPlot([Spines(respROINrs).DSIvect]', 'color', coc_prop(13,:))
    ylim([0 1])
    ylabel('Direction selectivity index (vect)')
    xticklabels('all spines')
    subplot(1,2,2)
    distributionPlot([Spines(respInputsNr).DSIvect]', 'color', coc_prop(3,:))
    ylim([0 1])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'36_DSIVectDistribution.png'))

    %Fig 360: Boxplot plot of DSIvect
    figure(360)
    subplot(1,2,1)
    boxplot([Spines(respROINrs).DSIvect], 'color',coc_prop(15,:))
    ylim([0 1])
    hold on
    scatter(repmat(1:1,size([Spines(respROINrs).DSIvect],1),1), [Spines(respROINrs).DSIvect], 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Direction selectivity index (vect)')
    xticklabels('all spines')
    box off
    subplot(1,2,2)
    boxplot([Spines(respInputsNr).DSIvect], 'color',coc_prop(3,:))
    hold on
    scatter(repmat(1:1,size([Spines(respInputsNr).DSIvect],1),1), [Spines(respInputsNr).DSIvect], 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylim([0 1])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'360_DSIvectDataPoints.png'))

    %Fig 37: Boxplot plot & distribution plot of bandwidth
    figure(37) 
    subplot(1,2,1)
    distributionPlot([Spines(respROINrs).Bandwidth]', 'color', coc_prop(13,:))
    ylim([0 90])
    ylabel('Bandwidth')
    xticklabels('all spines')
    subplot(1,2,2)
    distributionPlot([Spines(respInputsNr).Bandwidth]', 'color', coc_prop(3,:))
    ylim([0 90])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'37_Bandwidth.png'))

    %Fig 370: Boxplot plot of Bandwidth
    figure(370)
    subplot(1,2,1)
    boxplot([Spines(respROINrs).Bandwidth], 'color',coc_prop(15,:))
    ylim([0 90])
    hold on
    scatter(repmat(1:1,size([Spines(respROINrs).Bandwidth],1),1), [Spines(respROINrs).Bandwidth], 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Bandwidth')
    xticklabels('all spines')
    box off
    subplot(1,2,2)
    boxplot([Spines(respInputsNr).Bandwidth], 'color',coc_prop(3,:))
    hold on
    scatter(repmat(1:1,size([Spines(respInputsNr).Bandwidth],1),1), [Spines(respInputsNr).Bandwidth], 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylim([0 90])
    xticklabels('inputs')
    ax = gca;
    ax.YColor = 'none';
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'370_BandwidthDataPoints.png'))


    %--------------------------------------------------------------------------
    %40s:What is the aggregate of the responses
    %What is the summed responses
    %Fig 40: Synaptic aggregate vs. cell in orientation space based on curve
    figure(40)
    plot(anglesAligned, SomaRespNormShifted, 'Color', 'black', 'LineWidth', 5);
    hold on
    if ~isempty(curveSumInputShift)
        plot(anglesAligned, curveSumInputShift, 'Color', coc_prop(3,:), 'LineWidth', 2);
        hold on 
    end
    plot(anglesAligned, curveSumNoInputShift, 'Color', coc_prop(10,:), 'LineWidth', 2);
    hold on 
    plot(anglesAligned, curveSumFuncShift, 'Color', coc_prop(13,:), 'LineWidth', 2);
    xlabel('Angles in deg')
    xticks([-90 -45 0 45 90])
    ylabel('Normalized events')
    yticks([0 0.25 0.5 0.75 1])
    ylim([0 1]);
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'40_SummedInputsAligned.png'))

    %Fig 41: Synaptic aggregate vs. cell in orientation space based on event
    figure(41)
    plot(anglesAligned, SomaRespNormShifted, 'Color', 'black', 'LineWidth', 5);
    hold on
    if ~isempty(eventSumInputShift)
        plot(anglesAligned, eventSumInputShift, 'Color', coc_prop(3,:), 'LineWidth', 2);
        hold on 
    end
    plot(anglesAligned, eventSumFuncShift, 'Color', coc_prop(13,:), 'LineWidth', 2);
    xlabel('Angles in deg')
    xticks([-90 -45 0 45 90])
    ylabel('Normalized events')
    yticks([0 0.25 0.5 0.75 1])
    ylim([0 1]);
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'41_SummedEventsAligned.png'))

    %Fig 42: Synaptic aggregate vs. cell in direction space based on curve
    figure(42)
    polarPlot(anglesAll(1:end-1), SomaDirNorm, 1, [0 0 0], ':',1)
    polarPlot(anglesAll(1:end-1), curveSumDirFunc, 1, coc_prop(13,:), '-',0)
    if sum(~isnan(curveSumDirInput)) > 0
        polarPlot(anglesAll(1:end-1), curveSumDirInput, 1, coc_prop(3,:), '-',0)
    end
    saveas(gcf, fullfile(saveDir,'42_SummedInputsDir.png'))

    %Fig 43: Synaptic aggregate vs. cell in direction space based on event
    figure(43)
    polarPlot(anglesAll(1:end-1), SomaDirNorm, 1, [0 0 0], ':',1)
    polarPlot(anglesAll(1:end-1), eventSumDirFunc, 1, coc_prop(13,:), '-',0)
    if sum(~isnan(eventSumDirInput)) > 0
        polarPlot(anglesAll(1:end-1), eventSumDirInput, 1, coc_prop(3,:), '-',0)
    end
    saveas(gcf, fullfile(saveDir,'43_SummedEventsDir.png'))

    %--------------------------------------------------------------------------
    %50s & 60s: Local environment ori(50s) and dir (60s)

    getDeltaOri = @(x) x.funcData.deltaOri; %define this to be able to get the deltaOriValues
    deltaOriValues = cell2mat(arrayfun(getDeltaOri, Spines(oriSelectROIsNr), 'UniformOutput', false))';
    deltaOriInputValues = cell2mat(arrayfun(getDeltaOri, Spines(oriSelectInputLog), 'UniformOutput', false))';

    %Fig 50/60: local (dir) Dispersion
    localDisp_values = cell2mat(cellfun(@(x) x.localDispersion, num2cell(Spines(oriSelectROIsNr)), 'UniformOutput', false));
    localDisp_values = reshape(localDisp_values, 4, [])';
    localDispInput_values = cell2mat(cellfun(@(x) x.localDispersion, num2cell(Spines(oriSelectInputLog)), 'UniformOutput', false));
    localDispInput_values = reshape(localDispInput_values, 4, [])';
    localDirDisp_values = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(Spines(dirSelectROIsNr)), 'UniformOutput', false));
    localDirDisp_values = reshape(localDirDisp_values, 4, [])';
    localDirDispInput_values = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(Spines(dirSelectInputLog)), 'UniformOutput', false));
    localDirDispInput_values = reshape(localDirDispInput_values, 4, [])';

    figure(50)
    subplot(2,2,1)
    boxplot(localDisp_values,'color',coc_prop(15,:))
    ylim([0 45])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDisp_values,1),1), localDisp_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Local dispersion')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    boxplot(localDispInput_values,'color',coc_prop(3,:))
    ylim([0 45])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDispInput_values,1),1), localDispInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Local dispersion inputs')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,2)
    plot(deltaOriValues,localDisp_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 90])
    ylim([0 45])
    ylabel(sprintf('Local dispersion within \n(5 \\mum) of spine'))
    xlabel('delta Ori spine')
    box off
    subplot(2,2,4)
    plot(deltaOriInputValues,localDispInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 90])
    ylim([0 45])
    ylabel(sprintf('Local dispersion within \n(5 \\mum) of input'))
    xlabel('delta Ori inputs')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'50_Local dispersion inputs vs. Spines.png'))

    figure(60)
    subplot(2,2,1)
    boxplot(localDirDisp_values,'color',coc_prop(15,:))
    ylim([0 90])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDirDisp_values,1),1), localDirDisp_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Local dir dispersion spines')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    boxplot(localDirDispInput_values,'color',coc_prop(3,:))
    ylim([0 90])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDirDispInput_values,1),1), localDirDispInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Local dir dispersion inputs')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,2)
    plot([Spines(dirSelectROIsNr).prefDir],localDirDisp_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 360])
    ylim([0 90])
    ylabel(sprintf('Local dir dispersion \nwithin (5 \\mum) of spine'))
    xlabel('pref Dir spine')
    box off
    subplot(2,2,4)
    plot([Spines(dirSelectInputNr).prefDir],localDirDispInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 360])
    ylim([0 90])
    ylabel(sprintf('Local dir dispersion \nwithin (5 \\mum) of input'))
    xlabel('pref Dir inputs')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'60_Local dir dispersion inputs vs. Spines.png'))

    %Fig 51/61: HI (dir)
    HI_values = cell2mat(cellfun(@(x) x.HI, num2cell(Spines(oriSelectROIsNr)), 'UniformOutput', false));
    HI_values = reshape(HI_values, 4, [])';
    HIInput_values = cell2mat(cellfun(@(x) x.HI, num2cell(Spines(oriSelectInputNr)), 'UniformOutput', false));
    HIInput_values = reshape(HIInput_values, 4, [])';
    HIDir_values = cell2mat(cellfun(@(x) x.HIDir, num2cell(Spines(dirSelectROIsNr)), 'UniformOutput', false));
    HIDir_values = reshape(HIDir_values, 4, [])';
    HIDirInput_values = cell2mat(cellfun(@(x) x.HIDir, num2cell(Spines(dirSelectInputNr)), 'UniformOutput', false));
    HIDirInput_values = reshape(HIDirInput_values, 4, [])';

    figure(51)
    subplot(2,2,1)
    boxplot(HI_values,'color',coc_prop(15,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(HI_values,1),1), HI_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Homeogeneity index')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,2)
    plot(deltaOriValues,HI_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('HI of local environment \n(5 \\mum) of spine'))
    xlabel('delta Ori spine')
    box off
    subplot(2,2,3)
    boxplot(HIInput_values,'color',coc_prop(3,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(HIInput_values,1),1), HIInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Homeogeneity index')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot(deltaOriInputValues,HIInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('HI of local environment \n(5 \\mum) of input'))
    xlabel('delta Ori spine')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'51_Local HI to Spine vs Input.png'))

    figure(61)
    subplot(2,2,1)
    boxplot(HIDir_values,'color',coc_prop(15,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(HIDir_values,1),1), HIDir_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Local Homeogeneity index (dir)')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,2)
    plot([Spines(dirSelectROIsNr).prefDir],HIDir_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 360])
    ylim([0 1])
    ylabel(sprintf('HIdir of local environment \n(5 \\mum) of spine'))
    xlabel('Pref dir spine')
    box off
    subplot(2,2,3)
    boxplot(HIDirInput_values,'color',coc_prop(3,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(HIDirInput_values,1),1), HIDirInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('Local Homeogeneity index (dir)')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot([Spines(dirSelectInputNr).prefDir],HIDirInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 360])
    ylim([0 1])
    ylabel(sprintf('HIdir of local environment \n(5 \\mum) of input'))
    xlabel('pref Dir input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'61_Local HIdir to Spine vs Input.png'))

    %Fig 52/62: local (dir) dispersion vs. HI (dir)
    figure(52)
    subplot(1,2,1)
    plot(localDisp_values(:,2),HI_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 45])
    ylim([0 1])
    ylabel(sprintf('HI of local environment \n(5 \\mum) of spine'))
    xlabel(sprintf('circular Dispersion of local \nenvironment (5 \\mum) of spine'))
    box off
    subplot(1,2,2)
    plot(localDispInput_values(:,2),HIInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 45])
    ylim([0 1])
    ylabel(sprintf('HI of local environment \n(5 \\mum) of input'))
    xlabel(sprintf('circular Dispersion of local \nenvironment (5 \\mum) of input'))
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'52_Local HI vs local Dispersion spine vs. input.png'))

    figure(62)
    subplot(1,2,1)
    plot(localDirDisp_values(:,2),HIDir_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('HI of local environment \n(5 \\mum) of spine'))
    xlabel(sprintf('circular Dispersion of local \nenvironment (5 \\mum) of spine'))
    box off
    subplot(1,2,2)
    plot(localDirDispInput_values(:,2),HIDirInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('HI of local environment \n(5 \\mum) of input'))
    xlabel(sprintf('circular Dispersion of local \nenvironment (5 \\mum) of input'))
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'62_Local HIdir vs local dir Dispersion spine vs. input.png'))

    %Fig 53/63: local deltaOri/Dir
    localDeltaOri_values = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(Spines(oriSelectROIsNr)), 'UniformOutput', false));
    localDeltaOri_values = reshape(localDeltaOri_values, 4, [])';
    localDeltaOriInput_values = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(Spines(oriSelectInputNr)), 'UniformOutput', false));
    localDeltaOriInput_values = reshape(localDeltaOriInput_values, 4, [])';
    localDeltaDir_values = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(Spines(dirSelectROIsNr)), 'UniformOutput', false));
    localDeltaDir_values = reshape(localDeltaDir_values, 4, [])';
    localDeltaDirInput_values = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(Spines(dirSelectInputNr)), 'UniformOutput', false));
    localDeltaDirInput_values = reshape(localDeltaDirInput_values, 4, [])';

    figure(53)
    subplot(2,2,1)
    boxplot(localDeltaOri_values,'color',coc_prop(15,:))
    ylim([0 90])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaOri_values,1),1), localDeltaOri_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaOri of local environment to spine')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    plot([Spines(oriSelectROIsNr).prefOri],localDeltaOri_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 180])
    ylim([0 90])
    ylabel(sprintf('deltaOri of local environment \n(5 \\mum) to spine'))
    xlabel('pref Ori spine')
    box off
    subplot(2,2,2)
    boxplot(localDeltaOriInput_values,'color',coc_prop(3,:))
    ylim([0 90])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaOriInput_values,1),1), localDeltaOriInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaOri of local environment to input')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot([Spines(oriSelectInputLog).prefOri],localDeltaOriInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 180])
    ylim([0 90])
    ylabel(sprintf('deltaOri of local environment \n(5 \\mum) to input'))
    xlabel('pref Ori input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'53_Local deltaOri to Spine vs input.png'))

    figure(63)
    subplot(2,2,1)
    boxplot(localDeltaDir_values,'color',coc_prop(15,:))
    ylim([0 180])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaDir_values,1),1), localDeltaDir_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaDir of local environment to spine')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    plot([Spines(dirSelectROIsNr).prefDir],localDeltaDir_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 360])
    ylim([0 180])
    ylabel(sprintf('deltaDir of local environment \n(5 \\mum) to spine'))
    xlabel('pref Ori spine')
    box off
    subplot(2,2,2)
    boxplot(localDeltaDirInput_values,'color',coc_prop(3,:))
    ylim([0 180])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaDirInput_values,1),1), localDeltaDirInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaDir of local environment to input')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot([Spines(dirSelectInputLog).prefDir],localDeltaDirInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 360])
    ylim([0 180])
    ylabel(sprintf('deltaDir of local environment \n(5 \\mum) to input'))
    xlabel('pref Dir input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'63_Local deltaDir to Spine vs input.png'))

    %Fig 54/64: local deltaOriSoma/deltaDirSoma
    localDeltaOriSoma_values = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(Spines(oriSelectROIsNr)), 'UniformOutput', false));
    localDeltaOriSoma_values = reshape(localDeltaOriSoma_values, 4, [])';
    localDeltaOriSomaInput_values = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(Spines(oriSelectInputNr)), 'UniformOutput', false));
    localDeltaOriSomaInput_values = reshape(localDeltaOriSomaInput_values, 4, [])';
    localDeltaDirSoma_values = cell2mat(cellfun(@(x) x.localDeltaDirSoma, num2cell(Spines(dirSelectROIsNr)), 'UniformOutput', false));
    localDeltaDirSoma_values = reshape(localDeltaDirSoma_values, 4, [])';
    localDeltaDirSomaInput_values = cell2mat(cellfun(@(x) x.localDeltaDirSoma, num2cell(Spines(dirSelectInputNr)), 'UniformOutput', false));
    localDeltaDirSomaInput_values = reshape(localDeltaDirSomaInput_values, 4, [])';

    figure(54)
    subplot(2,2,1)
    boxplot(localDeltaOriSoma_values,'color',coc_prop(15,:))
    ylim([0 90])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaOriSoma_values,1),1), localDeltaOriSoma_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaOri of local environment to Soma')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    plot(deltaOriValues,localDeltaOriSoma_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 90])
    ylim([0 90])
    ylabel(sprintf('deltaOri of local environment \n(5 \\mum) to Soma'))
    xlabel('delta Ori spine')
    box off
    subplot(2,2,2)
    boxplot(localDeltaOriSomaInput_values,'color',coc_prop(3,:))
    ylim([0 90])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaOriSomaInput_values,1),1), localDeltaOriSomaInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaOri of local input environment to Soma')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot(deltaOriInputValues,localDeltaOriSomaInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 90])
    ylim([0 90])
    ylabel(sprintf('deltaOri of local input environment \n(5 \\mum) to Soma'))
    xlabel('delta Ori input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'54_Local deltaOri to Soma vs deltaOri Spine vs Input.png'))

    figure(64)
    subplot(2,2,1)
    boxplot(localDeltaDirSoma_values,'color',coc_prop(15,:))
    ylim([0 180])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaDirSoma_values,1),1), localDeltaDirSoma_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaDir of local environment \nto Soma')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    plot([Spines(dirSelectROIsNr).prefDir],localDeltaDirSoma_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 360])
    ylim([0 180])
    ylabel(sprintf('deltaDir of local environment \n(5 \\mum) to Soma'))
    xlabel('prefDir spine')
    box off
    subplot(2,2,2)
    boxplot(localDeltaDirSomaInput_values,'color',coc_prop(3,:))
    ylim([0 180])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDeltaDirSomaInput_values,1),1), localDeltaDirSomaInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
    ylabel('deltaDir of local input environment \nto Soma')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot([Spines(dirSelectInputNr).prefDir],localDeltaDirSomaInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 360])
    ylim([0 180])
    ylabel(sprintf('deltaDir of local input environment \n(5 \\mum) to Soma'))
    xlabel('prefDir input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'64_Local deltaDir to Soma vs prefDir Spine vs Input.png'))

    %Fig 55/65: nearest SimilarPrefOri/dir
    figure(55)
    subplot(1,2,1)
    plot([Spines(oriSelectROIsNr).prefOri],[Spines(oriSelectROIsNr).nearestSimilarPrefOri],'*', 'color', coc_prop(13,:));
    xlim([0 180])
    ylim([0 30])
    ylabel(sprintf('Nearest neighbor with \nsimilar preference in \\mum'))
    xlabel('pref Ori spine')
    box off
    subplot(1,2,2)
    plot([Spines(oriSelectInputNr).prefOri],[Spines(oriSelectInputNr).nearestSimilarPrefOri],'*', 'color', coc_prop(3,:));
    xlim([0 180])
    ylim([0 30])
    ylabel(sprintf('Nearest neighbor input with \nsimilar preference in \\mum'))
    xlabel('pref Ori input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'55_NearestNeighbor prefOri Spine vs. input.png'))

    randInputNearestNeighborAvgOri = NaN(10000,1);
    nonNanInputsOri = sum(~isnan([Spines(oriSelectInputNr).nearestSimilarPrefOri]));
    for r = 1:10000
        newOrder = randperm(length(oriSelectROIsNr));
        nearestNeighborRandOri = [Spines(oriSelectROIsNr(newOrder)).nearestSimilarPrefOri];
        %remove the nans
        nearestNeighborRandOri=rmmissing(nearestNeighborRandOri);
        randInputNearestNeighborAvgOri(r)=mean(nearestNeighborRandOri(1:nonNanInputsOri));
    end
    InputNearestNeighborOriAvg=nanmean([Spines(dirSelectInputNr).nearestSimilarPrefOri]);

    if nonNanInputsOri > 0
        figure
        h = cdfplot(randInputNearestNeighborAvgOri);
        xlabel('Mean nearest-neighbor with similar pref Ori distance in um')
        ylabel('Probablity')
        title('')
        xline(InputNearestNeighborOriAvg, '--r');
        legend('Random Data', 'Inputs')
        box off
        grid off
        h.Color = coc_prop(13,:);
        h.LineWidth = 3;
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDir,'55a_NearestNeighbor prefOri input to random.png'))
    end

    figure(65)
    subplot(1,2,1)
    plot([Spines(dirSelectROIsNr).prefDir],[Spines(dirSelectROIsNr).nearestSimilarPrefDir],'*', 'color', coc_prop(13,:));
    xlim([0 360])
    ylim([0 30])
    ylabel(sprintf('Nearest neighbor with \nsimilar preference in \\mum'))
    xlabel('pref Dir spine')
    box off
    subplot(1,2,2)
    plot([Spines(dirSelectInputNr).prefDir],[Spines(dirSelectInputNr).nearestSimilarPrefDir],'*', 'color', coc_prop(3,:));
    xlim([0 360])
    ylim([0 30])
    ylabel(sprintf('Nearest neighbor input with \nsimilar preference in \\mum'))
    xlabel('pref Dir input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'65_NearestNeighbor prefDir spine vs. input.png'))

    randInputNearestNeighborAvg = NaN(1000,1);
    nonNanInputs = sum(~isnan([Spines(dirSelectInputNr).nearestSimilarPrefDir]));
    for r = 1:10000
        newOrder = randperm(length(dirSelectROIsNr));
        nearestNeighborRand = [Spines(newOrder).nearestSimilarPrefDir];
        %remove the nans
        nearestNeighborRand=rmmissing(nearestNeighborRand);
        try
            randInputNearestNeighborAvg(r)=mean(nearestNeighborRand(1:nonNanInputs));
        catch
            randInputNearestNeighborAvg(r) = NaN;
        end
    end
    InputNearestNeighborAvg=nanmean([Spines(dirSelectInputNr).nearestSimilarPrefDir]);
    try
        figure
        h = cdfplot(randInputNearestNeighborAvg);
        xlabel('Mean nearest-neighbor with similar pref Dir distance in um')
        ylabel('Probablity')
        title('')
        xline(InputNearestNeighborAvg, '--r');
        legend('Random Data', 'Inputs')
        box off
        grid off
        h.Color = coc_prop(13,:);
        h.LineWidth = 3;
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDir,'65a_NearestNeighbor prefDir input to random.png'))
    catch
        disp('Not enough data for figure 65a')
    end
    %Fig 56/66: local OSI/DSI
    localOSI_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(Spines(respROINrs)), 'UniformOutput', false));
    localOSI_values = reshape(localOSI_values, 4, [])';
    localOSIInput_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(Spines(respInputsNr)), 'UniformOutput', false));
    localOSIInput_values = reshape(localOSIInput_values, 4, [])';
    localDSI_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(Spines(respROINrs)), 'UniformOutput', false));
    localDSI_values = reshape(localDSI_values, 4, [])';
    localDSIInput_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(Spines(respInputsNr)), 'UniformOutput', false));
    localDSIInput_values = reshape(localDSIInput_values, 4, [])';

    figure(56)
    subplot(2,2,1)
    boxplot(localOSI_values,'color',coc_prop(15,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localOSI_values,1),1), localOSI_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
    ylabel('local OSI')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    plot([Spines(respROINrs).OSI],localOSI_values(:,2),'*', 'color',coc_prop(13,:));
    xlim([0 1])
    ylim([0 1])
    ylabel(sprintf('local OSI (5 \\mum)'))
    xlabel('OSI spine')
    box off
    subplot(2,2,2)
    boxplot(localOSIInput_values,'color',coc_prop(3,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localOSIInput_values,1),1), localOSIInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
    ylabel('local OSI inputs')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot([Spines(respInputsNr).OSI],localOSIInput_values(:,2),'*', 'color',coc_prop(3,:));
    xlim([0 1])
    ylim([0 1])
    ylabel(sprintf('local OSI (5 \\mum)'))
    xlabel('OSI input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'56_LocalOSI spine vs. input.png'))

    figure(66)
    subplot(2,2,1)
    boxplot(localDSI_values,'color',coc_prop(15,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDSI_values,1),1), localDSI_values, 'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
    ylabel('local DSI')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,3)
    plot([Spines(respROINrs).DSI],localDSI_values(:,2),'*', 'color',coc_prop(13,:));
    xlim([0 1])
    ylim([0 1])
    ylabel(sprintf('local DSI (5 \\mum)'))
    xlabel('DSI spine')
    box off
    subplot(2,2,2)
    boxplot(localDSIInput_values,'color',coc_prop(3,:))
    ylim([0 1])
    xlim([0 5])
    hold on
    scatter(repmat(1:4,size(localDSIInput_values,1),1), localDSIInput_values, 'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
    ylabel('local DSI inputs')
    xlabel(sprintf('Radius size in \\mum'))
    xticklabels({'2.5', '5', '7.5', '10'})
    box off
    subplot(2,2,4)
    plot([Spines(respInputsNr).DSI],localDSIInput_values(:,2),'*', 'color',coc_prop(3,:));
    xlim([0 1])
    ylim([0 1])
    ylabel(sprintf('local DSI (5 \\mum)'))
    xlabel('DSI input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'66_LocalDSI spine vs. input.png'))

    %get those only from the selective ROIs
    localOSISel_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(Spines(oriSelectROIsNr)), 'UniformOutput', false));
    localOSISel_values = reshape(localOSISel_values, 4, [])';
    localOSISelInput_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(Spines(oriSelectInputNr)), 'UniformOutput', false));
    localOSISelInput_values = reshape(localOSISelInput_values, 4, [])';
    localDSISel_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(Spines(dirSelectROIsNr)), 'UniformOutput', false));
    localDSISel_values = reshape(localDSISel_values, 4, [])';
    localDSISelInput_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(Spines(dirSelectInputNr)), 'UniformOutput', false));
    localDSISelInput_values = reshape(localDSISelInput_values, 4, [])';

    %Fig 57/67: local OSI/DSI vs pref Ori/Dir
    figure(57)
    subplot(2,2,1)
    plot([Spines(oriSelectROIsNr).prefOri],localOSISel_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('local OSI (5 \\mum)'))
    xlabel('pref Ori spine')
    box off
    subplot(2,2,2)
    plot(deltaOriValues,localOSISel_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('local OSI (5 \\mum)'))
    xlabel('delta Ori spine')
    box off
    subplot(2,2,3)
    plot([Spines(oriSelectInputNr).prefOri],localOSISelInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('local OSI input (5 \\mum)'))
    xlabel('pref Ori input')
    box off
    subplot(2,2,4)
    plot(deltaOriInputValues,localOSISelInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 90])
    ylim([0 1])
    ylabel(sprintf('local OSI input (5 \\mum)'))
    xlabel('delta Ori input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'57_LocalOSI vs. prefOri and deltaOri spine vs. input.png'))

    figure(67)
    subplot(1,2,1)
    plot([Spines(dirSelectROIsNr).prefDir],localDSISel_values(:,2),'*', 'color', coc_prop(13,:));
    xlim([0 360])
    ylim([0 1])
    ylabel(sprintf('local DSI (5 \\mum)'))
    xlabel('pref Dir spine')
    box off
    subplot(1,2,2)
    plot([Spines(dirSelectInputNr).prefDir],localDSISelInput_values(:,2),'*', 'color', coc_prop(3,:));
    xlim([0 360])
    ylim([0 1])
    ylabel(sprintf('local DSI input (5 \\mum)'))
    xlabel('pref Dir input')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'67_LocalDSI vs. prefDir spine vs. input.png'))

    %--------------------------------------------------------------------------
    %70s: Properties of dendritic segments
    %What is the overal property of the dendritic segment? How homeogeneous is
    %it, how different from soma, how selective, ...

    %Fig 70: Branch circular dispersion ori/dir
    figure(70)
    subplot(1,4,1)
    hold on
    boxplot([Dendrites.circDispersion]', 'color', coc_prop(15,:))
    ylim([0 45])
    ylabel('Circular dispersion dendritic segment')
    xticklabels('all segments')
    box off
    scatter(repmat(1:1,size([Dendrites.circDispersion],1),1),[Dendrites.circDispersion],'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,2)
    hold on
    boxplot([Dendrites(inputBranches).circDispersion]', 'color', coc_prop(3,:))
    ylim([0 45])
    xticklabels('input segments')
    ax = gca;
    ax.YColor = 'none';
    scatter(repmat(1:1,size([Dendrites(inputBranches).circDispersion],1),1),[Dendrites(inputBranches).circDispersion],'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    box off
    subplot(1,4,3)
    hold on 
    boxplot([Dendrites.dirCircDispersion]', 'color', coc_prop(15,:))
    ylim([0 90])
    ylabel('Circular (dir) dispersion dendritic segment')
    xticklabels('all segments')
    box off
    scatter(repmat(1:1,size([Dendrites.dirCircDispersion],1),1),[Dendrites.dirCircDispersion],'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,4)
    hold on
    boxplot([Dendrites(inputBranches).dirCircDispersion]', 'color', coc_prop(3,:))
    ylim([0 90])
    xticklabels('input segments')
    ax = gca;
    ax.YColor = 'none';
    scatter(repmat(1:1,size([Dendrites(inputBranches).dirCircDispersion],1),1),[Dendrites(inputBranches).dirCircDispersion],'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'70_circularDispersion Inputbranches.png'))

    %Fig 71: Branch selectivity
    figure(71)
    subplot(1,4,1)
    boxplot([Dendrites.medianOSI]','color',coc_prop(15,:))
    ylim([0 1])
    xlim([0 2])
    hold on
    ylabel('median OSI')
    xticklabels('all segments')
    box off
    scatter(repmat(1:1,size([Dendrites.medianDSI],1),1),[Dendrites.medianOSI],'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,2)
    boxplot([Dendrites(inputBranches).medianOSI]','color',coc_prop(3,:))
    ylim([0 1])
    xlim([0 2])
    xticklabels('Input segments')
    hold on
    ax = gca;
    ax.YColor = 'none';
    box off
    scatter(repmat(1:1,size([Dendrites(inputBranches).medianOSI],1),1),[Dendrites(inputBranches).medianOSI],'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,3)
    boxplot([Dendrites.medianDSIvect]','color',coc_prop(15,:))
    ylim([0 1])
    hold on
    ylabel('median DSI')
    xticklabels('all segments')
    box off
    scatter(repmat(1:1,size([Dendrites.medianDSIvect],1),1),[Dendrites.medianDSIvect],'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,4)
    boxplot([Dendrites(inputBranches).medianDSIvect]','color',coc_prop(3,:));
    ylim([0 1])
    hold on
    xticklabels('Input segments')
    box off
    scatter(repmat(1:1,size([Dendrites(inputBranches).medianDSIvect],1),1),[Dendrites(inputBranches).medianDSIvect],'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    ax = gca;
    ax.YColor = 'none';
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'71_Branch selectivity Inputbranches.png'))

    %Fig 72: delta Ori/delta Dir
    figure(72)
    subplot(1,4,1)
    boxplot([Dendrites.deltaOri]','color',coc_prop(15,:))
    ylim([0 90])
    hold on
    ylabel('delta Ori')
    xticklabels('all segments')
    box off
    scatter(repmat(1:1,size([Dendrites.deltaOri],1),1),[Dendrites.deltaOri],'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,2)
    boxplot([Dendrites(inputBranches).deltaOri]','color',coc_prop(3,:))
    ylim([0 90])
    hold on
    xticklabels('input segments')
    box off
    scatter(repmat(1:1,size([Dendrites(inputBranches).deltaOri],1),1),[Dendrites(inputBranches).deltaOri],'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    ax = gca;
    ax.YColor = 'none';
    subplot(1,4,3)
    boxplot([Dendrites.deltaDir]','color',coc_prop(15,:))
    ylim([0 180])
    hold on
    ylabel('delta Dir')
    xticklabels('all segments')
    box off
    scatter(repmat(1:1,size([Dendrites.deltaDir],1),1),[Dendrites.deltaDir],'filled','MarkerFaceColor',coc_prop(13,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    subplot(1,4,4)
    boxplot([Dendrites(inputBranches).deltaDir]','color',coc_prop(3,:))
    ylim([0 180])
    hold on
    xticklabels('input segments')
    box off
    ax = gca;
    ax.YColor = 'none';
    scatter(repmat(1:1,size([Dendrites(inputBranches).deltaDir],1),1),[Dendrites(inputBranches).deltaDir],'filled','MarkerFaceColor',coc_prop(5,:), 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'72_DeltaOri_deltaDir Inputbranches.png'))

    %Fig 73: Branch circular dispersion vs. deltaOri
    figure(73)
    subplot(1,2,1)
    plot([Dendrites.deltaOri],[Dendrites.circDispersion], '*','color', coc_prop(13,:))
    xlim([0 90])
    xlabel('Delta Ori of dendritic segment to soma')
    ylim([0 90])
    ylabel('Dendritic segment circular dispersion')
    box off
    subplot(1,2,2)
    plot([Dendrites(inputBranches).deltaOri],[Dendrites(inputBranches).circDispersion], '*','color', coc_prop(3,:))
    xlim([0 90])
    xlabel('Delta Ori of dendritic input segment to soma')
    ylim([0 90])
    ylabel('Dendritic segment circular dispersion')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'73_DeltaOri_circDispersion Inputbranches.png'))


    %Fig 74: Branch circular dispersion dir vs. deltaDir
    figure(74) 
    subplot(1,2,1)
    plot([Dendrites.deltaDir],[Dendrites.dirCircDispersion], '*','color', coc_prop(13,:))
    xlim([0 180])
    xlabel('Delta dir of dendritic segment to soma')
    ylim([0 180])
    ylabel('Dendritic segment dir circular dispersion')
    box off
    subplot(1,2,2)
    plot([Dendrites(inputBranches).deltaDir],[Dendrites(inputBranches).dirCircDispersion], '*','color', coc_prop(3,:))
    xlim([0 180])
    xlabel('Delta dir of dendritic segment to soma')
    ylim([0 180])
    ylabel('Dendritic segment dir circular dispersion')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir,'74_DeltaDir_dirCircDispersion.png'))

    %--------------------------------------------------------------------------
    %80s: Pairwise distances
    %Fig 80: Pref Ori
    goodPairs = find([pwMeasures.goodPair] == 1);
    oriPairs = find([pwMeasures.oriSelect] == 1);
    dirPairs = find([pwMeasures.dirSelect] == 1);
    
    nonInputsPairs = find([pwMeasures.inputPair] == 0);
    mixedPairs = find([pwMeasures.inputPair] == 1);
    inputsPairs = find([pwMeasures.inputPair] == 2);
    
    goodNonInputsPairs = intersect(goodPairs, nonInputsPairs);
    goodMixedPairs = intersect(goodPairs, mixedPairs);
    goodInputPairs = intersect(goodPairs, inputsPairs);
    
    oriNonInputsPairs = intersect(oriPairs, nonInputsPairs);
    oriMixedPairs = intersect(oriPairs, mixedPairs);
    oriInputPairs = intersect(oriPairs, inputsPairs);
    
    dirNonInputsPairs = intersect(dirPairs, nonInputsPairs);
    dirMixedPairs = intersect(dirPairs, mixedPairs);
    dirInputPairs = intersect(dirPairs, inputsPairs);
    
    figure(80)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(oriNonInputsPairs).distance],[pwMeasures(oriNonInputsPairs).deltaOri],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(oriMixedPairs).distance],[pwMeasures(oriMixedPairs).deltaOri],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(oriInputPairs).distance],[pwMeasures(oriInputPairs).deltaOri],2);
    ylabel('deltaOri')
    ylim([0 90])
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'80_DeltaOri vs. Distance.png'))
    
    figure(81)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(dirNonInputsPairs).distance],[pwMeasures(dirNonInputsPairs).deltaDir],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(dirMixedPairs).distance],[pwMeasures(dirMixedPairs).deltaDir],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(dirInputPairs).distance],[pwMeasures(dirInputPairs).deltaDir],coc_prop(3,:));
    ylabel('deltaDir')
    ylim([0 180])
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'81_DeltaDir vs. Distance.png'))
    
    figure(82)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).deltaOSI],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedPairs).distance],[pwMeasures(goodMixedPairs).deltaOSI],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputPairs).distance],[pwMeasures(goodInputPairs).deltaOSI],coc_prop(3,:));
    ylabel('deltaOSI')
    ylim([0 1])
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'82_deltaOSI vs. Distance.png'))
    
    figure(83)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).deltaDSI],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedPairs).distance],[pwMeasures(goodMixedPairs).deltaDSI],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputPairs).distance],[pwMeasures(goodInputPairs).deltaDSI],coc_prop(3,:));
    ylabel('deltaDSI')
    ylim([0 1])
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'83_deltaDSI vs. Distance.png'))
    
    figure(84)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).deltaDSIvect],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedPairs).distance],[pwMeasures(goodMixedPairs).deltaDSIvect],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputPairs).distance],[pwMeasures(goodInputPairs).deltaDSIvect],coc_prop(3,:));
    ylabel('deltaDSIvect')
    ylim([0 1])
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'84_deltaDSIvect vs. Distance.png'))
    
    figure(85)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).curveCorr],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedPairs).distance],[pwMeasures(goodMixedPairs).curveCorr],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputPairs).distance],[pwMeasures(goodInputPairs).curveCorr],coc_prop(3,:));
    ylabel('Tuning curve correlation')
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'85_Tuning curve correlation vs. Distance.png'))
    
    figure(86)
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).trialCorr],coc_prop(13,:));
    hold on
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedPairs).distance],[pwMeasures(goodMixedPairs).trialCorr],coc_prop(10,:));
    plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputPairs).distance],[pwMeasures(goodInputPairs).trialCorr],coc_prop(3,:));
    ylabel('Trial-to-trial correlation')
    legend({'non-inputs'; 'non-input & inputs'; 'inputs'}, 'box','off')
    saveas(gcf, fullfile(saveDir,'86_Trial-to-trial correlation vs. Distance.png'))

end
%% 6.) Save
save([saveDirAna filesep '04_FunctionalInputReconstruction.mat'], 'Spines', 'Dendrites', 'Soma','pwMeasures', '-mat') 