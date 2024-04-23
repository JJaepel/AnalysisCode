function MultiCellInputAnalysis

%Combines all functional and anatomical data from all cells that were
%functionally characterized, matched to the confocal data and inputs
%analyzed

%0.) Read in xls file to see which cells can be added
%1.) Read in the files and combine the structures
%2.) Load the other data structure for all funtional data
%3.) Calculate
%4.) Plot

%% 0.) Define folders, colors and other standards for figures etc.
% folders
baseDir = 'Z:\Juliane\InputAnalysis\';
saveDir = [baseDir filesep 'Pooled Data\'];
saveDirImages = [saveDir filesep 'Inputs'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end
if ~exist(saveDirImages, 'dir')
    mkdir(saveDirImages)
end

%figures
% Defaults for Cell press. 1 col: 85mm, 1.5 col: 114mm, 2col:174mm
% Defaults for Nature press. 1 col: 89mm, 1.5 col: 136mm, 2col:183mm
fig.width= 5;               % fig.widthin cm
fig.height = fig.width;     % fig.height in cm (golden ratio default  (1 + n.sqrt(5)) / 2 1/GR)
fig.alw = 1;                % AxesLineWidth
fig.fsz = 14;               % Fontsize
fig.ln = 3;                 % lineWidth for e.g. averages
colorlevels = 12;
fig.cocA19 = cbrewer('seq', 'RdPu', colorlevels);
fig.cocV1 = cbrewer('seq', 'YlGn', colorlevels);
fig.cocAll = cbrewer('seq', 'Greys', colorlevels);

opengl software         % to get the axes on figures with transparency back!
set(0,'DefaultAxesTickDir', 'out')
set(0,'DefaultAxesFontSize', fig.fsz)
set(0,'DefaultAxesTickLength', [0.02 0.025]);

%% 1.) Load the pooled spine data

%load the data
temp = load([saveDir 'SpineData.mat']);
allSpines = temp.allSpines;
allBranches = temp.allBranches;
cellDetails = temp.cellDetails;
cellDetails = arrayfun(@(x) setfield(x, 'nrInputs', zeros(1,1)), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'apicalInputs', zeros(1,1)), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'basalInputs', zeros(1,1)), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'STEDDone', zeros(1,1)), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'apicalSTEDDone', zeros(1,1)), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'basalSTEDDone', zeros(1,1)), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'oriPref', []), cellDetails);
cellDetails = arrayfun(@(x) setfield(x, 'dirPref', []), cellDetails);

%add in a new field input
allSpines = arrayfun(@(x) setfield(x, 'input', zeros(1,1)), allSpines);
allSpines = arrayfun(@(x) setfield(x, 'STEDDone', zeros(1,1)), allSpines);
allSpines = arrayfun(@(x) setfield(x, 'inputType', zeros(1,1)), allSpines);

%add in the new fields for the Branches
allBranches = arrayfun(@(x) setfield(x, 'STEDDone', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'STEDTrace', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'STEDPerc', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'STEDLength', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'numInputs', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'InputDensity', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'percOfAllSpines', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'inputType', []), allBranches);
allBranches = arrayfun(@(x) setfield(x, 'ROIs', []), allBranches);

%start a Structure for pairwise measurements

for c = 1:length(cellDetails)
    if contains(cellDetails(c).input,  'Yes')
        %load the multimodaldata
        file = [baseDir char(cellDetails(c).animal) filesep char(cellDetails(c).cellName) '/E - Analysis/04_FunctionalInputReconstruction.mat'];
        if exist(file, 'file') == 2
            temp = load(file);
            InputSpines = temp.Spines;
            
            %find the inputspines
            inputs = find([InputSpines.Input] == 1);
            
            %classify the inputType based on the cellDetails
            inputType = cellDetails(c).inputType;
            switch inputType
                case 'A19'
                    iT = 1;
                case 'V1' 
                    iT = 2;
            end
            
            %find the different spines classes 
            apicalSpines =  find(cellfun(@(x) strcmp(x, 'apical'),{InputSpines.type}));
            basalSpines =  find(cellfun(@(x) strcmp(x, 'basal'),{InputSpines.type}));
            
            if ~isempty(inputs)
                cellDetails(c).nrInputs = length(inputs);
                cellDetails(c).apicalInputs = length(intersect(inputs, apicalSpines));
                cellDetails(c).basalInputs = length(intersect(inputs, basalSpines));
                %which of allSpines are coming from this cell?
                cellSpines = find([allSpines.CellNr] == cellDetails(c).cellNr);
                %go through all of them and mark those as inputs
                cellInputSpines = cellSpines(inputs);
                for ci = 1:length(cellInputSpines)
                    allSpines(cellInputSpines(ci)).input = iT;
                end
            else
                cellDetails(c).nrInputs = 0;
            end
        
            %add in whether they've been imaged with STED and the inputType
            cellSpines = find([allSpines.CellNr] == cellDetails(c).cellNr);
            for s = 1:length(cellSpines)
                %its the same list, so you can start at the first and add
%                 %them in
                allSpines(cellSpines(s)).STEDDone = InputSpines(s).STEDdone;
                allSpines(cellSpines(s)).inputType = iT;
            end
            
            %add in the total of STEDDone spines, separated by apical and
            %basal
            cellDetails(c).STEDDone = length(find([InputSpines.STEDdone] == 1));
            cellDetails(c).apicalSTEDDone = length(intersect(apicalSpines, find([InputSpines.STEDdone] == 1)));
            cellDetails(c).oriPref = temp.Soma.prefOri;
            cellDetails(c).dirPref = temp.Soma.prefDir;
            
            %get the pwMeasures
            temp2 = load(file, 'pwMeasures');
            temp2.pwMeasures = arrayfun(@(x) setfield(x, 'inputType', iT), temp2.pwMeasures);
            if c == 1
                pwMeasures = temp2.pwMeasures;
            else
                pwMeasures = [pwMeasures temp2.pwMeasures];
            end
        else
            cellDetails(c).nrInputs = 0;
        end
        
        %add in the input data  to the branches
        %find the branches for this cell
        InputBranches = temp.Dendrites;
        cellBranches = find([allBranches.CellNr] == cellDetails(c).cellNr);
        for cB = 1:length(cellBranches)
            allBranches(cellBranches(cB)).STEDDone = InputBranches(cB).STEDDone;
            allBranches(cellBranches(cB)).STEDTrace = InputBranches(cB).STEDTrace;
            allBranches(cellBranches(cB)).STEDPerc = InputBranches(cB).STEDPerc;
            allBranches(cellBranches(cB)).STEDLength = InputBranches(cB).STEDLength;
            allBranches(cellBranches(cB)).numInputs = InputBranches(cB).numInputs;
            allBranches(cellBranches(cB)).InputDensity = InputBranches(cB).InputDensity;
            allBranches(cellBranches(cB)).percOfAllSpines = InputBranches(cB).percOfAllSpines;
            allBranches(cellBranches(cB)).ROIs = InputBranches(cB).ROIs;
            allBranches(cellBranches(cB)).inputType = iT;
        end

    else
        cellDetails(c).nrInputs = 0;
    end
end

%% 2.) Selectors for easier data selection

%which ones are STEDdone
STEDROIs = find([allSpines.STEDDone] == 1);
A19STEDSpines = intersect(STEDROIs, find([allSpines.inputType] == 1));
V1STEDSpines = intersect(STEDROIs, find([allSpines.inputType] == 2));

%which ones are 2p spines
TwoPROIsNR = find([allSpines.TwoPMatch] ==1);

%which of those are responsive
goodROIs = find([allSpines.good] == 1);

%which of those are selective
oriTwoPROIs = find([allSpines.OSI] > 0.1); %all oriselectROIs, size = allfuncSpines
oriSelect = TwoPROIsNR(oriTwoPROIs);
oriGood = intersect(oriSelect, goodROIs);

dirTwoPROIs = find([allSpines.DSIvect] > 0.1); %all oriselectROIs, size = allfuncSpines
dirSelect = TwoPROIsNR(dirTwoPROIs);
dirGood = intersect(dirSelect, goodROIs);

%which ones are A19 spines
A19ROIs = find([allSpines.input] == 1);
V1ROIs = find([allSpines.input] == 2);

%intersect for all 2p inputs
A19Funct = intersect(TwoPROIsNR, A19ROIs);
V1Funct = intersect(TwoPROIsNR, V1ROIs);

%intersect for responsive inputs
goodA19 = intersect(goodROIs, A19ROIs);
goodV1 = intersect(goodROIs, V1ROIs);

%intersect for ori select input
oriGoodA19 = intersect(oriGood, A19ROIs);
oriGoodV1 = intersect(oriGood, V1ROIs);

%intersect for dir select input
dirGoodA19 = intersect(dirGood, A19ROIs);
dirGoodV1 = intersect(dirGood, V1ROIs);

%apical vs. basal
apicalROIs = find(cellfun(@(x) strcmp(x, 'apical'),{allSpines.type}));
basalROIs = find(cellfun(@(x) strcmp(x, 'basal'),{allSpines.type}));

%make sure they are STEDDone
apicalSTEDROIs = intersect(apicalROIs, STEDROIs);
basalSTEDROIs = intersect(basalROIs, STEDROIs);

apicalA19ROIs = intersect(apicalROIs, A19ROIs);
apicalV1ROIs = intersect(apicalROIs, V1ROIs);

basalA19ROIs = intersect(basalROIs, A19ROIs);
basalV1ROIs = intersect(basalROIs, V1ROIs);

apicalA19STEDSpines = intersect(A19STEDSpines, apicalSTEDROIs);
apicalV1STEDSpines = intersect(V1STEDSpines, apicalSTEDROIs);

basalA19STEDSpines = intersect(A19STEDSpines, basalSTEDROIs);
basalV1STEDSpines = intersect(V1STEDSpines, basalSTEDROIs);

%branches
A19Branches = find([allBranches.inputType] == 1);
V1Branches = find([allBranches.inputType] == 2);

apicalBranches = find(cellfun(@(x) strcmp(x, 'apical'),{allBranches.type}));
basalBranches = find(cellfun(@(x) strcmp(x, 'basal'),{allBranches.type}));

apicalA19Branches = intersect(A19Branches, apicalBranches);
apicalV1Branches = intersect(V1Branches, apicalBranches);

basalA19Branches = intersect(A19Branches, basalBranches);
basalV1Branches = intersect(V1Branches, basalBranches);


%% 3.) Plot overview
%Figure Overview
%100s - ALL STRUCTURAL
%100s - Plotting on top of constructed cell
%Fig 100: Plotting apical and basal of A19 inputs on top of respective dendrites
%Fig 101: Plotting apical and basal of V1 inputs on top of respective dendrites
%Fig 102: Plotting cell Nr of A19 inputs on top of respective dendrites
%Fig 103: Plotting cell Nr of V1 inputs on top of respective dendrites

%110s - Fraction of inputs

%120s - Distance to soma
%Fig 120-122: Distance from soma for all inputs, separated by apical and
%basal
%Fig 123-126: Same as Fig 120-126, but as cdfplot

%130s - Branch Order
%Fig 130-132: Histogram of branch order for all inputs, separated for 
%apical/basal
%Fig 133-135: Cdfplot of branch order for all inputs, separated for 
%apical/basal

%140s - Branch specific analysis
%Fig 140: Number of inputs per dendrites 
%Fig 141: Input density per dendrite for all/branches with more than one
%input
%Fig 142: Number of inputs vs. length of STED coverage

%150s - Clustering of inputs
%Fig 150: Distance to neares input on same branch
%Fig 151: Mean distance to nearest input
%Fig 152: Distribution of distance to nearest other input

%200s - ALL FUNCTIONAL
%200s - Plotting on top of constructed cell
%Fig 200: Plotting deltaOri of A19 inputs on top of respective dendrites
%Fig 201: Plotting deltaOri of V1 inputs on top of respective dendrites
%Fig 202: Plotting deltaDir of A19 inputs on top of respective dendrites
%Fig 203: Plotting deltaDir of V1 inputs on top of respective dendrites
%Fig 204: Plotting OSI of A19 inputs on top of respective dendrites
%Fig 205: Plotting OSI of V1 inputs on top of respective dendrites
%Fig 206: Plotting DSI of A19 inputs on top of respective dendrites
%Fig 207: Plotting DSI of V1 inputs on top of respective dendrites

%210s - Quick overview of general properties
%210 - Pie chart of fraction of responsive & ori-select
%211 - Distibution of delta Ori 
%212 - Distribution of OSI
%213 - Distribution of DSI
%214 - distribution of bandwith
%215 - Plotting a constructed cell with delta Ori

%220s & 30s: Local environment ori(220s) and dir (230s)
%Fig 220/230: local (dir) Dispersion
%Fig 221/231: HI (dir)
%Fig 222/232: local (dir) dispersion vs. HI (dir)
%Fig 223/233: local deltaOri/Dir
%Fig 224/234: local delteOriSoma/deltaDirSOma
%Fig 225/235: nearest SimilarPrefOri/dir
%Fig 226/236: local OSI/DSI

%240s: Properties of dendritic segments
%Fig 240: Branch circular dispersion ori/dir
%Fig 241: Branch selectivity
%Fig 242: delta Ori/delta Dir
%Fig 243: Branch circular dispersion vs. deltaOri
%Fig 244: Branch circular dispersion dir vs. deltaDir

%250s: Pairwise distances for A19 inputs
%Fig 250: Pref Ori
%Fig 251: Pref Dir
%Fig 252: OSI
%Fig 253: DSI
%Fig 254: DSIVect

%260s: Pairwise distances for V1 inputs
%Fig 260: Pref Ori
%Fig 261: Pref Dir
%Fig 262: OSI
%Fig 263: DSI
%Fig 264: DSIVect

%270s: Branch-wise analysis
%Fig 270: circular disperision
%Fig 271: dir circular dispersion
%Fig 272: deltaOri
%Fig 273: deltaDir
%Fig 274: medianOSI
%Fig 275: medianDSI
%Fig 276: medianDSIvect

%% 4.) Plots
%--------------------------------------------------------------------------
%% 100s STRUCTURE 
%Fig 100: Build a constructed cell A19 type
setFigProperties(100, fig)
constructInputCell(allSpines, allBranches, 'A19', 'struct')
saveas(gcf, fullfile(saveDirImages,'100_constructedInputCell_A19_type.png'))

%Fig 101: Build a constructed cell V1 type
setFigProperties(101, fig)
constructInputCell(allSpines, allBranches, 'V1', 'struct')
saveas(gcf, fullfile(saveDirImages,'101_constructedInputCell_V1_type.png'))

%Fig 102: Build a constructed cell A19 cell Nr
setFigProperties(102, fig)
constructInputCell(allSpines, allBranches, 'A19', 'CellNr')
saveas(gcf, fullfile(saveDirImages,'102_constructedInputCell_A19_cellNr.png'))

%Fig 103: Build a constructed cell V1 cell Nr
setFigProperties(103, fig)
constructInputCell(allSpines, allBranches, 'V1', 'CellNr')
saveas(gcf, fullfile(saveDirImages,'103_constructedInputCell_V1_cellNr.png'))

%% 110s - Fraction of spines
% Fig 110-112 - Fraction of A19/V1 spines of all STED spines, separated by
% apical & basal
setFigProperties(110, fig)
subplot(1,2,1)
bar(1,[length(A19ROIs)/length(A19STEDSpines)*100; (1-length(A19ROIs)/length(A19STEDSpines))*100], 'stacked');
barS = get(gca, 'Children');
set(barS(2), 'FaceColor', fig.cocA19(9,:));
set(barS(1), 'FaceColor', fig.cocAll(7,:));
text(1,length(A19ROIs)/length(A19STEDSpines)*100+5, [num2str(length(A19ROIs)/length(A19STEDSpines)*100,2) ' %'], 'HorizontalAlignment', 'center', 'Color', 'white');
xticklabels({'A19 Inputs'})
ylabel('Percentage of all spines')
box off
subplot(1,2,2)
bar(1,[length(V1ROIs)/length(V1STEDSpines)*100; (1-length(V1ROIs)/length(V1STEDSpines))*100], 'stacked');
barS = get(gca, 'Children');
set(barS(2), 'FaceColor', fig.cocV1(9,:));
set(barS(1), 'FaceColor', fig.cocAll(7,:));
text(1,length(V1ROIs)/length(V1STEDSpines)*100+5, [num2str(length(V1ROIs)/length(V1STEDSpines)*100,2) ' %'], 'HorizontalAlignment', 'center', 'Color', 'white');
xticklabels({'V1 Inputs'})
box off
set(gca, 'YColor', 'none');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '110_Percentage of inputs.png'))

setFigProperties(111, fig)
subplot(1,2,1)
bar(1,[length(apicalA19ROIs)/length(apicalA19STEDSpines)*100; (1-length(apicalA19ROIs)/length(apicalA19STEDSpines))*100], 'stacked');
barS = get(gca, 'Children');
set(barS(2), 'FaceColor', fig.cocA19(9,:));
set(barS(1), 'FaceColor', fig.cocAll(7,:));
text(1,length(apicalA19ROIs)/length(apicalA19STEDSpines)*100+5, [num2str(length(apicalA19ROIs)/length(apicalA19STEDSpines)*100,2) ' %'], 'HorizontalAlignment', 'center', 'Color', 'white');
xticklabels({'A19 Inputs'})
ylabel('Percentage of all apical spines')
box off
subplot(1,2,2)
bar(1,[length(apicalV1ROIs)/length(apicalV1STEDSpines)*100; (1-length(apicalV1ROIs)/length(apicalV1STEDSpines))*100], 'stacked');
barS = get(gca, 'Children');
set(barS(2), 'FaceColor', fig.cocV1(9,:));
set(barS(1), 'FaceColor', fig.cocAll(7,:));
text(1,length(apicalV1ROIs)/length(apicalV1STEDSpines)*100+5, [num2str(length(apicalV1ROIs)/length(apicalV1STEDSpines)*100,2) ' %'], 'HorizontalAlignment', 'center', 'Color', 'white');
xticklabels({'V1 Inputs'})
box off
set(gca, 'YColor', 'none');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '111_Percentage of apical inputs.png'))

setFigProperties(112, fig)
subplot(1,2,1)
bar(1,[length(basalA19ROIs)/length(basalA19STEDSpines)*100; (1-length(basalA19ROIs)/length(basalA19STEDSpines))*100], 'stacked');
barS = get(gca, 'Children');
set(barS(2), 'FaceColor', fig.cocA19(9,:));
set(barS(1), 'FaceColor', fig.cocAll(7,:));
text(1,length(basalA19ROIs)/length(basalA19STEDSpines)*100+5, [num2str(length(basalA19ROIs)/length(basalA19STEDSpines)*100,2) ' %'], 'HorizontalAlignment', 'center', 'Color', 'white');
xticklabels({'A19 Inputs'})
ylabel('Percentage of all basal spines')
box off
subplot(1,2,2)
bar(1,[length(basalV1ROIs)/length(basalV1STEDSpines)*100; (1-length(basalV1ROIs)/length(basalV1STEDSpines))*100], 'stacked');
barS = get(gca, 'Children');
set(barS(2), 'FaceColor', fig.cocV1(9,:));
set(barS(1), 'FaceColor', fig.cocAll(7,:));
text(1,length(basalV1ROIs)/length(basalV1STEDSpines)*100+5, [num2str(length(basalV1ROIs)/length(basalV1STEDSpines)*100,2) ' %'], 'HorizontalAlignment', 'center', 'Color', 'white');
xticklabels({'V1 Inputs'})
box off
set(gca, 'YColor', 'none');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '112_Percentage of basal inputs.png'))

% Fig 113 - 115: Fraction of A19/V1 spines for individual cells, separated
% by apical & basal
A19cells = intersect(find(cellfun(@(x) strcmp(x, 'A19'),{cellDetails.inputType})), find([cellDetails.nrInputs] > 0));
V1cells = intersect(find(cellfun(@(x) strcmp(x, 'V1'),{cellDetails.inputType})), find([cellDetails.nrInputs] > 0));

setFigProperties(113, fig)
subplot(1,2,1)
boxplot([cellDetails(A19cells).nrInputs]./[cellDetails(A19cells).STEDDone]*100, 'color', fig.cocA19(9,:))
hold on
xticklabels({'A19 Inputs'})
scatter(repmat(1:1, size(A19cells,1)), [cellDetails(A19cells).nrInputs]./[cellDetails(A19cells).STEDDone]*100, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9')
ylim([0 ceil(max([cellDetails.nrInputs]./[cellDetails.STEDDone]*100))])
ylabel('Percentage of all spines (per cell)')
box off
subplot(1,2,2)
boxplot([cellDetails(V1cells).nrInputs]./[cellDetails(V1cells).STEDDone]*100, 'color', fig.cocV1(9,:))
hold on
scatter(repmat(1:1, size(V1cells,1)), [cellDetails(V1cells).nrInputs]./[cellDetails(V1cells).STEDDone]*100, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9')
ylim([0 ceil(max([cellDetails.nrInputs]./[cellDetails.STEDDone]*100))])
xticklabels({'V1 Inputs'})
box off
set(gca, 'YColor', 'none');
saveas(gcf, fullfile(saveDirImages, '113_Cell based percentage of spines.png'))

setFigProperties(114, fig)
subplot(1,2,1)
boxplot([cellDetails(A19cells).apicalInputs]./[cellDetails(A19cells).apicalSTEDDone]*100, 'color', fig.cocA19(9,:))
hold on
xticklabels({'A19 Inputs'})
scatter(repmat(1:1, size(A19cells,1)), [cellDetails(A19cells).apicalInputs]./[cellDetails(A19cells).apicalSTEDDone]*100, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9')
ylim([0 ceil(max([cellDetails.apicalInputs]./[cellDetails.apicalSTEDDone]*100))])
ylabel('Percentage of all apical spines (per cell)')
box off
subplot(1,2,2)
boxplot([cellDetails(V1cells).apicalInputs]./[cellDetails(V1cells).STEDDone]*100, 'color', fig.cocV1(9,:))
hold on
scatter(repmat(1:1, size(V1cells,1)), [cellDetails(V1cells).apicalInputs]./[cellDetails(V1cells).apicalSTEDDone]*100, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9')
ylim([0 ceil(max([cellDetails.apicalInputs]./[cellDetails.apicalSTEDDone]*100))])
xticklabels({'V1 Inputs'})
box off
saveas(gcf, fullfile(saveDirImages, '114_Cell based percentage of apical spines.png'))

setFigProperties(115, fig)
subplot(1,2,1)
boxplot([cellDetails(A19cells).basalInputs]./[cellDetails(A19cells).basalSTEDDone]*100, 'color', fig.cocA19(9,:))
hold on
xticklabels({'A19 Inputs'})
scatter(repmat(1:1, size(A19cells,1)), [cellDetails(A19cells).basalInputs]./[cellDetails(A19cells).basalSTEDDone]*100, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9')
ylim([0 ceil(max([cellDetails.basalInputs]./[cellDetails.basalSTEDDone]*100))])
ylabel('Percentage of all basal spines (per cell)')
box off
subplot(1,2,2)
boxplot([cellDetails(V1cells).basalInputs]./[cellDetails(V1cells).basalSTEDDone]*100, 'color', fig.cocV1(9,:))
hold on
scatter(repmat(1:1, size(V1cells,1)), [cellDetails(V1cells).basalInputs]./[cellDetails(V1cells).basalSTEDDone]*100, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9')
ylim([0 ceil(max([cellDetails.basalInputs]./[cellDetails.basalSTEDDone]*100))])
xticklabels({'V1 Inputs'})
box off
set(gca, 'YColor', 'none');
saveas(gcf, fullfile(saveDirImages, '115_Cell based percentage of basal spines.png'))

% Fig 116 - 118: Fraction of A19/V1 spines for individual dendrites,
% separated by apical & basal
branchesWithSpines = find(cellfun(@(x) ~isempty(x),{allBranches.ROIs}));
getInputFraction = @(x) x.numInputs/sum([x.ROIs.STEDdone])*100; %define this to be able to get the deltaOriValues

%get data for branches
inputFractionA19Branches = cell2mat(arrayfun(getInputFraction, allBranches(intersect(branchesWithSpines, A19Branches)), 'UniformOutput', false))';
inputFractionA19Branches(isinf(inputFractionA19Branches)) = NaN;
inputFractionV1Branches = cell2mat(arrayfun(getInputFraction, allBranches(intersect(branchesWithSpines, V1Branches)), 'UniformOutput', false))';
inputFractionV1Branches(isinf(inputFractionV1Branches)) = NaN;
inputFractionApicalA19Branches = cell2mat(arrayfun(getInputFraction, allBranches(intersect(branchesWithSpines, apicalA19Branches)), 'UniformOutput', false))';
inputFractionApicalA19Branches(isinf(inputFractionApicalA19Branches)) = NaN;
inputFractionApicalV1Branches = cell2mat(arrayfun(getInputFraction, allBranches(intersect(branchesWithSpines, apicalV1Branches)), 'UniformOutput', false))';
inputFractionApicalV1Branches(isinf(inputFractionApicalV1Branches)) = NaN;
inputFractionBasalA19Branches = cell2mat(arrayfun(getInputFraction, allBranches(intersect(branchesWithSpines, basalA19Branches)), 'UniformOutput', false))';
inputFractionBasalA19Branches(isinf(inputFractionBasalA19Branches)) = NaN;
inputFractionBasalV1Branches = cell2mat(arrayfun(getInputFraction, allBranches(intersect(branchesWithSpines, basalV1Branches)), 'UniformOutput', false))';
inputFractionBasalV1Branches(isinf(inputFractionBasalV1Branches)) = NaN;

setFigProperties(116, fig)
subplot(1,2,1)
boxplot(inputFractionA19Branches, 'color', fig.cocA19(7,:))
hold on
xticklabels({'A19 Inputs'})
scatter(repmat(1:1, size(inputFractionA19Branches,1)), inputFractionA19Branches, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.7')
ylim([0 ceil(max([inputFractionA19Branches; inputFractionV1Branches]))])
ylabel('Percentage of all spines (per dendrite)')
box off
subplot(1,2,2)
boxplot(inputFractionV1Branches, 'color', fig.cocV1(7,:))
hold on
scatter(repmat(1:1, size(inputFractionV1Branches,1)), inputFractionV1Branches, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.7')
ylim([0 ceil(max([inputFractionA19Branches; inputFractionV1Branches]))])
xticklabels({'V1 Inputs'})
set(gca,'box','off','ycolor','w')
saveas(gcf, fullfile(saveDirImages, '116_Branch based percentage of spines.png'))

setFigProperties(117, fig)
subplot(1,2,1)
boxplot(inputFractionApicalA19Branches, 'color', fig.cocA19(7,:))
hold on
xticklabels({'A19 Inputs'})
scatter(repmat(1:1, size(inputFractionApicalA19Branches,1)), inputFractionApicalA19Branches, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.7')
ylim([0 ceil(max([inputFractionApicalA19Branches; inputFractionApicalV1Branches]))])
ylabel('Percentage of apical spines (per dendrite)')
box off
subplot(1,2,2)
boxplot(inputFractionApicalV1Branches, 'color', fig.cocV1(7,:))
hold on
scatter(repmat(1:1, size(inputFractionApicalV1Branches,1)), inputFractionApicalV1Branches, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.7')
ylim([0 ceil(max([inputFractionApicalA19Branches; inputFractionApicalV1Branches]))])
xticklabels({'V1 Inputs'})
set(gca,'box','off','ycolor','w')
saveas(gcf, fullfile(saveDirImages, '117_Branch based percentage of apical spines.png'))

setFigProperties(118, fig)
subplot(1,2,1)
boxplot(inputFractionBasalA19Branches, 'color', fig.cocA19(7,:))
hold on
xticklabels({'A19 Inputs'})
scatter(repmat(1:1, size(inputFractionBasalA19Branches,1)), inputFractionBasalA19Branches, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.7')
ylim([0 ceil(max([inputFractionBasalA19Branches; inputFractionBasalV1Branches]))])
ylabel('Percentage of basal spines (per dendrite)')
box off
subplot(1,2,2)
boxplot(inputFractionBasalV1Branches, 'color', fig.cocV1(7,:))
hold on
scatter(repmat(1:1, size(inputFractionBasalV1Branches,1)), inputFractionBasalV1Branches, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.7')
ylim([0 ceil(max([inputFractionBasalA19Branches; inputFractionBasalV1Branches]))])
xticklabels({'V1 Inputs'})
set(gca,'box','off','ycolor','w')
saveas(gcf, fullfile(saveDirImages, '118_Branch based percentage of basal spines.png'))

%% 120s - Distance from soma
%Fig 120-122: Distance from soma for all inputs, separated by apical and
%basal
%Fig 123-126: Same as Fig 110-116, but as cdfplot

%Fig 120-122: Distance from soma for all inputs, separated by apical and
%basal
setFigProperties(120, fig)
subplot(1,3,1)
distributionPlot([allSpines(STEDROIs).distToSoma]', 'color', fig.cocAll(7,:)); hold all
boxplot([allSpines(STEDROIs).distToSoma]')
box off
ylabel('Distance to soma')
xticklabels({'All spines'})
ylim([0 max([allSpines(STEDROIs).distToSoma])])
subplot(1,3,2)
distributionPlot([allSpines(A19ROIs).distToSoma]', 'color', fig.cocA19(9,:)); hold all
boxplot([allSpines(A19ROIs).distToSoma]')
ylim([0 max([allSpines(STEDROIs).distToSoma])])
set(gca,'box','off','ycolor','w')
xticklabels({'A19 Inputs'})
A19MeanDistSoma = NaN(length(A19STEDSpines),2);
A19MeanDistSoma(1:1:length(A19STEDSpines),1) = [allSpines(A19STEDSpines).distToSoma];
A19MeanDistSoma(1:1:length(A19ROIs),2) = [allSpines(A19ROIs).distToSoma];
[pA19MeanDistSoma,~] = kruskalwallis(A19MeanDistSoma, {'all', 'A19 inputs'}, 'off');
if pA19MeanDistSoma < 0.05
  text(1,max([allSpines(STEDROIs).distToSoma])+10,['p = ' num2str(pA19MeanDistSoma,2)], 'Color', fig.cocA19(9,:)) 
end
subplot(1,3,3)
distributionPlot([allSpines(V1ROIs).distToSoma]', 'color', fig.cocV1(9,:)); hold all
boxplot([allSpines(V1ROIs).distToSoma]')
xticklabels({'V1 Inputs'})
ylim([0 max([allSpines(STEDROIs).distToSoma])])
set(gca,'box','off','ycolor','w')
V1MeanDistSoma = NaN(length(V1STEDSpines),2);
V1MeanDistSoma(1:1:length(V1STEDSpines),1) = [allSpines(V1STEDSpines).distToSoma];
V1MeanDistSoma(1:1:length(V1ROIs),2) = [allSpines(V1ROIs).distToSoma];
[pV1MeanDistSoma,~] = kruskalwallis(V1MeanDistSoma, {'all', 'V1 inputs'}, 'off');
if pV1MeanDistSoma < 0.05
  text(1,max([allSpines(STEDROIs).distToSoma])+10,['p = ' num2str(pV1MeanDistSoma,2)], 'Color', fig.cocV1(9,:)) 
end
saveas(gcf, fullfile(saveDirImages, '120_DistanceFromBranchInputs.png'))

setFigProperties(121, fig)
subplot(1,3,1)
distributionPlot([allSpines(apicalSTEDROIs).distToSoma]', 'color', fig.cocAll(7,:)); hold all
boxplot([allSpines(apicalSTEDROIs).distToSoma]')
box off
ylabel('Distance to soma (apical spines)')
xticklabels({'All spines'})
ylim([0 max([allSpines(apicalSTEDROIs).distToSoma])])
subplot(1,3,2)
distributionPlot([allSpines(apicalA19ROIs).distToSoma]', 'color', fig.cocA19(9,:)); hold all
boxplot([allSpines(apicalA19ROIs).distToSoma]')
set(gca,'box','off','ycolor','w')
xticklabels({'A19 Inputs'})
ylim([0 max([allSpines(apicalSTEDROIs).distToSoma])])
A19MeanDistSomaApi = NaN(length(apicalA19STEDSpines),2);
A19MeanDistSomaApi(1:1:length(apicalA19STEDSpines),1) = [allSpines(apicalA19STEDSpines).distToSoma];
A19MeanDistSomaApi(1:1:length(apicalA19ROIs),2) = [allSpines(apicalA19ROIs).distToSoma];
[pA19MeanDistSomaApi,~] = kruskalwallis(A19MeanDistSomaApi, {'all', 'A19 inputs'}, 'off');
if pA19MeanDistSomaApi < 0.05
  text(1,max([allSpines(apicalSTEDROIs).distToSoma])+10,['p = ' num2str(pA19MeanDistSomaApi,2)], 'Color', fig.cocA19(9,:)) 
end
subplot(1,3,3)
distributionPlot([allSpines(apicalV1ROIs).distToSoma]', 'color', fig.cocV1(9,:)); hold all
boxplot([allSpines(apicalV1ROIs).distToSoma]')
title('V1 Inputs')
ylim([0 max([allSpines(apicalSTEDROIs).distToSoma])])
set(gca,'box','off','ycolor','w')
xticklabels({'V1 Inputs'})
V1MeanDistSomaApi = NaN(length(apicalV1STEDSpines),2);
V1MeanDistSomaApi(1:1:length(apicalV1STEDSpines),1) = [allSpines(apicalV1STEDSpines).distToSoma];
V1MeanDistSomaApi(1:1:length(apicalV1ROIs),2) = [allSpines(apicalV1ROIs).distToSoma];
[pV1MeanDistSomaApi,~] = kruskalwallis(V1MeanDistSomaApi, {'all', 'V1 inputs'}, 'off');
if pV1MeanDistSomaApi < 0.05
  text(1,max([allSpines(apicalSTEDROIs).distToSoma])+10,['p = ' num2str(pV1MeanDistSomaApi,2)], 'Color', fig.cocV1(9,:)) 
end
saveas(gcf, fullfile(saveDirImages, '121_DistanceFromBranchApicalInputs.png'))

setFigProperties(122, fig)
subplot(1,3,1)
distributionPlot([allSpines(basalSTEDROIs).distToSoma]', 'color', fig.cocAll(7,:)); hold all
boxplot([allSpines(basalSTEDROIs).distToSoma]')
box off
ylabel('Distance to soma (basal spines)')
ylim([0 max([allSpines(basalSTEDROIs).distToSoma])])
xticklabels({'All spines'})
subplot(1,3,2)
distributionPlot([allSpines(basalA19ROIs).distToSoma]', 'color', fig.cocA19(9,:)); hold all
boxplot([allSpines(basalA19ROIs).distToSoma]')
set(gca,'box','off','ycolor','w')
xticklabels({'A19 Inputs'})
ylim([0 max([allSpines(basalSTEDROIs).distToSoma])])
A19MeanDistSomaBas = NaN(length(basalA19STEDSpines),2);
A19MeanDistSomaBas(1:1:length(basalA19STEDSpines),1) = [allSpines(basalA19STEDSpines).distToSoma];
A19MeanDistSomaBas(1:1:length(basalA19ROIs),2) = [allSpines(basalA19ROIs).distToSoma];
[pA19MeanDistSomaBas,~] = kruskalwallis(A19MeanDistSomaBas, {'all', 'A19 inputs'}, 'off');
if pA19MeanDistSomaBas < 0.05
  text(1,max([allSpines(basalSTEDROIs).distToSoma])+10,['p = ' num2str(pA19MeanDistSomaBas,2)],'HorizontalAlignment', 'Center', 'Color', fig.cocA19(9,:)) 
end
subplot(1,3,3)
distributionPlot([allSpines(basalV1ROIs).distToSoma]', 'color', fig.cocV1(9,:)); hold all
boxplot([allSpines(basalV1ROIs).distToSoma]')
ylim([0 max([allSpines(basalSTEDROIs).distToSoma])])
set(gca,'box','off','ycolor','w')
xticklabels({'V1 Inputs'})
V1MeanDistSomaBas = NaN(length(basalV1STEDSpines),2);
V1MeanDistSomaBas(1:1:length(basalV1STEDSpines),1) = [allSpines(basalV1STEDSpines).distToSoma];
V1MeanDistSomaBas(1:1:length(basalV1ROIs),2) = [allSpines(basalV1ROIs).distToSoma];
[pV1MeanDistSomaBas,~] = kruskalwallis(V1MeanDistSomaBas, {'all', 'V1 inputs'}, 'off');
if pV1MeanDistSomaBas < 0.05
  text(1,max([allSpines(basalSTEDROIs).distToSoma])+10,['p = ' num2str(pV1MeanDistSomaBas,2)], 'HorizontalAlignment', 'Center', 'Color', fig.cocV1(9,:)) 
end
saveas(gcf, fullfile(saveDirImages, '122_DistanceFromBranchBasalInputs.png'))


%Fig 123-125: Same as Fig 120-122, but as cdfplot
setFigProperties(123, fig)
h(1,1) = cdfplot([allSpines(STEDROIs).distToSoma]); hold on
h(1,2) = cdfplot([allSpines(A19STEDSpines).distToSoma]);
h(1,3) = cdfplot([allSpines(V1STEDSpines).distToSoma]);
h(1,4) = cdfplot([allSpines(A19ROIs).distToSoma]);
h(1,5) = cdfplot([allSpines(V1ROIs).distToSoma]);

set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocAll(5,:), 'LineWidth', fig.ln);
set(h(1,3), 'Color', fig.cocAll(3,:), 'LineWidth', fig.ln);
set(h(1,4), 'Color', fig.cocA19(9,:), 'LineWidth', fig.ln);
set(h(1,5), 'Color', fig.cocV1(9,:), 'LineWidth', fig.ln);
legend('All Spines', 'A19 Cell Spines', 'V1 Cell Spines', 'A19 Inputs', 'V1 Inputs','Location', 'SouthEast')
legend('boxoff')
grid off
title('')
xlabel('Distance from start of branch')
ylabel('Cumulative fraction of all spines')
[hV1DistanceFromSoma,pV1DistanceFromSoma,~] = kstest2([allSpines(V1STEDSpines).distToSoma],[allSpines(V1ROIs).distToSoma]);
[hA19DistanceFromSoma,pA19DistanceFromSoma,~] = kstest2([allSpines(A19STEDSpines).distToSoma],[allSpines(A19ROIs).distToSoma]);
if hV1DistanceFromSoma
    text(5, 0.95,['p = ' num2str(pV1DistanceFromSoma,2)], 'Color', fig.cocV1(9,:));
end
if hA19DistanceFromSoma
    text(5, 0.90,['p = ' num2str(pA19DistanceFromSoma,2)], 'Color', fig.cocA19(9,:));
end
saveas(gcf, fullfile(saveDirImages, '123_CdfDistFromBranchCumulInputs.png'))

setFigProperties(124, fig)
h(1,1) = cdfplot([allSpines(apicalSTEDROIs).distToSoma]); hold on
h(1,2) = cdfplot([allSpines(apicalA19STEDSpines).distToSoma]);
h(1,3) = cdfplot([allSpines(apicalV1STEDSpines).distToSoma]);
h(1,4) = cdfplot([allSpines(apicalA19ROIs).distToSoma]);
h(1,5) = cdfplot([allSpines(apicalV1ROIs).distToSoma]);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocAll(5,:), 'LineWidth', fig.ln);
set(h(1,3), 'Color', fig.cocAll(3,:), 'LineWidth', fig.ln);
set(h(1,4), 'Color', fig.cocA19(9,:), 'LineWidth', fig.ln);
set(h(1,5), 'Color', fig.cocV1(9,:), 'LineWidth', fig.ln);
legend('All Spines', 'A19 Cell Spines', 'V1 Cell Spines', 'A19 Inputs', 'V1 Inputs','Location', 'SouthEast')
legend('boxoff')
grid off
title('')
xlabel('Distance from start of branch')
ylabel('Cumulative fraction of apical spines')
[hV1DistanceFromSomaApi,pV1DistanceFromSomaApi,~] = kstest2([allSpines(apicalV1STEDSpines).distToSoma],[allSpines(apicalV1ROIs).distToSoma]);
[hA19DistanceFromSomaApi,pA19DistanceFromSomaApi,~] = kstest2([allSpines(apicalA19STEDSpines).distToSoma],[allSpines(apicalA19ROIs).distToSoma]);
if hV1DistanceFromSomaApi
    text(5, 0.95,['p = ' num2str(pV1DistanceFromSomaApi,2)], 'Color', fig.cocV1(9,:));
end
if hA19DistanceFromSomaApi
    text(5, 0.90,['p = ' num2str(pA19DistanceFromSomaApi,2)], 'Color', fig.cocA19(9,:));
end
saveas(gcf, fullfile(saveDirImages, '124_CdfDistFromBranchCumulApicalInputs.png'))

setFigProperties(125, fig)
h(1,1) = cdfplot([allSpines(basalSTEDROIs).distToSoma]); hold on
h(1,2) = cdfplot([allSpines(basalA19STEDSpines).distToSoma]);
h(1,3) = cdfplot([allSpines(basalV1STEDSpines).distToSoma]);
h(1,4) = cdfplot([allSpines(basalA19ROIs).distToSoma]);
h(1,5) = cdfplot([allSpines(basalV1ROIs).distToSoma]);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocAll(5,:), 'LineWidth', fig.ln);
set(h(1,3), 'Color', fig.cocAll(3,:), 'LineWidth', fig.ln);
set(h(1,4), 'Color', fig.cocA19(9,:), 'LineWidth', fig.ln);
set(h(1,5), 'Color', fig.cocV1(9,:), 'LineWidth', fig.ln);
legend('All Spines', 'A19 Cell Spines', 'V1 Cell Spines', 'A19 Inputs', 'V1 Inputs','Location', 'SouthEast')
legend('boxoff')
grid off
title('')
xlabel('Distance from start of branch')
ylabel('Cumulative fraction of basal spines')
[hV1DistanceFromSomaBas,pV1DistanceFromSomaBas,~] = kstest2([allSpines(basalV1STEDSpines).distToSoma],[allSpines(basalV1ROIs).distToSoma]);
[hA19DistanceFromSomaBas,pA19DistanceFromSomaBas,~] = kstest2([allSpines(basalA19STEDSpines).distToSoma],[allSpines(basalA19ROIs).distToSoma]);
if hV1DistanceFromSomaBas
    text(5, 0.95,['p = ' num2str(pV1DistanceFromSomaBas,2)], 'Color', fig.cocV1(9,:));
end
if hA19DistanceFromSomaBas
    text(5, 0.90,['p = ' num2str(pA19DistanceFromSomaBas,2)], 'Color', fig.cocA19(9,:));
end
saveas(gcf, fullfile(saveDirImages, '125_CdfDistFromBranchCumulBasalInputs.png'))


%% 130s - Branch order
%Fig 130-132: Histogram of branch order for all inputs, separated for 
%apical/basal
%Fig 133-135: Cdfplot of branch order for all inputs, separated for 
%apical/basal

setFigProperties(130, fig)
binEdges = linspace(1,max([allSpines(STEDROIs).BranchOrder]), max([allSpines(STEDROIs).BranchOrder]));
subplot(1,3,1)
histogram([allSpines(STEDROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocAll(7,:))
title('All spines')
box off
subplot(1,3,2)
histogram([allSpines(A19ROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocA19(7,:))
title('A19 inputs')
box off
subplot(1,3,3)
histogram([allSpines(V1ROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocV1(7,:))
title('V1 inputs')
box off
xlabel('Branch order')
ylabel('Percentage of inputs')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '130_BranchOrderInputs.png'))

setFigProperties(131, fig)
binEdges = linspace(1,max([allSpines(apicalSTEDROIs).BranchOrder]), max([allSpines(apicalSTEDROIs).BranchOrder]));
subplot(1,3,1)
histogram([allSpines(apicalSTEDROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocAll(7,:))
title('All spines')
box off
subplot(1,3,2)
histogram([allSpines(apicalA19ROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocA19(7,:))
title('A19 inputs')
box off
subplot(1,3,3)
histogram([allSpines(apicalV1ROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocV1(7,:))
title('V1 inputs')
box off
xlabel('Branch order')
ylabel('Percentage of apical inputs')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '131_BranchOrderApicalInputs.png'))

setFigProperties(132, fig)
binEdges = linspace(1,max([allSpines(basalSTEDROIs).BranchOrder]), max([allSpines(basalSTEDROIs).BranchOrder]));
subplot(1,3,1)
histogram([allSpines(basalSTEDROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocAll(7,:))
title('All spines')
box off
subplot(1,3,2)
histogram([allSpines(basalA19ROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocA19(7,:))
title('A19 inputs')
box off
subplot(1,3,3)
histogram([allSpines(basalV1ROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', fig.cocV1(7,:))
title('V1 inputs')
box off
xlabel('Branch order')
ylabel('Percentage of basal inputs')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '132_BranchOrderBasalInputs.png'))

%Fig 133: Cdfplot of branch order for all inputs and separated for 
%apical/basal
setFigProperties(123, fig)
h(1,1) = cdfplot([allSpines(STEDROIs).BranchOrder]); hold on
h(1,2) = cdfplot([allSpines(A19ROIs).BranchOrder]);
h(1,3) = cdfplot([allSpines(V1ROIs).BranchOrder]);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocA19(7,:), 'LineWidth', fig.ln);
set(h(1,3), 'Color', fig.cocV1(7,:), 'LineWidth', fig.ln);
grid off
title('Branch Order (all)')
xticks([1 2 3 4 5])
xlabel('Distance from start of branch')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '133_BranchOrderCumulInputs.png'))

setFigProperties(134, fig)
h(1,1) = cdfplot([allSpines(apicalSTEDROIs).BranchOrder]); hold on
h(1,2) = cdfplot([allSpines(apicalA19ROIs).BranchOrder]);
h(1,3) = cdfplot([allSpines(apicalV1ROIs).BranchOrder]);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocA19(7,:), 'LineWidth', fig.ln);
set(h(1,3), 'Color', fig.cocV1(7,:), 'LineWidth', fig.ln);
grid off
title('Branch Order (apical)')
xticks([1 2 3 4 5])
xlabel('Distance from start of branch')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '134_BranchOrderCumulApicalInputs.png'))

setFigProperties(135, fig)
h(1,1) = cdfplot([allSpines(basalSTEDROIs).BranchOrder]); hold on
h(1,2) = cdfplot([allSpines(basalA19ROIs).BranchOrder]);
h(1,3) = cdfplot([allSpines(basalV1ROIs).BranchOrder]);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocA19(7,:), 'LineWidth', fig.ln);
set(h(1,3), 'Color', fig.cocV1(7,:), 'LineWidth', fig.ln);
grid off
title('Branch Order (basal)')
xticks([1 2 3 4 5])
xlabel('Distance from start of branch')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '135_BranchOrderCumulBasalInputs.png'))

%--------------------------------------------------------------------------
%% 140s:Branch specific analysis
%Fig 140: Number of inputs per dendrite
%Fig 141: Input density per dendrite 

BranchesA19Cells = intersect(find([allBranches.inputType] == 1), find([allBranches.STEDLength]>10));
BranchesV1Cells = intersect(find([allBranches.inputType] == 2), find([allBranches.STEDLength]>10));
BranchesA19Input = intersect( BranchesA19Cells, find([allBranches.numInputs] > 0));
BranchesV1Input = intersect( BranchesV1Cells, find([allBranches.numInputs] > 0));

%Fig 140: Number of inputs per dendrite
setFigProperties(130, fig)
subplot(1,2,1)
histogram(cell2mat({allBranches(BranchesA19Cells).numInputs}), 'FaceColor', fig.cocA19(7,:))
ylabel('Number of dendrites')
xlabel('Number of A19 inputs on dendrite')
box off
subplot(1,2,2)
histogram(cell2mat({allBranches(BranchesV1Cells).numInputs}), 'FaceColor', fig.cocV1(7,:))
ylabel('Number of dendrites')
xlabel('Number of A19 inputs on dendrite')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '140_Inputs per dendrite.png'))

%Fig 141: Input density per dendrite 
setFigProperties(141, fig)
subplot(1,4,1)
distributionPlot(cell2mat({allBranches(BranchesA19Input).InputDensity})', 'color', fig.cocA19(7,:)); hold all
boxplot(cell2mat({allBranches(BranchesA19Input).InputDensity})')
ylim([0 max([allBranches.InputDensity])])
ylabel('Input density (um per input)')
box off
subplot(1,4,2)
distributionPlot(cell2mat({allBranches(BranchesV1Input).InputDensity})', 'color', fig.cocV1(7,:)); hold all
boxplot(cell2mat({allBranches(BranchesV1Input).InputDensity})')
ylim([0 max([allBranches.InputDensity])])
box off

%find branches that have more than one input
multInputBranches = find(cell2mat(cellfun(@(x) x > 1, {allBranches.numInputs}, 'UniformOutput',false)));
BranchesA19MultipleInput = intersect(BranchesA19Input, multInputBranches);
BranchesV1MultipleInput = intersect(BranchesV1Input, multInputBranches);

subplot(1,4,3)
distributionPlot(cell2mat({allBranches(BranchesA19MultipleInput).InputDensity})', 'color', fig.cocA19(7,:)); hold all
boxplot(cell2mat({allBranches(BranchesA19MultipleInput).InputDensity}'))
ylabel('Input density (um per input) with more than 1 input')
ylim([0 max([allBranches.InputDensity])])
box off
subplot(1,4,4)
distributionPlot(cell2mat({allBranches(BranchesV1MultipleInput).InputDensity})', 'color', fig.cocV1(7,:)); hold all
boxplot(cell2mat({allBranches(BranchesV1MultipleInput).InputDensity}'))
ylim([0 max([allBranches.InputDensity])])
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '141_Inputs density.png'))

%Fig 132: Input number vs. STED length 
setFigProperties(142, fig)
subplot(1,2,1)
scatter(cell2mat({allBranches(BranchesA19Input).numInputs}), cell2mat({allBranches(BranchesA19Input).STEDLength}), 'filled', 'MarkerFaceColor',fig.cocA19(7,:));
xlabel('Number of A19 inputs on dendrite')
xlim([0 max([allBranches.numInputs])])
ylabel('Length of dendrite')
ylim([0 max([allBranches.STEDLength])])
subplot(1,2,2)
scatter(cell2mat({allBranches(BranchesV1Input).numInputs}), cell2mat({allBranches(BranchesV1Input).STEDLength}), 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
xlabel('Number of V1 inputs on dendrite')
xlim([0 max([allBranches.numInputs])])
ylim([0 max([allBranches.STEDLength])])
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '142_Number of inputs vs length.png'))

%--------------------------------------------------------------------------
%% 150s - Clustering of inputs
%calculations for A19
numMultiplesA19 = sum([allBranches(BranchesA19MultipleInput).numInputs]);  %how many inputs that have at least one other input on the same dendrite?
%make the arrays
DistToInputA19 = zeros(numMultiplesA19,1);
DistToInputSimA19 = zeros(numMultiplesA19, 10000);
DistCounA19 = 1;

for dendNr = 1:length(BranchesA19MultipleInput)
    inputROIsOnDendrite = find([allBranches(BranchesA19MultipleInput(dendNr)).ROIs.Input] == 1);
    %get distances from soma to measure distance
    allDist = [allBranches(BranchesA19MultipleInput(dendNr)).ROIs(inputROIsOnDendrite).distToSoma];
    for a = 1:length(allDist)
        distToSoma = allDist(a);
        diffDist= abs(allDist -distToSoma);
        diffDist(a) = 10000; %the distance to itself
        DistToInputA19(DistCounA19) = min(diffDist);
        DistCounA19 = DistCounA19 +1;
    end
end

for rep = 1:10000
    DistCounA19 = 1;
    for dendNr = 1:length(BranchesA19MultipleInput)
        numInputOnDendrite = length(find([allBranches(BranchesA19MultipleInput(dendNr)).ROIs.Input] == 1)); %how many inputs on the dendrite?
        allDist = ([allBranches(BranchesA19MultipleInput(dendNr)).ROIs.distToSoma]); %what are the distances to soma of all ROIs on the dendrite
        newlySortedDist = allDist(randperm(length(allBranches(BranchesA19MultipleInput(dendNr)).ROIs))); %let's shuffle them
        allDistRand = newlySortedDist(1:numInputOnDendrite); %the first two are the one of the input
        for a = 1:length(allDistRand)
            distToSomaRand = allDistRand(a);
            diffDistRand= abs(allDistRand -distToSomaRand);
            diffDistRand(a) = 10000;
            DistToInputSimA19(DistCounA19, rep) = min(diffDistRand);
            DistCounA19 = DistCounA19 +1;
        end
    end
end

%calculations for V1
numMultiplesV1 = sum([allBranches(BranchesV1MultipleInput).numInputs]);  %how many inputs that have at least one other input on the same dendrite?
%make the arrays
DistToInputV1 = zeros(numMultiplesV1,1);
DistToInputSimV1 = zeros(numMultiplesV1, 10000);
DistCounV1 = 1;

for dendNr = 1:length(BranchesV1MultipleInput)
    inputROIsOnDendrite = find([allBranches(BranchesV1MultipleInput(dendNr)).ROIs.Input] == 1);
    %get distances from soma to measure distance
    allDist = [allBranches(BranchesV1MultipleInput(dendNr)).ROIs(inputROIsOnDendrite).distToSoma];
    for a = 1:length(allDist)
        distToSoma = allDist(a);
        diffDist= abs(allDist -distToSoma);
        diffDist(a) = 10000; %the distance to itself
        DistToInputV1(DistCounV1) = min(diffDist);
        DistCounV1 = DistCounV1 +1;
    end
end

for rep = 1:10000
    DistCounV1 = 1;
    for dendNr = 1:length(BranchesV1MultipleInput)
        numInputOnDendrite = length(find([allBranches(BranchesV1MultipleInput(dendNr)).ROIs.Input] == 1)); %how many inputs on the dendrite?
        allDist = ([allBranches(BranchesV1MultipleInput(dendNr)).ROIs.distToSoma]); %what are the distances to soma of all ROIs on the dendrite
        newlySortedDist = allDist(randperm(length(allBranches(BranchesV1MultipleInput(dendNr)).ROIs))); %let's shuffle them
        allDistRand = newlySortedDist(1:numInputOnDendrite); %the first two are the one of the input
        for a = 1:length(allDistRand)
            distToSomaRand = allDistRand(a);
            diffDistRand= abs(allDistRand -distToSomaRand);
            diffDistRand(a) = 10000;
            DistToInputSimV1(DistCounV1, rep) = min(diffDistRand);
            DistCounV1 = DistCounV1 +1;
        end
    end
end

%Fig 150: Distance to nearest input on same branch
setFigProperties(150, fig)
subplot(1,2,1)
distributionPlot(DistToInputA19, 'color',fig.cocA19(7,:)); hold all
boxplot(DistToInputA19)
ylabel('Distance to nearest A19 input on same branch')
ylim([0 max(DistToInputA19)])
subplot(1,2,2)
distributionPlot(DistToInputV1, 'color',fig.cocV1(7,:)); hold all
boxplot(DistToInputV1)
ylabel('Distance to nearest V1 input on same branch')
ylim([0 max(DistToInputV1)])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '150_DistanceToNearestInput.png'))

%Fig 151: Mean distance to nearest input
setFigProperties(141, fig)
subplot(1,2,1)
histogram(mean(DistToInputSimA19,1), 200,'FaceColor',fig.cocAll(7,:))
xlabel('Mean nearest-neighbor distance in um')
xline(mean(DistToInputA19), '--', 'Color', fig.cocA19(7,:), 'LineWidth', fig.ln);
percA19 = find(sort(mean(DistToInputSimA19,1)) >= mean(DistToInputA19), 1, 'first') / numel(mean(DistToInputSimA19,1)) * 100;
disp(['Mean distance to nearest other A19 input: ' num2str(percA19,3) ' % of simulated data has a lower mean than the real data'])
legend('Random Data', 'A19')
subplot(1,2,2)
histogram(mean(DistToInputSimV1,1), 200,'FaceColor',fig.cocAll(7,:))
xlabel('Mean nearest-neighbor distance in um')
xline(mean(DistToInputV1), '--', 'Color', fig.cocV1(7,:), 'LineWidth', fig.ln);
percV1 = find(sort(mean(DistToInputSimV1,1)) >= mean(DistToInputV1), 1, 'first') / numel(mean(DistToInputSimV1,1)) * 100;
disp(['Mean distance to nearest other V1 input: ' num2str(percV1,3) ' % of simulated data has a lower mean than the real data'])
legend('Random Data', 'V1')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '151_MeanDistanceToNearestInput.png'))

%Fig 152: Distribution of distance to nearest other input
setFigProperties(152, fig)
subplot(2,1,1)
h(1,1) = cdfplot(DistToInputSimA19(:)); hold on
h(1,2) = cdfplot(DistToInputA19);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocA19(9,:), 'LineWidth', fig.ln);
legend('Random data', 'Real data','Location', 'SouthEast')
grid off
title('A19 Inputs')
xlabel('Distance to nearest other input on same dendrite')
ylabel('')
[hA19NearestNeighborInput,pA19NearestNeighborInput,~] = kstest2(DistToInputSimA19(:),DistToInputA19);
if hA19NearestNeighborInput
    disp(['Significance level for nearest Input neighbor for A19: ' num2str(pA19NearestNeighborInput,3)])
end

subplot(2,1,2)
h(1,1) = cdfplot(DistToInputSimV1(:)); hold on
h(1,2) = cdfplot(DistToInputV1);
set(h(1,1), 'Color', fig.cocAll(7,:), 'LineWidth', fig.ln);
set(h(1,2), 'Color', fig.cocV1(9,:), 'LineWidth', fig.ln);
legend('Random data', 'Real data','Location', 'SouthEast')
grid off
title('V1 Inputs')
xlabel('Distance to nearest other input on same dendrite')
ylabel('')
[hV1NearestNeighborInput,pV1NearestNeighborInput,~] = kstest2(DistToInputSimA19(:),DistToInputA19);
if hV1NearestNeighborInput
    disp(['Significance level for nearest Input neighbor for V1: ' num2str(pV1NearestNeighborInput,3)])
end
saveas(gcf, fullfile(saveDirImages, '152_DistributionDistanceToNearestInput.png'))

%--------------------------------------------------------------------------
%% FUNCTION

%% 200s - Plotting on top of constructed cell
%Fig 200: Build a constructed cell A19 delta Ori
setFigProperties(200, fig)
constructInputCell(allSpines, allBranches, 'A19', 'deltaOri')
saveas(gcf, fullfile(saveDirImages,'200_constructedInputCell_A19_deltaOri.png'))

%Fig 201: Build a constructed cell V1 delta Ori
setFigProperties(201, fig)
constructInputCell(allSpines, allBranches, 'V1', 'deltaOri')
saveas(gcf, fullfile(saveDirImages,'201_constructedInputCell_V1_deltaOri.png'))

%Fig 202: Build a constructed cell A19 delta Dir
setFigProperties(202, fig)
constructInputCell(allSpines, allBranches, 'A19', 'deltaDir')
saveas(gcf, fullfile(saveDirImages,'202_constructedInputCell_A19_deltaDir.png'))

%Fig 203: Build a constructed cell V1 delta Ori
setFigProperties(203, fig)
constructInputCell(allSpines, allBranches, 'V1', 'deltaDir')
saveas(gcf, fullfile(saveDirImages,'203_constructedInputCell_V1_deltaDir.png'))

%Fig 204: Build a constructed cell A19 OSI
setFigProperties(204, fig)
constructInputCell(allSpines, allBranches, 'A19', 'OSI')
saveas(gcf, fullfile(saveDirImages,'204_constructedInputCell_A19_OSI.png'))

%Fig 205: Build a constructed cell V1 OSI
setFigProperties(205, fig)
constructInputCell(allSpines, allBranches, 'V1', 'OSI')
saveas(gcf, fullfile(saveDirImages,'205_constructedInputCell_V1_OSI.png'))

%Fig 206: Build a constructed cell A19 OSI
setFigProperties(206, fig)
constructInputCell(allSpines, allBranches, 'A19', 'DSI')
saveas(gcf, fullfile(saveDirImages,'206_constructedInputCell_A19_DSI.png'))

%Fig 207: Build a constructed cell V1 OSI
setFigProperties(207, fig)
constructInputCell(allSpines, allBranches, 'V1', 'DSI')
saveas(gcf, fullfile(saveDirImages,'207_constructedInputCell_V1_DSI.png'))

%% 210s - Quick overview of general properties

%210 - Pie chart of fraction of responsive & ori-select
setFigProperties(210, fig)
subplot(1,3,1)
h = pie([length(oriGood)/length(TwoPROIsNR), length(goodROIs)/length(TwoPROIsNR)-length(oriGood)/length(TwoPROIsNR),1-length(goodROIs)/length(TwoPROIsNR)]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', fig.cocAll(9,:));
set(hp(2), 'FaceColor', fig.cocAll(7,:));
set(hp(3), 'FaceColor', fig.cocAll(4,:));
title('All spines')
legend({'Ori-select', 'Responsive', 'Non-responsive'}, 'Location', 'southoutside')
legend('boxoff')
subplot(1,3,2)
h = pie([length(oriGoodA19)/length(A19Funct), length(goodA19)/length(A19Funct)-length(oriGoodA19)/length(A19Funct),1-length(goodA19)/length(A19Funct)]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', fig.cocA19(9,:));
set(hp(2), 'FaceColor', fig.cocA19(7,:));
set(hp(3), 'FaceColor', fig.cocA19(4,:));
title('A19 inputs')
legend({'Ori-select', 'Responsive', 'Non-responsive'}, 'Location', 'southoutside')
legend('boxoff')
subplot(1,3,3)
h = pie([length(oriGoodV1)/length(V1Funct), length(goodV1)/length(V1Funct)-length(oriGoodV1)/length(V1Funct),1-length(goodV1)/length(V1Funct)]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', fig.cocV1(9,:));
set(hp(2), 'FaceColor', fig.cocV1(7,:));
set(hp(3), 'FaceColor', fig.cocV1(4,:));
title('V1 inputs')
legend({'Ori-select', 'Responsive', 'Non-responsive'}, 'Location', 'southoutside')
legend('boxoff')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'210_Fraction of responsive and ori-select.png'))

% Chi2 test for independendce
% Step 0: Get numbers for the individual fractions
noInput = [repmat('N', length(TwoPROIsNR)-length(goodROIs), 1); repmat('G', length(goodROIs)-length(oriGood),1); repmat('O', length(oriGood),1)];
A19Input = [repmat('N', length(A19Funct)-length(goodA19), 1); repmat('G', length(goodA19)-length(oriGoodA19),1); repmat('O', length(oriGoodA19),1)];
V1Input = [repmat('N', length(V1Funct)-length(goodV1), 1); repmat('G', length(goodV1)-length(oriGoodV1),1); repmat('O', length(oriGoodV1),1)];

% Step 1: do all vs. A19 
typeOfInputAllA19 = [repmat('N', length(TwoPROIsNR), 1); repmat('A', length(A19Funct),1)];
typeOfRespAllA19 = [noInput; A19Input];
[~,~,pvalAllA19] = crosstab(typeOfInputAllA19,typeOfRespAllA19);

% Step 2:  all vs. v1
typeOfInputAllV1 = [repmat('N', length(TwoPROIsNR), 1); repmat('V', length(V1Funct),1)];
typeOfRespAllV1 = [noInput; V1Input];
[~,~,pvalAllV1] = crosstab(typeOfInputAllV1,typeOfRespAllV1);

% Step 3: a19 vs. v1
typeOfInputA19V1 = [repmat('A', length(A19Funct),1); repmat('V', length(V1Funct),1)];
typeOfRespA19V1 = [A19Input;V1Input];
[~,~,pvalA19V1] = crosstab(typeOfInputA19V1,typeOfRespA19V1);

% Step 4: Do Bonferoni correction for 2 comparisons
numPopulations = 3;
pvalAllA19 = pvalAllA19/numPopulations;
pvalAllV1 = pvalAllV1/numPopulations;

%5) Distribution of responsiveness/ori select
pvalA19V1 = pvalA19V1/numPopulations;

% Step 5: Display result
disp(['P value for responder type All vs. A19: ' num2str(pvalAllA19)]);
disp(['P value for responder type All vs. V1: ' num2str(pvalAllV1)]);
disp(['P value for responder type A19 vs. V1: ' num2str(pvalA19V1)]);

%211 - Distibution of delta Ori 
getDeltaOri = @(x) x.funcData.deltaOri;
goodDeltaOri = cell2mat(arrayfun(getDeltaOri, allSpines(oriGood), 'UniformOutput', false))';
A19DeltaOri = cell2mat(arrayfun(getDeltaOri, allSpines(oriGoodA19), 'UniformOutput', false))';
V1DeltaOri = cell2mat(arrayfun(getDeltaOri, allSpines(oriGoodV1), 'UniformOutput', false))';
edges = linspace(0, 90, 10);

setFigProperties(211, fig)
plot(histcounts(V1DeltaOri,edges,'Normalization', 'probability'), 'color', fig.cocV1(9,:), 'LineWidth', fig.ln)
hold on
plot(histcounts(A19DeltaOri,edges,'Normalization', 'probability'), 'color', fig.cocA19(9,:), 'LineWidth', fig.ln)
plot(histcounts(goodDeltaOri,edges,'Normalization', 'probability'), 'color',fig.cocAll(9,:), 'LineWidth', fig.ln)
xticks(linspace(0, 9, 10))
xticklabels(edges)
xlabel('deltaOri soma - spine')
ylabel('Probability')
legend(['all Spines', 'V1 inputs', 'A19 inputs'])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'211_Delta Ori probability.png'))

goodDeltaOriA19Cells = cell2mat(arrayfun(getDeltaOri, allSpines(intersect(oriGood,find([allSpines.inputType] == 1))), 'UniformOutput', false))';
goodDeltaOriV1Cells = cell2mat(arrayfun(getDeltaOri, allSpines(intersect(oriGood,find([allSpines.inputType] == 2))), 'UniformOutput', false))';
edges = linspace(0, 90, 10);

setFigProperties(2110, fig)
plot(histcounts(goodDeltaOri,edges,'Normalization', 'probability'), 'color', fig.cocAll(9,:), 'LineWidth', fig.ln)
hold on
plot(histcounts(goodDeltaOriA19Cells,edges,'Normalization', 'probability'), 'color', fig.cocAll(5,:), 'LineWidth', fig.ln)
plot(histcounts(goodDeltaOriV1Cells,edges,'Normalization', 'probability'), 'color',fig.cocAll(7,:), 'LineWidth', fig.ln)
plot(histcounts(A19DeltaOri,edges,'Normalization', 'probability'), 'color', fig.cocA19(9,:), 'LineWidth', fig.ln)
plot(histcounts(V1DeltaOri,edges,'Normalization', 'probability'), 'color', fig.cocV1(9,:), 'LineWidth', fig.ln)
xticks(linspace(0, 9, 10))
xticklabels(edges)
xlabel('deltaOri soma - spine')
ylabel('Probability')
%legend(['all Spines', 'V1 input cell spines', 'A19 input cell spines'])
set(gcf, 'color', 'w');
box off
saveas(gcf, fullfile(saveDirImages,'2110_Delta Ori probability_with cell control.png'))

% Chi2 test for independendce of the first datapoint 
% Step 0: Get counts first and make vectors
AllCounts = histcounts(goodDeltaOri,edges,'Normalization', 'count');
A19Counts = histcounts(A19DeltaOri,edges,'Normalization', 'count');
V1Counts = histcounts(V1DeltaOri,edges,'Normalization', 'count');
AllFirstVsAll = [repmat('1',AllCounts(1),1); repmat('0',sum(AllCounts)-AllCounts(1),1)];
A19FirstVsAll = [repmat('1',A19Counts(1),1); repmat('0',sum(A19Counts)-A19Counts(1),1)];
V1FirstVsAll = [repmat('1',V1Counts(1),1); repmat('0',sum(V1Counts)-V1Counts(1),1)];

% Step 1: do all vs. A19 
typeOfInputAllA19 = [repmat('N', sum(AllCounts), 1); repmat('A', sum(A19Counts),1)];
typeOfOriDeltaAllA19 = [AllFirstVsAll; A19FirstVsAll];
[~,~,pvalDeltaOriAllA19] = crosstab(typeOfInputAllA19,typeOfOriDeltaAllA19);

% Step 2:  all vs. v1
typeOfInputAllV1 = [repmat('N', sum(AllCounts), 1); repmat('V', sum(V1Counts),1)];
typeOfOriDeltaAllV1 = [AllFirstVsAll; V1FirstVsAll];
[~,~,pvalDeltaOriAllV1] = crosstab(typeOfInputAllV1,typeOfOriDeltaAllV1);

% Step 3: a19 vs. v1
typeOfInputA19V1 = [repmat('A', sum(A19Counts),1); repmat('V', sum(V1Counts),1)];
typeOfOriDeltaA19V1 = [A19FirstVsAll;V1FirstVsAll];
[~,~,pvalDeltaOriA19V1] = crosstab(typeOfInputA19V1,typeOfOriDeltaA19V1);

% Step 4: Do Bonferoni correction for 2 comparisons
numPopulations = 3;
pvalDeltaOriAllA19 = pvalDeltaOriAllA19/numPopulations;
pvalDeltaOriAllV1 = pvalDeltaOriAllV1/numPopulations;
pvalDeltaOriA19V1 = pvalDeltaOriA19V1/numPopulations;

% Step 5: Display result
disp(['P value deltaOri All vs. A19: ' num2str(pvalDeltaOriAllA19)]);
disp(['P value deltaOri All vs. V1: ' num2str(pvalDeltaOriAllV1)]);
disp(['P value deltaOri A19 vs. V1: ' num2str(pvalDeltaOriA19V1)]);

%212 - Distribution of deltaDir
getDeltaDir = @(x) x.funcData.deltaDir;
goodDeltaDir = cell2mat(arrayfun(getDeltaDir, allSpines(dirGood), 'UniformOutput', false))';
A19DeltaDir = cell2mat(arrayfun(getDeltaDir, allSpines(dirGoodA19), 'UniformOutput', false))';
V1DeltaDir = cell2mat(arrayfun(getDeltaDir, allSpines(dirGoodV1), 'UniformOutput', false))';
edges = linspace(0, 180, 20);

setFigProperties(212, fig)
plot(histcounts(V1DeltaDir,edges,'Normalization', 'probability'), 'color', fig.cocV1(9,:), 'LineWidth', fig.ln)
hold on
plot(histcounts(A19DeltaDir,edges,'Normalization', 'probability'), 'color', fig.cocA19(9,:), 'LineWidth', fig.ln)
plot(histcounts(goodDeltaDir,edges,'Normalization', 'probability'), 'color',fig.cocAll(9,:), 'LineWidth', fig.ln)
xticks(linspace(0, 180, 9))
xticklabels(edges)
xlabel('deltaDir soma - spine')
ylabel('Probability')
legend(['all Spines', 'V1 inputs', 'A19 inputs'])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'212_Delta Dir probability.png'))

%213 - Distribution of OSI
getOSI = @(x) x.OSI;
goodOSI = cell2mat(arrayfun(getOSI, allSpines(goodROIs), 'UniformOutput', false))';
A19OSI = cell2mat(arrayfun(getOSI, allSpines(goodA19), 'UniformOutput', false))';
V1OSI = cell2mat(arrayfun(getOSI, allSpines(goodV1), 'UniformOutput', false))';

setFigProperties(213, fig)
subplot(1,3,1)
distributionPlot(goodOSI,'color', fig.cocAll(9,:)); hold on
boxplot(goodOSI','Label', {'All spines'})
ylim([0 1])
ylabel('OSI')
set(gca,'Box','off');
subplot(1,3,2)
distributionPlot(A19OSI,'color', fig.cocA19(9,:)); hold on
boxplot(A19OSI','Label', {'A19 Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
subplot(1,3,3)
distributionPlot(V1OSI,'color', fig.cocV1(9,:)); hold on
boxplot(V1OSI','Label', {'V1 Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'213_OSI distribution.png'))

%let's do statistics
maxNGood = max([length(goodV1) length(goodA19) length(goodROIs)]);
OSIall = NaN(maxNGood,3);
OSIall(1:1:length(goodROIs),1) = goodOSI;
OSIall(1:1:length(goodA19),2) = A19OSI;
OSIall(1:1:length(goodV1),3) = V1OSI;
group = {'all', 'A19 inputs', 'V1 inputs'};
[pOSI,~,statsOSI] = kruskalwallis(OSIall, group, 'off');

if pOSI < 0.05
    OSIStatsMulti = multcompare(statsOSI, 'CType','bonferroni', 'display', 'off');
    OSI_AllvsA19 = OSIStatsMulti(1,6);
    OSI_AllvsV1 = OSIStatsMulti(2,6);
    OSI_A19vsV1 = OSIStatsMulti(3,6);
    disp(['OSI: All vs. A19: p = ' num2str(OSI_AllvsA19) ', All vs. V1: p = ' num2str(OSI_AllvsV1) ', A19 vs. V1: p = ' num2str(OSI_A19vsV1)]);
    close gcf
else
    disp('No significant difference in OSI distribution')
end

%214 - Distribution of DSI
getDSI = @(x) x.DSIvect;
goodDSI = cell2mat(arrayfun(getDSI, allSpines(goodROIs), 'UniformOutput', false))';
A19DSI = cell2mat(arrayfun(getDSI, allSpines(goodA19), 'UniformOutput', false))';
V1DSI = cell2mat(arrayfun(getDSI, allSpines(goodV1), 'UniformOutput', false))';

setFigProperties(214, fig)
subplot(1,3,1)
distributionPlot(goodDSI,'color', fig.cocAll(9,:)); hold on
boxplot(goodDSI','Label', {'All spines'})
ylim([0 1])
ylabel('DSI')
set(gca,'Box','off');
subplot(1,3,2)
distributionPlot(A19DSI,'color', fig.cocA19(9,:)); hold on
boxplot(A19DSI','Label', {'A19 Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
subplot(1,3,3)
distributionPlot(V1DSI,'color', fig.cocV1(9,:)); hold on
boxplot(V1DSI','Label', {'Recurrent'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'214_DSI distribution.png'))

%let's do statistics
DSIall = NaN(maxNGood,3);
DSIall(1:1:length(goodROIs),1) = goodDSI;
DSIall(1:1:length(goodA19),2) = A19DSI;
DSIall(1:1:length(goodV1),3) = V1DSI;
[pDSI,~,statsDSI] = kruskalwallis(DSIall, group, 'off');

if pDSI < 0.05
    DSIStatsMulti = multcompare(statsDSI, 'CType','bonferroni', 'display', 'off');
    DSI_AllvsA19 = DSIStatsMulti(1,6);
    DSI_AllvsV1 = DSIStatsMulti(2,6);
    DSI_A19vsV1 = DSIStatsMulti(3,6);
    disp(['DSI: All vs. A19: p = ' num2str(DSI_AllvsA19) ', All vs. V1: p = ' num2str(DSI_AllvsV1) ', A19 vs. V1: p = ' num2str(DSI_A19vsV1)]);
    close gcf
else
    disp('No significant difference in DSI distribution')
end

%215 - distribution of bandwith
getBW = @(x) x.Bandwidth;
goodBW = cell2mat(arrayfun(getBW, allSpines(goodROIs), 'UniformOutput', false))';
A19BW = cell2mat(arrayfun(getBW, allSpines(goodA19), 'UniformOutput', false))';
V1BW = cell2mat(arrayfun(getBW, allSpines(goodV1), 'UniformOutput', false))';

setFigProperties(215, fig)
subplot(1,3,1)
distributionPlot(goodBW,'color', fig.cocAll(9,:)); hold on
boxplot(goodBW','Label', {'All spines'})
ylim([0 180])
ylabel('Bandwidth')
set(gca,'Box','off');
subplot(1,3,2)
distributionPlot(A19BW,'color', fig.cocA19(9,:)); hold on
boxplot(A19BW','Label', {'A19 Inputs'})
ylim([0 180])
set(gca,'box','off','ycolor','w')
subplot(1,3,3)
distributionPlot(V1BW,'color', fig.cocV1(9,:)); hold on
boxplot(V1BW','Label', {'Recurrent'})
ylim([0 180])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'215_Bandwidth distribution.png'))

%let's do statistics
BWall = NaN(maxNGood,3);
BWall(1:1:length(goodROIs),1) = goodBW;
BWall(1:1:length(goodA19),2) = A19BW;
BWall(1:1:length(goodV1),3) = V1BW;
[pBW,~,statsBW] = kruskalwallis(BWall, group, 'off');

if pBW < 0.05
    BWStatsMulti = multcompare(statsBW, 'CType','bonferroni', 'display', 'off');
    BW_AllvsA19 = BWStatsMulti(1,6);
    BW_AllvsV1 = BWStatsMulti(2,6);
    BW_A19vsV1 = BWStatsMulti(3,6);
    disp(['BW: All vs. A19: p = ' num2str(BW_AllvsA19) ', All vs. V1: p = ' num2str(BW_AllvsV1) ', A19 vs. V1: p = ' num2str(BW_A19vsV1)]);
    close gcf
else
    disp('No significant difference in BW distribution')
end

%--------------------------------------------------------------------------
%% 220s & 30s: Local environment ori(220s) and dir (230s)
%Fig 220-222: local Dispersion
localDispValues = cell2mat(cellfun(@(x) x.localDispersion, num2cell(allSpines(oriGood)), 'UniformOutput', false));
localDispValues = reshape(localDispValues, 4, [])';
localDispValuesA19 = cell2mat(cellfun(@(x) x.localDispersion, num2cell(allSpines(oriGoodA19)), 'UniformOutput', false));
localDispValuesA19 = reshape(localDispValuesA19, 4, [])';
localDispValuesV1 = cell2mat(cellfun(@(x) x.localDispersion, num2cell(allSpines(oriGoodV1)), 'UniformOutput', false));
localDispValuesV1 = reshape(localDispValuesV1, 4, [])';
localDispValuesLowDeltaOri = cell2mat(cellfun(@(x) x.localDispersion, num2cell(allSpines(oriGood(goodDeltaOri < 20))), 'UniformOutput', false));
localDispValuesLowDeltaOri = reshape(localDispValuesLowDeltaOri, 4, [])';
localDispValuesHighDeltaOri = cell2mat(cellfun(@(x) x.localDispersion, num2cell(allSpines(oriGood(goodDeltaOri > 20))), 'UniformOutput', false));
localDispValuesHighDeltaOri = reshape(localDispValuesHighDeltaOri, 4, [])';

localDeltaOriValues = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(allSpines(oriGood)), 'UniformOutput', false));
localDeltaOriValues = reshape(localDeltaOriValues, 4, [])';
localDeltaOriValuesA19 = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(allSpines(oriGoodA19)), 'UniformOutput', false));
localDeltaOriValuesA19 = reshape(localDeltaOriValuesA19, 4, [])';
localDeltaOriValuesV1 = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(allSpines(oriGoodV1)), 'UniformOutput', false));
localDeltaOriValuesV1 = reshape(localDeltaOriValuesV1, 4, [])';
localDeltaOriValuesLowDeltaOri = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(allSpines(oriGood(goodDeltaOri < 20))), 'UniformOutput', false));
localDeltaOriValuesLowDeltaOri = reshape(localDeltaOriValuesLowDeltaOri, 4, [])';
localDeltaOriValuesHighDeltaOri = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(allSpines(oriGood(goodDeltaOri > 20))), 'UniformOutput', false));
localDeltaOriValuesHighDeltaOri = reshape(localDeltaOriValuesHighDeltaOri, 4, [])';

%shuffle the identity of inputs within a branch to see if it changes the
%localDeltaOri and localDisp

[localDeltaOriA19Random, localDeltaOriSomaA19Random, localDispValuesA19Random] = randomizeInputPositionForlocalMeasurements(allSpines, allBranches, cellDetails, 'A19');
[localDeltaOriV1Random, localDeltaOriSomaV1Random, localDispValuesV1Random] = randomizeInputPositionForlocalMeasurements(allSpines, allBranches, cellDetails, 'V1');

setFigProperties(220, fig)
subplot(1,5,1)
boxplot(localDispValues(:,2)','Label', {'All spines'},'color',fig.cocAll(7,:))
ylim([0 45])
hold on
scatter(repmat(1:1,size(localDispValues,1),1), localDispValues(:,2), 'filled','MarkerFaceColor',fig.cocAll(7,:), 'MarkerFaceAlpha',0.05','jitter','on','jitterAmount',0.15)
box off
ylabel('Local dispersion')
subplot(1,5,2)
boxplot(localDispValuesA19(:,2)','Label', {'A19 Inputs'},'color',fig.cocA19(7,:))
ylim([0 45])
hold on
scatter(repmat(1:4,size(localDispValuesA19,1),1), localDispValuesA19(:,2), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,3)
boxplot(localDispValuesV1(:,2)','Label', {'V1 Inputs'},'color',fig.cocV1(7,:))
ylim([0 45])
hold on
scatter(repmat(1:1,size(localDispValuesV1,1),1), localDispValuesV1(:,2), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,4)
boxplot(localDispValuesLowDeltaOri(:,1),'color',fig.cocAll(9,:)); hold on
scatter(repmat(1:1,size(localDispValuesLowDeltaOri(:,1),1),1), localDispValuesLowDeltaOri(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
xticklabels({'deltaOri < 20'})
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,5)
boxplot(localDispValuesHighDeltaOri(:,1),'color',fig.cocAll(9,:)); hold on
scatter(repmat(1:1,size(localDispValuesHighDeltaOri(:,1),1),1), localDispValuesHighDeltaOri(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
xticklabels({'deltaOri > 20'})
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
saveas(gcf, fullfile(saveDirImages,'220_Local dispersion_comparison.png'))

setFigProperties(2201, fig)
subplot(1,5,1)
boxplot(localDispValues(:,2)','Label', {'All spines'},'color',fig.cocAll(7,:))
ylim([0 45])
hold on
scatter(repmat(1:1,size(localDispValues,1),1), localDispValues(:,2), 'filled','MarkerFaceColor',fig.cocAll(7,:), 'MarkerFaceAlpha',0.05','jitter','on','jitterAmount',0.15)
box off
ylabel('Local dispersion')
subplot(1,5,2)
boxplot(localDispValuesA19(:,2)','Label', {'A19 Inputs'},'color',fig.cocA19(7,:))
ylim([0 45])
hold on
scatter(repmat(1:4,size(localDispValuesA19,1),1), localDispValuesA19(:,2), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,4)
boxplot(localDispValuesV1(:,2)','Label', {'V1 Inputs'},'color',fig.cocV1(7,:))
ylim([0 45])
hold on
scatter(repmat(1:1,size(localDispValuesV1,1),1), localDispValuesV1(:,2), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,3)
boxplot(localDispValuesA19Random,'color',fig.cocAll(9,:)); hold on
xticklabels({'A19 random'})
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,5)
boxplot(localDispValuesV1Random,'color',fig.cocAll(9,:)); hold on
xticklabels({'V1 random'})
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
saveas(gcf, fullfile(saveDirImages,'2201_Local dispersion_randomComparison.png'))

setFigProperties(221, fig);
subplot(1,3,1)
boxplot(localDispValues,'color',fig.cocAll(7,:))
ylim([0 45])
hold on
scatter(repmat(1:4,size(localDispValues,1),1), localDispValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('Local dispersion')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDispValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 45])
hold on
scatter(repmat(1:4,size(localDispValuesA19,1),1), localDispValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDispValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 45])
hold on
scatter(repmat(1:4,size(localDispValuesV1,1),1), localDispValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'221_All local dispersion.png'))

%Fig 222: local Dispersion vs deltaOri
setFigProperties(222, fig);
scatter(goodDeltaOri,localDispValues(:,1),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaOri,localDispValuesA19(:,1),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaOri,localDispValuesV1(:,1),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
xlim([0 90])
ylim([0 45])
ylabel(sprintf('Local dispersion within (5 \\mum) of spine'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'222_Local dispersion vs delta Ori.png'))

%Fig 223: HI 
HIValues = cell2mat(cellfun(@(x) x.HI, num2cell(allSpines(oriGood)), 'UniformOutput', false));
HIValues = reshape(HIValues, 4, [])';
HIValuesA19 = cell2mat(cellfun(@(x) x.HI, num2cell(allSpines(oriGoodA19)), 'UniformOutput', false));
HIValuesA19 = reshape(HIValuesA19, 4, [])';
HIValuesV1 = cell2mat(cellfun(@(x) x.HI, num2cell(allSpines(oriGoodV1)), 'UniformOutput', false));
HIValuesV1= reshape(HIValuesV1, 4, [])';

setFigProperties(223, fig)
scatter(goodDeltaOri,HIValues(:,2),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaOri,HIValuesA19(:,2),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaOri,HIValuesV1(:,2),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
ylim([0 1])
xlim([0 90])
ylabel(sprintf('HI of local environment (5 \\mum) of spine'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'223_HI vs delta Ori.png'))

%Fig 224: local deltaOri
setFigProperties(224, fig)
subplot(1,3,1)
boxplot(localDeltaOriValues,'color',fig.cocAll(7,:))
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDeltaOriValues,1),1), localDeltaOriValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('deltaOri of local environment')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDeltaOriValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDeltaOriValuesA19,1),1), localDeltaOriValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:),'jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDeltaOriValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDeltaOriValuesV1,1),1), localDeltaOriValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'224_Local deltaOri.png'))

setFigProperties(2240, fig)
subplot(1,5,1)
boxplot(localDeltaOriValues(:,1),'color',fig.cocAll(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValues(:,1),1),1), localDeltaOriValues(:,1), 'filled','MarkerFaceColor',fig.cocAll(7,:), 'MarkerFaceAlpha',0.05','jitter','on','jitterAmount',0.15)
box off
xticklabels({'All inputs'})
ylabel('deltaOri of local environment')
ylim([0 90])
subplot(1,5,2)
boxplot(localDeltaOriValuesA19(:,1),'color',fig.cocA19(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValuesA19(:,1),1),1), localDeltaOriValuesA19(:,1), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'A19 inputs'})
ylim([0 90])
subplot(1,5,3)
boxplot(localDeltaOriValuesV1(:,1),'color',fig.cocV1(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValuesV1(:,1),1),1), localDeltaOriValuesV1(:,1), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'V1 inputs'})
ylim([0 90])
subplot(1,5,4)
boxplot(localDeltaOriValuesLowDeltaOri(:,1),'color',fig.cocAll(9,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValuesLowDeltaOri(:,1),1),1), localDeltaOriValuesLowDeltaOri(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
xticklabels({'deltaOri < 20'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,5)
boxplot(localDeltaOriValuesHighDeltaOri(:,1),'color',fig.cocAll(9,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValuesHighDeltaOri(:,1),1),1), localDeltaOriValuesHighDeltaOri(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
xticklabels({'deltaOri > 20'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
saveas(gcf, fullfile(saveDirImages,'2240_Local deltaOri_2.5um radius comparison.png'))

setFigProperties(2241, fig)
subplot(1,5,1)
boxplot(localDeltaOriValues(:,1),'color',fig.cocAll(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValues(:,1),1),1), localDeltaOriValues(:,1), 'filled','MarkerFaceColor',fig.cocAll(7,:), 'MarkerFaceAlpha',0.05','jitter','on','jitterAmount',0.15)
box off
xticklabels({'All inputs'})
ylabel('deltaOri of local environment')
ylim([0 90])
subplot(1,5,2)
boxplot(localDeltaOriValuesA19(:,1),'color',fig.cocA19(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValuesA19(:,1),1),1), localDeltaOriValuesA19(:,1), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'A19 inputs'})
ylim([0 90])
subplot(1,5,3)
boxplot(localDeltaOriA19Random,'color',fig.cocAll(7,:)); hold on
%scatter(repmat(1:1,size(localDeltaOriA19Random,2),1), localDeltaOriA19Random, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.05', 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'A19 random'})
ylim([0 90])
subplot(1,5,4)
boxplot(localDeltaOriValuesV1(:,1),'color',fig.cocAll(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriValuesV1(:,1),1),1), localDeltaOriValuesV1(:,1), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'V1 inputs'})
ylim([0 90])
subplot(1,5,5)
boxplot(localDeltaOriV1Random,'color',fig.cocAll(7,:)); hold on
%scatter(repmat(1:1,size(localDeltaOriV1Random,1),1), localDeltaOriV1Random, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.05', 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'V1 random'})
ylim([0 90])
saveas(gcf, fullfile(saveDirImages,'2240_Local deltaOri_compareRandom.png'))

%Fig 225: Local DeltaOriSoma
localDeltaOriSomaValues = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(oriGood)), 'UniformOutput', false));
localDeltaOriSomaValues = reshape(localDeltaOriSomaValues, 4, [])';
localDeltaOriSomaValuesA19 = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(oriGoodA19)), 'UniformOutput', false));
localDeltaOriSomaValuesA19 = reshape(localDeltaOriSomaValuesA19, 4, [])';
localDeltaOriSomaValuesV1 = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(oriGoodV1)), 'UniformOutput', false));
localDeltaOriSomaValuesV1 = reshape(localDeltaOriSomaValuesV1, 4, [])';
localDeltaOriSomaValuesLowDeltaOri = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(oriGood(goodDeltaOri < 20))), 'UniformOutput', false));
localDeltaOriSomaValuesLowDeltaOri = reshape(localDeltaOriSomaValuesLowDeltaOri, 4, [])';
localDeltaOriSomaValuesHighDeltaOri = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(oriGood(goodDeltaOri > 20))), 'UniformOutput', false));
localDeltaOriSomaValuesHighDeltaOri = reshape(localDeltaOriSomaValuesHighDeltaOri, 4, [])';

setFigProperties(225, fig)
subplot(1,3,1)
boxplot(localDeltaOriSomaValues,'color',fig.cocAll(7,:))
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDeltaOriSomaValues,1),1), localDeltaOriSomaValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('deltaOriSoma of local environment')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDeltaOriSomaValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDeltaOriSomaValuesA19,1),1), localDeltaOriSomaValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDeltaOriSomaValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDeltaOriSomaValuesV1,1),1), localDeltaOriSomaValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'225_Local deltaOriSoma.png'))

setFigProperties(2250, fig)
subplot(1,5,1)
boxplot(localDeltaOriSomaValues(:,1),'color',fig.cocAll(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValues(:,1),1),1), localDeltaOriSomaValues(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
box off
xticklabels({'All inputs'})
ylabel('deltaOriSoma of local environment')
ylim([0 90])
subplot(1,5,2)
boxplot(localDeltaOriSomaValuesA19(:,1),'color',fig.cocA19(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValuesA19(:,1),1),1), localDeltaOriSomaValuesA19(:,1), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'jitter','on','jitterAmount',0.15)
xticklabels({'A19 inputs'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,3)
boxplot(localDeltaOriSomaValuesV1(:,1),'color',fig.cocV1(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValuesV1(:,1),1),1), localDeltaOriSomaValuesV1(:,1), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'V1 inputs'})
ylim([0 90])
subplot(1,5,4)
boxplot(localDeltaOriSomaValuesLowDeltaOri(:,1),'color',fig.cocAll(9,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValuesLowDeltaOri(:,1),1),1), localDeltaOriSomaValuesLowDeltaOri(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
xticklabels({'deltaOri < 20'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,5)
boxplot(localDeltaOriSomaValuesHighDeltaOri(:,1),'color',fig.cocAll(9,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValuesHighDeltaOri(:,1),1),1), localDeltaOriSomaValuesHighDeltaOri(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
xticklabels({'deltaOri > 20'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
saveas(gcf, fullfile(saveDirImages,'2250_Local deltaOriSoma_2.5um radius comparison.png'))

setFigProperties(2251, fig)
subplot(1,5,1)
boxplot(localDeltaOriSomaValues(:,1),'color',fig.cocAll(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValues(:,1),1),1), localDeltaOriSomaValues(:,1), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.07','jitter','on','jitterAmount',0.15)
box off
xticklabels({'All inputs'})
ylabel('deltaOriSoma of local environment')
ylim([0 90])
subplot(1,5,2)
boxplot(localDeltaOriSomaValuesA19(:,1),'color',fig.cocA19(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValuesA19(:,1),1),1), localDeltaOriSomaValuesA19(:,1), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'jitter','on','jitterAmount',0.15)
xticklabels({'A19 inputs'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
subplot(1,5,3)
boxplot(localDeltaOriSomaA19Random,'color',fig.cocAll(9,:)); hold on
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'random A19'})
ylim([0 90])
subplot(1,5,4)
boxplot(localDeltaOriSomaValuesV1(:,1),'color',fig.cocV1(7,:)); hold on
scatter(repmat(1:1,size(localDeltaOriSomaValuesV1(:,1),1),1), localDeltaOriSomaValuesV1(:,1), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'jitter','on','jitterAmount',0.15)
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels({'V1 inputs'})
ylim([0 90])
subplot(1,5,5)
boxplot(localDeltaOriSomaV1Random,'color',fig.cocAll(9,:)); hold on
xticklabels({'random V1'})
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
saveas(gcf, fullfile(saveDirImages,'2250_Local deltaOriSoma_compareRandom.png'))

%Fig 222: local DeltaOriSoma vs deltaOri
setFigProperties(2251, fig);
scatter(goodDeltaOri,localDeltaOriSomaValues(:,1),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaOri,localDeltaOriSomaValuesA19(:,1),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaOri,localDeltaOriSomaValuesV1(:,1),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
xlim([0 90])
ylim([0 45])
ylabel(sprintf('Local delta Ori soma'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'222_Local DeltaOriSoma vs delta Ori.png'))

%Fig 226-227: nearest stimilar prefOri (vs. deltaOri)
setFigProperties(226, fig)
subplot(1,3,1)
boxplot([allSpines(oriGood).nearestSimilarPrefOri],'Label', {'All spines'}, 'color',fig.cocAll(7,:))
ylim([0 30])
hold on
scatter(repmat(1:1,size([allSpines(oriGood).nearestSimilarPrefOri],2),1), [allSpines(oriGood).nearestSimilarPrefOri], 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
ylabel('Nearest neighbor with similar pref')
subplot(1,3,2)
boxplot([allSpines(oriGoodA19).nearestSimilarPrefOri],'Label', {'A19 inputs'}, 'color',fig.cocA19(7,:))
ylim([0 30])
hold on
scatter(repmat(1:1,size([allSpines(oriGoodA19).nearestSimilarPrefOri],2),1), [allSpines(oriGoodA19).nearestSimilarPrefOri], 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
subplot(1,3,3)
boxplot([allSpines(oriGoodV1).nearestSimilarPrefOri],'Label', {'V1 inputs'}, 'color',fig.cocV1(7,:))
ylim([0 30])
hold on
scatter(repmat(1:1,size([allSpines(oriGoodV1).nearestSimilarPrefOri],2),1), [allSpines(oriGoodV1).nearestSimilarPrefOri], 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'226_Nearest neighbor same pref.png'))

setFigProperties(227, fig)
scatter(goodDeltaOri,[allSpines(oriGood).nearestSimilarPrefOri],'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaOri,[allSpines(oriGoodA19).nearestSimilarPrefOri],'o', 'filled','MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaOri,[allSpines(oriGoodV1).nearestSimilarPrefOri],'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
xlim([0 90])
ylim([0 25])
ylabel(sprintf('Nearest neighbor with similar preference in \\mum'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'227_Nearest neighbor similar pref vs. deltaOri.png'))

%Fig 228-229: local OSI
localOSIValues = cell2mat(cellfun(@(x) x.localOSI, num2cell(allSpines(goodROIs)), 'UniformOutput', false));
localOSIValues = reshape(localOSIValues, 4, [])';
localOSIValuesA19 = cell2mat(cellfun(@(x) x.localOSI, num2cell(allSpines(goodA19)), 'UniformOutput', false));
localOSIValuesA19 = reshape(localOSIValuesA19, 4, [])';
localOSIValuesV1 = cell2mat(cellfun(@(x) x.localOSI, num2cell(allSpines(goodV1)), 'UniformOutput', false));
localOSIValuesV1 = reshape(localOSIValuesV1, 4, [])';

setFigProperties(228, fig)
subplot(1,3,1)
boxplot(localOSIValues,'color',fig.cocAll(7,:))
ylim([0 1])
hold on
scatter(repmat(1:4,size(localOSIValues,1),1), localOSIValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('OSI of local environment')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localOSIValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 1])
hold on
scatter(repmat(1:4,size(localOSIValuesA19,1),1), localOSIValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localOSIValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 1])
hold on
scatter(repmat(1:4,size(localOSIValuesV1,1),1), localOSIValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'228_Local OSI.png'))

setFigProperties(229, fig)
scatter(goodOSI,localOSIValues(:,2),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19OSI,localOSIValuesA19(:,2),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1OSI,localOSIValuesV1(:,2),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
ylim([0 1])
xlim([0 1])
ylabel('local OSI')
xlabel('OSI of spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'229_Spine OSI vs local OSI.png'))


%% 230s: Local environment dir (230s)
% Fig 230-232: local Dir Dispersion
localDirDispValues = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(allSpines(dirGood)), 'UniformOutput', false));
localDirDispValues = reshape(localDirDispValues, 4, [])';
localDirDispValuesA19 = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(allSpines(dirGoodA19)), 'UniformOutput', false));
localDirDispValuesA19 = reshape(localDirDispValuesA19, 4, [])';
localDirDispValuesV1 = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(allSpines(dirGoodV1)), 'UniformOutput', false));
localDirDispValuesV1 = reshape(localDirDispValuesV1, 4, [])';

setFigProperties(230, fig)
subplot(1,3,1)
boxplot(localDirDispValues(:,2)','Label', {'All spines'})
ylim([0 90])
hold on
scatter(repmat(1:1,size(localDirDispValues,1),1), localDirDispValues(:,2), 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
box off
ylabel('Local dispersion dir')
subplot(1,3,2)
boxplot(localDirDispValuesA19(:,2)','Label', {'A19 Inputs'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDirDispValuesA19,1),1), localDirDispValuesA19(:,2), 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
subplot(1,3,3)
boxplot(localDirDispValuesV1(:,2)','Label', {'V1 Inputs'})
ylim([0 90])
hold on
scatter(repmat(1:1,size(localDirDispValuesV1,1),1), localDirDispValuesV1(:,2), 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'230_Local dir dispersion.png'))

setFigProperties(231, fig);
subplot(1,3,1)
boxplot(localDirDispValues,'color',fig.cocAll(7,:))
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDirDispValues,1),1), localDirDispValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('Local dispersion dir')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDirDispValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDirDispValuesA19,1),1), localDirDispValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDirDispValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 90])
hold on
scatter(repmat(1:4,size(localDirDispValuesV1,1),1), localDirDispValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'231_All local dispersion.png'))

%Fig 232: local dir Dispersion vs deltaDir
setFigProperties(232, fig);
scatter(goodDeltaDir,localDirDispValues(:,2),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaDir,localDirDispValuesA19(:,2),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaDir,localDirDispValuesV1(:,2),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
xlim([0 180])
ylim([0 90])
ylabel(sprintf('Local dir dispersion within (5 \\mum) of spine'))
xlabel('delta Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'232_Local dir dispersion vs delta Dir.png'))

%Fig 233: HI dir
HIDirValues = cell2mat(cellfun(@(x) x.HIDir, num2cell(allSpines(dirGood)), 'UniformOutput', false));
HIDirValues = reshape(HIDirValues, 4, [])';
HIDirValuesA19 = cell2mat(cellfun(@(x) x.HIDir, num2cell(allSpines(dirGoodA19)), 'UniformOutput', false));
HIDirValuesA19 = reshape(HIDirValuesA19, 4, [])';
HIDirValuesV1 = cell2mat(cellfun(@(x) x.HIDir, num2cell(allSpines(dirGoodV1)), 'UniformOutput', false));
HIDirValuesV1 = reshape(HIDirValuesV1, 4, [])';

setFigProperties(233, fig)
scatter(goodDeltaDir,HIDirValues(:,2),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaDir,HIDirValuesA19(:,2),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaDir,HIDirValuesV1(:,2),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
ylim([0 1])
xlim([0 180])
ylabel(sprintf('HI (dir) of local environment (5 \\mum) of spine'))
xlabel('delta Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'233_HI dir vs delta Dir.png'))

%Fig 234: local deltaDir
localDeltaDirValues = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(allSpines(dirGood)), 'UniformOutput', false));
localDeltaDirValues = reshape(localDeltaDirValues, 4, [])';
localDeltaDirValuesA19 = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(allSpines(dirGoodA19)), 'UniformOutput', false));
localDeltaDirValuesA19 = reshape(localDeltaDirValuesA19, 4, [])';
localDeltaDirValuesV1 = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(allSpines(dirGoodV1)), 'UniformOutput', false));
localDeltaDirValuesV1 = reshape(localDeltaDirValuesV1, 4, [])';

setFigProperties(234, fig)
subplot(1,3,1)
boxplot(localDeltaDirValues,'color',fig.cocAll(7,:))
ylim([0 180])
hold on
scatter(repmat(1:4,size(localDeltaDirValues,1),1), localDeltaDirValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('deltaDir of local environment')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDeltaDirValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 180])
hold on
scatter(repmat(1:4,size(localDeltaDirValuesA19,1),1), localDeltaDirValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDeltaDirValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 180])
hold on
scatter(repmat(1:4,size(localDeltaDirValuesV1,1),1), localDeltaDirValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'234_Local deltaDir.png'))

%Fig 235: local deltaDirSOma
localDeltaDirSomaValues = cell2mat(cellfun(@(x) x.localDeltaDirSoma, num2cell(allSpines(dirGood)), 'UniformOutput', false));
localDeltaDirSomaValues = reshape(localDeltaDirSomaValues, 4, [])';
localDeltaDirSomaValuesA19 = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(dirGoodA19)), 'UniformOutput', false));
localDeltaDirSomaValuesA19 = reshape(localDeltaDirSomaValuesA19, 4, [])';
localDeltaDirSomaValuesV1 = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(dirGoodV1)), 'UniformOutput', false));
localDeltaDirSomaValuesV1 = reshape(localDeltaDirSomaValuesV1, 4, [])';

setFigProperties(235, fig)
subplot(1,3,1)
boxplot(localDeltaDirSomaValues,'color',fig.cocAll(7,:))
ylim([0 180])
hold on
scatter(repmat(1:4,size(localDeltaDirSomaValues,1),1), localDeltaDirSomaValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('deltaDirSoma of local environment')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDeltaDirSomaValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 180])
hold on
scatter(repmat(1:4,size(localDeltaDirSomaValuesA19,1),1), localDeltaDirSomaValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDeltaDirSomaValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 180])
hold on
scatter(repmat(1:4,size(localDeltaDirSomaValuesV1,1),1), localDeltaDirSomaValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'235_Local deltaDirSoma.png'))

%Fig 236-237: nearest SimilarPrefDir
setFigProperties(236, fig)
subplot(1,3,1)
boxplot([allSpines(dirGood).nearestSimilarPrefDir],'Label', {'All spines'}, 'color',fig.cocAll(7,:))
ylim([0 30])
hold on
scatter(repmat(1:1,size([allSpines(dirGood).nearestSimilarPrefDir],2),1), [allSpines(dirGood).nearestSimilarPrefDir], 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
ylabel('Nearest neighbor with similar pref dir')
subplot(1,3,2)
boxplot([allSpines(dirGoodA19).nearestSimilarPrefDir],'Label', {'A19 inputs'}, 'color',fig.cocA19(7,:))
ylim([0 30])
hold on
scatter(repmat(1:1,size([allSpines(dirGoodA19).nearestSimilarPrefDir],2),1), [allSpines(dirGoodA19).nearestSimilarPrefDir], 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
subplot(1,3,3)
boxplot([allSpines(dirGoodV1).nearestSimilarPrefDir],'Label', {'V1 inputs'}, 'color',fig.cocV1(7,:))
ylim([0 30])
hold on
scatter(repmat(1:1,size([allSpines(dirGoodV1).nearestSimilarPrefDir],2),1), [allSpines(dirGoodV1).nearestSimilarPrefDir], 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'236_Nearest neighbor same pref dir.png'))

setFigProperties(237, fig)
scatter(goodDeltaDir,[allSpines(dirGood).nearestSimilarPrefDir],'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DeltaDir,[allSpines(dirGoodA19).nearestSimilarPrefDir],'o', 'filled','MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DeltaDir,[allSpines(dirGoodV1).nearestSimilarPrefDir],'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
xlim([0 90])
ylim([0 25])
ylabel(sprintf('Nearest neighbor with similar preference in \\mum'))
xlabel('delta Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'227_Nearest neighbor similar pref dir vs. deltaDir.png'))

%Fig 238-239: local DSI
localDSIValues = cell2mat(cellfun(@(x) x.localDSI, num2cell(allSpines(goodROIs)), 'UniformOutput', false));
localDSIValues = reshape(localDSIValues, 4, [])';
localDSIValuesA19 = cell2mat(cellfun(@(x) x.localDSI, num2cell(allSpines(goodA19)), 'UniformOutput', false));
localDSIValuesA19 = reshape(localDSIValuesA19, 4, [])';
localDSIValuesV1 = cell2mat(cellfun(@(x) x.localDSI, num2cell(allSpines(goodV1)), 'UniformOutput', false));
localDSIValuesV1 = reshape(localDSIValuesV1, 4, [])';

setFigProperties(228, fig)
subplot(1,3,1)
boxplot(localDSIValues,'color',fig.cocAll(7,:))
ylim([0 1])
hold on
scatter(repmat(1:4,size(localDSIValues,1),1), localDSIValues, 'filled','MarkerFaceColor',fig.cocAll(9,:), 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
box off
legend('All spines')
legend boxoff
ylabel('DSI of local environment')
xticklabels({'2.5', '5', '7.5', '10'})
subplot(1,3,2)
boxplot(localDSIValuesA19,'color',fig.cocA19(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 1])
hold on
scatter(repmat(1:4,size(localDSIValuesA19,1),1), localDSIValuesA19, 'filled','MarkerFaceColor',fig.cocA19(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('A19 input')
legend boxoff
subplot(1,3,3)
boxplot(localDSIValuesV1,'color',fig.cocV1(7,:))
xticklabels({'2.5', '5', '7.5', '10'})
ylim([0 1])
hold on
scatter(repmat(1:4,size(localDSIValuesV1,1),1), localDSIValuesV1, 'filled','MarkerFaceColor',fig.cocV1(9,:), 'MarkerFaceAlpha',0.9','jitter','on','jitterAmount',0.15)
box off
legend('V1 input')
legend boxoff
commonLabel = xlabel('Radius size in um');
commonLabel.Position = [-fig.width/2,-fig.height/2, 0];
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'238_Local DSI.png'))

setFigProperties(239, fig)
scatter(goodDSI,localDSIValues(:,2),'o', 'MarkerEdgeColor', fig.cocAll(7,:));
hold on
scatter(A19DSI,localDSIValuesA19(:,2),'o','filled', 'MarkerFaceColor', fig.cocA19(7,:));
scatter(V1DSI,localDSIValuesV1(:,2),'o', 'filled', 'MarkerFaceColor', fig.cocV1(7,:));
ylim([0 1])
xlim([0 1])
ylabel('local DSI')
xlabel('DSI of spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'239_Spine DSI vs local DSI.png'))

%% 240s: Pairwise distances

goodPairs = find([pwMeasures.goodPair] == 1);
oriPairs = find([pwMeasures.oriSelect] == 1);
dirPairs = find([pwMeasures.dirSelect] == 1);

nonInputsPairs = find([pwMeasures.inputPair] == 0);
mixedA19Pairs = intersect(find([pwMeasures.inputPair] == 1), find([pwMeasures.inputType] == 1));
mixedV1Pairs = intersect(find([pwMeasures.inputPair] == 1), find([pwMeasures.inputType] == 2));
inputsA19Pairs = intersect(find([pwMeasures.inputPair] == 2),find([pwMeasures.inputType] == 1));
inputsV1Pairs = intersect(find([pwMeasures.inputPair] == 2),find([pwMeasures.inputType] == 2));

goodNonInputsPairs = intersect(goodPairs, nonInputsPairs);
goodMixedA19Pairs = intersect(goodPairs, mixedA19Pairs);
goodMixedV1Pairs = intersect(goodPairs, mixedV1Pairs);
goodInputA19Pairs = intersect(goodPairs, inputsA19Pairs);
goodInputV1Pairs = intersect(goodPairs, inputsV1Pairs);

oriNonInputsPairs = intersect(oriPairs, nonInputsPairs);
oriMixedA19Pairs = intersect(oriPairs, mixedA19Pairs);
oriMixedV1Pairs = intersect(oriPairs, mixedV1Pairs);
oriInputA19Pairs = intersect(oriPairs, inputsA19Pairs);
oriInputV1Pairs = intersect(oriPairs, inputsV1Pairs);

dirNonInputsPairs = intersect(dirPairs, nonInputsPairs);
dirMixedA19Pairs = intersect(dirPairs, mixedA19Pairs);
dirMixedV1Pairs = intersect(dirPairs, mixedV1Pairs);
dirInputA19Pairs = intersect(dirPairs, inputsA19Pairs);
dirInputV1Pairs = intersect(dirPairs, inputsV1Pairs);

% Fig 240: delta Pref Ori
setFigProperties(240, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(oriNonInputsPairs).distance],[pwMeasures(oriNonInputsPairs).deltaOri],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(oriMixedA19Pairs).distance],[pwMeasures(oriMixedA19Pairs).deltaOri],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(oriInputA19Pairs).distance],[pwMeasures(oriInputA19Pairs).deltaOri],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(oriMixedV1Pairs).distance],[pwMeasures(oriMixedV1Pairs).deltaOri],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(oriInputV1Pairs).distance],[pwMeasures(oriInputV1Pairs).deltaOri],fig.cocV1(9,:));
ylabel('deltaOri')
ylim([0 90])
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'240_DeltaOri vs. Distance.png'))

%Fig 241: delta Pref Dir
setFigProperties(241, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(dirNonInputsPairs).distance],[pwMeasures(dirNonInputsPairs).deltaDir],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(dirMixedA19Pairs).distance],[pwMeasures(dirMixedA19Pairs).deltaDir],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(dirInputA19Pairs).distance],[pwMeasures(dirInputA19Pairs).deltaDir],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(dirMixedV1Pairs).distance],[pwMeasures(dirMixedV1Pairs).deltaDir],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(dirInputV1Pairs).distance],[pwMeasures(dirInputV1Pairs).deltaDir],fig.cocV1(9,:));
ylabel('deltaDir')
ylim([0 180])
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'241_DeltaDir vs. Distance.png'))

%Fig 242: deltaOSI
setFigProperties(242, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).deltaOSI],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedA19Pairs).distance],[pwMeasures(goodMixedA19Pairs).deltaOSI],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputA19Pairs).distance],[pwMeasures(goodInputA19Pairs).deltaOSI],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedV1Pairs).distance],[pwMeasures(goodMixedV1Pairs).deltaOSI],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputV1Pairs).distance],[pwMeasures(goodInputV1Pairs).deltaOSI],fig.cocV1(9,:));
ylabel('deltaOSI')
ylim([0 .5])
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'242_DeltaOSI vs. Distance.png'))

%Fig 243: deltaDSI
setFigProperties(243, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).deltaDSI],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedA19Pairs).distance],[pwMeasures(goodMixedA19Pairs).deltaDSI],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputA19Pairs).distance],[pwMeasures(goodInputA19Pairs).deltaDSI],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedV1Pairs).distance],[pwMeasures(goodMixedV1Pairs).deltaDSI],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputV1Pairs).distance],[pwMeasures(goodInputV1Pairs).deltaDSI],fig.cocV1(9,:));
ylabel('deltaDSI')
ylim([0 .7])
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'243_DeltaDSI vs. Distance.png'))

%Fig 244: deltaDSIvect
setFigProperties(244, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).deltaDSIvect],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedA19Pairs).distance],[pwMeasures(goodMixedA19Pairs).deltaDSIvect],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputA19Pairs).distance],[pwMeasures(goodInputA19Pairs).deltaDSIvect],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedV1Pairs).distance],[pwMeasures(goodMixedV1Pairs).deltaDSIvect],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputV1Pairs).distance],[pwMeasures(goodInputV1Pairs).deltaDSIvect],fig.cocV1(9,:));
ylabel('deltaDSIvect')
ylim([0 .5])
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'244_DeltaDSIvect vs. Distance.png'))

%Fig 245: Tuning curve correlation
setFigProperties(245, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).curveCorr],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedA19Pairs).distance],[pwMeasures(goodMixedA19Pairs).curveCorr],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputA19Pairs).distance],[pwMeasures(goodInputA19Pairs).curveCorr],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedV1Pairs).distance],[pwMeasures(goodMixedV1Pairs).curveCorr],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputV1Pairs).distance],[pwMeasures(goodInputV1Pairs).curveCorr],fig.cocV1(9,:));
ylabel('Tuning curve correlation')
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'245_Tuning curve correlation vs. Distance.png'))

%Fig 246: Trial-to-trial correlation
setFigProperties(246, fig)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodNonInputsPairs).distance],[pwMeasures(goodNonInputsPairs).trialCorr],fig.cocAll(7,:));
hold on
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedA19Pairs).distance],[pwMeasures(goodMixedA19Pairs).trialCorr],fig.cocA19(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputA19Pairs).distance],[pwMeasures(goodInputA19Pairs).trialCorr],fig.cocA19(9,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodMixedV1Pairs).distance],[pwMeasures(goodMixedV1Pairs).trialCorr],fig.cocV1(7,:));
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodInputV1Pairs).distance],[pwMeasures(goodInputV1Pairs).trialCorr],fig.cocV1(9,:));
ylabel('Trial-to-trial correlation')
legend({'non-inputs'; 'non-input & input A19'; 'inputs A19'; 'non-input & input V1'; 'inputs V1'}, 'box','off')
saveas(gcf, fullfile(saveDir,'246_Trial-to-trial correlation vs. Distance.png'))

%% 270s: Branch-wide analysis

%Fig 270: circular Dispersion
setFigProperties(270, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.circDispersion],2)), [allBranches.circDispersion], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.circDispersion],'color',fig.cocAll(7,:))
box off
ylim([0 45])
ylabel('Branch Dispersion')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).circDispersion],2)), [allBranches(A19Branches).circDispersion], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).circDispersion],'color',fig.cocA19(7,:))
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).circDispersion],2)), [allBranches(V1Branches).circDispersion], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).circDispersion],'color',fig.cocV1(7,:))
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'270_Branch dispersion.png'))

%Fig 271: circular Dispersion
setFigProperties(271, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.dirCircDispersion],2)), [allBranches.dirCircDispersion], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.dirCircDispersion],'color',fig.cocAll(7,:))
box off
ylim([0 90])
ylabel('Branch dir Dispersion')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).dirCircDispersion],2)), [allBranches(A19Branches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).dirCircDispersion],'color',fig.cocA19(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).dirCircDispersion],2)), [allBranches(V1Branches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).dirCircDispersion],'color',fig.cocV1(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'271_Branch dir dispersion.png'))

%Fig 272: deltaOri
setFigProperties(272, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.deltaOri],2)), [allBranches.deltaOri], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.deltaOri],'color',fig.cocAll(7,:))
box off
ylim([0 90])
ylabel('deltaOri')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).deltaOri],2)), [allBranches(A19Branches).deltaOri], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).deltaOri],'color',fig.cocA19(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).deltaOri],2)), [allBranches(V1Branches).deltaOri], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).deltaOri],'color',fig.cocV1(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'272_deltaOri.png'))

%Fig 273: deltaDir
setFigProperties(273, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.deltaDir],2)), [allBranches.deltaDir], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.deltaDir],'color',fig.cocAll(7,:))
box off
ylim([0 180])
ylabel('deltaDir')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).deltaDir],2)), [allBranches(A19Branches).deltaDir], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).deltaDir],'color',fig.cocA19(7,:))
ylim([0 180])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).deltaDir],2)), [allBranches(V1Branches).deltaDir], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).deltaDir],'color',fig.cocV1(7,:))
ylim([0 180])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'273_deltaDir.png'))

%Fig 274: medianOSI
setFigProperties(274, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.medianOSI],2)), [allBranches.medianOSI], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.medianOSI],'color',fig.cocAll(7,:))
box off
ylim([0 0.5])
ylabel('medianOSI')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).medianOSI],2)), [allBranches(A19Branches).medianOSI], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).medianOSI],'color',fig.cocA19(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).medianOSI],2)), [allBranches(V1Branches).medianOSI], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).medianOSI],'color',fig.cocV1(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'274_medianOSI.png'))

%Fig 275: medianOSI
setFigProperties(275, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.medianDSI],2)), [allBranches.medianDSI], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.medianDSI],'color',fig.cocAll(7,:))
box off
ylim([0 1])
ylabel('medianDSI')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).medianDSI],2)), [allBranches(A19Branches).medianDSI], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).medianDSI],'color',fig.cocA19(7,:))
ylim([0 1])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).medianDSI],2)), [allBranches(V1Branches).medianDSI], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).medianDSI],'color',fig.cocV1(7,:))
ylim([0 1])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'275_medianDSI.png'))

%Fig 276: medianDSIvect
setFigProperties(276, fig)
subplot(1,3,1)
scatter(repmat(1:1,size([allBranches.medianDSIvect],2)), [allBranches.medianDSIvect], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches.medianDSIvect],'color',fig.cocAll(7,:))
box off
ylim([0 0.5])
ylabel('medianDSIvect')
xticklabels('All branches')
subplot(1,3,2)
hold on
scatter(repmat(1:1,size([allBranches(A19Branches).medianDSIvect],2)), [allBranches(A19Branches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(A19Branches).medianDSIvect],'color',fig.cocA19(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,3,3)
hold on
scatter(repmat(1:1,size([allBranches(V1Branches).medianDSIvect],2)), [allBranches(V1Branches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(V1Branches).medianDSIvect],'color',fig.cocV1(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'276_medianDSIvect.png'))

%% 280s: Branch-wise analysis for apical

%Fig 280: circular Dispersion
setFigProperties(280, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).circDispersion],2)), [allBranches(apicalBranches).circDispersion], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).circDispersion],'color',fig.cocAll(7,:))
box off
ylim([0 45])
ylabel('Branch Dispersion')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).circDispersion],2)), [allBranches(apicalA19Branches).circDispersion], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).circDispersion],'color',fig.cocA19(7,:))
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).circDispersion],2)), [allBranches(apicalV1Branches).circDispersion], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).circDispersion],'color',fig.cocV1(7,:))
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).circDispersion],2)), [allBranches(basalBranches).circDispersion], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).circDispersion],'color',fig.cocAll(7,:))
box off
ylim([0 45])
ylabel('Branch Dispersion')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).circDispersion],2)), [allBranches(basalA19Branches).circDispersion], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).circDispersion],'color',fig.cocA19(7,:))
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).circDispersion],2)), [allBranches(basalV1Branches).circDispersion], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).circDispersion],'color',fig.cocV1(7,:))
ylim([0 45])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'280_apical basal Branch dispersion.png'))

%Fig 281: circular Dispersion
setFigProperties(281, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).dirCircDispersion],2)), [allBranches(apicalBranches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).dirCircDispersion],'color',fig.cocAll(7,:))
box off
ylim([0 90])
ylabel('Branch dir Dispersion')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).dirCircDispersion],2)), [allBranches(apicalA19Branches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).dirCircDispersion],'color',fig.cocA19(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).dirCircDispersion],2)), [allBranches(apicalV1Branches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).dirCircDispersion],'color',fig.cocV1(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).dirCircDispersion],2)), [allBranches(basalBranches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).dirCircDispersion],'color',fig.cocAll(7,:))
box off
ylim([0 90])
ylabel('Branch dir Dispersion')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).dirCircDispersion],2)), [allBranches(basalA19Branches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).dirCircDispersion],'color',fig.cocA19(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).dirCircDispersion],2)), [allBranches(basalV1Branches).dirCircDispersion], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).dirCircDispersion],'color',fig.cocV1(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'281_apical basal Branch dir dispersion.png'))

%Fig 282: deltaOri
setFigProperties(282, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).deltaOri],2)), [allBranches(apicalBranches).deltaOri], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).deltaOri],'color',fig.cocAll(7,:))
box off
ylim([0 90])
ylabel('deltaOri')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).deltaOri],2)), [allBranches(apicalA19Branches).deltaOri], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).deltaOri],'color',fig.cocA19(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).deltaOri],2)), [allBranches(apicalV1Branches).deltaOri], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).deltaOri],'color',fig.cocV1(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).deltaOri],2)), [allBranches(basalBranches).deltaOri], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).deltaOri],'color',fig.cocAll(7,:))
box off
ylim([0 90])
ylabel('deltaOri')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).deltaOri],2)), [allBranches(basalA19Branches).deltaOri], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).deltaOri],'color',fig.cocA19(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).deltaOri],2)), [allBranches(basalV1Branches).deltaOri], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).deltaOri],'color',fig.cocV1(7,:))
ylim([0 90])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'282_apicalbasalDeltaOri.png'))

%Fig 283: deltaDir
setFigProperties(283, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).deltaDir],2)), [allBranches(apicalBranches).deltaDir], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).deltaDir],'color',fig.cocAll(7,:))
box off
ylim([0 180])
ylabel('deltaDir')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).deltaDir],2)), [allBranches(apicalA19Branches).deltaDir], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).deltaDir],'color',fig.cocA19(7,:))
ylim([0 180])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).deltaDir],2)), [allBranches(apicalV1Branches).deltaDir], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).deltaDir],'color',fig.cocV1(7,:))
ylim([0 180])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).deltaDir],2)), [allBranches(basalBranches).deltaDir], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).deltaDir],'color',fig.cocAll(7,:))
box off
ylim([0 180])
ylabel('deltaDir')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).deltaDir],2)), [allBranches(basalA19Branches).deltaDir], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).deltaDir],'color',fig.cocA19(7,:))
ylim([0 180])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).deltaDir],2)), [allBranches(basalV1Branches).deltaDir], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).deltaDir],'color',fig.cocV1(7,:))
ylim([0 180])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'283_apicalbasaldeltaDir.png'))

%Fig 284: medianOSI
setFigProperties(284, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).medianOSI],2)), [allBranches(apicalBranches).medianOSI], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).medianOSI],'color',fig.cocAll(7,:))
box off
ylim([0 0.5])
ylabel('medianOSI')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).medianOSI],2)), [allBranches(apicalA19Branches).medianOSI], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).medianOSI],'color',fig.cocA19(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).medianOSI],2)), [allBranches(apicalV1Branches).medianOSI], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).medianOSI],'color',fig.cocV1(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).medianOSI],2)), [allBranches(basalBranches).medianOSI], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).medianOSI],'color',fig.cocAll(7,:))
box off
ylim([0 0.5])
ylabel('medianOSI')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).medianOSI],2)), [allBranches(basalA19Branches).medianOSI], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).medianOSI],'color',fig.cocA19(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).medianOSI],2)), [allBranches(basalV1Branches).medianOSI], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).medianOSI],'color',fig.cocV1(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'284_apicalbasalmedianOSI.png'))

%Fig 285: medianOSI
setFigProperties(285, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).medianDSI],2)), [allBranches(apicalBranches).medianDSI], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).medianDSI],'color',fig.cocAll(7,:))
box off
ylim([0 1])
ylabel('medianDSI')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).medianDSI],2)), [allBranches(apicalA19Branches).medianDSI], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).medianDSI],'color',fig.cocA19(7,:))
ylim([0 1])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).medianDSI],2)), [allBranches(apicalV1Branches).medianDSI], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).medianDSI],'color',fig.cocV1(7,:))
ylim([0 1])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).medianDSI],2)), [allBranches(basalBranches).medianDSI], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).medianDSI],'color',fig.cocAll(7,:))
box off
ylim([0 1])
ylabel('medianDSI')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).medianDSI],2)), [allBranches(basalA19Branches).medianDSI], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).medianDSI],'color',fig.cocA19(7,:))
ylim([0 1])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).medianDSI],2)), [allBranches(basalV1Branches).medianDSI], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).medianDSI],'color',fig.cocV1(7,:))
ylim([0 1])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'285_apicalbasalmedianDSI.png'))

%Fig 286: medianDSIvect
setFigProperties(286, fig)
subplot(1,6,1)
scatter(repmat(1:1,size([allBranches(apicalBranches).medianDSIvect],2)), [allBranches(apicalBranches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(apicalBranches).medianDSIvect],'color',fig.cocAll(7,:))
box off
ylim([0 0.5])
ylabel('medianDSIvect')
xticklabels('apical branches')
subplot(1,6,2)
hold on
scatter(repmat(1:1,size([allBranches(apicalA19Branches).medianDSIvect],2)), [allBranches(apicalA19Branches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalA19Branches).medianDSIvect],'color',fig.cocA19(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,3)
hold on
scatter(repmat(1:1,size([allBranches(apicalV1Branches).medianDSIvect],2)), [allBranches(apicalV1Branches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(apicalV1Branches).medianDSIvect],'color',fig.cocV1(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
subplot(1,6,4)
scatter(repmat(1:1,size([allBranches(basalBranches).medianDSIvect],2)), [allBranches(basalBranches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocAll(3,:), 'MarkerFaceAlpha',0.1')
hold on
boxplot([allBranches(basalBranches).medianDSIvect],'color',fig.cocAll(7,:))
box off
ylim([0 0.5])
ylabel('medianDSIvect')
xticklabels('basal branches')
subplot(1,6,5)
hold on
scatter(repmat(1:1,size([allBranches(basalA19Branches).medianDSIvect],2)), [allBranches(basalA19Branches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocA19(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalA19Branches).medianDSIvect],'color',fig.cocA19(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('A19 branches')
subplot(1,6,6)
hold on
scatter(repmat(1:1,size([allBranches(basalV1Branches).medianDSIvect],2)), [allBranches(basalV1Branches).medianDSIvect], 'filled','MarkerFaceColor',fig.cocV1(5,:), 'MarkerFaceAlpha',0.1')
boxplot([allBranches(basalV1Branches).medianDSIvect],'color',fig.cocV1(7,:))
ylim([0 0.5])
set(gca, 'YColor', 'none', 'box', 'off')
xticklabels('V1 branches')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'286_apicalbasalmedianDSIvect.png'))

%show them all on the screen
%tilefigs
%% 5.) Summary
disp(['Number of inputs A19: ' num2str(length(A19ROIs), '%d')])
disp(['Number of 2p matched inputs A19: ' num2str(length(A19Funct), '%d')])
disp(['Number of responsive A19 inputs: ' num2str(length(goodA19), '%d')])
disp(['Number of ori-selective A19 spines: ' num2str(length(oriGoodA19), '%d')])

disp(['Number of inputs V1: ' num2str(length(V1ROIs), '%d')])
disp(['Number of 2p matched inputs V1: ' num2str(length(V1Funct), '%d')])
disp(['Number of responsive V1 inputs: ' num2str(length(goodV1), '%d')])
disp(['Number of ori-selective V1 spines: ' num2str(length(oriGoodV1), '%d')])

%% 6.) Save
%write to text file
fid = fopen([saveDir filesep 'V1InputsSummary.txt'], 'w');
formatSpec = ['Number of V1 inputs: %d \n'...
    'Number of 2p matched V1 inputs: %d \n' ...
    'Number of responsive V1 inputs: %d \n'...
    'Number of ori-selective V1 spines:: %d'];
fprintf(fid, formatSpec, length(V1ROIs), length(V1Funct),length(goodV1), length(oriGoodV1));
fclose(fid);

fid = fopen([saveDir filesep 'A19InputsSummary.txt'], 'w');
formatSpec = ['Number of A19 inputs: %d \n'...
    'Number of 2p matched A19 inputs: %d \n' ...
    'Number of responsive A19 inputs: %d \n'...
    'Number of ori-selective A19 spines:: %d'];
fprintf(fid, formatSpec, length(A19ROIs), length(A19Funct),length(goodA19), length(oriGoodA19));
fclose(fid);


