function AnatomicalInputAnalysis(animal, varargin)

%Loads data of the confocal reconstruction of the cell, the excel sheet 
%containing the numbers of the confocal ROIs that are STED inputs and the 
%information of the STED coverage for the cell to combine all of it, plot
%the presentation of the (covered) cells and inputs, some analysis and the
%plotting of the analysis

%Inputs:
% - animal: which animal are we analysing?
% - cell: if there are multiple cells per animal, which one are we looking
% at?

%STEPS:
% - 0.) Define and if necessary, add the save dir
% - 1.) Load the confocal representation, the input data and the coverage
% - 2.) Add input data to ROIs 
% - 3.) Define certain ROI and dendrite categories as logical indicators
% - 4.) Make tracing of all the areas that were done with STED
% - 5.) Analysis:
%   - a) Coverage and inputs: How many inputs per dendrite? % of all spines
%   covered
%   - b.) Distance between inputs that are on the same dendrite
%   - c.) Cell wide: Distance from soma, branch order
% - 6.) Plot
% - 7.) Cell summary: STED coverage, input summary
% - 8.) Save

%Outputs:
% - Folder 03_AnatomicalInputAnalysis with Plots (see details at the
% beginning of Step 6)
% - 03_InputCellReconstruction.mat with variables:
%       - InputData: containing all data of a single ROI (previously: 
%       confData.ROIs), including whether it is an ROI
%       - cell: summary of the cell, including sted coverage, input
%       density, percentages, ...
%       - inputROINR: Numbers of all the ROIs that are inputs
% - AnatomyInputsSummary.txt: printed version of cell summary

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023

if nargin > 1
    cellName = varargin{1};
else
    cellName = [];
end

%% 0.) %define and if necessary, add the confocal and save dir 
baseDir = 'Z:\Juliane\InputAnalysis\';
cellDir = [baseDir char(animal) filesep char(cellName)];
confDir = [cellDir '\B - Confocal\'];
coverageDir = [cellDir '\D - Alignments\STED coverage'];

saveDirAna = [cellDir '\E - Analysis\'];
saveDir = ([saveDirAna filesep '03_AnatomicalInputAnalysis']);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%% 1.) Load the confocal representation and the matching

%a) Confocal image
%test if there is one or multiple slices
sliceFolders = dir([confDir 'SL*']);
isFolder=([sliceFolders.isdir]); %make sure those are folders, and not the .jpgs or .tif files
if ~isempty(sliceFolders(isFolder))
    for slice = 1:length(sliceFolders)
        imgFile = dir([confDir filesep sliceFolders(slice).name '\*.jpg']);
        temp = imread([confDir filesep sliceFolders(slice).name filesep imgFile.name]);
        template{slice} = imadjust(temp);
    end
    template = template{2}; %assumes that the cell body is always in the second slice if there are multiple slices
else
    imgFile = dir([confDir '\*.jpg']);
    template = imread([confDir filesep imgFile.name]);
    template = imadjust(template);
end

%b) Confocal reconstruction
temp = load([saveDirAna filesep '01_MorphologyCellReconstruction.mat']);
confData = temp.confData;
cell = temp.cell;

%c) Xls sheet with inpus
file = 'STED Inputs.xlsx';
if ~isempty(cellName)
    sheetName = [animal ' - ' cellName];
else
    sheetName = animal;
end
[~, ~, xls_all]=xlsread([baseDir file], sheetName);
inputMatrix = cell2mat(xls_all);
confData.ROIs = arrayfun(@(x) setfield(x, 'input', zeros(1,1)), confData.ROIs);

%d) STED coverage data
confData.ROIs = arrayfun(@(x) setfield(x, 'STEDdone', zeros(1,1)), confData.ROIs);
confData.allDendrites = arrayfun(@(x) setfield(x, 'STEDdone', []), confData.allDendrites);
coverFiles = dir([coverageDir filesep 'STEDCoverage*_Branch*.mat']);
coverageData = zeros(length(confData.allDendrites),3);
dc = 1;
for f=1:length(coverFiles) 
    %find the branch nr based on the file name - should end with
    %xx_dd-mmm-yyyy.mat with xx being the branch Nr
    brNr = str2double(coverFiles(f).name(end-17:end-16));
    
    %load the summary table
    new = load([coverFiles(f).folder filesep coverFiles(f).name]);
    toAdd = size(new.coverData,1);
    if toAdd > 1
        coverageData(dc:dc+toAdd-1,:) = new.coverData;
    else
        coverageData(dc,:) = new.coverData;
    end
    dc = dc+toAdd;
    
    %find all the dendrites/ROIs that are on this branch = fileNumber
    dendritesOnBranch = find([confData.allDendrites.Branch] == brNr);
    ROIsonConfBranch = find([confData.ROIs.Branch] == brNr); 
    %add the STED coverage info 
    for d = 1:length(new.DendDataT) %for every dendritic segment of this branch
          if ~isempty(new.DendDataT(d).ROIs) %for the ROIs
              STEDROIsOnBranch = [new.DendDataT(d).ROIs(find([new.DendDataT(d).ROIs.STEDdone])).ROINrOnBranch]; %which of the ROIs on that segment were done with STED
              for s = 1:length(STEDROIsOnBranch) %for every STED imaged ROI
                  ROIrow = intersect(ROIsonConfBranch, find([confData.ROIs.ROINrOnBranch]== STEDROIsOnBranch(s))); %find the one that is on the branch and has the same ROIonBranchNR
                  try
                    confData.ROIs(ROIrow).STEDdone = 1; %set it to 1
                  catch
                      disp(['Check ROI ' STEDROIsOnBranch(s) ' for an equivalent confocal ROI'])
                  end
              end
          end
          %add the STED coverage info to the tracing
          
          %MULTIPLE SLICE: Mulitple times the same dendrite, e. g. SL 02
          %and SL03 both have dendrite nr 1 -> correct!!!
          
          dendriteOnSlice = find([confData.allDendrites.dendNr] == new.DendDataT(d).dendNr);
          try
              confData.allDendrites(intersect(dendriteOnSlice, dendritesOnBranch)).STEDdone = new.DendDataT(d).STEDdone;
          catch
              disp('No STED done for this branch')
          end
    end
end
    
%% 2.) Add input data to ROIs 
headers = inputMatrix(1,:);
BranchNr = unique([confData.allDendrites.Branch]);
for d = 1:size(confData.Branches,2)
    %inputs
    denCol = find(headers == BranchNr(d));
    inputsOnBranch = inputMatrix(~isnan(inputMatrix(:,denCol)),denCol);
    if size(inputsOnBranch,1) > 1
        inputsOnBranch = inputsOnBranch(2:end); %remove header
        %find the right ROI
        ROIsOnBranch = find([confData.ROIs.Branch] == BranchNr(d));
        for input = 1:length(inputsOnBranch)
            ROINrMatches = find([confData.ROIs.ROINrOnBranch] == inputsOnBranch(input));
            ROIInputNr = intersect(ROIsOnBranch, ROINrMatches);
            %mark it down as 1
            try
                confData.ROIs(ROIInputNr).input = 1;
            catch
                disp(['No confocal ROI defined for branch ' num2str(BranchNr(d)) ', input ' num2str(inputsOnBranch(input))])
            end
        end
    end
end

%% 3.) Define certain ROI and dendrite categories as logical indicators
%apical vs. basal
apicalROIs = cellfun(@(x) strcmp(x, 'apical'),{confData.ROIs.type}); 
basalROIs = cellfun(@(x) strcmp(x, 'basal'),{confData.ROIs.type});

%inputs
inputROIs =  cellfun(@(x) (x==1),{confData.ROIs.input});
inputROINr = find(inputROIs ==1);

%combineThem
apicalInputROIs = and(apicalROIs,inputROIs);
basalInputROIs = and(basalROIs,inputROIs);

%spines imaged with STED
STEDDone = cellfun(@(x) (x == 1), {confData.ROIs.STEDdone});
apicalSTED = and(apicalROIs,STEDDone);
basalSTED = and(basalROIs,STEDDone);

%dendrites
%apicalDend = cellfun(@(x) strcmp(x, 'apical'),{confData.allDendrites.type}); 
%basalDend = cellfun(@(x) strcmp(x, 'basal'),{confData.allDendrites.type});

%% 4.) Make tracing of all the areas that were done with STED

confData.allDendrites = arrayfun(@(x) setfield(x, 'STEDTrace', []), confData.allDendrites);
for dend = 1:length(confData.allDendrites)
    %multiply tracing with whether that part is done with STED
    if ~isempty(confData.allDendrites(dend).STEDdone)
        rawTracing = confData.allDendrites(dend).normCoord .* confData.allDendrites(dend).STEDdone;
        %find which ones are zero 
        nonZeroTracing = ~(rawTracing(:,1) == 0); %assumes that if x is exactly zero than it can be removed

        %save only the ones that are nonZero to the STEDtrace
        confData.allDendrites(dend).STEDTrace = [rawTracing(nonZeroTracing,1) rawTracing(nonZeroTracing,2)];
    else
        confData.allDendrites(dend).STEDTrace = [];
    end
end

%% 5.Analysis
% a) Dendrite analysis: Coverage and inputs
for den = 1:size(confData.allDendrites,2)
    %coverage
    confData.allDendrites(den).STED = logical(coverageData(den,2));
    confData.allDendrites(den).STEDPerc = coverageData(den,3);
    confData.allDendrites(den).STEDLength = confData.allDendrites(den).STEDPerc*confData.allDendrites(den).length/100;
    
    %how many inputs per dendrite? dendrite density? Perc of all spines on
    %that dendrite (relative to STED coverage)
    ROIsOnDendrites = find([confData.ROIs.Dendrite] == den);
    inputsOnBranch = intersect(inputROINr, ROIsOnDendrites);
    if ~isempty(inputsOnBranch)
        confData.allDendrites(den).numInputs = length(inputsOnBranch);
        confData.allDendrites(den).InputDensity = confData.allDendrites(den).STEDLength/length(inputsOnBranch);
        confData.allDendrites(den).percOfAllSpines = 100 * length(inputsOnBranch)/(confData.allDendrites(den).numInputs *confData.allDendrites(den).STEDPerc);
    else
        confData.allDendrites(den).numInputs = 0;
        confData.allDendrites(den).InputDensity = NaN;
        confData.allDendrites(den).percOfAllSpines = 0;
    end
end

% b) Distance between inputs
allInputsDendrite = [confData.ROIs(inputROINr).Dendrite]; %which dendrites are the inputs on
dendritesWithInputs = unique(allInputsDendrite); %which dendrites have inputs
countsPerDendrite = histc(allInputsDendrite(:), dendritesWithInputs); %how many inputs per dendrites
dendritesWithMultiples = dendritesWithInputs(find(countsPerDendrite > 1)); %which dendrites have multiple inputs
numMultiples = sum(countsPerDendrite(countsPerDendrite>1)); %how many inputs are those?

%make the arrays
DistToInput = zeros(numMultiples,1);
DistToInputSim = zeros(numMultiples, 10000);
DistCoun = 1;
for dendNr = 1:length(dendritesWithMultiples)
    mulInp = find(allInputsDendrite == dendritesWithMultiples(dendNr)); %find all the inputs on the dendrite
    ROIsONdendrites = find([confData.ROIs.Dendrite] == dendritesWithMultiples(dendNr)); %find all ROIs on the dendrite
    %first measure real distances
    allDist = [confData.ROIs(inputROINr(mulInp)).distToSoma];
    for a = 1:length(allDist)
        distToSoma = allDist(a);
        diffDist= abs(allDist -distToSoma);
        diffDist(a) = 10000; %the distance to itself
        DistToInput(DistCoun) = min(diffDist);
        DistCoun = DistCoun +1;
    end
end

for rep = 1:10000
    DistCoun = 1;
    for dendNr = 1:length(dendritesWithMultiples)
        mulInp = find(allInputsDendrite == dendritesWithMultiples(dendNr)); %find all the inputs on the dendrite
        ROIsONdendrites = find([confData.ROIs.Dendrite] == dendritesWithMultiples(dendNr)); %find all ROIs on the dendrite
        %now let's repeat it by randomly selecting which spines on that
        %dendrite are the input spines
        newlySorted = ROIsONdendrites((randperm(length(ROIsONdendrites))));
        newInputs = newlySorted(1:length(mulInp));
        allDistRand = [confData.ROIs(newInputs).distToSoma];
        for a = 1:length(allDistRand)
            distToSomaRand = allDistRand(a);
            diffDistRand= abs(allDistRand -distToSomaRand);
            diffDistRand(a) = 10000;
            DistToInputSim(DistCoun, rep) = min(diffDistRand);
            DistCoun = DistCoun +1;
        end
    end
end

% c) Cell wide
DistFromSoma = cell2mat({confData.ROIs.distToSoma});
BranchOrder = cell2mat({confData.ROIs.BranchOrder});

InputDensity = cell2mat({confData.allDendrites.InputDensity});
STEDcoverDend = [confData.allDendrites.STEDLength];
suffSTEDcover = find(STEDcoverDend>10);

%% 6.) Plot
%Figure Overview
%10s - ROIs on top of tracing
%Fig 10: Input ROIs on top of tracing
%Fig 11: Apical input ROIs on top of tracing
%Fig 12: Basal input ROIs on top of tracing
%Fig 13: Both inputs separated by color on top of whole cell and STED
%covered cell tracing

%20s - Cell-wide analysis
%Fig 20: Distance from soma for all inputs and separated for
%apical/basal
%Fig 21: Same as Fig 20, but as cdfplot
%Fig 22: Histogram of branch order for all inputs and separated for 
%apical/basal
%Fig 23: Cdfplot of branch order for all inputs and separated for 
%apical/basal

%30s - Branch specific analysis
%Fig 30: Number of inputs per dendrite 
%Fig 31: Input density per dendrite for all/branches with more than one
%input
%Fig 32: Number of inputs vs. length of STED coverage

%40s - Clustering of inputs
%Fig 40: Distance to neares input on same branch
%Fig 41: Mean distance to nearest input

%--------------------------------------------------------------------------
%10s: Where are the ROIs
plotIDsOnCellTracing(confData.allDendrites,confData.ROIs, inputROIs, 10, 'blue',0)
saveas(gcf, fullfile(saveDir, '10_CellTracing_AllInputs.png'))
plotIDsOnCellTracing(confData.allDendrites,confData.ROIs, apicalInputROIs, 11, 'magenta',0)
saveas(gcf, fullfile(saveDir, '11_CellTracing_AllApicalInputs.png'))
plotIDsOnCellTracing(confData.allDendrites,confData.ROIs, basalInputROIs, 12, 'cyan',0)
saveas(gcf, fullfile(saveDir, '12_CellTracing_AllBasalInputs.png'))

%with STED tracing
plotIDsOnSTEDTracing(confData.allDendrites, confData.ROIs, inputROIs, 16)
saveas(gcf, fullfile(saveDir, '13_CellSTEDTracing_AllInputsCategorized.png'))

%--------------------------------------------------------------------------
%20s: Cell wide analysis
%Fig 20: Distance from soma for all inputs and separated for
%apical/basal
if sum(inputROIs > 1)
    figure(20)
    subplot(1,3,1)
    distributionPlot(DistFromSoma(inputROIs)', 'color', 'blue'); hold all
    boxplot(DistFromSoma(inputROIs)')
    title('Inputs')
    ylim([0 max(DistFromSoma)])
    subplot(1,3,2)
    try
        distributionPlot(DistFromSoma(apicalInputROIs)', 'color', 'red'); hold all
        boxplot(DistFromSoma(apicalInputROIs)')
    catch
        disp('Not enough apical spines')
    end
    title('Apical Spines')
    ylim([0 max(DistFromSoma)])
    axis off
    subplot(1,3,3)
    try
    distributionPlot(DistFromSoma(basalInputROIs)', 'color', 'green'); hold all
    boxplot(DistFromSoma(basalInputROIs)')
    catch
        disp('not enought basal spines')
    end
    title('basal Spines')
    ylim([0 max(DistFromSoma)])
    axis off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '20_DistanceFromBranchInputs.png'))

    %Fig 21: Distance from soma for all inputs and separated for
    %apical/basal as cdfplot
    figure(21);
    h(1,1) = cdfplot(DistFromSoma(inputROIs)); hold on
    try
        h(1,2) = cdfplot(DistFromSoma(apicalInputROIs));
        h(1,3) = cdfplot(DistFromSoma(basalInputROIs));
        set(h(1,1), 'Color', 'blue', 'LineWidth', 3);
        set(h(1,2), 'Color', 'red', 'LineWidth', 3);
        set(h(1,3), 'Color', 'green', 'LineWidth', 3);
        legend('All Inputs', 'All Apical Inpus', 'All Basal Inputs', 'Location', 'SouthEast')
    catch
        h(1,2) = cdfplot(DistFromSoma(basalInputROIs));
        set(h(1,1), 'Color', 'blue', 'LineWidth', 3);
        set(h(1,2), 'Color', 'green', 'LineWidth', 3);
        legend('All Inputs', 'All Basal Inputs', 'Location', 'SouthEast')
    end
    grid off
    title('')
    xlabel('Distance from start of branch')
    ylabel('Fraction of spines')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '21_DistFromBranchCumulInputs.png'))

    %Fig 22: Histogram of branch order for all inputs and separated for 
    %apical/basal
    figure(22)
    binEdges = linspace(1,max(BranchOrder), max(BranchOrder));
    subplot(1,3,1)
    histogram(BranchOrder(inputROIs),binEdges, 'Normalization','probability', 'FaceColor', 'blue')
    title('Inputs')
    box off
    subplot(1,3,2)
    histogram(BranchOrder(apicalInputROIs),binEdges, 'Normalization','probability', 'FaceColor', 'red')
    title('Apical inputs')
    box off
    subplot(1,3,3)
    histogram(BranchOrder(basalInputROIs),binEdges, 'Normalization','probability', 'FaceColor', 'green')
    title('Basal inputs')
    box off
    xlabel('Branch order')
    ylabel('Percentage of inputs')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '22_BranchOrderInputs.png'))

    %Fig 23: Cdfplot of branch order for all inputs and separated for 
    %apical/basal
    figure(23)
    try
        h(1,1) = cdfplot(BranchOrder(inputROIs)); hold on
        h(1,2) = cdfplot(BranchOrder(apicalInputROIs));
        h(1,3) = cdfplot(BranchOrder(basalInputROIs));
        set(h(1,1), 'Color', 'blue', 'LineWidth', 3);
        set(h(1,2), 'Color', 'red', 'LineWidth', 3);
        set(h(1,3), 'Color', 'green', 'LineWidth', 3);
    catch
        h(1,1) = cdfplot(BranchOrder(inputROIs)); hold on
        h(1,2) = cdfplot(BranchOrder(basalInputROIs));
        set(h(1,1), 'Color', 'blue', 'LineWidth', 3);
        set(h(1,2), 'Color', 'green', 'LineWidth', 3);
    end
    grid off
    title('Branch Order')
    xticks([1 2 3 4 5])
    xlabel('Distance from start of branch')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '23_BranchOrderCumulInputs.png'))

    %--------------------------------------------------------------------------
    %30s:Branch specific analysis
    %Fig 30: Number of inputs per dendrite
    figure(30)
    histogram(cell2mat({confData.allDendrites(suffSTEDcover).numInputs}), 'FaceColor', 'blue')
    ylabel('Number of dendrites')
    xlabel('Number of inputs on dendrite')
    title('Inputs per dendrite')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '30_Inputs per dendrite.png'))

    %Fig 31: Input density per dendrite 
    figure(31)
    subplot(1,2,1)
    distributionPlot(cell2mat({confData.allDendrites.InputDensity})', 'color', 'blue'); hold all
    boxplot(cell2mat({confData.allDendrites.InputDensity})')
    ylabel('Input dendity (Inputs/um)')
    title('Input density')
    box off

    %find branches that have more than one input
    multInputBranches = cell2mat(cellfun(@(x) x > 1, {confData.allDendrites.numInputs}, 'UniformOutput',false));
    subplot(1,2,2)
    distributionPlot(InputDensity(multInputBranches)', 'color', 'blue'); hold all
    boxplot(InputDensity(multInputBranches)')
    ylabel('Input dendity (Inputs/um) with more than 1 input')
    title('Input density')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '31_Inputs density.png'))

    %Fig 32: Number of inputs vs. length of STED coverage
    figure(32)
    allNumInputs = cell2mat({confData.allDendrites(suffSTEDcover).numInputs});
    allLength = cell2mat({confData.allDendrites(suffSTEDcover).length});
    scatter(allNumInputs,allLength)
    ylabel('Length of dendrite')
    xlabel('Number of inputs on dendrite')
    box off
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '32_Number of inputs vs length.png'))

    %--------------------------------------------------------------------------
    %40s - Clustering of inputs
    %Fig 40: Distance to neares input on same branch
    figure(40)
    distributionPlot(DistToInput, 'color', 'blue'); hold all
    boxplot(DistToInput)
    ylabel('Distance to nearest input on same branch')
    ylim([0 max(DistToInput)])
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '40_DistanceToNearestInput.png'))

    %Fig 41: Mean distance to nearest input
    figure(41)
    histogram(mean(DistToInputSim,1), 200,'FaceColor','Black')
    xlabel('Mean nearest-neighbor distance in um')
    xline(mean(DistToInput), '--r');
    legend('Random Data', 'Real data')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDir, '41_MeanDistanceToNearestInput.png'))
end

%% 7.) Cell summary
cell.totalSTEDcoverage = sum(cell2mat({confData.allDendrites.STEDLength}));
cell.PercSTEDcoverage = cell.totalSTEDcoverage/cell.totalLength*100;
cell.totalInputs = sum(cell2mat({confData.allDendrites.numInputs}));
cell.InputDensity =  cell.totalInputs/cell.totalSTEDcoverage;
cell.DistanceBetweenInputs = cell.totalSTEDcoverage/cell.totalInputs;
cell.PercOfInput = 100 * cell.totalInputs/(cell.PercSTEDcoverage*cell.totalSpines);
cell.PercApicalInputs = sum(apicalInputROIs)/sum(apicalSTED);
cell.PercBasalInputs = sum(basalInputROIs)/sum(basalSTED);

disp(['total cell length imaged with STED: ' num2str(cell.totalSTEDcoverage, '%.1f') ' ' char(181) 'm'])
disp(['Percentage of cell imaged with STED: ' num2str(cell.PercSTEDcoverage, '%.1f')])
disp(['total inputs: ' num2str(cell.totalInputs)])
disp(['Average input density: ' num2str(cell.InputDensity, '%.2f') ' Inputs per ' char(181) 'm'])
disp(['Average distance between inputs: ' num2str(cell.DistanceBetweenInputs, '%.2f') ' ' char(181) 'm'])
disp(['Percentage of input to the cell: ' num2str(cell.PercOfInput, '%.4f') ' %'])
disp(['Percentage of input to apical dendrites: ' num2str(cell.PercApicalInputs, '%.3f') ' %'])
disp(['Percentage of input to basal dendrites: ' num2str(cell.PercBasalInputs, '%.3f') ' %'])

%% 8.) Save
InputData = confData.ROIs;
STEDTracing = confData.allDendrites;
save([saveDirAna filesep '03_InputCellReconstruction.mat'], 'cell','InputData', 'inputROINr', 'STEDTracing','-mat') 

%write to text file
fid = fopen([saveDir filesep 'AnatomyInputsSummary.txt'], 'w');
formatSpec = ['total cell length imaged with STED: %.1f ' char(181) 'm\n'...
    'Percentage of cell imaged with STED: %.1f % \n'...
    'total number of inputs: %d \n'...
    'average input density: %.2f inputs per ' char(181) 'm\n'...
    'Average distance between inputs: %.2f' char(181) 'm',...
    'Percentage of input to the cell: %.3f %',...
    'Percentage of input to the apical dendrites: %.3f %',...
    'Percentage of input to the basal dendrites: %.3f %'];
fprintf(fid, formatSpec, cell.totalSTEDcoverage, cell.PercSTEDcoverage, cell.totalInputs, cell.InputDensity,cell.DistanceBetweenInputs, cell.PercOfInput, cell.PercApicalInputs,cell.PercBasalInputs);
fclose(fid);