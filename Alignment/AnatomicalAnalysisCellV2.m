function AnatomicalAnalysisCellV2(animal, varargin)

%Loads data of the confocal reconstruction of the cell, aligning if
%necessary several slides, plotting the reconstruction of the cell as well
%as doing some basic analysis of the structure and plotting those

%Inputs:
% - animal: which animal are we analysing?
% - cell: if there are multiple cells per animal, which one are we looking
% at?

%STEPS:
%0.) Define and if necessary, add the confocal and save dir
%1.) Load the confocal representation
%2.) Define certain ROI categories as logical indicators (apical/basal)
%3.) Align ROIs
%4.) Specify groups of dendrites for quicker access for plotting
%5.) Analysis
%   - a) Cell-wide: Distance from soma, Branch order
%   - b) Dendrite-specific: hoch many spines, density
%6.) Plot 
%7.) Cell summary: Total length, total spines, spine density
%8.) Save

%Outputs:
% - Folder 01_MorphologyCellReconstruction with Plots (see details at the
% beginning of Step 6
% - 01_MorphologyCellReconstruction.mat with variables:
%       - confData: containing all data of a single ROI, dendrites, and
%       branches
%       - template: confocal projection of one of the slices
%       - cell: cell summary
% - CellSummary.txt: printed version of cell summary

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.1: March, 2024

if nargin > 1
    cell = varargin{1};
else
    cell = [];
end

%% 0.) Define and if necessary, add the save dir

confDir = ['Z:\Juliane\InputAnalysis\' char(animal) filesep char(cell) '\B - Confocal\'];
saveDirAna = ['Z:\Juliane\InputAnalysis\' char(animal) filesep char(cell) '\E - Analysis\'];
if ~exist(saveDirAna,'dir')
    mkdir(saveDirAna)
end
saveDir = ([saveDirAna filesep '01_Morphology analysis']);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%% 1.) Load data
%test if there is one or multiple slices
sliceFolders = dir([confDir 'SL*']);
isFolder=([sliceFolders.isdir]); %make sure those are folders, and not the .jpgs or .tif files
if ~isempty(sliceFolders(isFolder))
    %load confData from all slices
    confData = struct('ROIs', [], 'Branches', [], 'allDendrites',[], 'Soma', []);
    numBranches = [];
    somaPos = [];
    for s = 1:length(sliceFolders)
        try
            temp = load([confDir filesep sliceFolders(s).name filesep 'cellReconstruction.mat']);
            %add in the information from which slice the ROI/dendrite is
            temp.ROIs = arrayfun(@(x) setfield(x, 'SL', s), temp.ROIs);
            temp.Branches= arrayfun(@(x) setfield(x, 'SL', s), temp.Branches);
            temp.allDendrites= arrayfun(@(x) setfield(x, 'SL', s), temp.allDendrites);

            %save to the new structure containing all ROIs
            confData.ROIs = [confData.ROIs temp.ROIs];
            confData.Branches = [confData.Branches temp.Branches];
            confData.allDendrites = [confData.allDendrites temp.allDendrites];
            numBranches = [numBranches length(temp.Branches)];
            
            somaPos = [somaPos; temp.Soma.xPos temp.Soma.yPos];
            clear temp
        catch
            disp(['No reconstruction for ' sliceFolders(s).name])
        end
    end
    
    %now get the image and soma info from the one that has the most
    %branches
    try
        [~,somaSlice] = max(numBranches);
        somaInfo = load([confDir filesep sliceFolders(somaSlice).name filesep 'cellReconstruction.mat']);
    catch
        somaSlice = somaSlice + 1;
        somaInfo = load([confDir filesep sliceFolders(somaSlice).name filesep 'cellReconstruction.mat']);
    end
    confData.Soma = somaInfo.Soma; clear somaInfo;

else
    %confocal reconstruction
    confData = load([confDir filesep 'cellReconstruction.mat']);
    confData.ROIs = arrayfun(@(x) setfield(x, 'SL', 1), confData.ROIs);
    confData.Branches= arrayfun(@(x) setfield(x, 'SL', 1), confData.Branches);
    confData.allDendrites= arrayfun(@(x) setfield(x, 'SL', 1), confData.allDendrites);
end

%% 2.) Calculate shifts for every branch
%for every branch, get the start of the branch and align it to the soma, so
%that the branch starts at 0/0 and goes from there
confData.Branches= arrayfun(@(x) setfield(x, 'BranchNr', []), confData.Branches);
confData.Branches= arrayfun(@(x) setfield(x, 'BranchStart', []), confData.Branches);
for b = 1:length(confData.Branches)
    confData.Branches(b).BranchNr = confData.Branches(b).dendrites(1).Branch;
    confData.Branches(b).BranchStart = confData.Branches(b).dendrites(1).pixelCoord(1,:);
end

%% 3.) Align branches and ROIs and plot on top of tracings
%shift dendrites according to the branch start 
for d = 1:length(confData.allDendrites)
    rightBranch = find([confData.Branches.BranchNr] == confData.allDendrites(d).Branch);
    confData.allDendrites(d).normCoord(:,1) = confData.allDendrites(d).pixelCoord(:,1)-confData.Branches(rightBranch).BranchStart(1);
    confData.allDendrites(d).normCoord(:,2) = confData.allDendrites(d).pixelCoord(:,2)-confData.Branches(rightBranch).BranchStart(2);
end

%now shift the slice ROIs to align with the branch start
for r = 1:length(confData.ROIs)
    rightBranch = find([confData.Branches.BranchNr] == confData.ROIs(r).Branch);
    confData.ROIs(r).xPosNorm = confData.ROIs(r).xPos -confData.Branches(rightBranch).BranchStart(1);
    confData.ROIs(r).yPosNorm = confData.ROIs(r).yPos -confData.Branches(rightBranch).BranchStart(2);
end

%% 4.) Specify groups of ROIs and dendrites for quicker access for plotting
%apical vs. basal ROIs
apicalROIs = cellfun(@(x) strcmp(x, 'apical'),{confData.ROIs.type}); 
basalROIs = cellfun(@(x) strcmp(x, 'basal'),{confData.ROIs.type});

%apical vs. basal dendrites
apicalDend = cellfun(@(x) strcmp(x, 'apical'),{confData.allDendrites.type}); 
basalDend = cellfun(@(x) strcmp(x, 'basal'),{confData.allDendrites.type}); 

%% 5.) Analysis
%a) Cell wide analysis
DistFromSoma = cell2mat({confData.ROIs.distToSoma});
BranchOrder = cell2mat({confData.ROIs.BranchOrder});

%b) Dendrite-specific analysis: how many spines? density?
for den = 1:length(confData.allDendrites)
    if confData.allDendrites(den).length > 10
        ROIsOnDendrites = find([confData.ROIs.Dendrite] == den);
        if ~isempty(ROIsOnDendrites)
            confData.allDendrites(den).numSpines = length(ROIsOnDendrites);
            confData.allDendrites(den).SpineDensity = length(ROIsOnDendrites)/confData.allDendrites(den).length;
        else
            confData.allDendrites(den).numSpines = 0;
            confData.allDendrites(den).SpineDensity = NaN;
        end
    else
        confData.allDendrites(den).numSpines = NaN;
        confData.allDendrites(den).SpineDensity = NaN;
    end
end

numSpines = cell2mat({confData.allDendrites.numSpines});
SpineDensity = cell2mat({confData.allDendrites.SpineDensity});

%% 6.) Plot
%Figure Overview
%10s - ROIs on top of tracing
%Fig 10: All ROIs on top of tracing color coded by apical vs. basal
%Fig 11: All ROIs on top of tracing, separated by slice

%20s - Cell wide analysis
%Fig 20: Distribution of distance from soma for all and separated for
%apical/basal
%Fig 21: Same as Fig 20, but as cdfplot
%Fig 22: Histogram of branch order for all and separated for apical/basal
%Fig 23: same as Fig 22, but as cdfplot

%30s - Branch specific analysis
%Fig 30: Number of spines per dendrite for all and separated for
%apical/basal
%Fig 31: Spine density per dendrite for all and separated for apical/basal

%--------------------------------------------------------------------------
%10s: Where are the ROIs
%plot them by type (apical/basal)
plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, ones(length(confData.ROIs),1),confData.Soma, 'type', 0, jet(2), 10, 0)
saveas(gcf, fullfile(saveDir, '10_CellTracing separated by type.png'))

%plot them by slice
if ~isempty(sliceFolders(isFolder))
    plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, ones(length(confData.ROIs),1),confData.Soma, 'SL', 0, jet(s), 11, 0)
    saveas(gcf, fullfile(saveDir, '11_CellTracing separted by SL.png'))
end

%plot them by branch
plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, ones(length(confData.ROIs),1),confData.Soma, 'Branch', 0, cbrewer('qual', 'Paired', 12), 12, 0)
saveas(gcf, fullfile(saveDir, '12_CellTracing separated by branch.png'))

%--------------------------------------------------------------------------
%20s: Cell wide analysis
% Distance from soma, branch order

%Fig 20: Distribution of distance from soma for all and separated for
%apical/basal
figure(20)
subplot(1,3,1)
distributionPlot(DistFromSoma', 'color', 'black'); hold all
boxplot(DistFromSoma')
title('All spines')
ylabel('Distance from branch start in um')
box off
ylim([0 max(DistFromSoma)])
if sum(apicalROIs)>0
    subplot(1,3,2)
    distributionPlot(DistFromSoma(apicalROIs)', 'color', 'red'); hold all
    boxplot(DistFromSoma(apicalROIs)')
    title('Apical Spines')
    ylim([0 max(DistFromSoma)])
    axis off
end
if sum(basalROIs)>0
    subplot(1,3,3)
    distributionPlot(DistFromSoma(basalROIs)', 'color', 'green'); hold all
    boxplot(DistFromSoma(basalROIs)')
    title('basal Spines')
    ylim([0 max(DistFromSoma)])
    axis off
end
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '20_SpineDistanceFromBranch.png'))

%Fig 21: Cdfplot of distance from soma for all and separated for
%apical/basal
figure(21);
h(1,1) = cdfplot(DistFromSoma); hold on
set(h(1,1), 'Color', 'black', 'LineWidth', 3);
if sum(apicalROIs)>0
    h(1,2) = cdfplot(DistFromSoma(apicalROIs));
    set(h(1,2), 'Color', 'red', 'LineWidth', 3);
end
if sum(basalROIs)>0
    h(1,3) = cdfplot(DistFromSoma(basalROIs));
    set(h(1,3), 'Color', 'green', 'LineWidth', 3);
end
grid off
title('')
xlabel('Distance from start of branch')
ylabel('Fraction of spines')
set(gcf, 'color', 'w');
legend('All spines', 'Apical Spines', 'Basal Spines', 'Location', 'SouthEast')
saveas(gcf, fullfile(saveDir, '21_DistFromBranchCumul.png'))

%Fig 22: Histogram of branch order for all and separated for apical/basal
figure(22)
binEdges = linspace(1,max(BranchOrder), max(BranchOrder));
subplot(1,3,1)
histogram(BranchOrder,binEdges, 'Normalization','probability', 'FaceColor', 'white')
title('All spines')
box off
if sum(apicalROIs)>0
    subplot(1,3,2)
    histogram(BranchOrder(apicalROIs),binEdges, 'Normalization','probability', 'FaceColor', 'red')
    title('Apical spines')
    box off
end
if sum(basalROIs)>0
    subplot(1,3,3)
    histogram(BranchOrder(basalROIs),binEdges, 'Normalization','probability', 'FaceColor', 'green')
    title('Basal spines')
    box off
end
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '22_BranchOrder.png'))

%Fig 22: Cdfplot of branch order for all and separated for apical/basal
h(1,1) = cdfplot(BranchOrder); hold on
set(h(1,1), 'Color', 'black', 'LineWidth', 3);

if sum(apicalROIs)>0
    h(1,2) = cdfplot(BranchOrder(apicalROIs));
    set(h(1,2), 'Color', 'red', 'LineWidth', 3);
end
if sum(basalROIs)>0
    h(1,3) = cdfplot(BranchOrder(basalROIs));
    set(h(1,3), 'Color', 'green', 'LineWidth', 3);
end
grid off
title('Branch Order')
xticks([1 2 3 4 5])
xlabel('Distance from start of branch')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '23_BranchOrderCumul.png'))

%--------------------------------------------------------------------------
%30s: Branch specific analysis
% Number of spines and spine density per dendrite

%Fig 30: Number of spines per dendrite for all and separated for
%apical/basal
figure(30)
subplot(1,3,1)
distributionPlot(numSpines', 'color', 'black'); hold all
boxplot(numSpines')
title('All dendrites')
ylim([0 150])
ylabel('Spines per dendrite')
box off
if sum(apicalROIs)>0
    subplot(1,3,2)
    distributionPlot(numSpines(apicalDend)', 'color', 'red'); hold all
    boxplot(numSpines(apicalDend)')
    title('Apical dendrites')
    ylim([0 150])
    box off
    axis off
end
if sum(basalROIs)>0
    subplot(1,3,3)
    distributionPlot(numSpines(basalDend)', 'color', 'green'); hold all
    boxplot(numSpines(basalDend)')
    ylim([0 150])
    title('Basal dendrites')
    box off
    axis off
end
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '30_Spines per dendrite.png'))

%Fig 31: Spine density per dendrite for all and separated for apical/basal
figure(31)
subplot(1,3,1)
distributionPlot(SpineDensity', 'color', 'black'); hold all
boxplot(SpineDensity')
title('All dendrites')
ylim([0 10])
ylabel('Spine density')
box off
if sum(apicalROIs)>0
    subplot(1,3,2)
    distributionPlot(SpineDensity(apicalDend)', 'color', 'red'); hold all
    boxplot(SpineDensity(apicalDend)')
    title('Apical dendrites')
    ylim([0 10])
    box off
    axis off
end
if sum(basalROIs)>0
    subplot(1,3,3)
    distributionPlot(SpineDensity(basalDend)', 'color', 'green'); hold all
    boxplot(SpineDensity(basalDend)')
    ylim([0 10])
    title('Basal dendrites')
    box off
    axis off
end
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '31_Spine density.png'))

%% 7.) Cell summary
cell = struct;
cell.totalLength = sum(cell2mat({confData.allDendrites.length}), 'omitnan');
cell.totalSpines = sum(cell2mat({confData.allDendrites.numSpines}), 'omitnan');
cell.SpineDensity = cell.totalSpines/cell.totalLength;

disp(['total dendrite length: ' num2str(cell.totalLength, '%.1f') ' ' char(181) 'm'])
disp(['total number of spines: ' num2str(cell.totalSpines)])
disp(['Average spine density: ' num2str(cell.SpineDensity, '%.2f') ' Spines per ' char(181) 'm'])

%% 8.) Save
save([saveDirAna filesep '01_MorphologyCellReconstruction.mat'], 'confData','cell','-mat') 

%write to text file
fid = fopen([saveDir filesep 'CellSummary.txt'], 'w');
formatSpec = ['total dendrite length: %.1f ' char(181) 'm\ntotal number of spines: %d \naverage spine density: %.2f spines per ' char(181) 'm'];
fprintf(fid, formatSpec, cell.totalLength, cell.totalSpines, cell.SpineDensity);
fclose(fid);

close all
