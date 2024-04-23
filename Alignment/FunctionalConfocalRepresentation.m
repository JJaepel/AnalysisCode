function FunctionalConfocalRepresentation(animal, varargin)

%Loads data of the confocal reconstruction of the cell, all 2p experimental
%data for all ROIs as well as the matching files to combine all the data,
%plotting the representation of a cell as well some basic analysis of the
%whole cell's functional data

%Inputs:
% - animal: which animal are we analysing?
% - cell: if there are multiple cells per animal, which one are we looking
% at?

%STEPS:
%0.) Define and if necessary, add the confocal and save dir
%1.) Load the confocal representation and the matching
%2.) Remove duplicates from multiple 2p experiments
%3.) Add 2p data to other confocal information
%4.) Specify groups of ROIs for quicker access for plotting
%5.) Futher analysis of the functional data
%   - a) Synaptic aggregate for ori/dir-selective ROIs
%   - b) Spine specific characteristic of environment (local ori/dir)
%   - c) Branch specific characteristic of dendritic segments
%   - d) Pairwise measurement of spines on the same dendritic segment
%6.) Plot
%7.) Save

%Outputs:
% - Folder 02_FunctionalCellReconstruction with Plots (see details at the
% beginning of Step 6)
% - 02_FunctionalCellReconstruction.mat with variables:
%       - confData: containing all data of a single ROI, dendrites, and
%       branches
%       - Grouping: contains the structure of all grouped 2p and confIDs
%       - Selectors: groups of ROIs with specific characteristics
%       - pwMeasures: pairwise measurements of all ROIs on the same
%       dendritic segment

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023


if nargin > 1
    cell = varargin{1};
else
    cell = [];
end

%switchboard
compareFirstLastAmount = 20;

%% 0.) %define and if necessary, add the confocal and save dir 
baseDir = ['Z:\Juliane\InputAnalysis\' char(animal) filesep char(cell)];
confDir = [baseDir '\B - Confocal\'];
funcDir = [baseDir '\A - 2p Imaging\'];

saveDirAna = [baseDir '\E - Analysis\'];
saveDir = ([saveDirAna filesep '02_FunctionalCellReconstruction']);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%% 1.) Load the confocal representation and the matching
%load the morphology reconstruction
load([saveDirAna filesep '01_MorphologyCellReconstruction.mat'], 'confData')

%load the groupings
%test if there is one or multiple slices
sliceFolders = dir([confDir 'SL*']);
isFolder=([sliceFolders.isdir]); %make sure those are folders, and not the .jpgs or .tif files
if ~isempty(sliceFolders(isFolder))
    GroupingFiles = [];
    for s = 1:length(sliceFolders)
        try
            tempGF = dir([baseDir filesep 'D - Alignments' filesep sliceFolders(s).name filesep 'Alignment 2p\Groups*.mat']);
            GroupingFiles = [GroupingFiles; tempGF];
        catch
            disp([sliceFolders(s).name ' not aligned yet'])
        end
    end
else
    GroupingFiles = dir([baseDir '\D - Alignments\Alignment 2p\Groups*.mat']);
end

%rewrite into grouping structure, consisting of ConfID, TwoPID, whether it
%is a good spine, maxResp and funcData (the last three needed for removing
%duplicates)
Grouping = struct;
counter = 1;
for b = 1:size(confData.Branches,2)
    allDendrites = unique([confData.allDendrites.Branch]);
    dendriteNr = allDendrites(b);
    %find the relevant matching files
    filesDend = cellfun(@(x) contains(x, ['Dendrite' sprintf('%02d',dendriteNr)]),{GroupingFiles.name}); %which are the files for that dendrite?
    for f = 1:length(filesDend)
        if filesDend(f) == 1
           GroupData = load([GroupingFiles(f).folder filesep GroupingFiles(f).name]); 
           %assumptions: 1) we never have more than 199 ROIs in a 2p experiment
           %-> convert the 2p ids into xxyy with xx = expNr*2 and yy = ROI nr
           %2) we never have more than 999 ROIs in a confocal dendrite and
           %not more than 99 dendrites
           %-> convert the confocal ids into xxyyy with xx = dendrite nr
           %and yy = ROI nr
           twoPnumber= strfind(GroupingFiles(f).name, 't000');
           baseTwoP = str2double(GroupingFiles(f).name(twoPnumber+1:twoPnumber+5))*200;
           
           %load the functional data for that exp
           expData = load([funcDir num2str(GroupingFiles(f).name(twoPnumber:twoPnumber+5)) filesep 'ROIsAna.mat']);
           
           for g = 1:length(GroupData.MatchData)
               Grouping(counter).ConfID = str2double(GroupData.MatchData(g).ConfocalID) + dendriteNr*1000;
               twoPID = str2double(GroupData.MatchData(g).TwoPID);
               Grouping(counter).TwoPID = twoPID + baseTwoP;
               Grouping(counter).good = expData.ce(twoPID).good;
               Grouping(counter).maxResp = max(expData.ce(twoPID).meanResp);
               Grouping(counter).funcData = expData.ce(twoPID);
               counter = counter+1;
           end
           
        end
    end
end

%% 2.) Remove duplicates (look which one responsiv and then which one has bigger amplitude), and sort 
[uniqueConfIDs, ~, uidx] = unique(cellfun(@(x) x, {Grouping.ConfID})); %find the unique ids
duplicates = find(accumarray(uidx, 1) > 1); %find the once that occur more often than once 

for d = 1:length(duplicates)
    %get the confIDs
    dConfID = uniqueConfIDs(duplicates(d));
    %find which groups they are
    dGroups = find(cellfun(@(x) x,{Grouping.ConfID})== dConfID);
    %those groups - are they good spines?
    goodGroups = find([Grouping(dGroups).good] == 1);
    if isempty(goodGroups) %if there are no good groups, take the one with the higer resp
        [~,higherAmpGroup] = max([Grouping(dGroups).maxResp]);
        selectedGroup = dGroups(higherAmpGroup);
    else
        if length(goodGroups) > 1  %if there is more than one good groups, take the one with the higer resp
            [~,higerAmpGroup] = max([Grouping(goodGroups).maxResp]);
            selectedGroup = dGroups(higerAmpGroup);
            
        else %take the only good group
            selectedGroup = dGroups(goodGroups);
        end
    end
    %remove the non-selectedGroups
    Grouping(dGroups(dGroups ~= selectedGroup)) = [];
end

[~,idx]=sort([Grouping.ConfID]);
Grouping=Grouping(idx);

%% 3.) Add 2p data to other confocal information
%load confocal reconstruction & cell summary
%confData = load([confDir filesep 'cellReconstruction.mat']);
TwoPFile = dir([funcDir filesep animal '_CellSummary*.mat']);
TwoPData = load([TwoPFile(1).folder filesep TwoPFile(1).name]);

%let's check if there are multiple cells from that animal
if size(TwoPData.completeCell,2) > 1
   if ~isempty(cell) %do we want to look at a specific cell and we have multiple cells to look at?
       cellInfo = animalParser;
       cellAnimalInfo = find(~cellfun(@isempty,cellfun(@(x) find(contains(x, animal)),cellInfo.animal,'UniformOutput',false))); %which ones contain the animal
       cellNameInfo = find(strcmp(cellInfo.cellName, cell), 1); %which ones contain the animal
       i = intersect(cellAnimalInfo, cellNameInfo);
       cellInd =cellInfo.cellNr{i};
   else %otherwise, let's just take the one that has the most spines
       try 
           spineNrs(1) = size(TwoPData.completeCell{1}.SpineData,2);
       catch 
           spineNrs(1) = 0;
       end
       try 
           spineNrs(2)  = size(TwoPData.completeCell{2}.SpineData,2);
       catch 
           spineNrs(2)  = 0;
       end
       [~, cellInd] = max(spineNrs);
   end
else
    cellInd = 1; %otherwise, set all of those to one
end

%add rows to quickly check if there is a 2p macth
confData.ROIs = arrayfun(@(x) setfield(x, 'TwoPMatch', zeros(1,1)), confData.ROIs);
confData.ROIs = arrayfun(@(x) setfield(x, 'good', zeros(1,1)), confData.ROIs);
confData.ROIs = arrayfun(@(x) setfield(x, 'expNr', NaN(1,1)), confData.ROIs);
%add info into it
dendriteROIs = []; %make this variable to check if there were any dendrites matched to the confocal ROIs
for g = 1:length(Grouping)
    %find the confocal ROI
    Branch = floor(Grouping(g).ConfID/1000);
    ROIonBranch = Grouping(g).ConfID-Branch*1000;
    sameBranch = find([confData.ROIs.Branch] == Branch);
    sameROINr = find([confData.ROIs.ROINrOnBranch] == ROIonBranch);
    row = intersect(sameBranch, sameROINr);
    if ~isempty(row)
        confData.ROIs(row).ConfID = Grouping(g).ConfID;

        %find the 2pID
        sameROI = find([TwoPData.completeCell{cellInd}.SpineData.uniqueID] == Grouping(g).TwoPID);

        %add the functional info
        confData.ROIs(row).TwoPMatch = 1;
        confData.ROIs(row).TwoPID = Grouping(g).TwoPID;
        confData.ROIs(row).good = Grouping(g).good;
        try
            %specific measurements
            confData.ROIs(row).cycRes = TwoPData.completeCell{cellInd}.SpineData(sameROI).cycRes;
            confData.ROIs(row).meanResp = TwoPData.completeCell{cellInd}.SpineData(sameROI).meanResp;
            confData.ROIs(row).prefOri = TwoPData.completeCell{cellInd}.SpineData(sameROI).prefOri;
            confData.ROIs(row).prefDir = TwoPData.completeCell{cellInd}.SpineData(sameROI).prefDir;
            confData.ROIs(row).OSI = TwoPData.completeCell{cellInd}.SpineData(sameROI).OSI;
            confData.ROIs(row).DSI = TwoPData.completeCell{cellInd}.SpineData(sameROI).DSI;
            confData.ROIs(row).DSIvect = TwoPData.completeCell{cellInd}.SpineData(sameROI).DSIvect;
            confData.ROIs(row).Bandwidth = TwoPData.completeCell{cellInd}.SpineData(sameROI).Bandwidth;

            %make sure everything else that might be interesting is saved as well
            confData.ROIs(row).funcData = TwoPData.completeCell{cellInd}.SpineData(sameROI);

            %add expNr
            confData.ROIs(row).expNr = TwoPData.completeCell{cellInd}.SpineData(sameROI).expNr;
        catch
            disp(['No functional data for group ' num2str(g) '- please check!'])
            confData.ROIs(row).TwoPMatch = 0;
            if isempty(dendriteROIs)
                fid = fopen([saveDir filesep 'DendriteROIs.txt'], 'w');
                formatSpec = 'ROIs that are dendrites in 2P: \nROI Nr: %d  \n';
                fprintf(fid, formatSpec, g');
                fclose(fid);
            else
                fid = fopen([saveDir filesep 'DendriteROIs.txt'], 'a');
                formatSpec = 'ROI Nr: %d  \n';
                fprintf(fid, formatSpec, g');
                fclose(fid);
            end
            dendriteROIs = [dendriteROIs g];
        end
    else 
        disp('No confocal ROI found')
    end
end

%add in cell info
confData.Soma.cyc = TwoPData.completeCell{cellInd}.cyc;
confData.Soma.meanResp = TwoPData.completeCell{cellInd}.meanResp;
confData.Soma.peakFit = TwoPData.completeCell{cellInd}.peakFit;
confData.Soma.prefOri = TwoPData.completeCell{cellInd}.prefOri;
confData.Soma.prefDir = TwoPData.completeCell{cellInd}.prefDir;
confData.Soma.OSI = TwoPData.completeCell{cellInd}.OSI;
confData.Soma.DSI = TwoPData.completeCell{cellInd}.DSI;
confData.Soma.DSIvect = TwoPData.completeCell{cellInd}.DSIvect;
confData.Soma.Bandwidth = TwoPData.completeCell{cellInd}.Bandwidth;

%% 4.) Specify functional groups 
TwoPROIs = ([confData.ROIs.TwoPMatch] == 1); %all the ROIs with a functional match, size = allSpines
TwoPROIsNR = find([confData.ROIs.TwoPMatch] ==1);
respROIs = ([confData.ROIs.good] == 1); %all ROIs with a good one,size = allSpines
respROINrs = find([confData.ROIs.good] == 1);

oriSelect = find([confData.ROIs.OSI] > 0.1); %all oriselectROIs, size = allfuncSpines
oriSelectROIsNr = TwoPROIsNR(oriSelect);%all oriselectROIs, size = allSpines
oriSelectROIsNr = intersect(oriSelectROIsNr,respROINrs); %make sure they are also good ones
dirSelect = find([confData.ROIs.DSI] > 0.1); %all dirSelectROIs, size = allFuncSpines
dirSelectROIsNr = TwoPROIsNR(dirSelect);%all oriselectROIs, size = allSpines
dirSelectROIsNr = intersect(dirSelectROIsNr,respROINrs); %make sure they are also good ones

oriSelectROIs = zeros(1,length(confData.ROIs));
oriSelectROIs(oriSelectROIsNr) = 1; %logical vector
dirSelectROIs = zeros(1,length(confData.ROIs));
dirSelectROIs(dirSelectROIsNr) = 1; %logical vector

%Find the first 100 vs. last 100 responsive spines
TwoPIDrespROIs = [confData.ROIs(oriSelectROIsNr).TwoPID];
[~, sortID] = sort(TwoPIDrespROIs);

if size(oriSelectROIsNr,2) > compareFirstLastAmount*2+1
    firstROIs = oriSelectROIsNr(sortID(1:compareFirstLastAmount));
    logicalFirstRespROIs = zeros(1,length(confData.ROIs));
    logicalFirstRespROIs(firstROIs) = 1;
    lastROIs = oriSelectROIsNr(sortID(end-compareFirstLastAmount:end));
    logicalLastRespROIs = zeros(1,length(confData.ROIs));
    logicalLastRespROIs(lastROIs) = 1;
end

%% 5.) Further analysis
% 5.a) Synaptic aggregate
angles = linspace(0,180,9);
angles = angles(1:end-1);
anglesAll = linspace(0,360,17);
anglesAll = anglesAll(1:end-1);

%get the soma tuning
SomaResp = (confData.Soma.meanResp(1:8)+confData.Soma.meanResp(9:16))/2; %fold in ori space
SomaRespNorm = (SomaResp)/(max(SomaResp)); %normalize
SomaDirNorm = confData.Soma.meanResp/(max(confData.Soma.meanResp));

%calculate the aggregate for the ori-selective/dir-selective ROIs
[~, curveSumFunc, ~] = calcSumSpineResponses(confData.ROIs, oriSelectROIsNr);
[~, curveSumFuncDir, ~] = calcSumSpineResponsesDir(confData.ROIs, oriSelectROIsNr);

if size(oriSelectROIsNr,2) > compareFirstLastAmount*2+1
    %lets compare the first 100 vs. the last 100
    [~, curveSumFirst, ~] = calcSumSpineResponses(confData.ROIs, firstROIs);
    [~, curveSumLast, ~] = calcSumSpineResponses(confData.ROIs, lastROIs);
    
    %how would it look like if we randomly pick 100?
    curveSumRandFirst = NaN(8,1000);
    curveSumRandLast = NaN(8,1000);
    for r = 1:1000
        newSort = randperm(length(sortID));
        firstRandROIs = oriSelectROIsNr(newSort(1:compareFirstLastAmount));
        lastRandROIs = oriSelectROIsNr(newSort(end-compareFirstLastAmount:end));
        [~, curveSumRandFirst(:,r), ~] = calcSumSpineResponses(confData.ROIs, firstRandROIs);
        [~, curveSumRandLast(:,r), ~] = calcSumSpineResponses(confData.ROIs, lastRandROIs);
    end
end

% 5.b) Spine specific characteristics of the environment
%Homeogeneity index, local dispersion, deltaOri vs. dist, nearest neighbor
%with similar pref, how different is environment to soma, meanOSI
confData.ROIs = calculateLocalOriMeasurements(confData.ROIs, confData.Soma);
confData.ROIs = calculateLocalDirMeasurements(confData.ROIs,confData.Soma);

% 5.c) Branch specific characteristics the dendritic segments
confData.allDendrites = calculateBranchMeasurements(confData.allDendrites, confData.ROIs, confData.Soma);

% 5.d) Pairwise measurements
pwMeasures= calculatePairwiseParameters(confData.ROIs);

%% 6.) Plot
%Figure Overview
%10s - ROIs on top of tracing
%Fig 10: Which confocal ones were characterized?
%Fig 11: Which confocal ones were responsive?

%20s - Properties of ROIs on top of tracing
%Fig 20: What is the ori pref of the responsive spines?
%Fig 21: What is the dir pref of the responsive spines?
%Fig 22: What is the OSI of the responsive spines?
%Fig 23: What is the delta Ori of the responsive spines?
%Fig 24: What is the DSI of the responsive spines?
%Fig 25: What is the bandwidth of the responsive spines?
%Fig 26: What is the ori pref for the first/last 100 resp spines?

%30s - Synaptic aggregate & histogram of prefOri/prefDir or as polarPlots
%Fig 30: Synaptic aggregate vs. cell in orientation space
%Fig 31: Synaptic aggregate of first 100 vs. last 100 vs. cell
%Fig 32: Synaptic aggregat of random first 100 vs. random last 100 vs. cell
%Fig 33: Distribution of prefOri as histogram
%Fig 34: Distribution of prefOri as polarPlot
%Fig 35: Synaptic aggregate vs. cell in direction space
%Fig 36: Distribution of prefDir as histogram
%Fig 37: Distribution of prefDir as polarPlot
%Fig 38: Distribution of bandwidth

%40s & 50s: Local environment ori(50s) and dir (60s)
%Fig 40/50: local (dir) Dispersion
%Fig 41/51: HI (dir)
%Fig 42/52: local (dir) dispersion vs. HI (dir)
%Fig 43/53: local deltaOri/Dir
%Fig 44/54: local delteOriSoma/deltaDirSOma
%Fig 45/55: nearest SimilarPrefOri/dir
%Fig 46/56: local OSI/DSI

%60s: Properties of dendritic segments
%Fig 60: Branch circular dispersion ori/dir
%Fig 61: Branch selectivity
%Fig 62: delta Ori/delta Dir
%Fig 63: Branch circular dispersion vs. deltaOri
%Fig 64: Branch circular dispersion dir vs. deltaDir

%70s: Pairwise distances
%Fig 70: Pref Ori
%Fig 71: Pref Dir
%Fig 72: OSI
%Fig 73: DSI
%Fig 74: DSIVect
%Fig 75: Tuning correlation
%Fig 76: Trial-to-trial correlation

%% Start plotting
%--------------------------------------------------------------------------
%10s: Where are the ROIs
plotIDsonTracing(confData.allDendrites, confData.ROIs, TwoPROIs, 10, 'blue')
saveas(gcf, fullfile(saveDir, '10_CellImage_TwoPROIs.png'))

plotIDsonTracing(confData.allDendrites, confData.ROIs, respROIs, 11, 'red')
saveas(gcf, fullfile(saveDir, '11_CellImage_RespROIs.png'))

%--------------------------------------------------------------------------
%20s: ROIs on top of tracing with cell function
%What is the ori pref, delta ori, OSI
plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, oriSelectROIs, confData.Soma, 'prefOri', 1, hsv(180), 20, 0)
saveas(gcf, fullfile(saveDir, '20_CellTracing_PrefOri.png'))

plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, dirSelectROIs, confData.Soma, 'prefDir', 1, hsv(360), 21, 0)
saveas(gcf, fullfile(saveDir, '21_CellTracing_PrefDir.png'))

plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, respROIs, confData.Soma, 'OSI', 0, inferno(100), 22, 0)
saveas(gcf, fullfile(saveDir, '22_CellTracing_OSI.png'))

plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, oriSelectROIs, confData.Soma,  'deltaOri', 0, inferno(9000), 23, 0)
saveas(gcf, fullfile(saveDir, '23_CellTracing_deltaOri.png'))

plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, respROIs, confData.Soma, 'DSI', 0, viridis(100), 24, 0)
saveas(gcf, fullfile(saveDir, '24_CellTracing_DSI.png'))

plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, respROIs, confData.Soma, 'Bandwidth', 0, jet(9000), 25, 0)
saveas(gcf, fullfile(saveDir, '25_CellTracing_Bandwidth.png'))

if size(oriSelectROIsNr,2) > compareFirstLastAmount*2+1
    plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, logicalFirstRespROIs, confData.Soma, 'prefOri', 1, hsv(180), 26, 0)
    saveas(gcf, fullfile(saveDir, '26_CellTracing_PrefOriFirst.png'))

    plotFieldsOnCellTracing(confData.allDendrites, confData.ROIs, logicalLastRespROIs, confData.Soma, 'prefOri', 1, hsv(180), 27, 0)
    saveas(gcf, fullfile(saveDir, '27_CellTracing_PrefOriLast.png'))
end

%--------------------------------------------------------------------------
%30s:What is the aggregate of the responses
%What is the summed responses, also in comparison to start of experiment
%and random selected subset

%Fig 30: Synaptic aggregate vs. cell in orientation space
figure(30)
plot(angles, SomaRespNorm, 'Color', 'black', 'LineWidth', 5);
hold on
plot(angles, curveSumFunc, 'Color', 'red', 'LineWidth', 2);
ylabel('normalized Response')
xlabel('Angle in deg')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '30_SummedResp.png'))

if size(oriSelectROIsNr,2) > compareFirstLastAmount*2+1
    %Fig 31: Synaptic aggregate of first 100 vs. last 100 vs. cell
    figure(31)
    plot(angles, SomaRespNorm, 'Color', 'black', 'LineWidth', 5);
    hold on
    plot(angles, curveSumFirst, 'Color', 'green', 'LineWidth', 2);
    hold on
    plot(angles, curveSumLast, 'Color', 'blue', 'LineWidth', 2);
    set(gcf, 'color', 'w');
    ylabel('normalized Response')
    xlabel('Angle in deg')
    saveas(gcf, fullfile(saveDir, '31_SummedResponsesFirstVsLast.png'))
end

%Fig 32: Synaptic aggregat of random first 100 vs. random last 100 vs. cell
figure(32)
plot(angles, SomaRespNorm, 'Color', 'black', 'LineWidth', 5);
hold on
if size(oriSelectROIsNr,2) > compareFirstLastAmount*2+1
    plot(angles, mean(curveSumRandFirst,2), 'Color', 'green', 'LineWidth', 2);
    hold on
    plot(angles, mean(curveSumRandLast,2), 'Color', 'blue', 'LineWidth', 2);
end
ylabel('normalized Response')
xlabel('Angle in deg')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '32_SummedResponsesRandFirstVsLast.png'))

%Fig 33: Distribution of prefOri as histogram
[counts, ~] = histcounts([confData.ROIs(oriSelectROIsNr).prefOri], [0 linspace(22.5/2,180-22.5/2,8) 180]); %do counts that have the oris as center of the bins
countsAdjusted = [counts(1)+counts(end) counts(2:end-1) counts(1)+counts(end)];
countsNormalized = countsAdjusted/max(countsAdjusted);
figure(33)
plot(angles, SomaRespNorm, 'Color', 'black', 'LineWidth', 5);
hold on
plot(angles, countsNormalized(1:end-1), 'Color', 'red', 'LineWidth', 2);
set(gcf, 'color', 'w');
xlabel('Pref ori in deg')
ylabel('Probability')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '33_PrefOriAll.png'))

%Fig 34: Distribution of prefOri as polarPlot
figure(34)
polarPlotOri([angles 180], [SomaRespNorm SomaRespNorm(1)], [], 'black')
hold on
polarPlotOri([angles 180], countsNormalized, [], 'red')
saveas(gcf, fullfile(saveDir, '34_PrefOriAll_PolarPlot.png'))

%Fig 35: Synaptic aggregate vs. cell in direction space
figure(35)
plot(anglesAll, SomaDirNorm, 'Color', 'black', 'LineWidth', 5);
hold on
plot(anglesAll, curveSumFuncDir, 'Color', 'red', 'LineWidth', 2);
ylabel('normalized Response')
xlabel('Angle in deg (dir)')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '35_SummedRespDir.png'))

%Fig 36: Distribution of prefDir as histogram
[countsDir, ~] = histcounts([confData.ROIs(dirSelectROIsNr).prefDir], [0 linspace(22.5/2,360-22.5/2,16) 360]); %do counts that have the oris as center of the bins
countsAdjustedDir = [countsDir(1)+countsDir(end) countsDir(2:end-1) countsDir(1)+countsDir(end)];
countsNormalizedDir = countsAdjustedDir/max(countsAdjustedDir);
figure(36)
plot(anglesAll, SomaDirNorm, 'Color', 'black', 'LineWidth', 5);
hold on
plot(anglesAll, countsNormalizedDir(1:end-1), 'Color', 'red', 'LineWidth', 2);
set(gcf, 'color', 'w');
xlabel('Pref dir in deg')
ylabel('Probability')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, '36_PrefDirAll.png'))

%Fig 37: Distribution of prefDir as polarPlot
figure(37)
polarPlot([anglesAll 360], [SomaDirNorm SomaDirNorm(1)], [], 'black','-',0)
hold on
polarPlot([anglesAll 360], countsNormalizedDir, [], 'red')
saveas(gcf, fullfile(saveDir, '37_PrefDirAll_PolarPlot.png'))

%Fig 38: Distribution of bandwidth
figure(38)
scatter(ones(size([confData.ROIs(oriSelectROIsNr).Bandwidth],1)), [confData.ROIs(oriSelectROIsNr).Bandwidth], 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
hold on
boxplot([confData.ROIs(oriSelectROIsNr).Bandwidth],'color','k')
ylabel('Bandwidth in deg')
xticklabels({'All spines'})
ylim([0 90])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'38_Bandwidth.png'))

%--------------------------------------------------------------------------
%40/50s: Properties of the local environment (50s: Ori, 60s: Dir)
getDeltaOri = @(x) x.funcData.deltaOri; %define this to be able to get the deltaOriValues
deltaOriValues = cell2mat(arrayfun(getDeltaOri, confData.ROIs(oriSelectROIsNr), 'UniformOutput', false))';

%Fig 40/50: local (dir)Dispersion
localDisp_values = cell2mat(cellfun(@(x) x.localDispersion, num2cell(confData.ROIs(oriSelectROIsNr)), 'UniformOutput', false));
localDisp_values = reshape(localDisp_values, 4, [])';
localDirDisp_values = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(confData.ROIs(dirSelectROIsNr)), 'UniformOutput', false));
localDirDisp_values = reshape(localDirDisp_values, 4, [])';

figure(40)
subplot(1,2,1)
boxplot(localDisp_values,'color','k')
ylim([0 45])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDisp_values,1),1), localDisp_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('Local dispersion')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot(deltaOriValues,localDisp_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 45])
ylabel(sprintf('Local dispersion within (5 \\mum) of spine'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'40_Local dispersion to Spine.png'))

figure(50)
subplot(1,2,1)
boxplot(localDirDisp_values,'color','k')
ylim([0 90])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDirDisp_values,1),1), localDirDisp_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('Local dir dispersion')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(dirSelectROIsNr).prefDir],localDirDisp_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 45])
ylabel(sprintf('Local dir dispersion within (5 \\mum) of spine'))
xlabel('pref Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'50_local dir dispersion to Spine.png'))

%Fig 41/51: HI
HI_values = cell2mat(cellfun(@(x) x.HI, num2cell(confData.ROIs(oriSelectROIsNr)), 'UniformOutput', false));
HI_values = reshape(HI_values, 4, [])';
HIDir_values = cell2mat(cellfun(@(x) x.HIDir, num2cell(confData.ROIs(dirSelectROIsNr)), 'UniformOutput', false));
HIDir_values = reshape(HIDir_values, 4, [])';

figure(41)
subplot(1,2,1)
boxplot(HI_values,'color','k')
ylim([0 1])
xlim([0 5])
hold on
scatter(repmat(1:4,size(HI_values,1),1), HI_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('Homeogeneity index')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot(deltaOriValues,HI_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 1])
ylabel(sprintf('HI of local environment (5 \\mum) of spine'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'41_Local HI to Spine.png'))

figure(51)
subplot(1,2,1)
boxplot(HIDir_values,'color','k')
ylim([0 1])
xlim([0 5])
hold on
scatter(repmat(1:4,size(HIDir_values,1),1), HIDir_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('Homeogeneity index')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(dirSelectROIsNr).prefDir],HIDir_values(:,2),'*', 'color', 'g');
xlim([0 360])
ylim([0 1])
ylabel(sprintf('HI of local environment (5 \\mum) of spine'))
xlabel('prefDir')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'51_Local HIdir to Spine.png'))

%Fig. 42/52 Local dispersion vs. HI
figure(42)
plot(localDisp_values(:,2),HI_values(:,2),'*', 'color', 'g');
xlim([0 45])
ylim([0 1])
ylabel(sprintf('HI of local environment (5 \\mum) of spine'))
xlabel(sprintf('circular Dispersion of local environment (5 \\mum) of spine'))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'42_Local HI vs local Dispersion.png'))

figure(52)
plot(localDirDisp_values(:,2),HIDir_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 1])
ylabel(sprintf('HIdir of local environment (5 \\mum) of spine'))
xlabel(sprintf('circular dir Dispersion of local environment (5 \\mum) of spine'))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'52_Local HIdir vs local dir Dispersion.png'))

%Fig 43/53: local deltaOri
localDeltaOri_values = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(confData.ROIs(oriSelectROIsNr)), 'UniformOutput', false));
localDeltaOri_values = reshape(localDeltaOri_values, 4, [])';
localDeltaDir_values = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(confData.ROIs(dirSelectROIsNr)), 'UniformOutput', false));
localDeltaDir_values = reshape(localDeltaDir_values, 4, [])';

figure(43)
subplot(1,2,1)
boxplot(localDeltaOri_values,'color','k')
ylim([0 90])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDeltaOri_values,1),1), localDeltaOri_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('deltaOri of local environment to spine')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(oriSelectROIsNr).prefOri],localDeltaOri_values(:,2),'*', 'color', 'g');
xlim([0 180])
ylim([0 90])
ylabel(sprintf('deltaOri of local environment (5 \\mum) to spine'))
xlabel('pref Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'43_Local deltaOri to Spine.png'))

figure(53)
subplot(1,2,1)
boxplot(localDeltaDir_values,'color','k')
ylim([0 180])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDeltaDir_values,1),1), localDeltaDir_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('deltaDir of local environment to spine')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(dirSelectROIsNr).prefDir],localDeltaDir_values(:,2),'*', 'color', 'g');
xlim([0 360])
ylim([0 180])
ylabel(sprintf('deltaDir of local environment (5 \\mum) to spine'))
xlabel('delta Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'53_Local deltaDir to Spine.png'))


%Fig 44/54: lodal delteOriSoma
localDeltaOriSoma_values = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(confData.ROIs(oriSelectROIsNr)), 'UniformOutput', false));
localDeltaOriSoma_values = reshape(localDeltaOriSoma_values, 4, [])';
localDeltaDirSoma_values = cell2mat(cellfun(@(x) x.localDeltaDirSoma, num2cell(confData.ROIs(dirSelectROIsNr)), 'UniformOutput', false));
localDeltaDirSoma_values = reshape(localDeltaDirSoma_values, 4, [])';

figure(44)
subplot(1,2,1)
boxplot(localDeltaOriSoma_values,'color','k')
ylim([0 90])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDeltaOriSoma_values,1),1), localDeltaOriSoma_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15)
ylabel('deltaOri of local environment to Soma')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot(deltaOriValues,localDeltaOriSoma_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 90])
ylabel(sprintf('deltaOri of local environment (5 \\mum) to Soma'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'44_Local deltaOri to Soma vs deltaOri Spine.png'))

figure(54)
subplot(1,2,1)
boxplot(localDeltaDirSoma_values,'color','k')
ylim([0 180])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDeltaDirSoma_values,1),1), localDeltaDirSoma_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
ylabel('deltaDir of local environment to Soma')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(dirSelectROIsNr).prefDir],localDeltaDirSoma_values(:,2),'*', 'color', 'g');
xlim([0 180])
ylim([0 180])
ylabel(sprintf('deltaDir of local environment (5 \\mum) to Soma'))
xlabel('pref Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'54_Local deltaDir to Soma vs pref Dir Spine.png'))

%Fig 45/55: nearest SimilarPrefOri/dir vs its own prefOri
figure(45)
plot([confData.ROIs(oriSelectROIsNr).prefOri],[confData.ROIs(oriSelectROIsNr).nearestSimilarPrefOri],'*', 'color', 'g');
xlim([0 180])
ylim([0 25])
ylabel(sprintf('Nearest neighbor with similar preference in \\mum'))
xlabel('pref Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'45_NearestNeighbor prefOri.png'))

figure(55)
plot([confData.ROIs(dirSelectROIsNr).prefDir],[confData.ROIs(dirSelectROIsNr).nearestSimilarPrefDir],'*', 'color', 'g');
xlim([0 360])
ylim([0 25])
ylabel(sprintf('Nearest neighbor with similar preference in \\mum'))
xlabel('pref Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'55_NearestNeighbor prefDir.png'))

%Fig 46/56: local OSI/DSI
%get all localOSI and DSI values & reshape them
localOSI_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(confData.ROIs(respROINrs)), 'UniformOutput', false));
localOSI_values = reshape(localOSI_values, 4, [])';

localDSI_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(confData.ROIs(respROINrs)), 'UniformOutput', false));
localDSI_values = reshape(localDSI_values, 4, [])';

figure(46)
subplot(1,2,1)
boxplot(localOSI_values,'color','k')
ylim([0 1])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localOSI_values,1),1), localOSI_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
ylabel('local OSI')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(respROINrs).OSI],localOSI_values(:,2),'*', 'color', 'g');
xlim([0 1])
ylim([0 1])
ylabel(sprintf('local OSI (5 \\mum)'))
xlabel('OSI spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'46_LocalOSI.png'))

figure(56)
subplot(1,2,1)
boxplot(localDSI_values,'color','k')
ylim([0 1])
xlim([0 5])
hold on
scatter(repmat(1:4,size(localDSI_values,1),1), localDSI_values, 'filled','MarkerFaceColor','green', 'MarkerFaceAlpha',0.1','jitter','on','jitterAmount',0.15)
ylabel('local DSI')
xlabel(sprintf('Radius size in \\mum'))
xticklabels({'2.5', '5', '7.5', '10'})
box off
subplot(1,2,2)
plot([confData.ROIs(respROINrs).DSI],localDSI_values(:,2),'*', 'color', 'g');
xlim([0 1])
ylim([0 1])
ylabel(sprintf('local DSI (5 \\mum)'))
xlabel('DSI spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'56_LocalDSI.png'))

%get those only from the selective ROIs
localOSISel_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(confData.ROIs(oriSelectROIsNr)), 'UniformOutput', false));
localOSISel_values = reshape(localOSISel_values, 4, [])';

localDSISel_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(confData.ROIs(dirSelectROIsNr)), 'UniformOutput', false));
localDSISel_values = reshape(localDSISel_values, 4, [])';

%Fig 47/57: local OSI/DSI vs pref Ori/Dir
figure(47)
subplot(1,2,1)
plot([confData.ROIs(oriSelectROIsNr).prefOri],localOSISel_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 1])
ylabel(sprintf('local OSI (5 \\mum)'))
xlabel('pref Ori spine')
box off
subplot(1,2,2)
plot(deltaOriValues,localOSISel_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 1])
ylabel(sprintf('local OSI (5 \\mum)'))
xlabel('delta Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'47_LocalOSI vs. prefOri and deltaOri.png'))

figure(57)
plot([confData.ROIs(dirSelectROIsNr).prefDir],localDSISel_values(:,2),'*', 'color', 'g');
xlim([0 360])
ylim([0 1])
ylabel(sprintf('local DSI (5 \\mum)'))
xlabel('prefDir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'57_LocalDSI vs prefOri.png'))

%--------------------------------------------------------------------------
%60s: Properties of dendritic segments
%What is the overal property of the dendritic segment? How homeogeneous is
%it, how different from soma, how selective, ...

%Fig 60: Branch dispersion
figure(60)
subplot(1,2,1)
boxplot([confData.allDendrites.circDispersion]','color','b')
ylim([0 90])
hold on
ylabel('Dendritic segment circular dispersion')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.circDispersion],1),1),[confData.allDendrites.circDispersion],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,2,2)
boxplot([confData.allDendrites.dirCircDispersion]','color','b')
ylim([0 180])
hold on
ylabel('Dendritic segment dir circular dispersion')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.dirCircDispersion],1),1),[confData.allDendrites.dirCircDispersion],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'60_circularDispersion.png'))

%Fig 61: Branch selectivity
figure(61)
subplot(1,3,1)
boxplot([confData.allDendrites.medianOSI]','color','b')
ylim([0 1])
hold on
ylabel('median OSI')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.medianOSI],1),1),[confData.allDendrites.medianOSI],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,3,2)
boxplot([confData.allDendrites.medianDSI]','color','b')
ylim([0 1])
hold on
ylabel('median DSI')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.medianDSI],1),1),[confData.allDendrites.medianDSI],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,3,3)
boxplot([confData.allDendrites.medianDSIvect]','color','b')
ylim([0 1])
hold on
ylabel('median DSIvect')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.medianDSIvect],1),1),[confData.allDendrites.medianDSIvect],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'61_Branch selectivity.png'))

%Fig 62: delta Ori/delta Dir
figure(62)
subplot(1,2,1)
boxplot([confData.allDendrites.deltaOri]','color','b')
ylim([0 90])
hold on
ylabel('delta Ori')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.deltaOri],1),1),[confData.allDendrites.deltaOri],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,2,2)
boxplot([confData.allDendrites.deltaDir]','color','b')
ylim([0 180])
hold on
ylabel('delta Dir')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([confData.allDendrites.deltaDir],1),1),[confData.allDendrites.deltaDir],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'62_DeltaOri_deltaDir.png'))

%Fig 63: Branch circular dispersion vs. deltaOri
figure(63) 
plot([confData.allDendrites.deltaOri],[confData.allDendrites.circDispersion], '*','color', 'blue')
xlim([0 90])
xlabel('Delta Ori of dendritic segment to soma')
ylim([0 90])
ylabel('Dendritic segment circular dispersion')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'63_DeltaOri_circDispersion.png'))

%Fig 64: Branch circular dispersion dir vs. deltaDir
figure(64) 
plot([confData.allDendrites.deltaDir],[confData.allDendrites.dirCircDispersion], '*','color', 'blue')
xlim([0 180])
xlabel('Delta dir of dendritic segment to soma')
ylim([0 180])
ylabel('Dendritic segment dir circular dispersion')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'64_DeltaDir_dirCircDispersion.png'))

%70s: Pairwise distances
%Fig 70: Pref Ori
oriPairs = find([pwMeasures.oriSelect] == 1);

figure(70)
plotBinnedPairWisePropertyVsDistance([pwMeasures(oriPairs).distance],[pwMeasures(oriPairs).deltaOri],[0.5, 0.5, 0.5]);
xlabel('distance')
ylabel('deltaOri')
saveas(gcf, fullfile(saveDir,'70_DeltaOri vs. Distance.png'))

%Fig 71: Pref Dir
dirPairs = find([pwMeasures.dirSelect] == 1);

figure(71)
plotBinnedPairWisePropertyVsDistance([pwMeasures(dirPairs).distance],[pwMeasures(dirPairs).deltaDir],[0.5, 0.5, 0.5]);
ylabel('deltaDir')
saveas(gcf, fullfile(saveDir,'71_DeltaDir vs. Distance.png'))

%Fig 72: OSI
goodPairs = find([pwMeasures.goodPair] == 1);

figure(72)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodPairs).distance],[pwMeasures(goodPairs).deltaOSI],[0.5, 0.5, 0.5]);
ylabel('deltaOSI')
ylim([0 0.25])
saveas(gcf, fullfile(saveDir,'72_DeltaOSI vs. Distance.png'))

%Fig 73: DSI
figure(73)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodPairs).distance],[pwMeasures(goodPairs).deltaDSI],[0.5, 0.5, 0.5]);
ylabel('deltaDSI')
saveas(gcf, fullfile(saveDir,'73_DeltaDSI vs. Distance.png'))

%Fig 74: DSIVect
figure(74)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodPairs).distance],[pwMeasures(goodPairs).deltaDSIvect],[0.5, 0.5, 0.5]);
ylabel('deltaDSIvect')
saveas(gcf, fullfile(saveDir,'74_DeltaDSIvect vs. Distance.png'))

%Fig 75: TuningCurve correlation
figure(75)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodPairs).distance],[pwMeasures(goodPairs).curveCorr],[0.5, 0.5, 0.5]);
ylabel('Tuning curve correlation')
saveas(gcf, fullfile(saveDir,'75_Tuning curve correlation vs. Distance.png'))

%Fig 76: TuningCurve correlation
figure(76)
plotBinnedPairWisePropertyVsDistance([pwMeasures(goodPairs).distance],[pwMeasures(goodPairs).trialCorr],[0.5, 0.5, 0.5]);
ylabel('Trial-to-trial correlation')
saveas(gcf, fullfile(saveDir,'76_Trial-to-trial correlation vs. Distance.png'))

%% Step 7: Save

Selectors = struct;
%logical 0/1
Selectors.TwoPROIs = TwoPROIs;
Selectors.respROIs = respROIs;
Selectors.oriSelectROIs = oriSelectROIs;
Selectors.dirSelectROIs = dirSelectROIs;

%ROI Nrs
Selectors.TwoPROIsNR = TwoPROIsNR;
Selectors.respROINrs = respROINrs;
Selectors.oriSelectROIsNr = oriSelectROIsNr;
Selectors.dirSelectROIsNr = dirSelectROIsNr;

save([saveDirAna filesep '02_FunctionalCellReconstruction.mat'], 'confData','Grouping','Selectors', 'pwMeasures','-mat') 
