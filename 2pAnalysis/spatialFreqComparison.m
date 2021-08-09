%% switch bord for variables
plotROIs = 0;
reloadData = 0;
server = 0;
SFSelective = 0.33;

%select colors for plots
cocCellsV1 = cbrewer('seq', 'Blues', 5);
cocCellsV3 = cbrewer('seq', 'Greens', 5);
cocAxonsV3 = cbrewer('seq', 'PuRd', 5);
fsz = 12;
%shufflenum = 100;

%% set folders and experiment files

if server == 1
    baseDir = 'Z:\Juliane\Data\';
    fileDir = 'Z:\Juliane\Organization\Animals\';
else
    baseDir = 'F:\Data\';
    fileDir = 'F:\Organization\Animals\';
end

axonFile = [fileDir '2pExpByStimulusAxon.xlsx'];
cellFile = [fileDir '2pExpByStimulus.xlsx'];
analysisDir = [baseDir 'ImageAnalysis\'];
saveDir = [analysisDir filesep 'SpatialFreqGrating' filesep];
if ~exist(saveDir)
    mkdir(saveDir)
end

%% read which experiments to use

%load xls files with experiment infos
[~, axonFileTxt, axonFileAll] = xlsread(axonFile, 'driftingGrating_ori_sf');
[~, cellFileTxt, cellFileAll] = xlsread(cellFile, 'driftingGrating_ori_sf');

%define the relevant colums
flagCol     = find(contains(cellFileTxt(1,:), 'Flag'), 1);
animalCol   = find(contains(cellFileTxt(1,:), 'animal'),1);
expCol      = find(contains(cellFileTxt(1,:), 'expNumber'), 1);
regionCol   = find(contains(cellFileTxt(1,:), 'region'), 1);

expInfoCells = struct;

%save information from cell file into the structure
k = 1;
for rowCell = 2:size(cellFileAll)
    flag = cellFileAll(rowCell, flagCol);
    if flag{1} == 0
        expInfoCells.animal{k} = cellFileAll(rowCell, animalCol);
        expInfoCells.region{k} = cellFileAll(rowCell, regionCol);
        exp = cellFileAll(rowCell, expCol);
        if exp{1} > 9
            if exp{1} > 99
                expInfoCells.expID{k} = ['t00' num2str(exp{1})];
            else
                expInfoCells.expID{k} = ['t000' num2str(exp{1})];
            end
        else
            expInfoCells.expID{k} = ['t0000' num2str(exp{1})];
        end
        k = k + 1;
    end
end

%save information from axon file into the structure (assuming same columns
%as in the cell file
for rowAxon = 2:size(axonFileAll)
    flag = axonFileAll(rowAxon, flagCol);
    if flag{1} == 0
        expInfoCells.animal{k} = axonFileAll(rowAxon, animalCol);
        expInfoCells.region{k} = axonFileAll(rowAxon, regionCol);
        exp = axonFileAll(rowAxon, expCol);
        if exp{1} > 9
            if exp{1} > 99
                expInfoCells.expID{k} = ['t00' num2str(exp{1})];
            else
                expInfoCells.expID{k} = ['t000' num2str(exp{1})];
            end
        else
            expInfoCells.expID{k} = ['t0000' num2str(exp{1})];
        end
        k = k + 1;
    end
end

%% load all data into several master files
cellsV1 = 1;
cellsV3 = 1;
axonsV3 = 1;

for il = 1:length(expInfoCells.animal)
    datapath = [analysisDir char(expInfoCells.animal{il}) filesep char(expInfoCells.expID{il}) filesep];
    region = expInfoCells.region{il};
    switch char(region{1})
        case 'V1/V2'
            masterCellsV1{cellsV1} = load(fullfile(datapath, 's1_ori_sf_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
            masterCellsV1{cellsV1}.metadata.animal = expInfoCells.animal{il};
            masterCellsV1{cellsV1}.metadata.expID = expInfoCells.expID{il};
            masterCellsV1{cellsV1}.responsive = find([masterCellsV1{cellsV1}.analysis.dff.roi.isResponseSignificant] == 1);
            masterCellsV1{cellsV1}.SFselective = find([masterCellsV1{cellsV1}.analysis.dff.roi.SFSI] > SFSelective); 
            cellsV1 = cellsV1 + 1;
        case 'V3'
            masterCellsV3{cellsV3} = load(fullfile(datapath, 's1_ori_sf_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
            masterCellsV3{cellsV3}.metadata.animal = expInfoCells.animal{il};
            masterCellsV3{cellsV3}.metadata.expID = expInfoCells.expID{il};
            masterCellsV3{cellsV3}.responsive = find([masterCellsV3{cellsV3}.analysis.dff.roi.isResponseSignificant] == 1);
            masterCellsV3{cellsV3}.SFselective = find([masterCellsV3{cellsV3}.analysis.dff.roi.SFSI] > SFSelective); 
            cellsV3 = cellsV3 + 1;
        case 'feedback'
            masterAxonsV3{axonsV3} = load(fullfile(datapath, 's1_ori_sf_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
            masterAxonsV3{axonsV3}.metadata.animal = expInfoCells.animal{il};
            masterAxonsV3{axonsV3}.metadata.expID = expInfoCells.expID{il};
            masterAxonsV3{axonsV3}.responsive = find([masterAxonsV3{axonsV3}.analysis.dff.roi.isResponseSignificant] == 1);
            masterAxonsV3{axonsV3}.SFselective = find([masterAxonsV3{axonsV3}.analysis.dff.roi.SFSI] > SFSelective); 
            axonsV3 = axonsV3 + 1;
    end
end

%% data consolidation
%prelocate variables
OSI_cellsV1     = [];
DSI_cellsV1     = [];
sFPref_cellsV1  = [];

OSI_cellsV3     = [];
DSI_cellsV3     = [];
sFPref_cellsV3  = [];

OSI_axonsV3     = [];
DSI_axonsV3     = [];
sFPref_axonsV3  = [];

%go through all files for cells in V1, V3 and feedback axons
for experiment = 1:length(masterCellsV1)
    responsive =  masterCellsV1{experiment}.responsive;
    SFselective = masterCellsV1{experiment}.SFselective;
    respSFselect = intersect(responsive, SFselective);
    OSI_cellsV1 = [OSI_cellsV1 masterCellsV1{experiment}.analysis.dff.roi(responsive).OSIFit];
    DSI_cellsV1 = [DSI_cellsV1 masterCellsV1{experiment}.analysis.dff.roi(responsive).DSI];
    for resp = 1:length(respSFselect)
        sFPref_cellsV1 = [sFPref_cellsV1 str2double(masterCellsV1{experiment}.analysis.dff.roi(respSFselect(resp)).prefSf)];
    end
end

for experiment = 1:length(masterCellsV3)
    responsive =  masterCellsV3{experiment}.responsive;
    SFselective = masterCellsV3{experiment}.SFselective;
    respSFselect = intersect(responsive, SFselective);
    OSI_cellsV3 = [OSI_cellsV3 masterCellsV3{experiment}.analysis.dff.roi(responsive).OSIFit];
    DSI_cellsV3 = [DSI_cellsV3 masterCellsV3{experiment}.analysis.dff.roi(responsive).DSI];
    for resp = 1:length(respSFselect)
        sFPref_cellsV3 = [sFPref_cellsV3 str2double(masterCellsV3{experiment}.analysis.dff.roi(respSFselect(resp)).prefSf)];
    end
end

for experiment = 1:length(masterAxonsV3)
    responsive =  masterAxonsV3{experiment}.responsive;
    SFselective = masterAxonsV3{experiment}.SFselective;
    respSFselect = intersect(responsive, SFselective);
    OSI_axonsV3 = [OSI_axonsV3 masterAxonsV3{experiment}.analysis.dff.roi(responsive).OSIFit];
    DSI_axonsV3 = [DSI_axonsV3 masterAxonsV3{experiment}.analysis.dff.roi(responsive).DSI];
    for resp = 1:length(respSFselect)
        sFPref_axonsV3 = [sFPref_axonsV3 str2double(masterAxonsV3{experiment}.analysis.dff.roi(respSFselect(resp)).prefSf)];
    end
end

%% plot OSI and DSI distribution
binNum = 10;

figure(1)
subplot(2,3,1)
plot(nanmedian(OSI_cellsV1), 1.05*max(histcounts(OSI_cellsV1, binNum)), 'v', 'MarkerSize', 8', 'MarkerEdgeColor', cocCellsV1(4,:), 'MarkerFaceColor', cocCellsV1(3,:)); hold on
text(nanmedian(OSI_cellsV1), 1.1*max(histcounts(OSI_cellsV1, binNum)), num2str(round(100*nanmedian(OSI_cellsV1))/100), 'HorizontalAlignment', 'center', 'Color', cocCellsV1(4,:), 'FontSize', fsz)
histogram(OSI_cellsV1, binNum, 'FaceColor', cocCellsV1(3,:), 'EdgeColor', cocCellsV1(4,:));
ylabel('Cells')
xlabel('OSI')
title('V1')
set(gca, 'Box', 'off')

subplot(2,3,2)
plot(nanmedian(OSI_cellsV3), 1.05*max(histcounts(OSI_cellsV3, binNum)), 'v', 'MarkerSize', 8', 'MarkerEdgeColor', cocCellsV3(4,:), 'MarkerFaceColor', cocCellsV3(3,:)); hold on
text(nanmedian(OSI_cellsV3), 1.1*max(histcounts(OSI_cellsV3, binNum)), num2str(round(100*nanmedian(OSI_cellsV3))/100), 'HorizontalAlignment', 'center', 'Color', cocCellsV3(4,:), 'FontSize', fsz)
histogram(OSI_cellsV3, binNum, 'FaceColor', cocCellsV3(3,:), 'EdgeColor', cocCellsV3(4,:));
ylabel('Cells')
xlabel('OSI')
xlim([0 1])
title('V3')
set(gca, 'Box', 'off')

subplot(2,3,3)
plot(nanmedian(OSI_axonsV3), 1.05*max(histcounts(OSI_axonsV3,binNum)), 'v', 'MarkerSize', 8', 'MarkerEdgeColor', cocAxonsV3(4,:), 'MarkerFaceColor', cocAxonsV3(3,:)); hold on
text(nanmedian(OSI_axonsV3), 1.1*max(histcounts(OSI_axonsV3,binNum)), num2str(round(100*nanmedian(OSI_axonsV3))/100), 'HorizontalAlignment', 'center', 'Color', cocAxonsV3(4,:), 'FontSize', fsz)
histogram(OSI_axonsV3, binNum, 'FaceColor', cocAxonsV3(3,:), 'EdgeColor', cocAxonsV3(4,:));
ylabel('Boutons')
xlabel('OSI')
xlim([0 1])
title('Axons from V3')
set(gca, 'Box', 'off')

subplot(2,3,4)
plot(nanmedian(DSI_cellsV1), 1.05*max(histcounts(DSI_cellsV1,binNum)), 'v', 'MarkerSize', 8', 'MarkerEdgeColor', cocCellsV1(4,:), 'MarkerFaceColor', cocCellsV1(3,:)); hold on
text(nanmedian(DSI_cellsV1), 1.1*max(histcounts(DSI_cellsV1,binNum)), num2str(round(100*nanmedian(DSI_cellsV1))/100), 'HorizontalAlignment', 'center', 'Color', cocCellsV1(4,:), 'FontSize', fsz)
histogram(DSI_cellsV1, binNum, 'FaceColor', cocCellsV1(3,:), 'EdgeColor', cocCellsV1(4,:));
ylabel('Cells')
xlabel('DSI')
xlim([0 1])
title('V1')
set(gca, 'Box', 'off')

subplot(2,3,5)
plot(nanmedian(DSI_cellsV3), 1.05*max(histcounts(DSI_cellsV3,binNum)), 'v', 'MarkerSize', 8', 'MarkerEdgeColor', cocCellsV3(4,:), 'MarkerFaceColor', cocCellsV3(3,:)); hold on
text(nanmedian(DSI_cellsV3), 1.1*max(histcounts(DSI_cellsV3,binNum)), num2str(round(100*nanmedian(DSI_cellsV3))/100), 'HorizontalAlignment', 'center', 'Color', cocCellsV3(4,:), 'FontSize', fsz)
histogram(DSI_cellsV3, binNum, 'FaceColor', cocCellsV3(3,:), 'EdgeColor', cocCellsV3(4,:));
ylabel('Cells')
xlabel('DSI')
xlim([0 1])
title('V3')
set(gca, 'Box', 'off')

subplot(2,3,6)
plot(nanmedian(DSI_axonsV3), 1.05*max(histcounts(DSI_axonsV3,binNum)), 'v', 'MarkerSize', 8', 'MarkerEdgeColor', cocAxonsV3(4,:), 'MarkerFaceColor', cocAxonsV3(3,:)); hold on
text(nanmedian(DSI_axonsV3), 1.1*max(histcounts(DSI_axonsV3,binNum)), num2str(round(100*nanmedian(DSI_axonsV3))/100), 'HorizontalAlignment', 'center', 'Color', cocAxonsV3(4,:), 'FontSize', fsz)
histogram(DSI_axonsV3, binNum, 'FaceColor', cocAxonsV3(3,:), 'EdgeColor', cocAxonsV3(4,:));
ylabel('Boutons')
xlabel('DSI')
xlim([0 1])
title('Axons from V3')
set(gca, 'Box', 'off')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, 'OSI_DSI_distribution.png'))

%% plot SF preference distribution
figure
subplot(1,3,1)
histogram(sFPref_cellsV1, 'FaceColor', cocCellsV1(3,:), 'EdgeColor', cocCellsV1(4,:));
ylabel('Cells')
xlabel('Optimal SF')
title('V1')
set(gca, 'Box', 'off')

subplot(1,3,2)
histogram(sFPref_cellsV3, 'FaceColor', cocCellsV3(3,:), 'EdgeColor', cocCellsV3(4,:));
ylabel('Cells')
xlabel('Optimal SF')
title('V3')
set(gca, 'Box', 'off')

subplot(1,3,3)
histogram(sFPref_axonsV3, 'FaceColor', cocAxonsV3(3,:), 'EdgeColor', cocAxonsV3(4,:));
ylabel('Boutons')
xlabel('Optimal SF')
title('Axons V3')
set(gca, 'Box', 'off')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, 'SF_distribution.png'))

%% divide SF pref into three categories: <=0.03, 0.02-0.6 and >0.6, plot
sFPrefCat_cellsV1 = ones(length(sFPref_cellsV1),1);
sFPrefCat_cellsV3 = ones(length(sFPref_cellsV3),1);
sFPrefCat_axonsV3 = ones(length(sFPref_axonsV3),1);

sFPrefCat_cellsV1(find([sFPref_cellsV1] <= 0.03)) = 0;
sFPrefCat_cellsV1(find([sFPref_cellsV1] > 0.06)) = 2;
sFPrefCat_cellsV1 = sFPrefCat_cellsV1+1;

sFPrefCat_cellsV3(find([sFPref_cellsV3] <= 0.03)) = 0;
sFPrefCat_cellsV3(find([sFPref_cellsV3] > 0.06)) = 2;
sFPrefCat_cellsV3 = sFPrefCat_cellsV3+1;

sFPrefCat_axonsV3(find([sFPref_axonsV3] <= 0.03)) = 0;
sFPrefCat_axonsV3(find([sFPref_axonsV3] > 0.06)) = 2;
sFPrefCat_axonsV3 = sFPrefCat_axonsV3+1;

%plot as a pie chart
binEdges = [0 1.1 2.1 3.1];

figure
subplot(1,3,1)
h = pie(histcounts(sFPrefCat_cellsV1,binEdges));
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocCellsV1(3,:));
set(hp(2), 'FaceColor', cocCellsV1(4,:));
set(hp(3), 'FaceColor', cocCellsV1(5,:));
title('cells V1')
legend({'<=0.03 cpd', '0.04-0.06 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,2)
h = pie(histcounts(sFPrefCat_cellsV3,binEdges));
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocCellsV3(3,:));
set(hp(2), 'FaceColor', cocCellsV3(4,:));
set(hp(3), 'FaceColor', cocCellsV3(5,:));
title('cells V3')
legend({'<=0.03 cpd', '0.04-0.06 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,3)
h = pie(histcounts(sFPrefCat_axonsV3,binEdges));
hp = findobj(h, 'Type', 'patch');
if size(hp,1) == 2
    emptyPos = find(histcounts(sFPrefCat_axonsV3,binEdges) == 0);
    if emptyPos == 1
       set(hp(1), 'FaceColor', cocAxonsV3(4,:));
       set(hp(2), 'FaceColor', cocAxonsV3(5,:));
       legend({'0.04-0.06 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
    elseif emptyPos == 2
       set(hp(1), 'FaceColor', cocAxonsV3(3,:));
       set(hp(3), 'FaceColor', cocAxonsV3(5,:))
       legend({'<=0.03 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
    else
       set(hp(1), 'FaceColor', cocAxonsV3(3,:));
       set(hp(2), 'FaceColor', cocAxonsV3(4,:));
       legend({'<=0.03 cpd', '0.04-0.06 cpd'}, 'Location', 'southoutside')
    end
else
    set(hp(1), 'FaceColor', cocAxonsV3(3,:));
    set(hp(2), 'FaceColor', cocAxonsV3(4,:));
    set(hp(3), 'FaceColor', cocAxonsV3(5,:));
    legend({'<=0.03 cpd', '0.04-0.06 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
end
title('axons from V3')
legend('boxoff')

set(gcf, 'color', 'w');
mtit('Spatial Frequency Preference')
saveas(gcf, fullfile(saveDir, 'spatialFreqPref.png'))

%% do the same, but only for cells
%plot as a pie chart
figure
subplot(1,2,1)
h = pie(histcounts(sFPrefCat_cellsV1,binEdges));
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocCellsV1(3,:));
set(hp(2), 'FaceColor', cocCellsV1(4,:));
set(hp(3), 'FaceColor', cocCellsV1(5,:));
title('cells V1')
legend({'<=0.03 cpd', '0.04-0.06 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,2,2)
h = pie(histcounts(sFPrefCat_cellsV3,binEdges));
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocCellsV3(3,:));
set(hp(2), 'FaceColor', cocCellsV3(4,:));
set(hp(3), 'FaceColor', cocCellsV3(5,:));
title('cells V3')
legend({'<=0.03 cpd', '0.04-0.06 cpd', '>0.06 cpd'}, 'Location', 'southoutside')
legend('boxoff')

set(gcf, 'color', 'w');
mtit('Spatial Frequency Preference')
saveas(gcf, fullfile(saveDir, 'spatialFreqPrefCells.png'))