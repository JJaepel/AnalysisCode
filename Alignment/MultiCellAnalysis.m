function MultiCellAnalysis

%Combines all functional and anatomical data from all cell that were
%functionally characterized and matched to the confocal data

%0.) Read in xls file to see which cells can be added
%1.) Read in the files and combine the structures
%2.) Specify functional groups 
%3.) Plot

%!!!NEXT STEP: THINK ABOUT HOW TO BEST POOL PWMEASURES

%% 0.) Read in xls file to see which cells can be added & define folders

baseDir = 'Z:\Juliane\InputAnalysis\';
saveDir = [baseDir filesep 'Pooled Data\'];
saveDirImages = [saveDir filesep 'All Spines'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end
if ~exist(saveDirImages, 'dir')
    mkdir(saveDirImages)
end

cellInfo = animalParser;
cellsInclude = find(cell2mat(cellfun(@(x) contains(x, 'Yes'),{cellInfo.funcFinished},'UniformOutput',false)) ==1);

%make sure that no cells are in there twice
animalNames = cellInfo.animal(cellsInclude);
indAnimals = unique(animalNames);

%% 1.) Read in the files and combine the structures
cellNr = 1;
allSpines = [];
allBranches = [];
pwMeasures = [];

%also save everything to be accessed by the next one to make sure that the
%right cells are connected
cellDetails = struct;

for f = 1:length(indAnimals)
    animal = indAnimals{f};
    entries = find(cell2mat((cellfun(@(x) contains(x, animal),{animalNames},'UniformOutput',false))) == 1);
    
    %are there multiple entries for that animal?
    if length(entries) > 1
        %is it multiple cells or just one cell across multiple slices?
        cellNames = cellInfo.cellName(cellsInclude(entries));
        cellsWithNames = sum(~cellfun(@isempty, cellNames));
        cellsInput = cellInfo.InputFinished(cellsInclude(entries));
        
        if cellsWithNames > 1
           for c = 1:length(cellNames)
              disp(['Loading cell Nr: ' num2str(cellNr) ', animal ' char(animal) ', ' cellNames{c}]) 
               
              tempDir = [baseDir char(animal) filesep cellNames{c} filesep 'E - Analysis' filesep];
              fileData = load([tempDir '02_FunctionalCellReconstruction.mat'],'confData', 'pwMeasures'); 
              fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'CellNr', cellNr),  fileData.confData.ROIs);
              fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'Animal', animal),  fileData.confData.ROIs);
              fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'CellName',cellNames{c}),  fileData.confData.ROIs);
              allSpines = [allSpines fileData.confData.ROIs];
              fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'CellNr', cellNr),  fileData.confData.allDendrites);
              fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'Animal', animal),  fileData.confData.allDendrites);
              fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'CellName',cellNames{c}),  fileData.confData.allDendrites);
              allBranches = [allBranches fileData.confData.allDendrites];
              
              %add info to cellDetails structure
              cellDetails(cellNr).cellNr = cellNr;
              cellDetails(cellNr).animal = animal;
              cellDetails(cellNr).cellName = cellNames{c};
              cellDetails(cellNr).input = cellsInput{c};
              cellDetails(cellNr).inputType = cell2mat(cellInfo.inputType(cellsInclude(entries(1))));
              
              cellNr = cellNr+1;              
           end
        else
            %load the data
            disp(['Loading cell Nr: ' num2str(cellNr) ', animal ' char(animal)]) 
            
            tempDir = [baseDir char(animal) filesep cell2mat(cellInfo.cellName(cellsInclude(entries(1)))) filesep 'E - Analysis' filesep];
            fileData = load([tempDir '02_FunctionalCellReconstruction.mat'],'confData', 'pwMeasures');
            fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'CellNr', cellNr),  fileData.confData.ROIs);
            fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'Animal', animal),  fileData.confData.ROIs);
            fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'CellName',cell2mat(cellInfo.cellName(cellsInclude(entries(1))))),  fileData.confData.ROIs);
            allSpines = [allSpines fileData.confData.ROIs];
            fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'CellNr', cellNr),  fileData.confData.allDendrites);
            fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'Animal', animal),  fileData.confData.allDendrites);
            fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'CellName',cell2mat(cellInfo.cellName(cellsInclude(entries(1))))),  fileData.confData.allDendrites);
            allBranches = [allBranches fileData.confData.allDendrites];
            
            %add info to cellDetails structure
            cellDetails(cellNr).cellNr = cellNr;
            cellDetails(cellNr).animal = animal;
            cellDetails(cellNr).cellName = cell2mat(cellInfo.cellName(cellsInclude(entries(1))));
            cellDetails(cellNr).input = cell2mat(cellInfo.InputFinished(cellsInclude(entries(1))));
            cellDetails(cellNr).inputType = cell2mat(cellInfo.inputType(cellsInclude(entries(1))));
            
            cellNr = cellNr+1;
        end
        
    else
        %load the data
        disp(['Loading cell Nr: ' num2str(cellNr) ', animal ' char(animal)]) 
        tempDir = [baseDir char(animal) filesep cell2mat(cellInfo.cellName(cellsInclude(entries))) filesep 'E - Analysis' filesep];
        fileData = load([tempDir '02_FunctionalCellReconstruction.mat'],'confData', 'pwMeasures');
        fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'CellNr', cellNr),  fileData.confData.ROIs);
        fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'Animal', animal),  fileData.confData.ROIs);
        fileData.confData.ROIs = arrayfun(@(x) setfield(x, 'CellName',[]),  fileData.confData.ROIs);
        allSpines = [allSpines fileData.confData.ROIs];
        fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'CellNr', cellNr),  fileData.confData.allDendrites);
        fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'Animal', animal),  fileData.confData.allDendrites);
        fileData.confData.allDendrites = arrayfun(@(x) setfield(x, 'CellName',[]),  fileData.confData.allDendrites);
        allBranches = [allBranches fileData.confData.allDendrites];
        
        %add info to cellDetails structure
        cellDetails(cellNr).cellNr = cellNr;
        cellDetails(cellNr).animal = animal;
        cellDetails(cellNr).cellName = [];
        cellDetails(cellNr).input = cellInfo.InputFinished(cellsInclude(entries));
        cellDetails(cellNr).inputType = cell2mat(cellInfo.inputType(cellsInclude(entries)));
        
        cellNr = cellNr+1;
    end
end

%% 2.) Specify functional groups 

TwoPROIs = ([allSpines.TwoPMatch] == 1); %all the ROIs with a functional match, size = allSpines
TwoPROIsNR = find([allSpines.TwoPMatch] ==1);
respROIs = ([allSpines.good] == 1); %all ROIs with a good one,size = allSpines
respROINrs = find([allSpines.good] == 1);

oriSelect = find([allSpines.OSI] > 0.1); %all oriselectROIs, size = allfuncSpines
oriSelectROIsNr = TwoPROIsNR(oriSelect);%all oriselectROIs, size = allSpines
oriSelectROIsNr = intersect(oriSelectROIsNr,respROINrs); %make sure they are also good ones
dirSelect = find([allSpines.DSI] > 0.1); %all dirSelectROIs, size = allFuncSpines
dirSelectROIsNr = TwoPROIsNR(dirSelect);%all oriselectROIs, size = allSpines
dirSelectROIsNr = intersect(dirSelectROIsNr,respROINrs); %make sure they are also good ones

oriSelectROIs = zeros(1,length(allSpines));
oriSelectROIs(oriSelectROIsNr) = 1; %logical vector
dirSelectROIs = zeros(1,length(allSpines));
dirSelectROIs(dirSelectROIsNr) = 1; %logical vector

%% 3.) Plot
%Figure Overview
%10s - Synaptic aggregate & histogram of deltaOri and deltaDir
%Fig 11: Synaptic aggregate per cell?
%Fig 12: Distribution of deltaOri as histogram
%Fig 13: Distribution of deltaOri as polarPlot
%Fig 14: Distribution of deltaDir as histogram
%Fig 15: Distribution of deltaDir as polarPlot
%Fig 16: Distribution of bandwidth

%20s & 30s: Local environment ori(20s) and dir (30s)
%Fig 20/30: local (dir) Dispersion
%Fig 21/31: HI (dir)
%Fig 22/32: local (dir) dispersion vs. HI (dir)
%Fig 23/33: local deltaOri/Dir
%Fig 24/34: local delteOriSoma/deltaDirSOma
%Fig 25/35: nearest SimilarPrefOri/dir
%Fig 26/36: local OSI/DSI

%40s: Properties of dendritic segments
%Fig 40: Branch circular dispersion ori/dir
%Fig 41: Branch selectivity
%Fig 42: delta Ori/delta Dir
%Fig 43: Branch circular dispersion vs. deltaOri
%Fig 44: Branch circular dispersion dir vs. deltaDir

%50s: Pairwise distances
%Fig 50: Pref Ori
%Fig 51: Pref Dir
%Fig 52: OSI
%Fig 53: DSI
%Fig 54: DSIVect

%60s - Cell wide analysis structure
%Fig 60: Distribution of distance from soma for all and separated for
%apical/basal
%Fig 61: Same as Fig 60, but as cdfplot
%Fig 62: Histogram of branch order for all and separated for apical/basal
%Fig 63: same as Fig 62, but as cdfplot

%70s - Branch specific analysis
%Fig 70: Number of spines per dendrite for all and separated for
%apical/basal
%Fig 71: Spine density per dendrite for all and separated for apical/basal
%Fig 72: Circular dispersion per dendrite for all and separated for apical/basal
%Fig 73: Dir circular dispersion per dendrite for all and separated for apical/basal
%Fig 74: DeltaOri per dendrite for all and separated for apical/basal
%Fig 75: DeltaDir per dendrite for all and separated for apical/basal
%Fig 76: MedianOSI per dendrite for all and separated for apical/basal
%Fig 77: MedianDSI per dendrite for all and separated for apical/basal
%Fig 78: MedianDSIvect per dendrite for all and separated for apical/basal

%--------------------------------------------------------------------------
%20/30s: Properties of the local environment (50s: Ori, 60s: Dir)
getDeltaOri = @(x) x.funcData.deltaOri; %define this to be able to get the deltaOriValues
deltaOriValues = cell2mat(arrayfun(getDeltaOri, allSpines(oriSelectROIsNr), 'UniformOutput', false))';

%Fig 20/30: local (dir)Dispersion
localDisp_values = cell2mat(cellfun(@(x) x.localDispersion, num2cell(allSpines(oriSelectROIsNr)), 'UniformOutput', false));
localDisp_values = reshape(localDisp_values, 4, [])';
localDirDisp_values = cell2mat(cellfun(@(x) x.localDirDispersion, num2cell(allSpines(dirSelectROIsNr)), 'UniformOutput', false));
localDirDisp_values = reshape(localDirDisp_values, 4, [])';

figure(20)
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
saveas(gcf, fullfile(saveDirImages,'20_Local dispersion to Spine.png'))

figure(30)
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
plot([allSpines(dirSelectROIsNr).prefDir],localDirDisp_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 45])
ylabel(sprintf('Local dir dispersion within (5 \\mum) of spine'))
xlabel('pref Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'30_local dir dispersion to Spine.png'))

%Fig 21/31: HI
HI_values = cell2mat(cellfun(@(x) x.HI, num2cell(allSpines(oriSelectROIsNr)), 'UniformOutput', false));
HI_values = reshape(HI_values, 4, [])';
HIDir_values = cell2mat(cellfun(@(x) x.HIDir, num2cell(allSpines(dirSelectROIsNr)), 'UniformOutput', false));
HIDir_values = reshape(HIDir_values, 4, [])';

figure(21)
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
saveas(gcf, fullfile(saveDirImages,'21_Local HI to Spine.png'))

figure(31)
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
plot([allSpines(dirSelectROIsNr).prefDir],HIDir_values(:,2),'*', 'color', 'g');
xlim([0 360])
ylim([0 1])
ylabel(sprintf('HI of local environment (5 \\mum) of spine'))
xlabel('prefDir')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'31_Local HIdir to Spine.png'))

%Fig. 22/32 Local dispersion vs. HI
figure(22)
plot(localDisp_values(:,2),HI_values(:,2),'*', 'color', 'g');
xlim([0 45])
ylim([0 1])
ylabel(sprintf('HI of local environment (5 \\mum) of spine'))
xlabel(sprintf('circular Dispersion of local environment (5 \\mum) of spine'))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'22_Local HI vs local Dispersion.png'))

figure(32)
plot(localDirDisp_values(:,2),HIDir_values(:,2),'*', 'color', 'g');
xlim([0 90])
ylim([0 1])
ylabel(sprintf('HIdir of local environment (5 \\mum) of spine'))
xlabel(sprintf('circular dir Dispersion of local environment (5 \\mum) of spine'))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'32_Local HIdir vs local dir Dispersion.png'))

%Fig 23/33: local deltaOri
localDeltaOri_values = cell2mat(cellfun(@(x) x.localDeltaOri, num2cell(allSpines(oriSelectROIsNr)), 'UniformOutput', false));
localDeltaOri_values = reshape(localDeltaOri_values, 4, [])';
localDeltaDir_values = cell2mat(cellfun(@(x) x.localDeltaDir, num2cell(allSpines(dirSelectROIsNr)), 'UniformOutput', false));
localDeltaDir_values = reshape(localDeltaDir_values, 4, [])';

figure(23)
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
plot([allSpines(oriSelectROIsNr).prefOri],localDeltaOri_values(:,2),'*', 'color', 'g');
xlim([0 180])
ylim([0 90])
ylabel(sprintf('deltaOri of local environment (5 \\mum) to spine'))
xlabel('pref Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'23_Local deltaOri to Spine.png'))

figure(33)
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
plot([allSpines(dirSelectROIsNr).prefDir],localDeltaDir_values(:,2),'*', 'color', 'g');
xlim([0 360])
ylim([0 180])
ylabel(sprintf('deltaDir of local environment (5 \\mum) to spine'))
xlabel('delta Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'33_Local deltaDir to Spine.png'))

%Fig 24/34: lodal delteOriSoma
localDeltaOriSoma_values = cell2mat(cellfun(@(x) x.localDeltaOriSoma, num2cell(allSpines(oriSelectROIsNr)), 'UniformOutput', false));
localDeltaOriSoma_values = reshape(localDeltaOriSoma_values, 4, [])';
localDeltaDirSoma_values = cell2mat(cellfun(@(x) x.localDeltaDirSoma, num2cell(allSpines(dirSelectROIsNr)), 'UniformOutput', false));
localDeltaDirSoma_values = reshape(localDeltaDirSoma_values, 4, [])';

figure(24)
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
saveas(gcf, fullfile(saveDirImages,'24_Local deltaOri to Soma vs deltaOri Spine.png'))

figure(34)
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
plot([allSpines(dirSelectROIsNr).prefDir],localDeltaDirSoma_values(:,2),'*', 'color', 'g');
xlim([0 180])
ylim([0 180])
ylabel(sprintf('deltaDir of local environment (5 \\mum) to Soma'))
xlabel('pref Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'34_Local deltaDir to Soma vs pref Dir Spine.png'))

%Fig 25/35: nearest SimilarPrefOri/dir vs its own prefOri
figure(25)
plot([allSpines(oriSelectROIsNr).prefOri],[allSpines(oriSelectROIsNr).nearestSimilarPrefOri],'*', 'color', 'g');
xlim([0 180])
ylim([0 25])
ylabel(sprintf('Nearest neighbor with similar preference in \\mum'))
xlabel('pref Ori spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'25_NearestNeighbor prefOri.png'))

figure(35)
plot([allSpines(dirSelectROIsNr).prefDir],[allSpines(dirSelectROIsNr).nearestSimilarPrefDir],'*', 'color', 'g');
xlim([0 360])
ylim([0 25])
ylabel(sprintf('Nearest neighbor with similar preference in \\mum'))
xlabel('pref Dir spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'35_NearestNeighbor prefDir.png'))

%Fig 26/36: local OSI/DSI
%get all localOSI and DSI values & reshape them
localOSI_values = cell2mat(cellfun(@(x) x.localOSI, num2cell(allSpines(respROINrs)), 'UniformOutput', false));
localOSI_values = reshape(localOSI_values, 4, [])';

localDSI_values = cell2mat(cellfun(@(x) x.localDSI, num2cell(allSpines(respROINrs)), 'UniformOutput', false));
localDSI_values = reshape(localDSI_values, 4, [])';

figure(26)
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
plot([allSpines(respROINrs).OSI],localOSI_values(:,2),'*', 'color', 'g');
xlim([0 1])
ylim([0 1])
ylabel(sprintf('local OSI (5 \\mum)'))
xlabel('OSI spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'26_LocalOSI.png'))

figure(36)
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
plot([allSpines(respROINrs).DSI],localDSI_values(:,2),'*', 'color', 'g');
xlim([0 1])
ylim([0 1])
ylabel(sprintf('local DSI (5 \\mum)'))
xlabel('DSI spine')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'36_LocalDSI.png'))

%--------------------------------------------------------------------------
%40s: Properties of dendritic segments
%What is the overal property of the dendritic segment? How homeogeneous is
%it, how different from soma, how selective, ...

%Fig 40: Branch dispersion
figure(40)
subplot(1,2,1)
boxplot([allBranches.circDispersion]','color','b')
ylim([0 90])
hold on
ylabel('Dendritic segment circular dispersion')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.circDispersion],1),1),[allBranches.circDispersion],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,2,2)
boxplot([allBranches.dirCircDispersion]','color','b')
ylim([0 180])
hold on
ylabel('Dendritic segment dir circular dispersion')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.dirCircDispersion],1),1),[allBranches.dirCircDispersion],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'40_circularDispersion.png'))

%Fig 41: Branch selectivity
figure(41)
subplot(1,3,1)
boxplot([allBranches.medianOSI]','color','b')
ylim([0 1])
hold on
ylabel('median OSI')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.medianOSI],1),1),[allBranches.medianOSI],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,3,2)
boxplot([allBranches.medianDSI]','color','b')
ylim([0 1])
hold on
ylabel('median DSI')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.medianDSI],1),1),[allBranches.medianDSI],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,3,3)
boxplot([allBranches.medianDSIvect]','color','b')
ylim([0 1])
hold on
ylabel('median DSIvect')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.medianDSIvect],1),1),[allBranches.medianDSIvect],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'41_Branch selectivity.png'))

%Fig 42: delta Ori/delta Dir
figure(42)
subplot(1,2,1)
boxplot([allBranches.deltaOri]','color','b')
ylim([0 90])
hold on
ylabel('delta Ori')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.deltaOri],1),1),[allBranches.deltaOri],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
subplot(1,2,2)
boxplot([allBranches.deltaDir]','color','b')
ylim([0 180])
hold on
ylabel('delta Dir')
set(gca,'xtick',[])
box off
scatter(repmat(1:1,size([allBranches.deltaDir],1),1),[allBranches.deltaDir],'filled','MarkerFaceColor','blue', 'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'42_DeltaOri_deltaDir.png'))

%Fig 43: Branch circular dispersion vs. deltaOri
figure(43) 
plot([allBranches.deltaOri],[allBranches.circDispersion], '*','color', 'blue')
xlim([0 90])
xlabel('Delta Ori of dendritic segment to soma')
ylim([0 90])
ylabel('Dendritic segment circular dispersion')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'43_DeltaOri_circDispersion.png'))

%Fig 44: Branch circular dispersion dir vs. deltaDir
figure(44) 
plot([allBranches.deltaDir],[allBranches.dirCircDispersion], '*','color', 'blue')
xlim([0 180])
xlabel('Delta dir of dendritic segment to soma')
ylim([0 180])
ylabel('Dendritic segment dir circular dispersion')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages,'44_DeltaDir_dirCircDispersion.png'))

% %--------------------------------------------------------------------------
% %50s: Pairwise distances
% %Fig 50: Pref Ori
% plotPairWisePropertyVsDistance(pwMeasures.distOriSelect,pwMeasures.deltaOri,80);
% ylabel('deltaOri')
% ylim([0 90])
% saveas(gcf, fullfile(saveDirImages,'50_DeltaOri vs. Distance.png'))
% 
% %Fig 51: Pref Dir
% plotPairWisePropertyVsDistance(pwMeasures.distDirSelect,pwMeasures.deltaDir,81);
% ylabel('deltaDir')
% saveas(gcf, fullfile(saveDirImages,'51_DeltaDir vs. Distance.png'))
% 
% %Fig 52: OSI
% plotPairWisePropertyVsDistance(pwMeasures.distResp,pwMeasures.deltaOSI,82);
% ylabel('deltaOSI')
% ylim([0 0.25])
% saveas(gcf, fullfile(saveDirImages,'52_DeltaOSI vs. Distance.png'))
% 
% %Fig 53: DSI
% plotPairWisePropertyVsDistance(pwMeasures.distResp,pwMeasures.deltaDSI,83);
% ylabel('deltaDSI')
% saveas(gcf, fullfile(saveDirImages,'53_DeltaDSI vs. Distance.png'))
% 
% %Fig 54: DSIVect
% plotPairWisePropertyVsDistance(pwMeasures.distResp,pwMeasures.deltaDSIvect,84);
% ylabel('deltaDSIvect')
% saveas(gcf, fullfile(saveDirImages,'54_DeltaDSIvect vs. Distance.png'))

%--------------------------------------------------------------------------
%60s - Cell wide analysis structure
%Fig 60: Distribution of distance from soma for all and separated for
%apical/basal
apicalROIs = cellfun(@(x) strcmp(x, 'apical'),{allSpines.type}); 
basalROIs = cellfun(@(x) strcmp(x, 'basal'),{allSpines.type});

figure(60)
subplot(1,3,1)
distributionPlot([allSpines.distToSoma]', 'color', 'black'); hold all
boxplot([allSpines.distToSoma]')
title('All spines')
ylabel('Distance from branch start in um')
box off
ylim([0 max([allSpines.distToSoma])])
subplot(1,3,2)
distributionPlot([allSpines(apicalROIs).distToSoma]', 'color', 'red'); hold all
boxplot([allSpines(apicalROIs).distToSoma]')
title('Apical Spines')
ylim([0 max([allSpines.distToSoma])])
axis off
subplot(1,3,3)
distributionPlot([allSpines(basalROIs).distToSoma]', 'color', 'green'); hold all
boxplot([allSpines(basalROIs).distToSoma]')
title('basal Spines')
ylim([0 max([allSpines.distToSoma])])
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '60_SpineDistanceFromBranch.png'))

%Fig 61: Same as Fig 60, but as cdfplot
figure(61);
h(1,1) = cdfplot([allSpines.distToSoma]); hold on
set(h(1,1), 'Color', 'black', 'LineWidth', 3);
h(1,2) = cdfplot([allSpines(apicalROIs).distToSoma]);
set(h(1,2), 'Color', 'red', 'LineWidth', 3);
h(1,3) = cdfplot([allSpines(basalROIs).distToSoma]);
set(h(1,3), 'Color', 'green', 'LineWidth', 3);
grid off
title('')
xlabel('Distance from start of branch')
ylabel('Fraction of spines')
set(gcf, 'color', 'w');
legend('All spines', 'Apical Spines', 'Basal Spines', 'Location', 'SouthEast')
saveas(gcf, fullfile(saveDirImages, '61_DistFromBranchCumul.png'))

%Fig 62: Histogram of branch order for all and separated for apical/basal
figure(22)
binEdges = linspace(1,max([allSpines.BranchOrder]), max([allSpines.BranchOrder]));
subplot(1,3,1)
histogram([allSpines.BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', 'white')
title('All spines')
box off
subplot(1,3,2)
histogram([allSpines(apicalROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', 'red')
title('Apical spines')
box off
subplot(1,3,3)
histogram([allSpines(basalROIs).BranchOrder],binEdges, 'Normalization','probability', 'FaceColor', 'green')
title('Basal spines')
box off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '62_BranchOrder.png'))

%Fig 63: same as Fig 62, but as cdfplot
figure(63);
h(1,1) = cdfplot([allSpines.BranchOrder]); hold on
set(h(1,1), 'Color', 'black', 'LineWidth', 3);
h(1,2) = cdfplot([allSpines(apicalROIs).BranchOrder]);
set(h(1,2), 'Color', 'red', 'LineWidth', 3);
h(1,3) = cdfplot([allSpines(basalROIs).BranchOrder]);
set(h(1,3), 'Color', 'green', 'LineWidth', 3);
grid off
title('')
xlabel('Branchorder')
ylabel('Fraction of spines')
set(gcf, 'color', 'w');
legend('All spines', 'Apical Spines', 'Basal Spines', 'Location', 'SouthEast')
saveas(gcf, fullfile(saveDirImages, '63_BranchOrderCumul.png'))

%70s - Branch specific analysis
%Fig 70: Number of spines per dendrite for all and separated for
%apical/basal
apicalBranches = cellfun(@(x) strcmp(x, 'apical'),{allBranches.type}); 
basalBranches = cellfun(@(x) strcmp(x, 'basal'),{allBranches.type});

figure(70)
subplot(1,3,1)
distributionPlot([allBranches.numSpines]', 'color', 'black'); hold all
boxplot([allBranches.numSpines]')
title('All dendrites')
ylim([0 150])
ylabel('Spines per dendrite')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).numSpines]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).numSpines]')
title('Apical dendrites')
ylim([0 150])
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).numSpines]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).numSpines]')
ylim([0 150])
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '70_Spines per dendrite.png'))

%Fig 71: Spine density per dendrite for all and separated for apical/basal
figure(71)
subplot(1,3,1)
distributionPlot([allBranches.SpineDensity]', 'color', 'black'); hold all
boxplot([allBranches.SpineDensity]')
title('All dendrites')
ylim([0 5])
ylabel('Spine density')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).SpineDensity]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).SpineDensity]')
title('Apical dendrites')
ylim([0 5])
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).SpineDensity]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).SpineDensity]')
ylim([0 5])
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '71_Spine density.png'))

%Fig 72: circular dispersion
figure(72)
subplot(1,3,1)
distributionPlot([allBranches.circDispersion]', 'color', 'black'); hold all
boxplot([allBranches.circDispersion]')
ylabel('Circular Dispersion')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).circDispersion]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).circDispersion]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).circDispersion]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).circDispersion]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '72_circDispersion.png'))

%Fig 73: directional circular dispersion
figure(73)
subplot(1,3,1)
distributionPlot([allBranches.dirCircDispersion]', 'color', 'black'); hold all
boxplot([allBranches.dirCircDispersion]')
ylabel('dir circular dispersion')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).dirCircDispersion]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).dirCircDispersion]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).dirCircDispersion]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).dirCircDispersion]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '73_dirCircDispersion.png'))

%Fig 74: deltaOri
figure(74)
subplot(1,3,1)
distributionPlot([allBranches.deltaOri]', 'color', 'black'); hold all
boxplot([allBranches.deltaOri]')
ylabel('deltaOri')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).deltaOri]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).deltaOri]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).deltaOri]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).deltaOri]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '74_deltaOri.png'))

%Fig 75: deltaDir
figure(75)
subplot(1,3,1)
distributionPlot([allBranches.deltaDir]', 'color', 'black'); hold all
boxplot([allBranches.deltaDir]')
ylabel('deltaOri')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).deltaDir]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).deltaDir]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).deltaDir]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).deltaDir]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '75_deltaDir.png'))

%Fig 76: medianOSI
figure(76)
subplot(1,3,1)
distributionPlot([allBranches.medianOSI]', 'color', 'black'); hold all
boxplot([allBranches.medianOSI]')
ylabel('medianOSI')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).medianOSI]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).medianOSI]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).medianOSI]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).medianOSI]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '76_medianOSI.png'))

%Fig 77: medianDSI
figure(77)
subplot(1,3,1)
distributionPlot([allBranches.medianDSI]', 'color', 'black'); hold all
boxplot([allBranches.medianDSI]')
ylabel('medianDSI')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).medianDSI]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).medianDSI]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).medianDSI]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).medianDSI]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '77_medianDSI.png'))

%Fig 78: medianDSIvect
figure(78)
subplot(1,3,1)
distributionPlot([allBranches.medianDSIvect]', 'color', 'black'); hold all
boxplot([allBranches.medianDSIvect]')
ylabel('medianDSIvect')
title('All dendrites')
box off
subplot(1,3,2)
distributionPlot([allBranches(apicalBranches).medianDSIvect]', 'color', 'red'); hold all
boxplot([allBranches(apicalBranches).medianDSIvect]')
title('Apical dendrites')
box off
axis off
subplot(1,3,3)
distributionPlot([allBranches(basalBranches).medianDSIvect]', 'color', 'green'); hold all
boxplot([allBranches(basalBranches).medianDSIvect]')
title('Basal dendrites')
box off
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirImages, '78_medianDSIvect.png'))

%% 6.) Summary

disp(['Number of cells: ' num2str(cellNr-1, '%d')])
disp(['Number of spines: ' num2str(length(allSpines), '%d')])
disp(['Number of 2p matched spines: ' num2str(length(TwoPROIsNR), '%d')])
disp(['Number of responsive spines: ' num2str(length(respROINrs), '%d')])
disp(['Number of orientation-selective spines: ' num2str(length(oriSelectROIsNr), '%d')])

%% 7.) Save
%write to text file
fid = fopen([saveDir filesep 'MultiCellSummary.txt'], 'w');
formatSpec = ['Number of cells: %d \n'...
    'Number of spines: %d \n'...
    'Numbrt og 2p matched spines: %d \n' ...
    'Number of responsive spines: %d \n'...
    'Number of orientation-selective spines: %d'];
fprintf(fid, formatSpec, cellNr-1, length(allSpines), length(TwoPROIsNR), length(respROINrs), length(oriSelectROIsNr));
fclose(fid);

%save variables
save([saveDir filesep 'SpineData.mat'], 'allSpines', 'allBranches', 'pwMeasures','cellDetails', '-mat') 

