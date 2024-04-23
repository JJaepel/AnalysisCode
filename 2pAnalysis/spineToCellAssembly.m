function spineToCellAssembly(cellInfo, i)

%This function links the spines back to their original cells on an animal 
%base, does calculations of their relationship and saves everything both in
%a table and a .mat file
% Step 1: List all experiments
% Step 2: Load all files and calc deltaOri
% Step 3: Construct data holders for each cell
% Step 4: Make plot 
% Step 5: Save
%
% Input:
%  - cellInfo: contains information about all the cells that we are
%  analyzing
%  - i: which cell are you intersted in
% (- xls sheet with exp infos)
% (-ROIsAna.mat for each experiment: contains all ROIs and their
% characterization)
%
%Output:
% - xls File: table with overview of all ROIs
% - completeCell.mat: File with variable completeCell
% - Figure that shows the correlation fo all dendrites to the soma as a
% quick overview on how diverse the cell is

%% 1.) List all experiments & define saveDirs on a cell basis
%where do we find the experiments?
filePath =  'Z:\Juliane\Organization\Animals\';
file = 'SpinePerAnimal.xlsx';
sheetName = 'bimodal';

%get all the information about the experiments
expInfo = findExpInfo([filePath file], sheetName);

%select experiments containing the cell
animalID = cellInfo.animal(i);
onlyAnimal= cellfun(@(x) find(contains(x, char(animalID))),expInfo.animal,'UniformOutput',false);
AnimalExp = find(~cellfun(@isempty,onlyAnimal)); 
allCellNrs = unique(cell2mat(cellInfo.cellNr(find(~cellfun('isempty', onlyAnimal)))));

%get the spine & cell experiments
spineExp = find(cellfun(@isempty,cellfun(@(x) find(contains(x, 'cells')),expInfo.region,'UniformOutput',false)));
spineExp = intersect(AnimalExp, spineExp);
cellExp = find(~cellfun(@isempty,cellfun(@(x) find(contains(x, 'cells')),expInfo.region,'UniformOutput',false)));
cellExp = intersect(AnimalExp, cellExp);

if iscellstr(cellInfo.cellName(i)) 
    saveDir = ['Z:\Juliane\InputAnalysis\' char(animalID) filesep char(cellInfo.cellName(i)) filesep 'A - 2p Imaging' filesep];
else
    saveDir = ['Z:\Juliane\InputAnalysis\' char(animalID) filesep 'A - 2p Imaging' filesep];
end
%% 2.) Load all files & calculate deltaOri

%cellData
if iscellstr(cellInfo.cellName(i)) 
    loadDirCell = ['Z:\Juliane\InputAnalysis\' char(animalID) filesep char(cellInfo.cellName(i)) filesep 'A - 2p Imaging' filesep char(expInfo.name{cellExp}) filesep];
else
    loadDirCell = ['Z:\Juliane\InputAnalysis\' char(animalID) filesep 'A - 2p Imaging' filesep char(expInfo.name{cellExp}) filesep];
end
load([loadDirCell 'ROIsAna.mat'], 'ce');
allCellData = ce; clear ce;

%SpineData & DendriteData
for s = 1:length(spineExp)
    %first check if that experiment is part of that cell
    
    
    if iscellstr(cellInfo.cellName(i)) 
        loadDirSpine = ['Z:\Juliane\InputAnalysis\' char(animalID) filesep char(cellInfo.cellName(i)) filesep 'A - 2p Imaging' filesep char(expInfo.name{spineExp(s)}) filesep];
    else
        loadDirSpine = ['Z:\Juliane\InputAnalysis\' char(animalID) filesep 'A - 2p Imaging' filesep char(expInfo.name{SpineExp(s)}) filesep];
    end
    load([loadDirSpine 'ROIsAna.mat'], 'ce');
    CellROINr = expInfo.cellROIInd(SpineExp(s));
    SpineData = ce; clear ce;
    for r = 1:length(SpineData)
        deltaOri = abs(SpineData(r).prefOri - allCellData(CellROINr).prefOri);
        if deltaOri >90
            SpineData(r).deltaOri = 180 - deltaOri;
        else
            SpineData(r).deltaOri = deltaOri;
        end
        SpineData(r).cellROINr = CellROINr;
    end
    allSpineData{s} = SpineData;
    allCellROIs(s) = CellROINr;
end


%% 3.) Construct data for this cell

%write cell info 
completeCell{c}.animalID = animalID;
completeCell{c}.expNr = expInfo.name{cellExp};
completeCell{c}.ROINr = cellNrs(c);
completeCell{c}.depth = allCellData(cellNrs(c)).depth;
completeCell{c}.denType = allCellData(cellNrs(c)).denType;
completeCell{c}.cyc = allCellData(cellNrs(c)).cyc;
completeCell{c}.meanResp = allCellData(cellNrs(c)).meanResp;
completeCell{c}.peakFit = allCellData(cellNrs(c)).peakFit;
completeCell{c}.prefOri = allCellData(cellNrs(c)).prefOri;
completeCell{c}.prefDir = allCellData(cellNrs(c)).prefDir;
completeCell{c}.OSI = allCellData(cellNrs(c)).OSI;
completeCell{c}.DSI = allCellData(cellNrs(c)).DSI;
completeCell{c}.DSIvect = allCellData(cellNrs(c)).DSIvect;
completeCell{c}.Bandwidth = allCellData(cellNrs(c)).bandwidth;

cellTable=table(string(completeCell{c}.denType), "na","na","na", ...
    str2double(completeCell{c}.depth), "na", completeCell{c}.prefOri,completeCell{c}.prefDir,completeCell{c}.OSI, ...
    completeCell{c}.DSI, completeCell{c}.DSIvect, completeCell{c}.Bandwidth, "na", "na");    
%get spineData & dendriteData
CellSpines = struct;
CellDendrites = struct;
spineCounter = 1;
dendCounter = 1;
for s = 1:size(allSpineData,2)
    if allSpineData{s}(1).cellROINr == cellNrs(c)
        for sp = 1:size(allSpineData{s},2)
            if allSpineData{s}(sp).spine

                baseTwoP = char(expInfo.name{SpineExp(s)});
                try
                   CellSpines(spineCounter).uniqueID = [char(expInfo.altName{SpineExp(s)}) '-' num2str(expInfo.dendROIInd(SpineExp(s)), '%02d') '-' num2str(sp,'%02d')];
                catch
                   CellSpines(spineCounter).uniqueID = str2double(baseTwoP(2:6))*100+sp;
                end
                CellSpines(spineCounter).ID = spineCounter;
                CellSpines(spineCounter).expNr =  allSpineData{s}(sp).file;
                CellSpines(spineCounter).dendriteID =  expInfo.dendROIInd(SpineExp(s));
                CellSpines(spineCounter).spineROINr =  sp;
                CellSpines(spineCounter).depth =  allSpineData{s}(sp).depth;
                CellSpines(spineCounter).denType =  allSpineData{s}(sp).denType;
                CellSpines(spineCounter).cycRes = allSpineData{s}(sp).cycRes;
                CellSpines(spineCounter).meanResp = allSpineData{s}(sp).meanResp;
                CellSpines(spineCounter).peakFit = allSpineData{s}(sp).peakFit;
                CellSpines(spineCounter).prefOri = allSpineData{s}(sp).prefOri;
                CellSpines(spineCounter).prefDir = allSpineData{s}(sp).prefDir;
                CellSpines(spineCounter).OSI = allSpineData{s}(sp).OSI;
                CellSpines(spineCounter).DSI = allSpineData{s}(sp).DSI;
                CellSpines(spineCounter).DSIvect = allSpineData{s}(sp).DSIvect;
                CellSpines(spineCounter).corr = allSpineData{s}(sp).corr;
                CellSpines(spineCounter).good = allSpineData{s}(sp).good;
                CellSpines(spineCounter).deltaOri = allSpineData{s}(sp).deltaOri;
                CellSpines(spineCounter).events = allSpineData{s}(sp).events;
                CellSpines(spineCounter).eventAmps = allSpineData{s}(sp).eventAmps;
                CellSpines(spineCounter).Bandwidth = allSpineData{s}(sp).bandwidth;

                spineTable = table(string(CellSpines(spineCounter).denType),CellSpines(spineCounter).dendriteID, CellSpines(spineCounter).spineROINr, string(CellSpines(spineCounter).uniqueID), ...
                    str2double(CellSpines(spineCounter).depth), CellSpines(spineCounter).deltaOri, CellSpines(spineCounter).prefOri,CellSpines(spineCounter).prefDir,CellSpines(spineCounter).OSI,...
                    CellSpines(spineCounter).DSI, CellSpines(spineCounter).DSIvect,CellSpines(spineCounter).Bandwidth, CellSpines(spineCounter).corr, CellSpines(spineCounter).good...
                    );
                cellTable = [cellTable; spineTable];                        
                spineCounter = spineCounter + 1;
            else
                baseTwoP = char(expInfo.name{SpineExp(s)});
                try
                    CellDendrites(dendCounter).uniqueID = [char(expInfo.altName{SpineExp(s)}) '-' num2str(expInfo.dendROIInd(SpineExp(s)), '%02d') '-' num2str(sp,'%02d')];
                catch
                    CellDendrites(dendCounter).uniqueID = str2double(baseTwoP(2:6))*100+sp;
                end 
                CellDendrites(dendCounter).ID = dendCounter;
                CellDendrites(dendCounter).expNr =  allSpineData{s}(sp).file;
                CellDendrites(dendCounter).dendriteID =  expInfo.dendROIInd(SpineExp(s));
                CellDendrites(dendCounter).depth =  allSpineData{s}(sp).depth;
                CellDendrites(dendCounter).denType =  allSpineData{s}(sp).denType;
                CellDendrites(dendCounter).cyc = allSpineData{s}(sp).cyc;
                CellDendrites(dendCounter).meanResp = allSpineData{s}(sp).meanResp;
                CellDendrites(dendCounter).good = evalDend(CellDendrites(dendCounter).cyc);
                CellDendrites(dendCounter).peakFit = allSpineData{s}(sp).peakFit;
                CellDendrites(dendCounter).prefOri = allSpineData{s}(sp).prefOri;
                CellDendrites(dendCounter).prefDir = allSpineData{s}(sp).prefDir;
                CellDendrites(dendCounter).OSI = allSpineData{s}(sp).OSI;
                CellDendrites(dendCounter).DSI = allSpineData{s}(sp).DSI;
                CellDendrites(dendCounter).DSIvect = allSpineData{s}(sp).DSIvect;
                CellDendrites(dendCounter).Bandwidth = allSpineData{s}(sp).bandwidth;
                CellDendrites(dendCounter).deltaOri = allSpineData{s}(sp).deltaOri;
                CellDendrites(dendCounter).events = allSpineData{s}(sp).events;
                CellDendrites(dendCounter).eventAmps = allSpineData{s}(sp).eventAmps;
                dendCounter = dendCounter + 1;
            end
        end
    end
end
completeCell{c}.SpineData = CellSpines;
completeCell{c}.DendriteData = CellDendrites;
cellTable.Properties.VariableNames = ["dendrite Type","dendrite ID", "spine ID", "unique ID", ...
   "depth","deltaOri", "preferred Orientation", "preferred Direction", "OSI", "DSI", "vectorized DSI", "Bandwidth",...
   "dendrite Correlation", "spine Evaluation"];
modality = expInfo.Modality{CellExp};
switch modality{1}
    case 'EM'
        saveName = [char(expInfo.altName{CellExp}) '_' char(animalID) '_CellSummary'];
    otherwise
        saveName = [char(animalID) '_CellSummary'];
end
writetable(cellTable,[saveDir saveName '.xls'],'Sheet',['Cell ' num2str(c)])
%% 4.) Check plots

%select the dendrites that are good and selective
goodDend = find([CellDendrites.good] == 1);
selectDend = find([CellDendrites.OSI] > 0.1);
goodSelectDend = intersect(goodDend, selectDend);

%get their deltaOri and uniqueID
deltaOriGSD = [CellDendrites(goodSelectDend).deltaOri];
uniqueGSD = [CellDendrites(goodSelectDend).uniqueID];

plot(uniqueGSD,deltaOriGSD, '*')
xlabel('Dendrite ID')
ylabel('Delta Ori to soma')
set(gcf, 'color', 'w');
box off
saveas(gcf, fullfile(saveDir,'delteOriDendritesToSoma.png'))

%% 5.) Save the data
save([saveDir saveName '.mat'], 'completeCell', '-mat') 
