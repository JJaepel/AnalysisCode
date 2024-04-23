function spineToCellComputations(animalID)

%This function links the spines back to their original cells on an animal 
%base, does calculations of their relationship and saves everything both in
%a table and a .mat file
% Step 1: List all experiments
% Step 2: Load all files per cell, calc deltaOri & construct data holder
% for each cell
% Step 3: Make plots
% Step 4: Write it to xls file & make quick summary
% Step 5: Save
%
% Input:
%  - animalID: Which animal are you interested in?
% (- xls sheet with exp infos)
% (- ROIsAna.mat for each experiment: contains all ROIs and their
% characterization)
%
%Output:
% - xls File: table with overview of all ROIs
% - completeCell.mat: File with variable completeCell
% - Figure that shows the correlation fo all dendrites to the soma as a
% quick overview on how diverse the cell is

%% 1.) List all experiments
%where do we find the experiments?
filePath =  'Z:\Juliane\Organization\Animals\';
file = 'SpinePerAnimal.xlsx';
sheetName = 'bimodal';

%get all the information about the experiments
expInfo = findExpInfo([filePath file], sheetName);

%select experiments containing the animal
onlyAnimal= cellfun(@(x) find(contains(x, char(animalID))),expInfo.animal,'UniformOutput',false);
AnimalExp = find(~cellfun(@isempty,onlyAnimal)); 
%AnimalExp = setdiff(ind,find(onlyOthers));

%get the cellID
onlyCells = cellfun(@(x) find(contains(x, 'cells')),expInfo.region,'UniformOutput',false);
onlySpines = find(cellfun(@isempty,onlyCells)); 
CellInd = find(~cellfun(@isempty,onlyCells)); 
CellExp = intersect(AnimalExp, CellInd);
SpineExp = intersect(AnimalExp, onlySpines);

%get the info about the experiment to know if it is one of the included
%animals, which changes where they are being saved
info = makeInfoStructExp(0, expInfo, CellExp(1));

%% 2.) Load all files, calculate deltaOri & construct data for the cell

%lets find out how many cells we have in this animal
cellNrs = unique([expInfo.cellROIInd(SpineExp)]);

for c = 1:length(cellNrs)
    if str2double(info.include)
        if length(cellNrs) > 1
            cellInfo = animalParser; %get the cellInfo
            cellAnimalInfo = find(~cellfun(@isempty,cellfun(@(x) find(contains(x, info.animalShort)),cellInfo.animal,'UniformOutput',false))); %which ones contain the animal
            cellNrInfo = find(cellfun(@(x) ~isempty(x) && x == c, cellInfo.cellNr)); %which ones have the right cellNr?
            i = intersect(cellAnimalInfo, cellNrInfo);
            if ~isempty(i)
                saveDir = ['Z:\Juliane\InputAnalysis\' char(info.animalShort) filesep char(cellInfo.cellName{i(1)}) filesep 'A - 2p Imaging' filesep];
                loadDirCell = [saveDir cellInfo.somaExpNr{i(1)} filesep];
                cellID = cellNrs(cellInfo.cellNr{i(1)});
            else
                continue
            end
        else
            saveDir = ['Z:\Juliane\InputAnalysis\' char(info.animalShort) filesep 'A - 2p Imaging' filesep];
            loadDirCell = cell2mat(strcat(saveDir, expInfo.name{CellExp(1)}, filesep));
            cellID = 1;
        end
    else
        saveDir = ['Z:\Juliane\Data\ImageAnalysis\' char(info.animal) filesep];
        loadDirCell = cell2mat(strcat(saveDir, expInfo.name{CellExp(1)}, filesep));
        cellID = cellNrs(c);
    end
    
    %load cellData
    load([loadDirCell 'ROIsAna.mat'], 'ce');
    allCellData = ce; clear ce;
    
    %write cell info 
    completeCell{c}.animalID = animalID;
    completeCell{c}.expNr = allCellData(cellID).file;
    completeCell{c}.ROINr = cellID;
    completeCell{c}.depth = allCellData(cellID).depth;
    completeCell{c}.denType = allCellData(cellID).denType;
    completeCell{c}.cyc = allCellData(cellID).cyc;
    completeCell{c}.meanResp = allCellData(cellID).meanResp;
    completeCell{c}.peakFit = allCellData(cellID).peakFit;
    completeCell{c}.prefOri = allCellData(cellID).prefOri;
    completeCell{c}.prefDir = allCellData(cellID).prefDir;
    completeCell{c}.OSI = allCellData(cellID).OSI;
    completeCell{c}.DSI = allCellData(cellID).DSI;
    completeCell{c}.DSIvect = allCellData(cellID).DSIvect;
    completeCell{c}.Bandwidth = allCellData(cellID).bandwidth;
    
    cellTable=table(string(completeCell{c}.denType), "na","na","na", ...
        str2double(completeCell{c}.depth), "na", completeCell{c}.prefOri,completeCell{c}.prefDir,completeCell{c}.OSI, ...
        completeCell{c}.DSI, completeCell{c}.DSIvect, completeCell{c}.Bandwidth, "na", "na"); 

    %SpineData & DendriteData - make structures
    CellSpines = struct;
    CellDendrites = struct;
    spineCounter = 1;
    dendCounter = 1;
    
    %go through all files
    for s = 1:length(SpineExp)
        if expInfo.cellROIInd(SpineExp(s)) == cellID %only if that file is part of that cell, load it
            loadDirSpine = cell2mat(strcat(saveDir, expInfo.name{SpineExp(s)}, filesep));
            load([loadDirSpine 'ROIsAna.mat'], 'ce');
            SpineData = ce; clear ce; 
            
            for r = 1:length(SpineData)
                %calculate the deltaOri
                deltaOri = abs(SpineData(r).prefOri - allCellData(cellID).prefOri);
                if deltaOri >90
                    SpineData(r).deltaOri = 180 - deltaOri;
                else
                    SpineData(r).deltaOri = deltaOri;
                end
                
                %calculate the deltaDir
                deltaDir = abs(SpineData(r).prefDir - allCellData(cellID).prefDir);
                if deltaOri >180
                    SpineData(r).deltaDir = 360 - deltaDir;
                else
                    SpineData(r).deltaDir = deltaDir;
                end
                
                %check if it is a spine or dendrite
                if SpineData(r).spine
                    %if it is a spine, add it to the spine structure
                    baseTwoP = char(expInfo.name{SpineExp(s)});
                    try
                       CellSpines(spineCounter).uniqueID = [char(expInfo.altName{SpineExp(s)}) '-' num2str(expInfo.dendROIInd(SpineExp(s)), '%02d') '-' num2str(sp,'%02d')];
                    catch
                       CellSpines(spineCounter).uniqueID = str2double(baseTwoP(2:6))*200+r;
                    end
                    %save the relevant info in the structure
                    CellSpines(spineCounter).ID = spineCounter;
                    CellSpines(spineCounter).expNr =  SpineData(r).file;
                    CellSpines(spineCounter).dendriteID =  expInfo.dendROIInd(SpineExp(s));
                    CellSpines(spineCounter).spineROINr = r;
                    CellSpines(spineCounter).depth =  SpineData(r).depth;
                    CellSpines(spineCounter).denType =  SpineData(r).denType;
                    CellSpines(spineCounter).cycRes = SpineData(r).cycRes;
                    CellSpines(spineCounter).meanResp = SpineData(r).meanResp;
                    CellSpines(spineCounter).peakFit = SpineData(r).peakFit;
                    CellSpines(spineCounter).prefOri = SpineData(r).prefOri;
                    CellSpines(spineCounter).prefDir = SpineData(r).prefDir;
                    CellSpines(spineCounter).OSI = SpineData(r).OSI;
                    CellSpines(spineCounter).DSI = SpineData(r).DSI;
                    CellSpines(spineCounter).DSIvect = SpineData(r).DSIvect;
                    CellSpines(spineCounter).corr = SpineData(r).corr;
                    CellSpines(spineCounter).good = SpineData(r).good;
                    CellSpines(spineCounter).deltaOri = SpineData(r).deltaOri;
                    CellSpines(spineCounter).deltaDir = SpineData(r).deltaDir;
                    CellSpines(spineCounter).events = SpineData(r).events;
                    CellSpines(spineCounter).eventAmps = SpineData(r).eventAmps;
                    CellSpines(spineCounter).Bandwidth = SpineData(r).bandwidth;
                    
                    %add it to the table
                    spineTable = table(string(CellSpines(spineCounter).denType),CellSpines(spineCounter).dendriteID, CellSpines(spineCounter).spineROINr, string(CellSpines(spineCounter).uniqueID), ...
                        str2double(CellSpines(spineCounter).depth), CellSpines(spineCounter).deltaOri, CellSpines(spineCounter).prefOri,CellSpines(spineCounter).prefDir,CellSpines(spineCounter).OSI,...
                        CellSpines(spineCounter).DSI, CellSpines(spineCounter).DSIvect,CellSpines(spineCounter).Bandwidth, CellSpines(spineCounter).corr, CellSpines(spineCounter).good...
                        );
                    cellTable = [cellTable; spineTable];
                    
                    %go up with the spineCounter
                    spineCounter = spineCounter + 1;
                
                else
                    %if it is a dendrite, add it to the dendrite
                    baseTwoP = char(expInfo.name{SpineExp(s)});
                    try
                        CellDendrites(dendCounter).uniqueID = [char(expInfo.altName{SpineExp(s)}) '-' num2str(expInfo.dendROIInd(SpineExp(s)), '%02d') '-' num2str(sp,'%02d')];
                    catch
                        CellDendrites(dendCounter).uniqueID = str2double(baseTwoP(2:6))*200+r;
                    end 
                    %save the relevant info in the structure
                    CellDendrites(dendCounter).ID = dendCounter;
                    CellDendrites(dendCounter).expNr =  SpineData(r).file;
                    CellDendrites(dendCounter).dendriteID =  expInfo.dendROIInd(SpineExp(s));
                    CellDendrites(dendCounter).depth =  SpineData(r).depth;
                    CellDendrites(dendCounter).denType =  SpineData(r).denType;
                    CellDendrites(dendCounter).cyc = SpineData(r).cyc;
                    CellDendrites(dendCounter).meanResp = SpineData(r).meanResp;
                    CellDendrites(dendCounter).good = evalDend(CellDendrites(dendCounter).cyc);
                    CellDendrites(dendCounter).peakFit = SpineData(r).peakFit;
                    CellDendrites(dendCounter).prefOri = SpineData(r).prefOri;
                    CellDendrites(dendCounter).prefDir = SpineData(r).prefDir;
                    CellDendrites(dendCounter).OSI = SpineData(r).OSI;
                    CellDendrites(dendCounter).DSI = SpineData(r).DSI;
                    CellDendrites(dendCounter).DSIvect = SpineData(r).DSIvect;
                    CellDendrites(dendCounter).Bandwidth = SpineData(r).bandwidth;
                    CellDendrites(dendCounter).deltaOri = SpineData(r).deltaOri;
                    CellDendrites(dendCounter).events = SpineData(r).events;
                    CellDendrites(dendCounter).eventAmps = SpineData(r).eventAmps;
                    
                    %go up with the spineCounter
                    dendCounter = dendCounter + 1;
                end
                
            end
        end
    end
    
    completeCell{c}.SpineData = CellSpines;
    completeCell{c}.DendriteData = CellDendrites;
    
    %% 3.) Check plots

    %select the dendrites that are good and selective
    goodDend = find([CellDendrites.good] == 1);
    selectDend = find([CellDendrites.OSI] > 0.1);
    goodSelectDend = intersect(goodDend, selectDend);

    %get their deltaOri and uniqueID
    deltaOriGSD = [CellDendrites(goodSelectDend).deltaOri];
    uniqueGSD = [CellDendrites(goodSelectDend).uniqueID];
    
    %X1: deltaOri vs. dendrite ID
    figure(c*10+1)
    plot(uniqueGSD,deltaOriGSD, '*')
    xlabel('Dendrite ID')
    ylabel('Delta Ori to soma')
    set(gcf, 'color', 'w');
    box off
    saveas(gcf, fullfile(saveDir,['delteOriDendritesToSoma_Cell' num2str(c) '.png']))
    
    %X2: pies for responsiveness, ori select & dir select
    figure(c*10+2)
    subplot(1,2,1)
    nonResp = length(find([completeCell{c}.SpineData.good] == 0))/length(completeCell{c}.SpineData);
    oriSelect = length(find([completeCell{c}.SpineData.good] == 1 & [completeCell{c}.SpineData.OSI] > 0.15))/length(completeCell{c}.SpineData);
    h = pie([oriSelect 1-nonResp-oriSelect nonResp]);
    hp = findobj(h, 'Type', 'patch');
    set(hp(1), 'FaceColor', 'red');
    set(hp(2), 'FaceColor', 'white');
    set(hp(3), 'FaceColor', 'black');
    legend({'Ori-select','Resp','Non-resp'}, 'Location', 'southoutside')
    legend('boxoff')

    subplot(1,2,2)
    dirSelect = length(find([completeCell{c}.SpineData.good] == 1 & [completeCell{c}.SpineData.DSIvect] > 0.15))/length(completeCell{c}.SpineData);
    h = pie([dirSelect 1-nonResp-dirSelect nonResp]);
    hp = findobj(h, 'Type', 'patch');
    set(hp(1), 'FaceColor', 'red');
    set(hp(2), 'FaceColor', 'white');
    set(hp(3), 'FaceColor', 'black');
    legend({'Dir-select','Resp','Non-resp'}, 'Location', 'southoutside')
    legend('boxoff')
    set(gcf, 'color', 'w')
    saveas(gcf, fullfile(saveDir,['PercentageSelectivity_Cell' num2str(c) '.png']))
    
    %X3: histogram or oriplot for pref ori & dir of all spines
    angles = linspace(0,180,9);

    oriSelectROIsNr = find([completeCell{c}.SpineData.good] == 1 & [completeCell{c}.SpineData.OSI] > 0.15);
    [counts, ~] = histcounts([completeCell{c}.SpineData(oriSelectROIsNr).prefOri], [0 linspace(22.5/2,180-22.5/2,8) 180]);
    countsAdjusted = [counts(1)+counts(end) counts(2:end-1) counts(1)+counts(end)];
    countsNormalized = countsAdjusted/max(countsAdjusted);

    SomaResp = (completeCell{c}.meanResp(1:8)+completeCell{c}.meanResp(9:16))/2;
    SomaRespNorm = (SomaResp)/(max(SomaResp));

    figure(c*10+3)
    polarPlotOri(angles, [SomaRespNorm SomaRespNorm(1)], [], 'black','-',0)
    hold on
    polarPlotOri(angles, countsNormalized, [], 'red')
    saveas(gcf, fullfile(saveDir,['OriTuning_Cell' num2str(c) '.png']))
    
    %% 4.) Write it to xls file & make quick summary
    %get the heads for the table
    cellTable.Properties.VariableNames = ["dendrite Type","dendrite ID", "spine ID", "unique ID", ...
       "depth","deltaOri", "preferred Orientation", "preferred Direction", "OSI", "DSI", "vectorized DSI", "Bandwidth",...
       "dendrite Correlation", "spine Evaluation"];
   
   %change the name depending on the modality
    modality = expInfo.Modality{CellExp};
    switch modality{1}
        case 'EM'
            saveName = cell2mat(strcat(expInfo.altName{1}, '_', char(animalID), '_CellSummary_Cell_', num2str(c)));
        otherwise
            saveName = [char(animalID) '_CellSummary_' num2str(c)];
    end
    
    %write the table
    writetable(cellTable,[saveDir saveName '.xls'],'Sheet',['Cell ' num2str(c)])
    clear cellTable
    
    %CellSummary with numbers    
    fid = fopen([saveDir filesep 'QuickSummary_Cell' num2str(c) '.txt'], 'w');
    formatSpec = ['total spines imaged: %d \n'...
        'responsive spines: %d \n'...
        'ori-selective spines: %d \n'...
        'dir-selective spines: %d \n'];
    fprintf(fid, formatSpec, length([completeCell{c}.SpineData]), length(find([completeCell{c}.SpineData.good] == 1)),oriSelect*length([completeCell{c}.SpineData]), dirSelect*length([completeCell{c}.SpineData]));
    fclose(fid);
    
    %% 5.) Save the data    
    %write the Matlab file
    save([saveDir saveName '.mat'], 'completeCell', '-mat') 
end
