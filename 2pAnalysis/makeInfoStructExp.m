function info = makeInfoStructExp(server, expInfo, i)

% This function combines the information from the expInfo according to the
% analysisParams, basically mostly a rewriting task
%
% Inputs
% - analysisParams: what kind of steps are we doing? where do we run?
% - expInfo: all the info from the excel sheet
% - i: which experiment are we looking at right now?
% 
% Output:
% - info: structure with all the relevant information for the following
% analysis that will be used all over in the coming steps

%% Step 0: Initialize the structure & get the essential 
info = struct;
info.animal = char(expInfo.animal{i});
info.name = char(expInfo.name{i});
info.sp2 = char(expInfo.sp2_id{i});
info.cellNr = num2str(expInfo.cellROIInd(i));
info.ROINr = num2str(expInfo.dendROIInd(i));
info.depth = num2str(expInfo.Depth(i));
info.denType = char(expInfo.denType{i});
info.Modality = char(expInfo.Modality{i});

%% Step 1: Based on the server and selected experiments, make a few decisions

%do we work from server or RAID?
if server
    info.drive = 'Z:\Juliane\Data\';
else
    info.drive = 'F:\Data\';
end

info.isCell = contains(expInfo.region{i}, 'cells');
info.include =  num2str(expInfo.include(i)); %is it part of input analysis?

%but if we are looking at experiment to include in the final experiments
%and we are not running it for the first time, then let's save this
%somewhere special
if str2double(info.include)
    %run the animal parser
    cellInfo = animalParser;
    info.animalShort = info.animal(1:5); %remove the date
    
    %do we have multiple cells for this animal?
    allCells = cellfun(@(x) find(contains(x, info.animalShort)),cellInfo.animal,'UniformOutput',false); %find all cells in the sheet
    allCellNrs = unique(cell2mat(cellInfo.cellNr(find(~cellfun('isempty', allCells))))); %get the cell numbers of those
    if length(allCellNrs) > 1 %if we have multiple ones, get the unique names
       info.allCellNames = unique(cellInfo.cellName(find(~cellfun('isempty', allCells))));
       
       %now let's check if we are looking at a cell and will need this file
       %multiple time or if we are looking at spines and just need to know
       %the correct folderName
       if info.isCell
           info.copyToOtherCells = 1; %we need to copy all the data later to the other folders, but for now let's save it in the first one
           info.cellName = info.allCellNames{allCellNrs(1)};
       else
           %we still need to know the cellName to copy the data correctly, so search for the cellROINr and get the name for it           
           %find all exp of this animal
           thisAnimalExp = cellfun(@(x) find(contains(x, char(info.animalShort))),expInfo.animal,'UniformOutput',false);
           %get all unique cellROIs and sort them
           cellROIs = sort(unique(expInfo.cellROIInd(find(~cellfun('isempty', thisAnimalExp)))));
           %make sure there is no zero in it, which would be the first
           %value as it is sorted
           if any(cellROIs == 0)
               cellROIs = cellROIs(2:end);
           end
           %at which position does the spine exp ROI meet the cellROIs?
           thisSpineCellROI = find(cellROIs == expInfo.cellROIInd(i));
           %take the name of this cellROI
           info.cellName = info.allCellNames{thisSpineCellROI};
           
           %no need to copy afterwards
           info.copyToOtherCells = 0;
       end
    else
        info.cellName = []; %no need to make an extra folder
        info.copyToOtherCells = 0; %no need to copy afterwards 
    end
    
    
    info.saveDir = ['Z:\Juliane\InputAnalysis\' char(info.animalShort) filesep info.cellName filesep 'A - 2p Imaging\' char(info.name) filesep];
    info.saveDir2 = [info.drive 'ImageAnalysis\' char(info.animal) filesep char(info.name) filesep];
    
    %let's save those additional saveDirs for copying data later
    if info.copyToOtherCells
        for i = 1:length(allCellNrs)
            info.addSaveDirs{i} = ['Z:\Juliane\InputAnalysis\' char(info.animalShort) filesep info.allCellNames{allCellNrs(i)} filesep 'A - 2p Imaging\' char(info.name) filesep];
        end
    end
    
else
    info.saveDir = [info.drive 'ImageAnalysis\' char(info.animal) filesep char(info.name) filesep];
    info.saveDir2 = info.saveDir;
    info.copyToOtherCells = 0;
end

%% Step 2: Add the TwoP & Sp2 dir, zoom and make the saveDirs

info.TwoPDir = [info.drive '2P_data\' char(info.animal) filesep char(info.name) '\Registered\'];
info.Sp2Dir = [info.drive 'Spike2Data\' char(info.animal) filesep char(info.sp2) filesep];

info.zoom = getZoomSpines([info.drive '2P_data\' char(info.animal) filesep char(info.name)]);

info.saveDirROI = [info.saveDir 'ROIs\'];
if ~exist(info.saveDir, 'dir')% make new file directory 
    mkdir(info.saveDir); 
end
if ~exist(info.saveDirROI, 'dir')% make new file directory 
    mkdir(info.saveDirROI); 
end

%% Step 3: Save a copy
save([info.saveDir 'expInfo.mat'], 'info', '-mat')
