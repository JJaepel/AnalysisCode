function UpdateInputAnalysis

%% STEP 1: Get all experiments and all necessary data
% need: Animal name, whether multiple cells, and cell name as counted per
% cell
cellinfo = animalParser;
animals = unique(cellinfo.animal);

%if some of the animals contain temp, remove those
flaggedAnimals = cellfun(@(str) contains(str, 'TEMP'), animals);
animals = animals(~flaggedAnimals);

%% STEP 2: Set up a loop to go through all single animal and single cell analysis functions

%single animal
% for f = 1:length(animals)
%     SpineImagingAnalysis([],animals{f})
%     disp(['Currently working on: ' char(animals{f})])
%     spineToCellComputations(animals{f})
%     close all
% end

%single cells
for c = 1:length(cellinfo.cellNr)
%     if contains(cellinfo.animal{c}, 'TEMP')
%         continue
%     end
    
%     disp(['Currently working on: ' char(cellinfo.animal{c})])
%     confocalTracingtoPlot(cellinfo, c)
%     close all
    
%     if ~isempty(cellinfo.cellName{c})
%        AnatomicalAnalysisCellV2(cellinfo.animal{c}, cellinfo.cellName{c})
%     else
%        AnatomicalAnalysisCellV2(cellinfo.animal{c})
%     end
%     close all
% 
%     %let's check if the functional alignment has been done
%     if contains('Yes', cellinfo.funcFinished{c})
%         if ~isempty(cellinfo.cellName{c})
%             FunctionalConfocalRepresentation(cellinfo.animal{c},cellinfo.cellName{c})
%         else
%             FunctionalConfocalRepresentation(cellinfo.animal{c})
%         end
%     end
%     close all
%     
%     %let's check if the input analysis has been done
%     if contains('Yes', cellinfo.InputFinished{c})
%         if ~isempty(cellinfo.cellName{c})
%            AnatomicalInputAnalysis(cellinfo.animal{c},cellinfo.cellName{c})
%            close all
%            MultiModalAnalysisNew(cellinfo.animal{c},cellinfo.cellName{c})
%         else
%            AnatomicalInputAnalysis(cellinfo.animal{c})
%            close all
%            MultiModalAnalysisNew(cellinfo.animal{c})
%         end
%     end
%     close all
end

%% STEP 3: Combine everything
MultiCellAnalysis
MultiCellInputAnalysis
