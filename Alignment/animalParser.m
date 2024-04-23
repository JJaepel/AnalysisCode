  function cellInfo = animalParser
%Reads in the information about the different cells, like animal ID, slice
%number, number of dendrites, input type and which dendrites are apical and
%saves all of that in an array for further computations


%% 1.) Read in file
filePath =  'Z:\Juliane\InputAnalysis\';
file = 'STED animals.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'cells');

%% 2.) Define the readout columns
animalCol = find(contains(xls_txt(1,:),'Animal'),1);
cellNrCol = find(contains(xls_txt(1,:),'CellNr'),1);
cellNameCol = find(contains(xls_txt(1,:),'CellName'),1);
sliceCol = find(contains(xls_txt(1,:),'Slice'),1);
dendriteNrCol = find(contains(xls_txt(1,:),'DendriteNr'),1);
apicalCol = find(contains(xls_txt(1,:),'ApicalDendrites'),1);
inputCol = find(contains(xls_txt(1,:),'InputType'),1);
expCol = find(contains(xls_txt(1,:),'SomaExp'),1);
flagCol = find(contains(xls_txt(1,:),'Flag'),1);
funcFinCol = find(contains(xls_txt(1,:),'FunctFinished'),1);
InputFinCol = find(contains(xls_txt(1,:),'Inputs Finished'),1);

%% 3.) Define the structure and read in the info

cellInfo = struct;
%go through all the rows & save the info for each cell
for i = 2:size(xls_txt,1)
    cellInfo.animal{i-1} = cell2mat(xls_all(i,animalCol)); %animal
    
    cellInfo.cellNr{i-1} = cell2mat(xls_all(i,cellNrCol)); %cellNr
    %if that is empty, remove the nan and replace it with empty
    if isnan(cellInfo.cellNr{i-1})
        cellInfo.cellNr{i-1} = [];
    end

    cellInfo.cellName{i-1}= cell2mat(xls_all(i,cellNameCol)); %cellName
    if isnan(cellInfo.cellName{i-1})
        cellInfo.cellName{i-1} = [];
    end
    
    cellInfo.slice{i-1}= cell2mat(xls_all(i,sliceCol)); %slice
    if isnan(cellInfo.slice{i-1})
        cellInfo.slice{i-1} = [];
    end
    
    cellInfo.dendriteNr{i-1}= cell2mat(xls_all(i,dendriteNrCol)); %number of dendrites
    
    cellInfo.apicalDend{i-1} = cell2mat(xls_all(i,apicalCol)); %apical dendrites
    if isnan(cellInfo.apicalDend{i-1})
        cellInfo.apicalDend{i-1} = [];
    end
    %if there is content, make sure to change its format
    if ~isempty(cellInfo.apicalDend{i-1}) %if there is content, 
        %if it is longer than 2, it is multiple dendrites, so split them
        if size(cellInfo.apicalDend{i-1},2) > 2
            %split it by the comma, read in the individual ones & save as
            %vector in the cell
            numbers = strsplit(cellInfo.apicalDend{i-1}, ', ');
            vector = str2double(numbers);
            cellInfo.apicalDend{i-1} = vector;
        end
    end

    cellInfo.inputType{i-1}= cell2mat(xls_all(i,inputCol)); %inputType
    
    cellInfo.somaExpNr{i-1}= cell2mat(xls_all(i,expCol)); %somaExpNr
    
    cellInfo.flag{i-1} = cell2mat(xls_all(i, flagCol)); %flag
    
    cellInfo.funcFinished{i-1}= cell2mat(xls_all(i,funcFinCol)); %finished functional analysis
    
    cellInfo.InputFinished{i-1}= cell2mat(xls_all(i,InputFinCol)); %finished functional analysis

end

