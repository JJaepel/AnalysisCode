function combineSWCfiles
%Combine swc files that might be separated due to not being able to fully
%trace it

%% Step 1: Select files
baseDir = 'Z:/Juliane/InputAnalysis/';

[startFile, path] = uigetfile('.swc', 'Select start file', baseDir);
endfile = uigetfile('.swc', 'Select end file', path);

%% Step 2: Read out file information and combine it to one

%load start file
b = loadfilelist([path startFile]);
headerInfo = b(1:6);
traceStart = cell2mat(cellfun(@str2num, b(4:end), 'UniformOutput', false)');

%load end file
c = loadfilelist([path endfile]);
traceEnd = cell2mat(cellfun(@str2num, c(4:end), 'UniformOutput', false)');

%check the gap between the two files
gapX = abs(traceStart(end,3)-traceEnd(1,3));
gapY = abs(traceStart(end,4)-traceEnd(1,4));
if gapX > 1.5 || gapY > 1.5
    %whaat is the size of the biggest gap
    maxGap = max([gapX,gapY]);
    
    %each gap needs to be smaller than 1.5, so divide that to get amount of
    %needed steps
    neededSteps = ceil(maxGap/1.5);
    
    %make a new array to add to the trace
    gapInfo = zeros(neededSteps,7);
    gapInfo(:,1) = linspace(length(traceStart)+1, length(traceStart)+7,neededSteps); %indices
    xBridge = linspace(traceStart(end,3),traceEnd(1,3), neededSteps+2); %x coordinates
    gapInfo(:,3) = xBridge(2:end-1);
    yBridge = linspace(traceStart(end,4),traceEnd(1,4), neededSteps+2); %y coordinates
    gapInfo(:,4) = yBridge(2:end-1);
    zBridge = linspace(traceStart(end,5),traceEnd(1,5), neededSteps+2); %z slices
    gapInfo(:,5) = round(zBridge(2:end-1));
    gapInfo(:,7) = linspace(length(traceStart), length(traceStart)+6,neededSteps); %parent
    
    %add the gap to the traceStar
    traceStart = [traceStart; gapInfo];
end

%modify the second file to match the first one
traceEnd(:,1) = traceEnd(:,1) + ones(length(traceEnd),1)*length(traceStart);
traceEnd(:,7) = traceEnd(:,7) + ones(length(traceEnd),1)*length(traceStart);
traceEnd(1,7) = traceEnd(1,7)+ 1;

%combine them & convert to cell
fullTrace = [traceStart; traceEnd];

%% Step 3: Save file as swc

%determine fileName as combination of both
fileName = [path startFile(1:end-4) '_' endfile(end-6:end-4) '.swc'];

%open it
fileID = fopen(fileName, 'w');

%write header data
formatSpecHeader = '%s\n';
for h = 1:length(headerInfo)
    fprintf(fileID, formatSpecHeader, headerInfo{h});
end

%write data
formatSpec = '%s %s %s %s %s %s %s\n';
for nodeIndex = 1:length(fullTrace)
     output{nodeIndex, 1} = num2str(fullTrace(nodeIndex, 1));   % index
     output{nodeIndex, 2} = num2str(fullTrace(nodeIndex, 2));   % type
     output{nodeIndex, 3} = num2str(fullTrace(nodeIndex, 3));   % x
     output{nodeIndex, 4} = num2str(fullTrace(nodeIndex, 4));   % y
     output{nodeIndex, 5} = num2str(fullTrace(nodeIndex, 5));   % z
     output{nodeIndex, 6} = num2str(fullTrace(nodeIndex, 6));   % radius
     output{nodeIndex, 7} = num2str(fullTrace(nodeIndex, 7));   % parent
     fprintf(fileID, formatSpec, output{nodeIndex, :});
end

%close the file
fclose(fileID);