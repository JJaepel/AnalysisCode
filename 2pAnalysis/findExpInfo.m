function expInfo = findExpInfo(filePath,varargin)
% parses the xls file for the correct experiment information to forward to
% the analysis scripts
%
% Input:
% - filePath: path and filename of corresponding xlsx file
% - (sheetName): name of corresponding sheet
%
% Output:
% - expInfo: structure with the experiment information ordered by row in
% the excel sheet

%% Step 1: Read in experiment file
if nargin > 1
    sheet = varargin{1}; %add sheetName
else 
    sheetNames = sheetnames(filePath);
    sheet = sheetNames(1); %otherwise assume first sheet
end

[~, xls_txt, xls_all]=xlsread(filePath, sheet);

%% Step 2: Find the corresponding columns

%same in all experiments
animalCol = find(contains(xls_txt(1,:),'animal'),1);
expCol = find(contains(xls_txt(1,:),'expNumber'),1);
spk2Col = find(contains(xls_txt(1,:),'spk2'),1);
volCol = find(contains(xls_txt(1,:),'volume'),1);
nameCol = find(contains(xls_txt(1,:),'name'),1);
runCol = find(contains(xls_txt(1,:),'Run'),1);
flagCol = find(contains(xls_txt(1,:),'Flag'),1); 
regCol = find(contains(xls_txt(1,:),'region'), 1);
ageCol = find(contains(xls_txt(1,:), 'age'),1);
try
    EOCol = find(contains(xls_txt(1,:), 'EO'),1);
catch 
    EOCol = NaN;
end

%special columns based on sheet
switch sheet
    case 'Retinotopy'
        specCol = find(contains(xls_txt(1,:),'special'),1);
        comCol = find(contains(xls_txt(1,:), 'comments'),1);
        stimCol = find(contains(xls_txt(1,:),'stimulus'),1);
    case 'bimodal'
        cellROICol = find(contains(xls_txt(1,:), 'cell ROI'),1);
        dendROICol = find(contains(xls_txt(1,:), 'dendrite ROI'),1);
        includeCol = find(contains(xls_txt(1,:), 'Include'),1);
        DepthCol = find(contains(xls_txt(1,:), 'Depth'),1);
        denTypeCol = find(contains(xls_txt(1,:), 'denType'),1);
        modCol = find(contains(xls_txt(1,:), 'mod'),1);
        altNameCol = find(contains(xls_txt(1,:), 'altName'),1);   
end

%% Step 3: Define the structure and read in the info

expInfo = struct;
%go through all the rows & save the info for each experiment

for i = 2:size(xls_txt,1)
    %columns there for all xls files
    expInfo.animal{i-1} = xls_all(i,animalCol); %animal
    
    exp = xls_all(i,expCol); %expNR
    if iscell(exp)
        exp = cell2mat(exp);
        if ischar(exp); exp = str2double(exp); end
    end
    if exp > 9
        if exp > 99
            expInfo.exp_id{i-1} = ['t00' num2str(exp)];
        else
            expInfo.exp_id{i-1} = ['t000' num2str(exp)];
        end
    else
        expInfo.exp_id{i-1} = ['t0000' num2str(exp)];
    end
    expInfo.exp_series{i-1} = ['tseries_' num2str(exp)];
    
    sp2 = xls_all(i,spk2Col); %sp2id
    if iscell(sp2)
        sp2 = cell2mat(sp2);
        if ischar(sp2); sp2 = str2double(sp2); end
    end
    if sp2 > 9
        if sp2>99
            expInfo.sp2_id{i-1} = ['t00' num2str(sp2)];
        else
            expInfo.sp2_id{i-1} = ['t000' num2str(sp2)];
        end
    else
        expInfo.sp2_id{i-1} = ['t0000' num2str(sp2)];
    end
    
    vol = xls_all(i,volCol); %volume
    if contains(vol{1}, 'yes')
        expInfo.vol{i-1} = 1;
    elseif contains(vol{1}, 'no')
        expInfo.vol{i-1} = 0;
    else
        expInfo.vol{i-1} = [];
    end
    
    expInfo.name{i-1} = xls_all(i,nameCol); %name
    
    runInd = xls_all{i,runCol}; %should you analyze this one? 
    if ischar(runInd); runInd = str2double(runInd); end
    expInfo.run(i-1) = runInd; 
    
    flagInd = xls_all{i,flagCol}; %is there something bad about this one?
    if ischar(flagInd); flagInd = str2double(flagInd); end
    expInfo.flag(i-1) = flagInd;
    
    regInd = xls_all{i, regCol}; %region: A19 vs. V1, cell vs. dendrite
    if ischar(regInd)
        regInd = str2double(regInd);
        if isnan(regInd); regInd = xls_all{i, regCol}; end
    end
    expInfo.region{i-1} = regInd;
    
    expInfo.age{i-1} = xls_all(i,ageCol); %age
    
    switch sheet
        case 'Retinotopy'
            specInd = xls_all{i,specCol}; %to account for some variations
            if ischar(specInd); specInd = str2double(specInd); end
            expInfo.special(i-1) = specInd;
            
            expInfo.stimulus{i-1} = xls_all{i, stimCol}; %which stim?
            
            expInfo.comments{i-1} = xls_all{i, comCol}; %sparseNoise file?
        case 'bimodal'
            cellROIInd = xls_all{i,cellROICol}; %what is the soma ROI?
            if ischar(cellROIInd); cellROIInd = str2double(cellROIInd); end
            expInfo.cellROIInd(i-1) = cellROIInd;
            
            dendROIInd = xls_all{i,dendROICol}; %what is the dendrite Nr?
            if ischar(dendROIInd); dendROIInd = str2double(dendROIInd); end
            expInfo.dendROIInd(i-1) = dendROIInd;
            
            depthInd = xls_all{i,DepthCol}; %what is the depth?
            if ischar(depthInd); depthInd = str2double(depthInd); end
            expInfo.Depth(i-1) = depthInd;
            
            expInfo.denType{i-1} = xls_all(i,denTypeCol); %apical or basal?
            
            expInfo.Modality{i-1} = xls_all(i,modCol); %STED or EM?
            
            expInfo.altName{i-1} = xls_all(i,altNameCol); %EM name?
            
            includeInd = xls_all{i,includeCol}; %include in inputanalysis?
            if ischar(includeInd); includeInd = str2double(includeInd); end
            expInfo.include(i-1) = includeInd;
    end
    
    if ~isnan(EOCol); expInfo.EO{i-1} = xls_all{i, EOCol}; end
end