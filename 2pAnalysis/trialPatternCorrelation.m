function [corrCoeff, corrTrialsMatched, corrTrialsOrtho] = trialPatternCorrelation(analysis, metadata, field, saveDirectory)
% Computs correlation of patterns of cells across trials independent of the
% stimulus
% Input:
% - metadata: structure containing StimParams information
% - analysis: structure containing for each ROI (analysis.roi) the
% stim response trace and whether they are significant responders
% - field: which kind of trace to use
%
% Ouput:
% - corrCoeff: Correlation data

%% for each trial, independent of stimulus, get the cell pattern
%start with all trials of one stimulus first, then the next stimulus, etc.,
%eg. stim1 trial 1, stim1 trial 2,..., stim 1 last trial, stim2 trial 1,...

numTrials = (metadata.StimParams.uniqStims-1) * metadata.StimParams.numTrials;
stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);

respCells = linspace(1,length(analysis.dff.roi), length(analysis.dff.roi)); %taking all cells
%respCells = find([analysis.(field).roi.isResponseSignificant] == 1); %taking only significant responders
numCells = length(respCells); 

cellPattern = zeros(numCells, numTrials);

%make a vector containing the stimulus identity for all trials
stimSeq = linspace(1, metadata.StimParams.uniqStims-1, metadata.StimParams.uniqStims-1);
stimIDs =  repelem(stimSeq,metadata.StimParams.numTrials);

for cell =1:length(respCells)
    cellNr = respCells(cell);
    Response = analysis.(field).roi(cellNr).stimResponseTrace(1:end-1, :, :);
    Response = median(Response(:,:,stimWindow),3); %average over the stimWindow
    Response = reshape(Response',[(metadata.StimParams.uniqStims-1) * metadata.StimParams.numTrials,1]); %reshape to have it all stims in one vector as desired
    
%     %zscore response
%     mu = nanmean(Response(:));
%     sd = nanstd(Response(:));
%     zScoreResponse = (Response - mu)/sd;
    
    %normalize Response
    xmin=min(Response(:));
    xmax=max(Response(:));
    normResponse = (Response-xmin)/(xmax-xmin);
    
    %save
    cellPattern(cell, :) = normResponse;
end

%order cells py orientation preference
prefOri = [analysis.(field).roi(respCells).preferredOrientation]; 
[OriSorted, OriOrder] = sort(prefOri);
cellPatternSorted = cellPattern(OriOrder',:);

figure
imagesc(cellPatternSorted)
c = gray;
c = flipud(c);
colormap(c)
colorbar
xlabel('Trial number')
ylabel('Cell number')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'TrialResponses.png'))

%% for each trial pair, do the correlation
corrCoeff = zeros(numTrials, numTrials);
for i = 1:numTrials
    for j = 1:numTrials
        corrCoeff(i,j) = corr(cellPattern(:,i),cellPattern(:,j));
    end
end

%% visualize the pattern in a corelogramm
figure
imagesc(corrCoeff)
colormap(redblue)
caxis([-1 1])
colorbar
set(gca, 'YDir', 'normal')
xlabel('Trial Number')
ylabel('Trial number')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'CorrelationTrialPatterns.png'))

%% calculate indices for matched pairs and orthogonal pairs
orthTrialShift = (metadata.StimParams.uniqStims-1)/4;
corrTrialsMatched = zeros((metadata.StimParams.uniqStims-1)*metadata.StimParams.numTrials,metadata.StimParams.numTrials);
corrTrialsOrthoLeft= zeros((metadata.StimParams.uniqStims-1)*metadata.StimParams.numTrials,metadata.StimParams.numTrials);
corrTrialsOrthoRight= zeros((metadata.StimParams.uniqStims-1)*metadata.StimParams.numTrials,metadata.StimParams.numTrials);

for stimNr = 1:(metadata.StimParams.uniqStims-1)
    %find all matched trials and do the correlation
    matchedTrials = find(stimIDs == stimNr);
    for i = 1:metadata.StimParams.numTrials
        for j = 1:metadata.StimParams.numTrials
            corrTempMatched(i,j) = corr(cellPattern(:,matchedTrials(i)),cellPattern(:,matchedTrials(j)));
        end
    end
    %add it to one big structure
    corrTrialsMatched(1+(stimNr-1)*metadata.StimParams.numTrials:stimNr*metadata.StimParams.numTrials,:)=corrTempMatched;

    %find all orthogonal trials
    orthPairLeft = stimNr + orthTrialShift;
    orthPairRight = stimNr - orthTrialShift;
    if orthPairLeft > (metadata.StimParams.uniqStims-1)
        orthPairLeft = orthPairLeft-(metadata.StimParams.uniqStims-1);
    end
    if orthPairRight <= 0
        orthPairRight = orthPairRight+(metadata.StimParams.uniqStims-1);
    end
    orthoTrialsLeft = find(stimIDs == orthPairLeft);
    orthoTrialsRight = find(stimIDs == orthPairRight);
    %do correlation for left orthogonal trials
    for i = 1:metadata.StimParams.numTrials
        for j = 1:metadata.StimParams.numTrials
            corrTempLeft(i,j) = corr(cellPattern(:,matchedTrials(i)),cellPattern(:,orthoTrialsLeft(j)));
        end
    end
    
    %do correlation for right orthogonal trials
    for i = 1:metadata.StimParams.numTrials
        for j = 1:metadata.StimParams.numTrials
            corrTempRight(i,j) = corr(cellPattern(:,matchedTrials(i)),cellPattern(:,orthoTrialsRight(j)));
        end
    end
    
    %add it to one big structure
    corrTrialsOrthoLeft(1+(stimNr-1)*metadata.StimParams.numTrials:stimNr*metadata.StimParams.numTrials,:)=corrTempLeft;
    corrTrialsOrthoRight(1+(stimNr-1)*metadata.StimParams.numTrials:stimNr*metadata.StimParams.numTrials,:)=corrTempRight;
end
%add all ortho values together
corrTrialsOrtho = [corrTrialsOrthoLeft; corrTrialsOrthoRight];

figure
subplot(1,2,1)
boxplot(corrTrialsMatched(:), 'labels', 'Matched')
ylim([-0.5 1])
subplot(1,2,2)
boxplot(corrTrialsOrtho(:), 'labels', 'Orthogonal')
ylim([-0.5 1])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'OrthoVsMatchedTrialPatterns.png'))
