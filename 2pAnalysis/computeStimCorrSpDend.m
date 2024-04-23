function [newCorr] = computeStimCorrSpDend(cc, ce)
% Computes the correlation between a spine and the dendrite that it is on
%
% Input:
% - cc: Spine number
% - ce: structure containing information about ROI, each row is one roi,
% containing anatomical information, raw & dff traces, chopped stimulus
% trace, cyc traces & traces after dendritic substraction
%
% Output:
% - newCorr: correlation between spine and its dendrite

newCorr = [];
% for each stim and trial
for stim = 1:size(ce(cc).cycRes,1)-1
    for trial = 1:size(ce(cc).cycRes,2)
        %% Step 1: Get the spine signal
        cycSp = squeeze(ce(cc).cycRes(stim,:,:));
        cycSp(cycSp==0) = NaN;
        
        %% Step 2: Find the first dendrite that has a higher number than the spine number = its dendrite
        dendrites = find([ce.dendrite]==1);
        dendrNr = find(dendrites > cc);
        cycDen = squeeze(ce(dendrites(dendrNr(1))).cyc(stim,:,:));
        
        %% Step 3: Save the correlation between the two in a variable
        newCorr = [newCorr;corr(cycSp(:),cycDen(:),'rows','pairwise')];
    end
end

%% Step 4: Get the mean correlation
newCorr = nanmean(newCorr);