function ce = computeDffSpinesNew(info, ce)

%This function takes the raw traces and converts it to the dff trace.
%
% Inputs:
% - info: contains information about the experiment, here important for
% knowing where to save the data
% - ce: structure with the info about the ROIs, so far containing
% anatomical data and raw trace
%
% Output:
% - ce: now updated to also include the dff trace as 

%Go through all ROIs after each other
for cc = 1:length(ce)
    %% Step 1: Get the raw trace from the structure
    raw = ce(cc).raw;
    
    %% Step 2: Calculate dff as 99pt medfilter after resample (1/4)
    dff = filterBaseline_dFcomp(resample(raw,1,4)); 
    
    %% Step 3: Remove nans and inifity numbers by setting them to 0
    dff(isnan(dff)) = 0;  
    dff(isinf(dff)) = 0;
    
    %% Step 4: Save to ce as dff trace
    ce(cc).dff = dff;
end

%% Step 5: Save all of them as ROIs.mat
save([info.saveDir 'ROIs.mat'], 'ce', '-mat') 