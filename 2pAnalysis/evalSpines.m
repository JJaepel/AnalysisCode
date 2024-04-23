function ce = evalSpines(ce, info)

% Evaluate a recorded spine for use in analyses based on set criterion.
% Currently setting 'good' structure field to '1' if passed criterion. 

% Criteria: 
% (1) >10% df/f for a stimulus
% (2) >1 SNR for stimulus
% (3) <0.4 residual correlation with dendrite
%
% Inputs:
% - ce: structure containing information about ROI, each row is one roi,
% containing anatomical information, raw & dff traces, chopped stimulus
% trace, cyc traces & traces after dendritic substraction
% - info: all information about the experiment, such as whether it was a
% cell and where to save things
%
% Outputs:
% - ce: updated with info whether it is a good spine or not
% - ROIsAna.mat: File with the stored ce

%% Step 1: Go through all ROIs after each 
for cc = 1:length(ce)
    %if it is a spine
    if ce(cc).spine
        %% Step 2: Calculate necessary parameters
        %Check (1) >10% df/f for a stimulus:
        resp = ce(cc).meanResp; % get mean Resp 
        
        %Check (2) >1 SNR for stimulus
        resperr = ce(cc).mpeakErr'; %get error
        blank = squeeze(ce(cc).cycRes(end,:,:));
        
        spont = mean(blank(:)); %calculate spontaneous signal and error for SNR calculation
        spont(spont < 0) = 0;
        sponterr = std(mean(blank,2))./sqrt(size(ce(cc).cycRes,2));
        
        snr = (resp - spont) ./ (resperr + sponterr);  %calculate SNR
        
        %Check (3) <0.4 residual correlation with dendrite
        residualcorr = computeStimCorrSpDend(cc, ce); % compute correlation with dendrite 
        
        %% Step 3: If all conditions are met, set good to 1, else 0
        if (residualcorr<.4) && sum((resp>0.10) & (snr>1))>0            
            ce(cc).good = 1;
        else            
            ce(cc).good = 0;
        end
    else
        ce(cc).good = nan;
    end
end

%% Step 4: Save in mat file
save([info.saveDir 'ROIsAna.mat'], 'ce', '-mat') 