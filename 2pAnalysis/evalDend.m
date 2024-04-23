function good = evalDend(cyc)

% PURPOSE:
% Evaluate a recorded dendrite for use in analyses based on set criterion.
% Currently setting 'good' structure field to '1' if passed criterion. 
%
% PARAMS:
% - cyc: reponse trace to stimulus, dimensionality stims (incl. blank) x
% trials x timeWindow
%
% Criteria: 
% (1) >10% df/f for a stimulus
% (2) >1 SNR for stimulus

%calc response and respError
resp = mean(cyc(1:end-1,:,:),3);
respMean = mean(resp,2);
respErr = std(resp,0,2)/sqrt(size(resp,1));

%errResp = std(peak(1:n,:),0,2)/sqrt(size(peak(1:n,:),2)); %SEM

%calculate spontaneous signal and error for SNR calculation
blank = squeeze(cyc(end,:,:));
spontMean = mean(blank(:));
spontMean(spontMean < 0) = 0;
spontErr = std(mean(blank,2))./sqrt(size(cyc,2));

%calculate SNR
SNR = (respMean - spontMean) ./ (respErr + spontErr);

%decide if good
if sum((respMean > 0.10) & (SNR > 1)) > 0
    good = 1;
else
    good = 0;
end 