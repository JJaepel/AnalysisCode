function trialCorr = trialToTrialCorrelation(RespSpine1, RespSpine2)

% Calculates trial-to-trial correlation for two spines
%
% Input:
% - RespSpine1: response of first spine, format: stims x trials x full
% response
% - RespSpine2: response of second spine format: stims x trials x full
% response

% Steps:
% 1.) Reshape responses
% 2.) Do correlation for each trial
% 3.) Average
% 
% Output: 
% - trialCorr: correlation between trials

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: March, 2024

%% Step 1: Reshape responses

%get peaks for each trial and stim
for ii = 1:size(RespSpine1,1)
    for jj = 1:size(RespSpine1,2)
        r = squeeze(RespSpine1(ii,jj,:));
        ff = fft(r)./length(r);
        peakSpine1(ii,jj) = ff(1)+2*abs(ff(2));
        if ii==size(RespSpine1,1)
            peakSpine1(ii,jj) = ff(1);
        end
    end
end

for ii = 1:size(RespSpine2,1)
    for jj = 1:size(RespSpine2,2)
        r = squeeze(RespSpine2(ii,jj,:));
        ff = fft(r)./length(r);
        peakSpine2(ii,jj) = ff(1)+2*abs(ff(2));
        if ii==size(RespSpine2,1)
            peakSpine2(ii,jj) = ff(1);
        end
    end
end

%% Step 2: Get correlation for each trial
for trial = 1:size(peakSpine1,2)
    singleTrialCorr(trial) = corr(peakSpine1(:,trial), peakSpine2(:,trial));
end

%% Step 3: Average across trials

trialCorr = mean(singleTrialCorr);