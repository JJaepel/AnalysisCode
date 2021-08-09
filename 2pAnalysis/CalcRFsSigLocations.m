function [thzRFtempONci]=CalcRFsSigLocations(ONdF,metadata, threshold, SamplingFreq, DisplayWindow)

%calculate parameters of WindowSize averaging
StimDur_corr=round(SamplingFreq * metadata.StimParams.stimDuration);

%Interpolates the RF grids such that it consists of 1*1 degree patches
[X, Y] = meshgrid(metadata.StimParams.minAzi:metadata.StimParams.stimSize(1):metadata.StimParams.maxAzi,metadata.StimParams.minElev:metadata.StimParams.stimSize(2):metadata.StimParams.maxElev);
Y = flip(Y);
[Xi, Yi] = meshgrid(metadata.StimParams.minAzi:1:metadata.StimParams.maxAzi, metadata.StimParams.minElev:1: metadata.StimParams.maxElev);
Yi = flip(Yi);
if isempty(X)
    [X, Y] = meshgrid(metadata.StimParams.minAzi:metadata.StimParams.stimSize(1):metadata.StimParams.maxAzi,metadata.StimParams.maxElev:metadata.StimParams.stimSize(2):metadata.StimParams.minElev);
    if min(metadata.StimParams.stimSize) < 2
        [Xi, Yi] = meshgrid(metadata.StimParams.minAzi:.5:metadata.StimParams.maxAzi, metadata.StimParams.maxElev:1: metadata.StimParams.minElev);
    else
        [Xi, Yi] = meshgrid(metadata.StimParams.minAzi:1:metadata.StimParams.maxAzi, metadata.StimParams.maxElev:1: metadata.StimParams.minElev);
    end
end

for knd = 1:(length(DisplayWindow)-StimDur_corr)
    RFtempON(:,:,knd) = mean(ONdF(knd:(StimDur_corr-1)+knd,:,:),1);
    RFtempON(:,:,knd)=smooth2a(RFtempON(:,:,knd),0,0);
    
    % interpolate
    RFtempONci(:,:,knd)= interp2(X,Y,RFtempON(:,:,knd),Xi,Yi,'cubic');
    RFtempONcitemp=RFtempONci(:,:,knd);
    
    % convert into z scores
    zRFtempONci(:,:,knd)=(RFtempONcitemp-nanmean(RFtempONcitemp(1:end)))/nanstd(RFtempONcitemp(1:end));
        
    % apply threshold equivalent to the sigma (of a normal distribution)
    thzRFtempONci(:,:,knd)=zRFtempONci(:,:,knd);
    thzRFtempONci(thzRFtempONci<threshold)=0;
end
