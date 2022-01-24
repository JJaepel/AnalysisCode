function [VI] = computeVI(response, numBootsTrap, type)
% Calculates fanoFactor based on the responses sorted by stimulus
%
% Input:
% - response: a numStims x Trials  array containing the response
% traces
% - type: whether to compare all pairs or just sequential
%
% Ouput:
% - Varibility index 
%
    if nargin < 3
        type = 'seq';
    end
    
    if nargin < 2
        numBootsTrap = 250;
    end
    
    numStims = size(response,1);
    numTrials = size(response,2);
    
    %resample trial responses within stimuli to recompute resampled tuning
    %curves
    resampledCurves = zeros(numStims, numBootsTrap);
    for i = 1:numBootsTrap
        grabbedTrials = randi([1 numTrials],1,numStims);
        for stim = 1:numStims
            resampledCurves(stim,i) = response(stim,grabbedTrials(stim));
        end
    end
    
    switch type
        case 'all'
            %calculate Pearson's correlation coefficient for all pairs
            corrCoeff = zeros(numBootsTrap, numBootsTrap-1);
            for i = 1:numBootsTrap
                for j= 2:numBootsTrap-1
                    corrCoeff(i,j) = corr(resampledCurves(:,i),resampledCurves(:,j));
                end
            end
            
            %caclulate variability index
            normalizationFactor = (numBootsTrap * numBootsTrap-1)/2; %divided by two as the order is not important and each pair represented twice

        case 'seq'
           %calculate Pearson's correlation coefficient for all pairs
           corrCoeff = zeros(numBootsTrap-1,1);
           for n = 2:numBootsTrap
               corrCoeff(n) = corr(resampledCurves(:,n),resampledCurves(:,n-1));
           end
           
           %caclulate variability index
           normalizationFactor = numBootsTrap-1;
    end
    
    %caclulate variability index
    VI = 1- sum(corrCoeff(:))/normalizationFactor;     
end
