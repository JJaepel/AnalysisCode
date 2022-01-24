function [fanoFactor] = computeFanoFactor(responses)
% Calculates fanoFactor based on the responses sorted by stimulus
%
% Input:
% - numStims: number of stimulations
% - responses: a numStims x Trials  array containing the response
%
% Ouput:
% - fanoFactor at prefered stimulus
%
    
    %look for preferred stimulus
    meanResponses = mean(responses,2);
    [~, prefIdx]= max(meanResponses);
    
    %for preferred stimulus calculate the variance
    stimResponses = responses(prefIdx, :);
    
    %caclulate fanoFactor
    fanoFactor = nanvar(stimResponses)/nanmean(stimResponses);       
end
