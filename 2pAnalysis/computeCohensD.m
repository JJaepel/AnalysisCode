function [cohensD] = computeCohensD(numStims, responses)
% Calculates cohen's D based on the responses sorted by stimulus
%
% Input:
% - numStims: number of stimulations
% - responses: a numStims x Trials  array containing the response
%
% Ouput:
% - cohensD
%
    
    meanResponses = mean(responses,2);

    [~, prefIdx]= max(meanResponses);
    
    orthIdx=[ (mod(numStims/4 + prefIdx -1, numStims) +1), (mod( prefIdx -1- numStims/4 , numStims) +1)];
        
    prefTrials = responses(prefIdx,:);
    orthTrials = horzcat(responses(orthIdx(1),:), responses(orthIdx(2),:));
    
    s1sq=sum((prefTrials(:)- mean(prefTrials(:))).^2)/(length(prefTrials(:)-1));
    s2sq=sum((orthTrials(:)- mean(orthTrials(:))).^2)/(length(orthTrials(:)-1));
    pooledSD= sqrt((length(prefTrials(:)-1) * s1sq + (length(orthTrials(:)-1) * s2sq) )/ (length(prefTrials(:))+length(orthTrials(:)) -2));
    
    cohensD = (mean(prefTrials,2) - mean(orthTrials,2)) / pooledSD;
       
end
