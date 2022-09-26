function SFVar = computeSFVar(SfResponses,Sfs)
% 
% returns the weighted spatial Frequency variance
%
% Input:
% - SfResponses: Array with mean response at preferred direction
% - Sfs: spatial frequencies used
%
% Ouput:
% - SFVar: weighted spatial Frequency variance
%
%

%what is the prefered SF and its value (prefSfResp)
[prefSfResp, prefInd]=max(SfResponses);
prefSf = Sfs(prefInd);

numSf = length(Sfs);
weight = zeros(numSf,1);
weigthedResp = zeros(numSf,1);

%for each sf, calculate:
for sf = 1:numSf
    %the weight that goes into it - the further from the prefered sf, the
    %lower the weight
    if Sfs(sf) == prefSf
        weight(sf) = 0;
    else
        weight(sf) = 1/abs((log(prefSf/Sfs(sf))));
    end
end
%normalize weights to make it more comparable
weight = weight/sum(weight);
 
for sf = 1:numSf
    %the weighted responses as a difference between that responses and the
    %maximum response
    weigthedResp(sf) = weight(sf)*(prefSfResp-SfResponses(sf))/SfResponses(sf);
end

%sum the weighted responses 
%-> the bigger %the difference, the higher the factor
SFVar = sum(weigthedResp);
