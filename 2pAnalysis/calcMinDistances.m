function PopAnalysis = calcMinDistances(PopAnalysis, analysisParams)
% Calculates homeogeneity index for a field of view for a variety of radiuses 
%
% Input:
%
%
% Ouput:
% - meanMinDistRe
%
for ROI = 1:length(PopAnalysis.(analysisParams.field).distROIs)
    validROIs = [PopAnalysis.(analysisParams.field).distROIs(ROI,:)> 0];
    minDistROIs(ROI) = min(PopAnalysis.(analysisParams.field).distROIs(ROI,validROIs));
end

PopAnalysis.(analysisParams.field).minDistROIs = minDistROIs;
PopAnalysis.(analysisParams.field).meanMinDistROIs = mean(minDistROIs);
PopAnalysis.(analysisParams.field).semMinDistROIs = std(minDistROIs)/length(minDistROIs);

for respROI = 1:length(PopAnalysis.(analysisParams.field).distRespROIs)
    validRespROIs = [PopAnalysis.(analysisParams.field).distRespROIs(respROI,:)> 0];
    minDistRespROIs(respROI) = min(PopAnalysis.(analysisParams.field).distRespROIs(respROI,validRespROIs));
end

PopAnalysis.(analysisParams.field).minDistRespROIs = minDistRespROIs;
PopAnalysis.(analysisParams.field).meanMinDistRespROIs = mean(minDistRespROIs);
PopAnalysis.(analysisParams.field).semMinDistRespROIs = std(minDistRespROIs)/length(minDistRespROIs);