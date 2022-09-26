function Plot2DResponse(metadata, analysis, field, roi,saveDir)

h=figure('units','normalized','outerposition',[0 0 1 1]);
numCon = metadata.StimParams.numCon;

if metadata.StimParams.numSf > 1
    conLabels = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
elseif metada.StimParams.numTf > 1
    conLabels = strsplit(metadata.StimParams.temporalFreq(2:end-1), ',');
end

respMatrix = analysis.(field).roi(roi).avgStimResponse(1:end-1);
respMatrix = reshape(respMatrix, metadata.StimParams.numDirections,numCon)';

imagesc(respMatrix)
xlabel('directions in degree')
xticks = 1:metadata.StimParams.numDirections;  
xlabels = metadata.StimParams.directions'; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);


yticks = 1:numCon;
set(gca, 'YTick', yticks, 'YTickLabel', conLabels);
if metadata.StimParams.numSf > 1
    ylabel('spatial frequency in cpd');
elseif metada.StimParams.numTf > 1
    ylabel('temporal frequency in Hz');
end

colorbar
colormap(jet)

set(gcf, 'Color', 'w')
saveas(gcf, fullfile(saveDir, ['ROI_Nr_' num2str(roi) '2D_AvgStimResp.png']))
