function Plot2DResponseTrace(metadata, analysis, field, roi,saveDir)

h=figure('units','normalized','outerposition',[0 0 1 1]);
numCon = metadata.StimParams.numCon;

if metadata.StimParams.numSf > 1
    conLabels = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
elseif metada.StimParams.numTf > 1
    conLabels = strsplit(metadata.StimParams.temporalFreq(2:end-1), ',');
end

respMatrix = analysis.(field).roi(roi).avgResponseTrace(1:end-1,:);
TraceLength = size(respMatrix,2);
respTraceMatrix = zeros(numCon,metadata.StimParams.numDirections*TraceLength);
for con=1:numCon
    temp = respMatrix((con-1)*metadata.StimParams.numDirections+1:con*metadata.StimParams.numDirections,:);
    respTraceMatrix(con,:) = reshape(temp',1,metadata.StimParams.numDirections*TraceLength)';
end


imagesc(respTraceMatrix)
hold all
for dir = 1:metadata.StimParams.numDirections
    %make line for stimulus separation
    line([dir*TraceLength dir*TraceLength], [0 4], 'Color', 'black', 'LineWidth', 1)
    %make line for stim on- and offset
    line([(dir-1)*TraceLength+analysis.(field).windowStart (dir-1)*TraceLength+analysis.(field).windowStart], [0 4], 'Color', 'black', 'LineStyle', '--')
    line([(dir-1)*TraceLength+analysis.(field).windowStop (dir-1)*TraceLength+analysis.(field).windowStop], [0 4], 'Color', 'black', 'LineStyle', '--')
end

xlabel('directions in degree')
xticks = 1:TraceLength:TraceLength*metadata.StimParams.numDirections;  
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

title(['ROI ' num2str(roi) ', prefSf: ' analysis.(field).roi(roi).prefSf ', SFSI: ' num2str(analysis.(field).roi(roi).SFSI) ', pref Dir:' num2str(analysis.(field).roi(roi).preferredDirection)])

set(gcf, 'Color', 'w')
saveas(gcf, fullfile(saveDir, ['ROI_Nr_' num2str(roi) '2D_AvgStimTrace.png']))
close gcf
