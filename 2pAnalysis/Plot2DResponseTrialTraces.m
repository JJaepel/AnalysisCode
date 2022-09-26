function Plot2DResponseTrialTraces(metadata, analysis, field, roi,saveDir)

h=figure(2000+roi);
set(h,'units','normalized','outerposition',[0 0 1 1]);
numCon = metadata.StimParams.numCon;

if metadata.StimParams.numSf > 1
    conLabels = metadata.StimParams.spatialFreq;
elseif metadata.StimParams.numTf > 1
    conLabels = metadata.StimParams.temporalFreq;
else
    conLabels = 'all con';
end

respMatrix = analysis.(field).roi(roi).stimResponseTrace(1:end-1,:,:);
numTrials = size(respMatrix,2);
TraceLength = size(respMatrix,3);
respTraceMatrix = zeros(numCon*numTrials,metadata.StimParams.numDirections*TraceLength);

ind = 1;
for con = 1:numCon
    temp = respMatrix((con-1)*metadata.StimParams.numDirections+1:con*metadata.StimParams.numDirections,:,:);
    temp = reshape(temp,numTrials*metadata.StimParams.numDirections, TraceLength);
    for trial=1:numTrials
        tempTrial = temp((trial-1)*metadata.StimParams.numDirections+1:trial*metadata.StimParams.numDirections,:);
        respTraceMatrix(ind,:) = reshape(tempTrial',1,metadata.StimParams.numDirections*TraceLength);
        ind = ind+1;    
    end
end

imagesc(respTraceMatrix)
hold all
for dir = 1:metadata.StimParams.numDirections
    %make line for stimulus separation
    line([dir*TraceLength dir*TraceLength], [0 numCon*numTrials+1], 'Color', 'black', 'LineWidth', 1)
    %make line for stim on- and offset
    line([(dir-1)*TraceLength+analysis.(field).stimStart (dir-1)*TraceLength+analysis.(field).stimStart], [0 numCon*numTrials+1], 'Color', 'black', 'LineStyle', '--')
    line([(dir-1)*TraceLength+analysis.(field).stimStop (dir-1)*TraceLength+analysis.(field).stimStop], [0 numCon*numTrials+1], 'Color', 'black', 'LineStyle', '--')
end
for con = 1:numCon
    line([0 TraceLength*metadata.StimParams.numDirections],[(con-1)*numTrials+0.5 (con-1)*numTrials+0.5], 'Color', 'black', 'LineWidth', 1)
end

xlabel('directions in degree')
xticks = 1:TraceLength:TraceLength*metadata.StimParams.numDirections;  
xlabels = metadata.StimParams.directions'; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);

yticks = numTrials/2:numTrials:numTrials*numCon;
set(gca, 'YTick', yticks, 'YTickLabel', conLabels);
if metadata.StimParams.numSf > 1
    ylabel('spatial frequency in cpd');
elseif metadata.StimParams.numTf > 1
    ylabel('temporal frequency in Hz');
end

colorbar
colormap(jet)

if metadata.StimParams.numSf > 1
    title(['ROI ' num2str(roi) ', prefSf: ' num2str(analysis.(field).roi(roi).prefSf) ', SFSI: ' num2str(analysis.(field).roi(roi).SFSI) ', Sf var:' num2str(analysis.(field).roi(roi).SFVar)])
    subtitle(['prefDir indiv. SF: ' num2str(analysis.(field).roi(roi).preferredDirection_allSf)]);
elseif metadata.StimParams.numTf > 1
    title(['ROI ' num2str(roi) ', prefTf: ' num2str(analysis.(field).roi(roi).prefTf) ', TFSI: ' num2str(analysis.(field).roi(roi).TFSI) ', pref Dir:' num2str(analysis.(field).roi(roi).preferredDirection)])
    subtitle(['prefDir indiv. TF: ' num2str(analysis.(field).roi(roi).preferredDirection_allTf)]);
else
    title(['ROI ' num2str(roi) ', pref Ori:' num2str(round(analysis.(field).roi(roi).preferredOrientation,0)) ', pref Dir:' num2str(round(analysis.(field).roi(roi).preferredDirection,0))])
    subtitle(['OSIFit: ' num2str(round(analysis.(field).roi(roi).OSIFit,2)) ', circVar: ' num2str(round(analysis.(field).roi(roi).OriCircVar,2)) ', DSI: ' num2str(round(analysis.(field).roi(roi).DSI,2))])
end

set(gcf, 'Color', 'w')
saveas(gcf, fullfile(saveDir, ['2D_TrialStimTrace_ROI_Nr_' num2str(roi) '.png']))
%close gcf
