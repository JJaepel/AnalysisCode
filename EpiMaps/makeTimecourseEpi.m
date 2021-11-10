function makeTimecourseEpi(metadata, analysis, analysisParams)

eval(sprintf('cmap = %s(%d);','hsv',metadata.StimParams.uniqStims-1));
LUT = cat(1,cmap, [0 0 0]);
numberOfFrames     = size(analysis.(analysisParams.field).roi.stimResponseTrace,1);

figure
for condition = 1:metadata.StimParams.uniqStims
    t = ((0:(numberOfFrames-1))/metadata.Imaging.rate)-analysisParams.pre;
    y = nanmedian(analysis.(analysisParams.field).roi.stimResponseTrace(:,:,condition,analysis.ROI(:)),4); % Takes median response of the pixel values. Reduces noise relative to mean
    yMean     = squeeze(mean(y,2));
    yStdError = std(y,[],2)/sqrt(size(y,2));
    shadedErrorBar(t,yMean,yStdError,'lineProps',{'-','LineWidth',3,'Color',LUT(condition,:)},'patchSaturation',0.4); hold on;
end
xlabel('Time (s)');
ylabel('Response amplitude');
axis square;
set(gca,'Box','off');
set(gcf, 'color', 'w');

% Add stimulus box
stimStart = 0;
stimStop  = metadata.StimParams.stimDuration;
yLimits = get(gca,'YLim');
rectangle('Position',[stimStart yLimits(2) stimStop-stimStart 0.025*range(yLimits)],'FaceColor','k')
saveas(gcf, fullfile(analysisParams.saveDirectory, 'Timecourse.png'))
close gcf