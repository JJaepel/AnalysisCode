function plotAllContoursSummary(analysisParams, analysis, metadata)

if length(metadata.StimParams.directions) == 16
    nRows = 4;
    nCol = 5;
    plotNrStims = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; 
    plotNrOriMap = 17;
    plotNrQuant = [18 19 20];
elseif length(metadata.StimParams.directions) == 8
    nRows = 4;
    nCol = 4;
    plotNrStims = [1 2 3 4 5 6 7 8];
    plotNrOriMap = 9;
    plotNrQuant = [10 11 12];
end

clipValue = prctile(analysis.(analysisParams.field).AvgRespMaps(:),[analysisParams.clippingPercentile 100-analysisParams.clippingPercentile]); 
figure('units','normalized','outerposition',[0 0 1 1])
for stim = 1:length(metadata.StimParams.directions)
    subplot(nRows, nCol,plotNrStims(stim))
    for trial=1:size(analysis.(analysisParams.field).zScore,4)
        for i = 1:size(analysis.(analysisParams.field).ContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).ContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).ContourXData{stim,trial}{i}, analysis.(analysisParams.field).ContourYData{stim,trial}{i}, 'color', analysisParams.conC(trial,:))
                hold all
            end
        end
    end
    set(gca,'ydir','reverse')
    axis off
    axis equal
    title([num2str(metadata.StimParams.directions(stim)) ' deg'])
end
subplot(nRows,nCol, plotNrOriMap)
cmap=hsv;
OriMap = analysis.(analysisParams.field).orientationMap;
OriMap(~analysis.maskBV(:)) = NaN;
imagesc(polarMapEpi(OriMap, [0 clipValue(2)]));
colormap(cmap), colorbar; axis image ;
axis off
axis equal
title('Orientation map')

subplot(nRows,nCol, plotNrQuant(1))
allOverlap = [analysis.(analysisParams.field).perOverlapA19(:); analysis.(analysisParams.field).perOverlapA19Rnd(:); analysis.(analysisParams.field).perOverlapV1(:); analysis.(analysisParams.field).perOverlapV1Rnd(:)];
boxHelp = [zeros(length(analysis.(analysisParams.field).perOverlapA19(:)), 1); ones(length(analysis.(analysisParams.field).perOverlapA19Rnd(:)), 1); 2*ones(length(analysis.(analysisParams.field).perOverlapV1(:)), 1); 3*ones(length(analysis.(analysisParams.field).perOverlapV1Rnd(:)), 1)];
boxplot(allOverlap, boxHelp, 'Labels',{'A19','A19 Random', 'V1', 'V1 Random'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),analysisParams.cocA19(3,:),'FaceAlpha',.75)
patch(get(h(2),'XData'),get(h(2),'YData'),analysisParams.cocA19(5,:),'FaceAlpha',.75)
patch(get(h(3),'XData'),get(h(3),'YData'),analysisParams.cocV1(3,:),'FaceAlpha',.75)
patch(get(h(4),'XData'),get(h(4),'YData'),analysisParams.cocV1(5,:),'FaceAlpha',.75)
box off
ylabel('percentage overlap')

subplot(nRows,nCol, plotNrQuant(2))
allCorr = [analysis.(analysisParams.field).corrMapsA19(:); analysis.(analysisParams.field).corrMapsA19Rnd(:); analysis.(analysisParams.field).corrMapsV1(:); analysis.(analysisParams.field).corrMapsV1Rnd(:)];
boxHelp = [zeros(length(analysis.(analysisParams.field).corrMapsA19(:)), 1); ones(length(analysis.(analysisParams.field).corrMapsA19Rnd(:)), 1); 2*ones(length(analysis.(analysisParams.field).corrMapsV1(:)), 1); 3*ones(length(analysis.(analysisParams.field).corrMapsV1Rnd(:)), 1)];
boxplot(allCorr, boxHelp, 'Labels',{'A19','A19 Random', 'V1', 'V1 Random'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),analysisParams.cocA19(3,:),'FaceAlpha',.75)
patch(get(h(2),'XData'),get(h(2),'YData'),analysisParams.cocA19(5,:),'FaceAlpha',.75)
patch(get(h(3),'XData'),get(h(3),'YData'),analysisParams.cocV1(3,:),'FaceAlpha',.75)
patch(get(h(4),'XData'),get(h(4),'YData'),analysisParams.cocV1(5,:),'FaceAlpha',.75)
box off
ylabel('correlation')

subplot(nRows,nCol, plotNrQuant(3))
allConF = [analysis.(analysisParams.field).convFactA19(:); analysis.(analysisParams.field).convFactV1(:)];
boxHelp = [zeros(length(analysis.(analysisParams.field).convFactA19(:)), 1); ones(length(analysis.(analysisParams.field).convFactV1(:)), 1)];
boxplot(allConF, boxHelp, 'Labels',{'A19','V1'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),analysisParams.cocA19(3,:),'FaceAlpha',.75)
patch(get(h(2),'XData'),get(h(2),'YData'),analysisParams.cocV1(3,:),'FaceAlpha',.75)
box off
ylabel('Convergence Factor')

set(gcf, 'color', 'w')
saveas(gcf, fullfile(analysisParams.saveDirectory, 'TrialCorrelationAllStims.png'))