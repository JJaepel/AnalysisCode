function plotPopulationMetrices(analysisParams, PopAnalysis, saveDirectory)
% Figure 301: Delta ori vs. distance
% Figure 302: homogeneity index
% Figure 303: Correlation s. distance
% Figure 304: Trial patterns
% Figure 305: Stimulus decoder

%% Figure 301: Delta ori vs. distance & homogeneity index
figure(301)
errorbar(PopAnalysis.(analysisParams.field).edges(2:end),PopAnalysis.(analysisParams.field).meanDeltaOri, PopAnalysis.(analysisParams.field).semDeltaOri, 'o-', 'Color', analysisParams.coc_prop(2,:), 'MarkerFaceColor', analysisParams.coc_prop(2,:))
hold all
errorbar(PopAnalysis.(analysisParams.field).edges(2:end),mean(PopAnalysis.(analysisParams.field).meanDeltaOriShuffle),mean(PopAnalysis.(analysisParams.field).semDeltaOriShuffle), 'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel('Distance in \mum')
ylabel('\DeltaOrientation preference (\circ)')
ylim([0 90])
try
    xlim([0 PopAnalysis.(analysisParams.field).edges(end)])
end
legend('Data', 'Shuffled')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '301_OSI_vs_distance.png'))

%% Figure 302: homogeneity index
HI = PopAnalysis.(analysisParams.field).ori_cells.HomeogeneityIndex;
figure(302)
subplot(1,3,1)
plot(nanmedian(HI(:,1)), max(histcounts(HI(:,1))),'v','MarkerSize', 8','MarkerEdgeColor',analysisParams.coc_prop(5,:),'MarkerFaceColor',analysisParams.coc_prop(5,:)); hold on
text(nanmedian(HI(:,1)), 1.05 * max(histcounts(HI(:,1))),num2str(round(100*nanmedian(HI(:,1)))/100),'HorizontalAlignment','center', 'Color', analysisParams.coc_prop(6,:), 'FontSize', 12')
histogram(HI(:,1),10,'FaceColor', analysisParams.coc_prop(5,:), 'EdgeColor', analysisParams.coc_prop(6,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (50 \mum)');
set(gca,'Box','off');
try
    subplot(1,3,2)
    hdl(1) = cdfplot(HI(:,1)); hold all
    hdl(2) = cdfplot(HI(:,2)); hold all
    hdl(3) = cdfplot(HI(:,3)); hold all
    hdl(4) = cdfplot(HI(:,4)); hold all
    hdl(5) = cdfplot(HI(:,5)); hold all
    set(hdl(1), 'Color', analysisParams.coc_prop(1,:), 'LineWidth', 2)
    set(hdl(2), 'Color', analysisParams.coc_prop(2,:), 'LineWidth', 2)
    set(hdl(3), 'Color', analysisParams.coc_prop(3,:), 'LineWidth', 2)
    set(hdl(4), 'Color', analysisParams.coc_prop(4,:), 'LineWidth', 2)
    set(hdl(5), 'Color', analysisParams.coc_prop(5,:), 'LineWidth', 2)
    grid off
    xlabel('Homogeneity Index')
    xlim([0 1])
    ylabel('Cumulative fraction of cells')
    legend('0-50', '50-100', '100-150', '150-200', '200-250', 'Location', 'NorthWest'); legend('boxoff')
    set(gca,'Box','off');
    title('')
end
subplot(1,3,3)
avg_HI = nanmean(HI);
sem_HI = nanstd(HI)/sqrt(length(HI));
radius = linspace(0,250,6);
try
    errorbar(radius(2:end), avg_HI(1:end-1), sem_HI(1:end-1),'o-', 'Color', analysisParams.coc_prop(6,:), 'MarkerFaceColor', analysisParams.coc_prop(6,:))
end
xlabel('Maximal radius in \mum')
ylabel('Homeogeneity Index')
xlim([0 250])
ylim([0 1])
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '302_HomeogeneityIndex.png'))

%% Figure 303: Correlation vs. distance
figure(303)
subplot(1,2,1)
boxplot(PopAnalysis.(analysisParams.field).shortDistanceCorr(:), 'labels', '< 100 um')
ylim([-1 1])
subplot(1,2,2)
boxplot(PopAnalysis.(analysisParams.field).longDistanceCorr(:), 'labels', '300-400 um')
ylim([-1 1])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '303_ResponseCorrelationDistance.png'))

%% Figure 304: Trial patterns
figure(304)
subplot(1,6,[1 2])
imagesc(PopAnalysis.(analysisParams.field).cellPatternSorted)
c = gray;
c = flipud(c);
colormap(c)
colorbar
xlabel('Trial number')
ylabel('Cell number')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'TrialResponses.png'))

subplot(1,6,[3 4])
imagesc(PopAnalysis.(analysisParams.field).corrCoeff)
colormap(redblue)
caxis([-1 1])
colorbar
set(gca, 'YDir', 'normal')
xlabel('Trial Number')
ylabel('Trial number')

subplot(1,6,5)
boxplot(PopAnalysis.(analysisParams.field).corrTrialsMatched(:), 'labels', 'Matched')
ylim([-0.5 1])
subplot(1,6,6)
boxplot(PopAnalysis.(analysisParams.field).corrTrialsOrtho(:), 'labels', 'Orthogonal')
ylim([-0.5 1])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '304_TrialPatterns.png'))

%% Figure 305: Stimulus decoder
if analysisParams.predictor
    figure(305)
    DecoderPerformanceAll = [PopAnalysis.(analysisParams.field).decoderData{1}, PopAnalysis.(analysisParams.field).decoderData{2}];
    groups = {'Data', 'Shuffle'};
    boxplot(DecoderPerformanceAll, 'labels', groups);
    ylim([0 1])
    ylabel('Fraction of correctly decoded trials');
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, '305_DecoderPerformance.png'))
end
