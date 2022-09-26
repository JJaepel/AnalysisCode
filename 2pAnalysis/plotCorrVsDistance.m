function plotCorrVsDistance(analysisParams, PopAnalysis, saveDirectory)

%% Figure 101:Corr vs. distance
figure(101)
errorbar(PopAnalysis.(analysisParams.field).edges(2:end),PopAnalysis.(analysisParams.field).meanCorr, PopAnalysis.(analysisParams.field).semCorr, 'o-', 'Color', analysisParams.coc_prop(2,:), 'MarkerFaceColor', analysisParams.coc_prop(2,:))
hold all
errorbar(PopAnalysis.(analysisParams.field).edges(2:end),mean(PopAnalysis.(analysisParams.field).meanCorrShuffle),mean(PopAnalysis.(analysisParams.field).semCorrShuffle), 'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel('Distance in \mum')
ylabel('Correlation')
%ylim([0 90])
xlim([0 PopAnalysis.(analysisParams.field).edges(end)])
legend('Data', 'Shuffled')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '101_Correlation_vs_distance.png'))

%% Figure 102: long vs. short distance correlation
figure(102)
subplot(1,2,1)
boxplot(PopAnalysis.(analysisParams.field).shortDistanceCorr(:), 'labels', '< 100 um')
ylim([0 1])
subplot(1,2,2)
boxplot(PopAnalysis.(analysisParams.field).longDistanceCorr(:), 'labels', '300-400 um')
ylim([0 1])
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '102_CorrelationDistanceLongvsShort.png'))