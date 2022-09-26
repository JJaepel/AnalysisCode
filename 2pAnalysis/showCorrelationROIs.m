function showCorrelationROIs(analysisParams, PopAnalysis, saveDirectory)

figure(100)
imagesc(PopAnalysis.(analysisParams.field).corrROIs);
%caxis(analysis.clippingRange); 
colormap(rwb(64)); 
axis image;
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '100_CorrelationTable.png'))