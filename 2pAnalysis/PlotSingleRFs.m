function PlotSingleRFs(SingleRF,Significance,DisplayWindow,metadata,ind,figSaveDir,fieldType, fieldX, fieldY)

% load colormaps
load('C:\Users\jaepelj\Documents\GIT\AnalysisCode\HelperFunctions\BlueToWhiteToRed2')
screen_size = get(0, 'ScreenSize');

%UsedScreenWidth = metadata.StimParams.numStimAzi * metadata.StimParams.stimSize(1);
%UsedScreenHeight = metadata.StimParams.numStimElev * metadata.StimParams.stimSize(2);

%create figure
f=figure(10000+ind);
set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%plots the subtracted ON OFF maps
subplot(2,3,1)
colormap(BlueToWhiteToRed2)
x=[metadata.StimParams.minAzi metadata.StimParams.maxAzi];
y=[metadata.StimParams.maxElev metadata.StimParams.minElev];
imagesc(x,y,SingleRF),

set(gca,'TickDir','out');
set(gca, 'YDir','normal')
set(gca,'TickLen',[0.01 0.01]);
xlabel('azimuth');
ylabel('elevation');
title('Single RF position map');

%plots the perimeter of the fields
subplot(2,3,2)
title('single RF perimeter')
hold on
if fieldType == 1 %ON field
    plot(metadata.StimParams.minAzi+fieldX,metadata.StimParams.maxElev-fieldY,'r.','MarkerSize',6)
else
    plot(metadata.StimParams.minAzi+fieldX,metadata.StimParams.maxElev -fieldY,'b.','MarkerSize',6)
end
set(gca,'YDir','normal')
axis tight
xlim([metadata.StimParams.minAzi metadata.StimParams.maxAzi]);
try
    ylim([metadata.StimParams.minElev metadata.StimParams.maxElev]);
catch
    ylim([metadata.StimParams.maxElev metadata.StimParams.minElev]);
end
xlabel('azimuth');
ylabel('elevation');

%plots the significance of the responses at each frame
subplot(2,6,7:12),
hold on,

if fieldType == 1 
    plot(abs(log10(Significance(1,:,ind))'),'r')
    plot(abs(log10(Significance(3,:,ind))'),'r:')
else
    plot(abs(log10(Significance(2,:,ind))'),'b')
    plot(abs(log10(Significance(4,:,ind))'),'b:')
end
plot([1 length(DisplayWindow)],[1.3 1.3],'k')
xlim([1 length(DisplayWindow)])
set(gca,'XTick',1:length(DisplayWindow))
set(gca,'XTickLabel',DisplayWindow(1:end))
ylabel('negative log10 of p')
xlabel('frame relative to stimulus onset')
legend('anova','Krusk.-W.','p=0.05')
title('Time course and significance over response window') 
saveas(f, fullfile(figSaveDir, ['RFs_lat, single RF, ROI' num2str(ind) '.png']))
pause(0.1);
close(f);
pause(0.1);