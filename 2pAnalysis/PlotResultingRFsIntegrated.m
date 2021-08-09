function PlotResultingRFsIntegrated(smthONOFF,ONx,ONy,OFFx,OFFy,Significance,DisplayWindow,metadata,ind, figSaveDir)

% load colormaps
load('C:\Users\jaepelj\Documents\GIT\AnalysisCode\HelperFunctions\BlueToWhiteToRed2')
screen_size = get(0, 'ScreenSize');

UsedScreenWidth = abs(metadata.StimParams.minAzi-metadata.StimParams.maxAzi);
UsedScreenHeight = abs(metadata.StimParams.maxElev-metadata.StimParams.minElev);

%create figure
f=figure(10000+ind);
set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%plots the subtracted ON OFF maps
subplot(2,3,1) 
colormap(BlueToWhiteToRed2)
x=[metadata.StimParams.minAzi metadata.StimParams.maxAzi];
y=[metadata.StimParams.maxElev metadata.StimParams.minElev];
imagesc(x,y,smthONOFF),

set(gca,'TickDir','out');
set(gca, 'YDir', 'normal');
set(gca,'TickLen',[0.01 0.01]);
xlabel('azimuth');
ylabel('elevation');
title('RF position map');

%plots the perimeter of the fields
subplot(2,3,2)
hold on
plot(metadata.StimParams.minAzi+ONx,metadata.StimParams.maxElev -ONy,'r.','MarkerSize',6)
plot(metadata.StimParams.minAzi+OFFx,metadata.StimParams.maxElev -OFFy,'b.','MarkerSize',6)
set(gca,'YDir','normal')
axis tight
xlim([metadata.StimParams.minAzi metadata.StimParams.maxAzi]);
ylim([metadata.StimParams.minElev metadata.StimParams.maxElev]);
xlabel('azimuth');
ylabel('elevation');
title('RF perimeter')

%plots the significance of the responses at each frame
subplot(2,6,7:12),hold on,
plot(abs(log10(Significance(1,:))'),'r')
plot(abs(log10(Significance(2,:))'),'b')
plot(abs(log10(Significance(3,:))'),'r:')
plot(abs(log10(Significance(4,:))'),'b:')
plot([1 length(DisplayWindow)],[1.3 1.3],'k')
xlim([1 length(DisplayWindow)])
set(gca,'XTick',1:length(DisplayWindow))
set(gca,'XTickLabel',DisplayWindow(1:end))
ylabel('negative log10 of p')
xlabel('frame relative to stimulus onset')
legend('anova','anova','Krusk.-W.','Krusk.-W.','p=0.05')
title('Time course and significance over response window') 
saveas(f, fullfile(figSaveDir, ['RFs_lat, both RF, ROI ',num2str(ind) '.png']))
pause(0.1);
close(f);