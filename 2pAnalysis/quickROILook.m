ROINR = 86;


cyc = ce(ROINR).cycRes;
meanResp = ce(ROINR).meanResp;
peakFit = ce(ROINR).peakFit;
n = 16;
angs = 0:pi/8:2*pi-(pi/8);
xs =   0:0.1:(2*pi);

figure
subplot(2,1,1)
plotcycRes2(cyc,1:n)
box off
set(gca,'TickDir','out')
set(gca,'Xtick',[])
ylabel 'NormResp'

subplot(2,1,2)
plot(rad2deg(angs),meanResp,'ok',rad2deg(xs),peakFit,'--k')
box off
ylabel 'Peak'
xticks([0 90 180 270]) 
set(gca,'TickDir','out')
ylim([0 max(max(meanResp), max(peakFit))])
xlim([0 360])
xlabel 'Direction in deg'
set(gcf, 'Color', 'w')