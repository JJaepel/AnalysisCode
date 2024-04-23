function plotPosSpines(ce, info)

figure
imshow(cat(3,info.template,info.template,info.template)/prctile(info.template(:),99.9));
axis image
hold on
for cc = 1:length(ce)
    xpos= ce(cc).xPos;
    ypos= ce(cc).yPos;
    plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor', 'blue');
    text(xpos, ypos, num2str(cc),'HorizontalAlignment','center', 'VerticalAlignment','middle','color', 'white') 
end
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(info.saveDir, '10_ROI_positions.png'))