function plotROIpositions(data, saveDirectory)

imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
axis image
hold on
for i = 1:length(data.roi)
    xpos= data.roi(i).xPos;
    ypos= data.roi(i).yPos;
    plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor', 'blue');
    text(xpos, ypos, num2str(i),'HorizontalAlignment','center', 'VerticalAlignment','middle','color', 'white') 
end
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '101_ROI_positions.png'))