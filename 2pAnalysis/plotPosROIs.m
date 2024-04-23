function plotPosROIs(ce, info)

%Plots the position of the ROIs with their ROINr on top of the template
%
% Inputs:
% - ce: structure containing information about ROI, each row is one roi,
% containing anatomical information as well as functional
% - info: structure containing all information about an experiment,
% including the template and the saveDir
%
% Output:
% - Figure with template and ROIs on top, saved in save dir

%% Step 1: Make a figure and plot the template
figure
imshow(cat(3,info.template,info.template,info.template)/prctile(info.template(:),99.9));
axis image
hold on

%% Step 2: For each ROI, get the position, plot the marker and the number of the ROI
for cc = 1:length(ce)
    xpos= ce(cc).xPos;
    ypos= ce(cc).yPos;
    plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor', 'blue');
    text(xpos, ypos, num2str(cc),'HorizontalAlignment','center', 'VerticalAlignment','middle','color', 'white') 
end

%% Step 3: Remove box, set background to white and save it
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(info.saveDir, '10_ROI_positions.png'))