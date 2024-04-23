function plotIDsonTracing(dendrites, ROIs, IDs, figNR, color)

%Plots the ID of certain ROIs in a predefined color on top of the
%dendrite tracing and all other ROIs in white

%Inputs:
% - dendrites: structure containing all information about the dendritic 
% segments that need to be plotted
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at, logical vector the size of the
% length of the ROIs
% (- figNR: What is the figure Nr of this figure? 1 by default)
% (- color: which color should they be plotted in? Red by default)

%Steps: 
% - 1.) Plot whole cell tracing in black
% - 2.) Plot soma in white
% - 3.) Plot ROIs of spines on top

%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: March, 2024

%% Step 1.) Plot whole cell tracing in black
figure(figNR) %make figure
set(gca, 'YDir','reverse')
hold on

%plot cell tracing
for allD = 1:length(dendrites)
    plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', 'black');
    hold on
end

%% Step 2.) Plot soma

plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', 'white');

%% Step 3.) Plot ROIs of spines on top
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    if IDs(r) == 1 %is it one of the specific ROIs?
        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', color); %plot in predefined color
    else
        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white'); %plot in white
    end
    hold on
end

axis off %remove axis
set(gcf, 'color', 'w'); %set background to white
