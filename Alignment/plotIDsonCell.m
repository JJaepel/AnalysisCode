function plotIDsonCell(template, ROIs, IDs, figNR, color)

%Plots the ID of certain ROIs in a predefined color on top of the
%projection of the confocal image and all other ROIs in white

%Inputs:
% - template: projection of the confocal image
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at, logical vector the size of the
% length of the ROIs
% (- figNR: What is the figure Nr of this figure? 1 by default)
% (- color: which color should they be plotted in? Red by default)

%Steps: 
% - 1.) Show image of cell in gray colors
% - 2.) Plot ROIs of spines on top

%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023

if nargin < 5
    color = 'red';
    if nargin < 4
        figNR = 1;
    end
end

%% Step 1.) Show image of cell in gray colors

figure(figNR) %make figure
imagesc(template) %show template
axis image %set axis to size of image
colormap('gray') %set colormap to gray
hold on

%% Step 2.) Plot ROIs of spines on top
for r = 1:length(ROIs)
    xpos= ROIs(r).xPos; %get the xPos of the ROI
    ypos= ROIs(r).yPos; %get the yPos of the ROI
    if IDs(r) == 1 %is it one of the specific ROIs?
        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', color); %plot in predefined color
    else
        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white'); %plot in white
    end
    hold on
end

axis off %remove axis
set(gcf, 'color', 'w'); %set background to white