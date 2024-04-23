function plotIDsOnCellTracing(dendrites, ROIs, IDs, figNR, color, plotOther)

%Plots the inputs in the color of choise on top of the whole cell tracing
%in black

%Inputs:
% - dendrites: structure containing all information about the dendritic 
% segments that need to be plotted
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at
% (- figNR: What is the figure Nr of this figure? 1 by default)
% (- color: Which color should they be plotted in? Red by default)
% (- plotOther: Should you also plot the other ROIs of this cell? No by
% default)

%Steps: 
% - 1.) Plot whole cell tracing in black
% - 2.) Plot ROIs of spines

%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.1: March, 2024

if nargin < 6
    plotOther = 1;
    if nargin < 5
        color = 'red';
        if nargin < 4
            figNR = 1;
        end
    end
end

%% Step 1: plot whole cell tracing in gray
figure(figNR) %make figure
set(gca, 'YDir','reverse')
hold on

%plot cell tracing
for allD = 1:length(dendrites)
    plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', 'black');
    hold on
end

%% Step 2: plot ROIs of spines
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    if IDs(r) == 1 %is it one of the specific ROIs?
        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', color);
    else
        if plotOther %if not one of those, plot in white if selected
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end

axis off %remove axis
set(gcf, 'color', 'w'); %set background to white