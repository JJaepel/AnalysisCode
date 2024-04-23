function plotFunctInputsOnCellTracing(dendrites, ROIs, IDs, figNr)

%Plots the 2p matched inputs on top of the whole cell tracing in grey as 
%well as the STEDcovered tracing in black, separated by those that are just
%there (asterix), those that are functionally responsive, but not
%ori-selective and those that are ori-selective

%Inputs:
% - dendrites: structure containing all information about the dendritic 
% segments that need to be plotted
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at
% (- figNR: What is the figure Nr of this figure? 1 by default)

%Steps: 
% - 1.) Plot whole cell tracing in gray
% - 2.) Plot STED tracing in black
% - 3.) Plot ROIs of spines in different ways, dependent on whether ROI is
% good and selective

%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023

if nargin < 4
    figNr = 1;
end

%% Step 1: plot whole cell tracing in gray
figure(figNr) %make figure
set(gca, 'YDir','reverse')
hold on

for allD = 1:length(dendrites)
    plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', [.7 .7 .7]);
    hold on
end

%% Step 2: plot STED tracing in black
for allSTED = 1:length(dendrites)
    if ~isempty(dendrites(allSTED).STEDTrace)
        plot(dendrites(allSTED).STEDTrace(:,1) ,dendrites(allSTED).STEDTrace(:,2), 'LineWidth',2, 'color', 'black');
        hold on
    end
end

%% Step 3: plot ROIs of spines
LUT = hsv(180);
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    
    if IDs(r) == 1 %is it one of the specific ROIs?
        prefOri = ROIs(r).funcData.prefOri;
        if ROIs(r).good %is it a functionally responsive ROI?
            if ROIs(r).OSI > 0.1 %is it ori-selective? then plot the oripref
                plot(xpos,ypos,'ok','MarkerSize',7,'MarkerFaceColor', LUT(1+floor(prefOri),:));
            else %else just a black dot
                plot(xpos, ypos, 'ok', 'MarkerSize', 7, 'MarkerFaceColor', 'black');
            end
            
        else %else just a black asterix11
            plot(xpos, ypos, '*r', 'MarkerSize', 7);
        end

    end
    hold on
end

axis off %remove axis
set(gcf, 'color', 'w'); %set background to white

%% Step 4: Add colorbar as legend and title
axis off %remove axis
colormap(LUT) 
colorbar %get the legend
%define the axis in colorrange
caxis([0 size(LUT,1)])
title('Inputs', 'FontSize', 14) %add the field as the title
