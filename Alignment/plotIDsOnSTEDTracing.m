function plotIDsOnSTEDTracing(dendrites, ROIs, IDs, figNr, plotOther)

%Plots the apical (green) and basal (red) inputs on top of the whole cell 
%tracing in grey as well as the STEDcovered tracing in black 

%Inputs:
% - dendrites: structure containing all information about the dendritic 
% segments that need to be plotted
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at
% (- figNR: What is the figure Nr of this figure? 1 by default)
% (- plotOther: Should you also plot the other ROIs of this cell? No by
% default)

%Steps: 
% - 1.) Plot whole cell tracing in gray
% - 2.) Plot STED tracing in black
% - 3.) Plot ROIs of spines

%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.1: March, 2024

if nargin < 5
    plotOther = 0;
    if nargin < 4
        figNr = 1;
    end
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
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    if IDs(r) == 1 %is it one of the specific ROIs?
        type = ROIs(r).type; %apical or basal?
        switch type
            case 'apical'
                plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'green');
            case 'basal'
                plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'red');
        end
    else
        if plotOther %if not one of those, plot in white if selected
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end

axis off %remove axis
set(gcf, 'color', 'w'); %set background to white