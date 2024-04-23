function plotCompareFieldsOnSTEDTracing(dendrites, ROIs, IDs, inputIDs, Soma, field, circ, LUT, figNR, plotOther)

%Plots the functional or anatomical properties of input ROIs on top of the 
%whole cell tracing in grey as well as the STEDcovered tracing in black and
%compares it to the cell tuning

%Inputs:
% - dendrites: structure containing all information about the dendritic 
% segments that need to be plotted
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at
% - inputIDs: which ROIs are the corresponding inputs
% - Soma: information about the soma
% - field: what property are we looking at?
% - circ: is it a circular measurement such as prefOri and prefDir?
% (- LUT: what is the color scheme? jet with 100 nuances by default)
% (- figNR: What is the figure Nr of this figure? 1 by default)
% (- plotOther: Should you also plot the other ROIs of this cell? No by
% default)

%Steps: 
% A) All ROIs
% - 1.) Plot whole cell tracing in black
% - 2.) Plot whole cell tracing in black
% - 3.) Plot soma and its property
% - 4.) Plot ROIs of spines
% - 5.) Add colorbar as legend and title
% B) Repeat for input ROIs


%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023

if nargin < 10
    plotOther = 0;
    if nargin < 9
        figNR = 1;
        if nargin < 8
            LUT = jet(100);
            if nargin < 7
                circ = 1;
            end
        end
    end
end

%% A) Start with all ROIs

% Step 1.) Plot whole cell tracing in gray
figure(figNR) %make figure
subplot(1,2,1)
set(gca, 'YDir','reverse')
hold on

%plot cell tracing
for allD = 1:length(dendrites)
    plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', [0.7 0.7 0.7]);
    hold on
end

% Step 2.) Plot STED tracing in balck
for allSTED = 1:length(dendrites)
    if ~isempty(dendrites(allSTED).STEDTrace)
        plot(dendrites(allSTED).STEDTrace(:,1) ,dendrites(allSTED).STEDTrace(:,2), 'LineWidth',2, 'color', 'black');
        hold on
    end
end

% Step 3.) Plot soma and its property
switch field 
    case 'deltaOri' %deltaOri is 0 by definition
        plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', LUT((1+0),:));
    otherwise
        if isfield(Soma,field)
            fieldContent = Soma.(field); %get the property for the soma
            if circ %decide if it is circular & plot
                plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', LUT(1+floor(fieldContent),:)); 
            else
                plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:));
            end
        end
end

% Step 4.) Plot ROIs and their property
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    if IDs(r) == 1 %is it one of the specific ROIs?
        fieldContent = ROIs(r).funcData.(field);
        if fieldContent > length(LUT)-1
            fieldContent = length(LUT)-1;
        end
        
        if circ %decide if circular
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT(1+floor(fieldContent),:));
        else 
            try
                plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:));
            end
        end
    else
        if plotOther %if not one of those, plot in white if selected
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end

% Step 5.) Add colorbar as legend and title
axis off %remove axis
colormap(LUT) 
colorbar %get the legend
%define the axis in colorrange
if circ 
    caxis([0 size(LUT,1)])
else
    caxis([0 size(LUT,1)/100])
end
title('All ROIs', 'FontSize', 12) %add the field as the title

%% B.) Repeat with input ROIs
subplot(1,2,2)
set(gca, 'YDir','reverse')
hold on

%plot cell tracing
for allD = 1:length(dendrites)
    plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', [0.7 0.7 0.7]);
    hold on
end

% Step 2.) Plot STED tracing in balck
for allSTED = 1:length(dendrites)
    if ~isempty(dendrites(allSTED).STEDTrace)
        plot(dendrites(allSTED).STEDTrace(:,1) ,dendrites(allSTED).STEDTrace(:,2), 'LineWidth',2, 'color', 'black');
        hold on
    end
end

% Step 3.) Plot soma and its property
switch field 
    case 'deltaOri' %deltaOri is 0 by definition
        plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', LUT((1+0),:));
    otherwise
        if isfield(Soma,field)
            fieldContent = Soma.(field); %get the property for the soma
            if circ %decide if it is circular & plot
                plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', LUT(1+floor(fieldContent),:)); 
            else
                plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:));
            end
        end
end

% Step 4.) Plot ROIs and their property
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    if inputIDs(r) == 1 %is it one of the specific ROIs?
        fieldContent = ROIs(r).funcData.(field);
        if circ %decide if circular
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT(1+floor(fieldContent),:));
        else 
            try
                plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:));
            end
        end
    else
        if plotOther %if not one of those, plot in white if selected
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end

% Step 5.) Add colorbar as legend and title
axis off %remove axis
colormap(LUT) 
colorbar %get the legend
%define the axis in colorrange
if circ 
    caxis([0 size(LUT,1)])
else
    caxis([0 size(LUT,1)/100])
end
title('Inputs', 'FontSize', 14) %add the field as the title

%% Finish them
sgtitle(field, 'FontSize', 20) %add the field as the title
set(gcf, 'color', 'w'); %set background to white


screenSize = get(groot, 'ScreenSize');
set(gcf, 'Position', [0.05*screenSize(3),0.05*screenSize(4),0.8*screenSize(3),0.4*screenSize(3)])