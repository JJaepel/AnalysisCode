function plotFieldsOnCellTracing(dendrites, ROIs, IDs,Soma, field, circ, LUT, figNR, plotOther)

%Plots the functional or anatomical properties of ROIs on top of the whole
%cell tracing in black

%Inputs:
% - dendrites: structure containing all information about the dendritic 
% segments that need to be plotted
% - ROIs: structure containing all information about the ROIs of the cell
% - IDs: which ROIs should be looked at
% - Soma: information about the soma
% - field: what property are we looking at?
% - circ: is it a circular measurement such as prefOri and prefDir?
% (- LUT: what is the color scheme? jet with 100 nuances by default)
% (- figNR: What is the figure Nr of this figure? 1 by default)
% (- plotOther: Should you also plot the other ROIs of this cell? No by
% default)

%Steps: 
% - 1.) Plot whole cell tracing in black
% - 2.) Plot soma and its property
% - 3.) Plot ROIs of spines
% - 4.) Add colorbar as legend and title

%Output:
% Figure showing the above

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: November, 2023

if nargin < 7
    plotOther = 1;
    if nargin < 6
        figNR = 1;
        if nargin < 5
            LUT = jet(100);
            if nargin < 4
                circ = 1;
            end
        end
    end
end

%% Step 1.) Plot whole cell tracing in black
figure(figNR) %make figure
set(gca, 'YDir','reverse')
hold on

%plot cell tracing
if contains(field, 'Branch')
   for allD = 1:length(dendrites)
        i = dendrites(allD).Branch;
        plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', LUT(i,:));
    end 
else
    for allD = 1:length(dendrites)
        plot(dendrites(allD).normCoord(:,1) ,dendrites(allD).normCoord(:,2), 'LineWidth',1.5, 'color', 'black');
    end
end

%% Step 2.) Plot soma and its property
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

%% Step 3.) Plot ROIs and their property
for r = 1:length(ROIs)
    xpos= ROIs(r).xPosNorm; %get the xPos of the ROI
    ypos= ROIs(r).yPosNorm; %get the yPos of the ROI
    if IDs(r) == 1 %is it one of the specific ROIs?
        switch field
            case 'SL' %except for when you look at the slice of an ROI, you can get the property from the functional data of the ROI
                fieldContent = ROIs(r).(field);
            case 'type'
                if contains(ROIs(r).(field), 'apical')
                    fieldContent = 1;
                else
                    fieldContent = 2;
                end
            case 'Branch'
                fieldContent = ROIs(r).(field);
            otherwise
                fieldContent = ROIs(r).funcData.(field);
        end
        if circ %decide if circular
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT(1+floor(fieldContent),:));
        else
            switch field 
                case 'SL'
                    plot(xpos, ypos, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', LUT(fieldContent,:));
                case 'type' %slice is much more limited, so different scaling
                    plot(xpos, ypos, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', LUT(fieldContent,:));
                case 'Branch'
                    i = ROIs(r).Branch;
                    plot(xpos,ypos, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', LUT(i,:))
                otherwise  
                    try
                        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:));
                    end
            end
        end
    else
        if plotOther %if not one of those, plot in white if selected
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end

%% Step 4.) Add colorbar as legend and title
axis off %remove axis
colormap(LUT) 
%define the axis in colorrange
if circ 
    colorbar %get the legend
    caxis([0 size(LUT,1)])
else
    if ~contains(field, 'SL') && ~contains(field, 'type') && ~contains(field, 'Branch')
        colorbar %get the legend
        caxis([0 size(LUT,1)/100])
    end
end
title(field, 'FontSize', 14) %add the field as the title
set(gcf, 'color', 'w'); %set background to white