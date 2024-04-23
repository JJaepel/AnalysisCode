function plotMetricesROIs(ce, info)
% Plots the distribution of OSI and DSI as well as the distribution of the
% pref ori and pref dir
%
% Inputs:
% - ce: structure containing information about ROI, each row is one roi,
% containing anatomical information as well as functional
% - info: structure containing all information about an experiment,
% including the template and the saveDir
%
% Output:
% - Figure 30: distribution OSI & DSI
% - Figure 31: direction pref Ori and pref Dir

%% Step 0: Set the color palette
cocOri = cbrewer('seq', 'Blues', 12);
cocDir = cbrewer('seq', 'Greens', 12);

%% Step 1: Select the ROIs that are good spines or all cells that are ori-selective
if info.isCell
    goodROIs = linspace(1,length(ce),length(ce)); %all cells
else
    goodROIs = find([ce.good]==1); %select for good spines
end


if ~isempty(goodROIs)
    %% Step 2: Plot OSI/DSI distributions
    figure(30)
    subplot(1,3, 1)
    distributionPlot([ce(goodROIs).OSI]','color', cocOri(6,:)); hold all
    boxplot([ce(goodROIs).OSI])
    title('OSI')
    ylim([0 1])
    
    subplot(1,3,2)
    distributionPlot([ce(goodROIs).DSIvect]','color', cocDir(4,:)); hold all
    boxplot([ce(goodROIs).DSIvect])
    title('DSI (vect)')
    ylim([0 1])
    
    subplot(1,3,3)
    distributionPlot([ce(goodROIs).DSI]','color', cocDir(6,:)); hold all
    boxplot([ce(goodROIs).DSI])
    title('DSI)')
    ylim([0 1])
    
    %set background color to white and save
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(info.saveDir, '30_GratingMetrices.png'))
    
    %% Step 3: Plot pref Ori/Dir distributions
    figure(31)
    subplot(1,2,1)
    distributionPlot([ce(goodROIs).prefOri]','color', cocOri(8,:)); hold all
    boxplot([ce(goodROIs).prefOri])
    title('prefered Orientation')
    ylim([0 180])
    
    subplot(1,2,2)
    distributionPlot([ce(goodROIs).prefDir]','color', cocDir(8,:)); hold all
    boxplot([ce(goodROIs).prefDir])
    title('prefered Direction')
    ylim([0 360])
    
    %set background color to white and save
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(info.saveDir, '31_DistrPref.png'))
end
