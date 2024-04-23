function plotPrefROIs(ce,info)

% Plots the prefered orientation or direction of spines or cells on top of
% the template as a quick way to see the distribution of preferences. 
%
% Inputs:
% - ce: structure containing information about ROI, each row is one roi,
% containing anatomical information as well as functional
% - info: structure containing all information about an experiment,
% including the template and the saveDir
%
% Output:
% - Figure 20: orientation preference on top of template
% - Figure 21: direction preference on top of template

%% Step 1: Select the ROIs that are good spines or all cells that are ori-selective
if info.isCell
    goodROIs = find([ce.OSI] > 0.15);
else
    goodROIs = find([ce.good]==1); %select for good spines
end

%only continue if those are not empty
if ~isempty(goodROIs)
    
    %% Step 2: Plot ori pref on template
    %Make a figure with the template and an orientation colorbar
    figure(20)
    imshow(cat(3,info.template,info.template,info.template)/prctile(info.template(:),99));
    colormap(hsv)
    LUT = hsv(180);
    title('Orientation preference map')
    caxis([0 180]); colorbar('Location', 'southoutside');

    axis image
    hold on
    
    % Plot ROIs and their prefOri on top
    for r = 1:length(goodROIs)
        cc = goodROIs(r);
        xpos=ce(cc).xPos;
        ypos= ce(cc).yPos;
        plot(xpos,ypos,'ok','MarkerSize',10','MarkerFaceColor',LUT(1+floor(ce(cc).prefOri),:));
    end
    
    %remove box, set background to white and save
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(info.saveDir, '20_OriPref.png'))
    
    %% Step 3: Plot dir pref on template
    %Make a figure with the template and an direction colorbar
    figure(21)
    imshow(cat(3,info.template,info.template,info.template)/prctile(info.template(:),99));
    colormap(hsv)
    LUT = hsv(360);
    title('Direction preference map')
    caxis([0 360]); colorbar('Location', 'southoutside');

    axis image
    hold on
    
    % Plot ROIs and their pref dir on top
    for r = 1:length(goodROIs)
        cc = goodROIs(r);
        xpos=ce(cc).xPos;
        ypos= ce(cc).yPos;
        plot(xpos,ypos,'ok','MarkerSize',10','MarkerFaceColor',LUT(1+floor(ce(cc).prefDir),:));
    end
    
    %remove box, set background to white and save
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(info.saveDir, '21_Dir.png'))
end

