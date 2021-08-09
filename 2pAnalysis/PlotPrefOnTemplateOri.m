function PlotPrefOnTemplateOri(analysis, data, type, field,template,rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,template,template,template)/prctile(template(:),99));
    colormap(hsv) 
    coc = cbrewer('qual', 'Set1', 8);
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    elseif type == 3
        title('Spatial frequency preference map')
        colorbar('Location', 'southoutside');
    elseif type == 4
        title('Temporal frequency preference map')
        colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        l = rois(i);
        xpos= data.roi(l).xPos;
        ypos= data.roi(l).yPos;
        try
            if type == 1
                plot(xpos,ypos,'ok','MarkerSize',10','MarkerFaceColor',LUT(1+floor(analysis.(field).roi(l).preferredOrientation),:));
            elseif type == 2
                plot(xpos, ypos, 'ok', 'MarkerSize', 10', 'MarkerFaceColor', LUT(1+floor(analysis.(field).roi(l).preferredDirection),:));
            elseif type == 3
                plot(xpos, ypos, 'ok', 'MarkerSize', 10', 'MarkerFacecolor', coc(analysis.(field).roi(l).prefSfStim,:));
            elseif type == 4
                plot(xpos, ypos, 'ok', 'MarkerSize', 10', 'MarkerFacecolor', coc(analysis.(field).roi(l).prefTfStim,:))
            end
        catch
        end
    end
end