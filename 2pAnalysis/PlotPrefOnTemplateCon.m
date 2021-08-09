function PlotPrefOnTemplateCon(pref, SI, data, type, template,rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,template,template,template)/prctile(template(:),99));
    colormap(hsv) 
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        xpos= data.roi(rois(i)).xPos;
        ypos= data.roi(rois(i)).yPos;
        sizeMarker = 10;
        try
            if type == 1
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor((pref(i))),:));
            elseif type == 2
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor((pref(i))),:));
            end
        catch
        end
    end
end