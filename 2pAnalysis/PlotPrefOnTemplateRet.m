function PlotPrefOnTemplateRet(analysis, data, metadata,type, field,template,rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    temp = cat(3,template,template,template)/prctile(template(:),99);
    if type == 1
        NumStims = metadata.StimParams.uniqStims-1;
        title('Patch map')
    elseif type == 2
        title('Elevation map')
        NumStims = metadata.StimParams.numElevation;
    elseif type == 3
        NumStims = metadata.StimParams.numAzimuth;
        title('Azimuth map')
    end
    LUT = jet(NumStims);
    %make a legend
    CM = colormap('jet');
    IX = round(linspace(1,64,NumStims));
    for cc = 1:length(IX)
        C{cc} = CM(IX(cc),:);
    end
    
    [yRes,xRes,~] = size(temp);
    xRes = xRes - 100;
    for s = 1:NumStims
        xRange = round(50+(((s*(xRes/NumStims)) - (xRes/(NumStims*2))):(s*(xRes/NumStims))));
        for c = 1:3
            temp(yRes:yRes+10,xRange,c) = C{s}(c);
        end
    end
    imshow(temp);
    
    %axis image
    hold on
    for i = 1:length(rois)
        l = rois(i);
        xpos= data.roi(l).xPos;
        ypos= data.roi(l).yPos;
        try
           if type == 1
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(l).prefPatch,:));
            elseif type == 2
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(l).prefElev,:));
            elseif type == 3
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(l).prefAzi,:));
           end
        catch
        end
    end
end