function showCorrelationStructure(corrTable,ROI, analysis, saveDirectory)
    % Show correlation table and an interactive clicker to inspect correlation patterns

    % Parse out variable input arguments and load default parameters where necessary
    if(~isfield(analysis,'clippingRange')), analysis.clippingRange = [-1 1];                    end
    if(~isfield(analysis,'LUT')),           analysis.LUT           = rwb(64);                   end
    if(~isfield(analysis,'sortingMatrix')), analysis.sortingMatrix = 1:sum(analysis.ROI(:));     end

    analysis.LUT = cat(1,[0,0,0],analysis.LUT);

    % Create a labeled binary mask (Storing the indexed position of each seed point)
    ROILabeled = zeros(size(ROI));
    ROILabeled(ROI) = 1:sum(ROI(:));
    index = round(size(corrTable,1)/2); % % Default starting indexed position
    [x,y]=find(ROILabeled==index); %Default starting x- and y- position

    % Show correlation table and correlation patterns
    figure; 
    % Show sorted correlation table
    subplot(1,3,1); 
        [~,sortingMatrix] = sort(analysis.sortingMatrix(analysis.maskBVActive(:)));
        sortedCorrTable = corrTable(sortingMatrix(:),:);
        sortedCorrTable = sortedCorrTable(:,sortingMatrix(:));
        imagesc(sortedCorrTable); 
        caxis(analysis.clippingRange); 
        colormap([analysis.LUT]); 
        axis image;
        colorbar; 
        
    %show rawImage 
    subplot(1,3,2)
    imagesc(analysis.rawFMeanImg)
    colormap('gray')
    axis image;
    axis off
    
    % Show an interactive clicker of correlation patterns
    subplot(1,3,3);
    while(1)
        try
            imagesc(recomputeImage(real(corrTable(index,:)),ROI)); 
            colormap(analysis.LUT); axis image; caxis(analysis.clippingRange);
            title(sprintf('x=%d,y=%d,index=%d',round(x),round(y),index)); 
            colorbar
            drawnow;
            % get interactive point to compute different correlation patterns
            [x,y] = ginput(1);
            index = ROILabeled(round(y),round(x));
            saveas(gcf, [saveDirectory, 'CorrelationPatterns' num2str(index) '.png'])
        catch
            sprintf('Error: Outside ROI. Breaking interactive plotter');
            break;
            %index = round(size(corrTable,1)/2);
            %pause(1);
            %continue;
        end
    end
end