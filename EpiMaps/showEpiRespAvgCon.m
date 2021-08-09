function showEpiRespAvgCon(analysis, metadata, trialAveragedMaps, orientationAveragedMaps, directionAveragedMaps, sfAveragedMaps, field, mapType,saveDirectory)
    clippingPercentile = 0.2;
    clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]); 
    spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
    for j = 1:3
        name = mapType{j};
        switch name
            case 'Orientation'
                mapSet    = orientationAveragedMaps;
                maxTheta  = 180;
                fieldName = 'orientationMap';
            case 'Direction'
                mapSet    = directionAveragedMaps;
                maxTheta  = 360;
                fieldName = 'directionMap';
            case 'SpatialFreq'
                mapSet    = sfAveragedMaps;
                maxTheta  = spatFreq_cell{size(sfAveragedMaps,3)};
                fieldName = 'sFMap';
        end
        nStims = size(mapSet,3);
        nRows  = ceil(sqrt(nStims+1));
        nCols  = nRows;

        % Show trial-averaged single-condition maps for orientation/direction
        h = makeFigureFullScreen(figure);
        set(h,'Name',name)
        thetas=(0:maxTheta/(nStims):maxTheta-1);
        for i = 1:nStims
            figure(h); subplot(nRows,nCols,i);
                imagesc(mapSet(:,:,i));
                colorbar; 
                colormap('gray');
                title(num2str(thetas(i)));
                switch name
            	case 'SpatialFreq'
                    title(num2str(spatFreq_cell{i}))
                end
                
                axis image; axis off;
                caxis(clipValue);
        end

        % Show polar map
        figure(h); subplot(nRows,nCols,nStims+1);
            imagesc(polarMapEpi(analysis.(field).roi.(fieldName), [0 clipValue(2)])); 
            axis image; axis off;
            try
                caxis([0 maxTheta]);
            catch
                switch name
                    case 'SpatialFreq'
                        caxis([0 str2num(maxTheta)]);
                end
            end
            cb=colorbar; colormap(cb, hsv);
            title(sprintf('%s map',name))
       saveas(gcf, fullfile(saveDirectory, [fieldName '.png']))
    end
    
    %analysis.(field).roi.stimResponseTrace = permute(analysis.(field).roi.stimResponseTrace, [5 4 3 1 2]);
end