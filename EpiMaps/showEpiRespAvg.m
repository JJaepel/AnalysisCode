function analysis = showEpiRespAvg(analysis, metadata, field, saveDirectory)
    analysis.(field).roi.stimResponseTrace = permute(analysis.(field).roi.stimResponseTrace, [4 5 3 2 1]);
    includedFrames = [round(metadata.Imaging.rate * metadata.StimParams.isi/2)+1:round(metadata.Imaging.rate * metadata.StimParams.isi/2)+ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration)];
    stimResponseTrace = mean(analysis.(field).roi.stimResponseTrace(:,:,1:(end-1),:,includedFrames),5);
    
    trialAveragedMaps = squeeze(mean(stimResponseTrace,4));
    trialAveragedMaps(isnan(trialAveragedMaps(:))) = 0;
    analysis.(field).trialAveragedMaps = trialAveragedMaps;
    
    analysis.(field).orientationMap = times(vectorSum(trialAveragedMaps,2,3), analysis.maskBV);
    analysis.(field).directionMap   = times(vectorSum(trialAveragedMaps,1,3), analysis.maskBV);

    clippingPercentile = 0.2;
    clipValue = prctile(trialAveragedMaps(:),[clippingPercentile 100-clippingPercentile]); 
    
    mapType = {'Orientation','Direction'};
    for j = 1:2 
        name = mapType{j};
        switch name
            case 'Orientation'
                nOri      = (size(trialAveragedMaps,3)/2);
                mapSet    = (trialAveragedMaps(:,:,1:nOri) ...
                            +trialAveragedMaps(:,:,(nOri+1):end))/2;
                maxTheta  = 180;
                fieldName = 'orientationMap';
            case 'Direction'
                mapSet    = trialAveragedMaps;
                maxTheta  = 360;
                fieldName = 'directionMap';
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
                stimMap = mapSet(:,:,i);
                stimMap(~analysis.maskBV(:)) = NaN;
                imagesc(stimMap);
                colorbar; 
                colormap('gray');
                title(num2str(thetas(i)))
                axis image; axis off;
                %caxis(clipValue);
        end

        % Show polar map
        figure(h); subplot(nRows,nCols,nStims+1);
        imagesc(polarMapEpi(analysis.(field).(fieldName), [0 clipValue(2)])); 
        axis image; axis off;
        caxis([0 maxTheta]);
        cb=colorbar; colormap(cb, hsv);
        title(sprintf('%s map',name))
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDirectory, [fieldName '.png']))
        close gcf
        
        % Show magnitude map
        figure
        subplot(1,2,1)
        imagesc(magMap(analysis.(field).(fieldName),[0 clipValue(2)]));
        colormap('gray');
        axis image; axis off;
        caxis(clipValue);
        title(sprintf('%s magnitude map',name))

        %show polar map next to magnitude map
        subplot(1,2,2)
        cmap=hsv;
        imagesc(polarMapEpi(analysis.(field).(fieldName), [0 clipValue(2)]));
        colormap(cmap), cbh=colorbar; axis image ;
        title(sprintf('%s hsv map',name))
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDirectory, [fieldName '_Magnitude_HSV.png']))
        close gcf
    end
    
    analysis.(field).roi.stimResponseTrace = permute(analysis.(field).roi.stimResponseTrace, [5 4 3 1 2]);
end