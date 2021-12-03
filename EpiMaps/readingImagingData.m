function [tifStack] = readingImagingData(EpiDirectory, spatialDownSamplingFactor)
    if nargin < 2
        spatialDownSamplingFactor = 1;
    end
    
    downsample = 1/spatialDownSamplingFactor;
    tifStack= [];
    tifFiles = dir([EpiDirectory filesep '*.tif']);
    numberOfFiles = size(tifFiles, 1);
    if(numberOfFiles ==0), error('The Image path does not contain any tifs'); end
    for currentFile = 1:length(tifFiles)
        disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);

        % Specify stack name
        filePath = [EpiDirectory tifFiles(currentFile).name];
        newlyLoaded = readTiffsEpi(filePath);
        
        scaledX = round(size(newlyLoaded,1)*downsample);
        scaledY = round(size(newlyLoaded,2)*downsample);
        
        %spatially downsample if specified
        newlyAdded = zeros(scaledX,scaledY,size(newlyLoaded,3),'uint16');
        if spatialDownSamplingFactor == 1
            newlyAdded = newlyLoaded;
        else
            for i = 1:size(newlyLoaded,3)
                newlyAdded(:,:,i) = imresize(newlyLoaded(:,:,i),[scaledX scaledY]);
            end
        end
        
        % Read images into tifStack
        tifStack = cat(3,tifStack,newlyAdded);
    end
end