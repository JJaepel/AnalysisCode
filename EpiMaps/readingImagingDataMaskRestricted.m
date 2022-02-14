function [tifStack] = readingImagingDataMaskRestricted(EpiDirectory, spatialDownSamplingFactor,restrict)
    if nargin < 2
        spatialDownSamplingFactor = 1;
    end
    if nargin < 3
        restrict = 1;
    end
    
    downsample = 1/spatialDownSamplingFactor;
    tifStack= [];
    tifFiles = dir([EpiDirectory filesep '*.tif']);
    numberOfFiles = size(tifFiles, 1);
    if(numberOfFiles ==0), error('The Image path does not contain any tifs'); end
    
    if restrict
        disp(['Loading Imaging Stack to make mask']);
        filePath = [EpiDirectory tifFiles(1).name];
        stack = readTiffsEpi(filePath);
        changeTime = std(double(stack),[],3);
        threshold = (changeTime-min(changeTime(:)))/(max(changeTime(:))-min(changeTime(:)));
        threshold(threshold<0.1)=0; threshold(threshold>0.1)=1;
        xVals = find(sum(threshold,2)>5); yVals = find(sum(threshold,1)>5);
        if min(xVals) < 6
            xRange = 1:max(xVals)+5;
        elseif max(xVals)+5 > size(threshold,2)
            xRange= min(xVals)-5:max(xVals);
        else
            xRange = min(xVals)-5:max(xVals)+5;
        end
        if min(yVals) < 6
            yRange = 1:max(yVals)+5;
        elseif max(yVals)+5 > size(threshold,1)
            xRange= min(yVals)-5:max(yVals);
        else
            yRange = min(yVals)-5:max(yVals)+5;
        end   
    end
    
    for currentFile = 1:length(tifFiles)
        disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);

        % Specify stack name
        filePath = [EpiDirectory tifFiles(currentFile).name];
        newlyLoaded = readTiffsEpi(filePath);
        newlyLoaded = newlyLoaded(xRange, yRange,:);
        
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