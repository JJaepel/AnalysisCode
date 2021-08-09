function [data] = LoadImagingData(EpiDirectory)
    filename = [EpiDirectory filesep '*.tif'];
    files = dir(filename);
    numberOfFiles = size(files, 1);
    if(numberOfFiles ==0), error('The Image path does not contain any tifs'); end
    framesPerFile = zeros(numberOfFiles,1);
    filesToCheck = [1 numberOfFiles];
    for n = filesToCheck
        fileName  = char(strcat(EpiDirectory,'\',files(n).name));
        imageInfo = imfinfo(fileName);
        framesPerFile(n)=numel(imageInfo);
    end
    framesPerFile(framesPerFile==0) = framesPerFile(1);
    numberOfFrames = sum(framesPerFile);
    data.rawF  = zeros([imageInfo(1).Height,imageInfo(1).Width,numberOfFrames],'uint16');
    data.ROI   =true( [imageInfo(1).Height,imageInfo(1).Width]); 
    
    frameCounter = 0;
    for n = 1:numberOfFiles
        fileName  = char(strcat(EpiDirectory,'\',files(n).name));
        data.rawF(:,:,frameCounter+[1:framesPerFile(n)])=readTiffsScaled(fileName,1,100);
        frameCounter=frameCounter+framesPerFile(n);
    end
    data.rawFMeanImg     = mean(data.rawF,3);
end