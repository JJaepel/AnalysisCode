function [tifStack] = readingImagingData(EpiDirectory)
    tifStack= [];
    tifFiles = dir([EpiDirectory filesep '*.tif']);
    numberOfFiles = size(tifFiles, 1);
    if(numberOfFiles ==0), error('The Image path does not contain any tifs'); end
    for currentFile = 1:length(tifFiles)
        disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);

        % Specify stack name
        filePath = [EpiDirectory tifFiles(currentFile).name];

        % Read images into tifStack
        tifStack = cat(3,tifStack,readTiffsEpi(filePath));
    end
end