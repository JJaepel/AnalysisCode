function tifstack = readTiffsEpi(filePath)
    tHeader = tic();
    disp(['Reading Image Stack - ' filePath]);
    tiffReader = ScanImageTiffReader(filePath);
    tifstack = tiffReader.data();
    tifstack = permute(tifstack,[2 1 3]);
    tiffReader.close();
    disp(['Finished Reading Image Stack - ' num2str(toc(tHeader)) ' seconds Elapsed']);
end