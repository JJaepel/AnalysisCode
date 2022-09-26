function removeOldFiles(folder)

filePattern = fullfile(folder, '*.png'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(folder, baseFileName);
    delete(fullFileName);
end