function [ContourXData, ContourYData, ContourMap, FilledMap] = makeContour(image, threshold) 

%input:
% - image: which data would you like to draw contours on
% - threshold: which threshold to use for drawing the contours
%
% output:
% - ContourXData: cell containing all x data for the contours
% - ContourYData: cell containing all y data for the contours
% - ContourMap: 2D map of the contours (of at least 5 pixels)
% - FilledMap: thresholded 2D map

[c, ~] = contour(image);
cdata = contourdata(c);
allLevels = unique([cdata.level]);
level = find(allLevels >= threshold);
try
    levelData = find([cdata.level] == allLevels(level(1)));
catch
    levelData = find([cdata.level] == allLevels(end));
end
ContourXData = {cdata(levelData).xdata};
ContourYData = {cdata(levelData).ydata};
ContourMap = zeros(size(image,1), size(image,2));
for con = 1:size(ContourXData,2)
    if size(ContourXData{con},1) > 5
        xValues = round(ContourXData{con});
        yValues = round(ContourYData{con});
        for dots = 1:length(xValues)
            ContourMap(yValues(dots)-1:yValues(dots)+1, xValues(dots)-1:xValues(dots)+1) = 1;
        end
    end
end

FilledMap = zeros(size(image,1), size(image,2));
FilledMap(image<threshold) = 0;
FilledMap(image>threshold) = 1;