function corrTable = computeCorrelationTable(imgStack,ROI)
    % Compute correlations of the input imaging stacking

    % Get image dimensions, reshape tif stack into a transposed 2d array ==> (x,y,n) to (x*y,n)
    imgStackSize = size(imgStack);
        x  = imgStackSize(1);
        y  = imgStackSize(2);
        t  = imgStackSize(3);
    
    %saturate imageStack to account for different maximal activity
    imgStackQuart = quantile(imgStack, 0.9,3);
    imgStack = imgStack ./ imgStackQuart;
    imgStack(imgStack>1) = 1;
    
    imgStack = reshape(imgStack,[x*y t]);
    imgStack = imgStack(ROI(:),:);

    % Compute correlations
    corrTable = corr(imgStack');
end