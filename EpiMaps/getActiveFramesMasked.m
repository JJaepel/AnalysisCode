function [activeFrames, maskBVActive, eventOnsets] = getActiveFramesMasked(imageStack,maskBV)
    %this function takes raw traces from the imaging and gives filtered
    %active frames
    
    imageStack = im2double(imageStack);
    imageStack = imageStack .* maskBV;
        imgStackSize = size(imageStack);
            x  = imgStackSize(1);
            y  = imgStackSize(2);
            t  = imgStackSize(3);

    imageStack = reshape(imageStack,[x*y t]);
    imageStack_raw = reshape(imageStack,[x y t]);
    
    all_frames = sum(imageStack,1);
    [peaks,locations] = findpeaks(smooth(all_frames(:)));
    
    peaks_std = std(peaks);
    peask_mean = mean(peaks);
    threshold = peask_mean + .05 * peaks_std;
    
    peak_thresholded = find(peaks>threshold);
    peak_thresholded = locations(peak_thresholded,1);
    activeFrames = imageStack_raw(:,:,peak_thresholded);
    activeFrames = double(activeFrames);
    
    %cut NaN border values
    [xRange,yRange] = find(squeeze(activeFrames(:,:,1)));
    xRange = sort(xRange);
    x_bottom = min(xRange(:)); x_ceil = max(xRange(:)); y_bottom = min(yRange(:)); y_ceil = max(yRange(:));
    
    %if x > --> zeroes to xdim; if y > --> zeroes to ydim;
    numb_x = x_ceil - x_bottom +1;
    numb_y = y_ceil - y_bottom +1;
    diff_x = numb_x - numb_y;
    activeFrames = activeFrames((x_bottom:x_ceil),(y_bottom:y_ceil+diff_x),:);
    
    maskBVActive = maskBV((x_bottom:x_ceil),(y_bottom:y_ceil+diff_x));
    
    eventOnsets = peak_thresholded;
end