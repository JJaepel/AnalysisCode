function img = recomputeImage(array,ROI)
    % recompute an image from a 1d array
    img = zeros(size(ROI));
    img(ROI) = array;
end