function [rgbMap] = polarMapEpi(z, clippingValue)
    
    selectedRange    = [0 1];

    %%1.) Magnitude Map
    magnitudeMap = abs(z);
    
    % Clip Magnitude Map
    typeOfClipping = length(clippingValue);
    switch typeOfClipping
        case 1 % denotes clipping by std around the mean 
            excludingNaNs = ~isnan(magnitudeMap(:));
            highClipVal = mean(magnitudeMap(excludingNaNs(:)))+clippingValue*std(magnitudeMap(excludingNaNs(:)));
            lowClipVal  = mean(magnitudeMap(excludingNaNs(:)))-clippingValue*std(magnitudeMap(excludingNaNs(:)));
        case 2 % denotes absolute threshold clipping
            highClipVal = clippingValue(2);
            lowClipVal  = clippingValue(1);
    end
    
    magnitudeMap(magnitudeMap > highClipVal) = highClipVal;
    magnitudeMap(magnitudeMap < lowClipVal ) = lowClipVal;
    
    % Normalize Magnitude Map Between 0 and 1
    switch typeOfClipping
        case 1
            offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));
        case 2
            offsetMagnitudeMap = magnitudeMap - min(min(lowClipVal));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(highClipVal));
    end
    
    %% 2.) Phase Response Map
    phaseMap = angle(z)*180/pi;
    normalizedPhaseMap = wrapTo360(phaseMap)./360;
    normalizedPhaseMap = normalizedPhaseMap - selectedRange(1);
    normalizedPhaseMap = normalizedPhaseMap ./ (selectedRange(2)-selectedRange(1));
    normalizedPhaseMap(normalizedPhaseMap>1) = 1;
    normalizedPhaseMap(normalizedPhaseMap<0) = 0;
    
    %% 3.) RGB Polar Response Map
    HueMap        = normalizedPhaseMap;
    SaturationMap = normalizedMagnitudeMap; 
    ValueMap      = ones(size(z)); % normally ignored unless responseImg is specified
    %rgbMap        = hsv2rgb(cat(3,HueMap,SaturationMap,ValueMap));
    
    LUT = hsv;
    rgbMagnitudeMap = ind2rgb(im2uint8(normalizedMagnitudeMap)*length(LUT),gray);
    for(currentChannel = 1:3)
        rgbMagnitudeMap(:,:,currentChannel) = double(normalizedMagnitudeMap);
    end  
    rgbPhaseMap = double(ind2rgb(uint8(normalizedPhaseMap*length(LUT)),LUT)); 
    rgbMap = rgbPhaseMap.*rgbMagnitudeMap;
end