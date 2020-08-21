function [rgbPolar] = polarMapNew(z, LUT, clipValue, selectedRange)
    if(nargin < 2), LUT = hsv;     end
    if(nargin < 3), clipValue = 3; end
    if(nargin < 4), selectedRange = [0 1]; end

    %%1.) Magnitude Map
    magnitudeMap = abs(z);
    
    % Clip Magnitude Map
    typeOfClipping = length(clipValue);
    switch typeOfClipping
        case 1 % denotes clipping by std around the mean
            highClipVal = 0*mean(magnitudeMap(:))+clipValue*std(magnitudeMap(:));
            lowClipVal  = mean(magnitudeMap(:))-clipValue*std(magnitudeMap(:));
        case 2 % denotes absolute threshold clipping
            highClipVal = clipValue(2);
            lowClipVal  = clipValue(1);
    end
    upperLimit = find(magnitudeMap > highClipVal );
    lowerLimit = find(magnitudeMap < lowClipVal );
    magnitudeMap(upperLimit) = highClipVal;
    magnitudeMap(lowerLimit) = lowClipVal;
    
    % Normalize Magnitude Map Between 0 and 1
    switch typeOfClipping
        case 1
            offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));
        case 2
            offsetMagnitudeMap = magnitudeMap - min(min(lowClipVal));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(highClipVal));
    end
    
    rgbMagnitudeMap = ind2rgb(im2uint8(normalizedMagnitudeMap)*length(LUT),gray);
    for currentChannel = 1:3
        rgbMagnitudeMap(:,:,currentChannel) = double(normalizedMagnitudeMap);
    end
    
    %% 2.) Phase Response Map
    phaseMap = angle(z)*180/pi;
    normalizedPhaseMap = wrapTo360(phaseMap)./360;
    normalizedPhaseMap = normalizedPhaseMap - selectedRange(1);
    normalizedPhaseMap = normalizedPhaseMap ./ (selectedRange(2)-selectedRange(1));
    normalizedPhaseMap(normalizedPhaseMap>1) = 1;
    normalizedPhaseMap(normalizedPhaseMap<0) = 0;
    rgbPhaseMap = double(ind2rgb(uint8(normalizedPhaseMap*length(LUT)),LUT));
    
    %% 3.) RGB Polar Response Map
    rgbPolar = rgbPhaseMap.*rgbMagnitudeMap;
end