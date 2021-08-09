function polarMap = polarMap(z,LUT,clipValue,discardMagnitude)
% function polarMap = polarMap(z,LUT,clipValue)
% Create a polar response map from a complex field matrix, the resultant 
% data array from the vector addition of single condition responses.
%
% LUT - optional parameter specifying the LUT to use
% Modified by David on 7/23/2011

if(nargin < 2), LUT = hsv;     end
if(nargin < 3), clipValue = 3; end
if(nargin < 4), discardMagnitude = false; end

%isIntrinsic = false;
%if(~isIntrinsic),LUT = LUT([(length(LUT)/2+1):length(LUT),1:length(LUT)/2],:); end % shifts LUT since this is calcium. 

% Construct Magnitude Response Map, Angular Response Map, Polar Response Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.) Magnitude Response Map
if(discardMagnitude)
    rgbMagnitudeMap = repmat(double(ones(size(z))),[1 1 3]);
else
    magnitudeMap = abs(z); 

    % Clip Magnitude Map
    typeOfClipping = length(clipValue);
    switch typeOfClipping
        case 1 % denotes clipping by std around the mean 
            highClipVal = mean(magnitudeMap(:))+clipValue*std(magnitudeMap(:));
            lowClipVal  = mean(magnitudeMap(:))-clipValue*std(magnitudeMap(:));
        case 2 % denotes absolute threshold clipping
            highClipVal = clipValue(2);
            lowClipVal  = clipValue(1);
    end
    upperLimit = find(magnitudeMap > highClipVal );
    lowerLimit = find(magnitudeMap < lowClipVal );
    magnitudeMap(upperLimit) = highClipVal;
    magnitudeMap(lowerLimit) = lowClipVal;


    % Normalize Magnitude Map Between 0 and 1 and convert to RGB
    switch typeOfClipping
        case 1
            offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));
        case 2
            offsetMagnitudeMap = magnitudeMap - min(min(lowClipVal));
            normalizedMagnitudeMap = offsetMagnitudeMap / max(max(highClipVal));
    end
    rgbMagnitudeMap = ind2rgb(im2uint8(normalizedMagnitudeMap)*length(LUT),gray);
    for(currentChannel = 1:3)
        rgbMagnitudeMap(:,:,currentChannel) = double(normalizedMagnitudeMap);
    end  
end

%% 2.) Phase Response Map
phaseMap = angle(z)*180/pi;
normalizedPhaseMap = wrapTo360(phaseMap)./360;
%normalizedPhaseMap=phaseMap;
rgbPhaseMap = double(ind2rgb(uint8(normalizedPhaseMap*length(LUT)),LUT)); 


%% 3.) RGB Polar Response Map
polarMap = rgbPhaseMap.*rgbMagnitudeMap;
