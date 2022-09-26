function analysis = findActivity(analysis, data, field, ROI, globalThreshold)

%% 1.) get exp trace from ROI

expTrace = data.roi(ROI).dff;

%% 2.) Find peaks and threshold them to be above the average
smoothTrace = smooth(smooth(expTrace,10),10);
[peaksValue,peakTime] = findpeaks(smoothTrace);
cellThreshold = mean(expTrace(:));
threshold = max(cellThreshold, globalThreshold);
peakThresholded = find(peaksValue>threshold);
peakTime = peakTime(peakThresholded);
peaksValue = peaksValue(peakThresholded);

%% 3.) Find event onset and offset
dydx = gradient(expTrace);
derivatives = zscore(dydx);
onsetThreshold = std(derivatives);
offsetThreshold = -0.5*std(derivatives);

for event = 1:length(peakTime)      
    %find onset
    %first look for the first value that is half the size of the peak, but
    %at least higher than one std
    halfMax = expTrace(peakTime(event))*0.5;
    halfMax = max(halfMax,std(expTrace(:))); 
    onsetHalfMax = peakTime(event)-1;
    try
        while  expTrace(onsetHalfMax) > halfMax
            onsetHalfMax = onsetHalfMax-1;
        end
        %then look for how long from there the derivativ is positive
        onsetTimeTotal = onsetHalfMax;
        while derivatives(onsetTimeTotal) > onsetThreshold
            onsetTimeTotal = onsetTimeTotal-1;
        end
        onsetTime(event) = onsetTimeTotal+1;
    catch
        onsetTime(event) = 1;
    end
    
    %find offset
    %first look for the first value that is half the size of the peak
    offsetHalfMax = peakTime(event)+1;
    try
        while  expTrace(offsetHalfMax) > halfMax
            offsetHalfMax = offsetHalfMax+1;
        end
        %then look for how long from there the derivative is above threshold
        offsetTimeTotal = offsetHalfMax;
        while derivatives(offsetTimeTotal) < offsetThreshold
            offsetTimeTotal = offsetTimeTotal+1;
        end
        offsetTime(event) = offsetTimeTotal-1;
    catch
        offsetTime(event) = length(expTrace);
    end
end

%% 4.) Separate multiple events
% some events are close by each other and share the same on- and offset -
% separate them into individual events
uniqOnset = unique(onsetTime);
multiOnsets = uniqOnset((1<histc(onsetTime,unique(onsetTime))));
for event = 1:length(multiOnsets)
    onsetEvent = find(onsetTime == multiOnsets(event));
    trace = expTrace(onsetTime(onsetEvent(1)):offsetTime(onsetEvent(1)));
    deriv = derivatives(onsetTime(onsetEvent(1)):offsetTime(onsetEvent(1)));
    startFrame = onsetTime(onsetEvent(1));
    
    localMinsDer = find(islocalmax(deriv));
    peaksFound = peakTime(onsetEvent)-startFrame;
    
    for doubleEvents = 1:length(peaksFound)-1
        %how many localMins are between the events?
        borders = find(localMinsDer > peaksFound(doubleEvents) & localMinsDer < peaksFound(doubleEvents+1));
        %if there is exactly one, take that as the border the closest one
        %as on/offsets, else use the closest ones
        if ~isempty(borders) %if you do not find borders, leave it as it is
            if length(borders) == 1
                offsetTime(onsetEvent(doubleEvents)) = localMinsDer(borders)+startFrame;
                onsetTime(onsetEvent(doubleEvents+1)) = localMinsDer(borders)+startFrame;
            else
                %closet border to peak left = offset of that peak
                [~, closBorder] = min(abs(peaksFound(doubleEvents)-localMinsDer(borders)));
                offsetTime(onsetEvent(doubleEvents)) = localMinsDer(borders(closBorder))+startFrame;
                %closest border to peak right = onset of that peak
                [~, closBorder] = min(abs(peaksFound(doubleEvents+1)-localMinsDer(borders)));
                onsetTime(onsetEvent(doubleEvents+1)) = localMinsDer(borders(closBorder))+startFrame;
            end
        end
    end
end

%% 4.) Write everything into the analysis structure
analysis.(field).roi(ROI).eventOnsets = onsetTime;
analysis.(field).roi(ROI).eventPeaks = peakTime;
analysis.(field).roi(ROI).eventOffsets = offsetTime;
analysis.(field).roi(ROI).duration = offsetTime - onsetTime;
analysis.(field).roi(ROI).meanPeak = mean(peakTime);
analysis.(field).roi(ROI).meanDuration = mean(analysis.(field).roi(ROI).duration);
analysis.(field).roi(ROI).numEvents = length(peakTime);

for i = 1:length(peakTime)
    startFrame = max(1,peakTime(i)-49);
    endFrame = min(size(expTrace,1),peakTime(i)+49);
    analysis.(field).roi(ROI).events(i).trace = expTrace(startFrame:endFrame);
    analysis.(field).roi(ROI).events(i).deriv = derivatives(startFrame:endFrame);
    analysis.(field).roi(ROI).events(i).peakAmp = expTrace(peakTime(i));
end