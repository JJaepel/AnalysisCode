function events = findEvents(expTrace,AreaName, verbose, saveDir)

if nargin < 2
    AreaName = '';
end
if nargin < 3
    verbose = 0;
end
if nargin < 4
    saveDir = cd;
end

%% 1.) Find peaks and threshold them to be above the average
[peaksValue,peakTime] = findpeaks(smooth(smooth(expTrace)));
threshold = mean(expTrace(:));
peakThresholded = find(peaksValue>threshold);
peakTime = peakTime(peakThresholded);

peakTimes = zeros(length(expTrace),1);
peakTimes(peakTime) = 1;

[peaksValue,peakTime] = findpeaks(smooth(smooth(expTrace)));
threshold = mean(expTrace(:));
peakThresholded = find(peaksValue>threshold);
peakTime = peakTime(peakThresholded);

peakTimes = zeros(length(expTrace),1);
peakTimes(peakTime) = 1;

%% 2.) Find event onset and offset
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
    while  expTrace(onsetHalfMax) > halfMax
        onsetHalfMax = onsetHalfMax-1;
    end
    %then look for how long from there the derivativ is positive
    onsetTimeTotal = onsetHalfMax;
    while derivatives(onsetTimeTotal) > onsetThreshold
        onsetTimeTotal = onsetTimeTotal-1;
    end
    onsetTime(event) = onsetTimeTotal+1;
    
    %find offset
    %first look for the first value that is half the size of the peak
    offsetHalfMax = peakTime(event)+1;
    while  expTrace(offsetHalfMax) > halfMax
        offsetHalfMax = offsetHalfMax+1;
    end
    %then look for how long from there the derivative is above threshold
    offsetTimeTotal = offsetHalfMax;
    while derivatives(offsetTimeTotal) < offsetThreshold
        offsetTimeTotal = offsetTimeTotal+1;
    end
    offsetTime(event) = offsetTimeTotal-1;
end

%% 3.) Separate multiple events
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

%% 4.) Write everything into a structure and, if wanted, check the output
events=struct;
for i=1:length(peakTime)
    events(i).onset = onsetTime(i);
    events(i).peakTime = peakTime(i);
    events(i).offset = offsetTime(i);
    
    startFrame = max(1,peakTime(i)-49);
    endFrame = min(size(expTrace,1),peakTime(i)+49);
    events(i).trace = expTrace(startFrame:endFrame);
    events(i).deriv = derivatives(startFrame:endFrame);
    
    if verbose
        figure
        subplot(2,1,1)
        ymin= min(events(i).trace); ymax = max(events(i).trace);
        window = (events(i).onset-startFrame+1):1:(events(i).offset-startFrame+1);
        patch([window fliplr(window)], [ymin*ones(1,length(window)) ymax*ones(1,length(window))], [1 .9 .9], 'LineStyle', 'none');
        hold all
        plot(events(i).trace)
        hold all
        plot(events(i).peakTime-startFrame+1,events(i).trace(events(i).peakTime-startFrame),'b*') 
        xlabel('frame')
        ylabel('normalized trace')
        title(['Event Nr. ' num2str(i)])
        subplot(2,1,2)
        plot(events(i).deriv)
        hold all
        plot(onsetTime(i)-startFrame+1,events(i).deriv(onsetTime(i)-startFrame+1),'r*')
        plot(offsetTime(i)-startFrame+1,events(i).deriv(offsetTime(i)-startFrame+1),'r*')
        xlabel('frame')
        ylabel('normalized derivative')
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDir, [AreaName ' Event Nr ' num2str(i)]))
    end
end
