function eventData = findLocalEvents(AreaData,eventData,AreaName,verbose, saveDir)

if nargin <3 
    AreaName = '';
end
if nargin < 4
    verbose = 0;
end
if nargin < 5
    saveDir = cd;
end

for event = 1:length(eventData)
    % 1.) Find local trace and derivatived based on event location
    locationMask = eventData(event).peakTimeArea;
    trace = getAreaTrace(AreaData,locationMask);
    startFrame = max(1,eventData(event).peakTime-49);
    endFrame = min(size(trace,1),eventData(event).peakTime+49);
    
    dydx = gradient(trace);
    derivatives = zscore(dydx);
    eventData(event).localTrace = trace(startFrame:endFrame);
    eventData(event).localDeriv = derivatives(startFrame:endFrame);
    
    % 2.) find local peak
    [peakAmp,peakPos] = findpeaks(smooth(smooth(eventData(event).localTrace)));
    threshold = mean(eventData(event).localTrace(:));
    peakThresholded = find(peakAmp>threshold);
    peakPos = peakPos(peakThresholded);
    peakAmp = peakAmp(peakThresholded);
    
    %if there are multiple peaks, find the one closest to the whole Area
    %trace peak
    if length(peakThresholded) > 1
        [~,closestPeak] = min(abs(peakPos-(eventData(event).peakTime-startFrame)));
        peakPos = peakPos(closestPeak);
        peakAmp = peakAmp(closestPeak);
    end
    eventData(event).localPeak = peakPos+startFrame;
        
    % 3.)find onset
    %first look for the first value that is half the size of the peak, but
    %at least higher than one std
    halfMax = peakAmp*0.5;
    halfMax = max(halfMax,std(trace(:))); 
    onsetHalfMax =  peakPos-1;
    while  eventData(event).localTrace(onsetHalfMax) > halfMax
        onsetHalfMax = onsetHalfMax-1;
    end
    
    %then look for how long from there the derivativ is positive
    onsetThreshold = std(derivatives);  
    onsetTimeTotal = onsetHalfMax;
    try
        while eventData(event).localDeriv(onsetTimeTotal) > onsetThreshold
            onsetTimeTotal = onsetTimeTotal-1;
        end
    catch
        onsetTimeTotal = 0;
    end
    onsetTimeTotal = onsetTimeTotal+1;
    eventData(event).localOnset = onsetTimeTotal+startFrame;
    
    % 4.) find offset
    %first look for the first value that is half the size of the peak
    offsetHalfMax = peakPos+1;
    while   eventData(event).localTrace(offsetHalfMax) > halfMax
        offsetHalfMax = offsetHalfMax+1;
    end
    %then look for how long from there the derivative is above threshold
    offsetThreshold = -0.5*std(derivatives);
    offsetTimeTotal = offsetHalfMax;
    while eventData(event).localDeriv(onsetTimeTotal)< offsetThreshold
        offsetTimeTotal = offsetTimeTotal+1;
    end
    offsetTimeTotal = offsetTimeTotal-1;
    eventData(event).localOffset = offsetTimeTotal+startFrame;
    
    if verbose
        figure
        subplot(2,1,1)
        ymin= min(eventData(event).localTrace); ymax = max(eventData(event).localTrace);
        window = (onsetTimeTotal+1):1:(offsetTimeTotal+1);
        patch([window fliplr(window)], [ymin*ones(1,length(window)) ymax*ones(1,length(window))], [1 .9 .9], 'LineStyle', 'none');
        hold all
        plot(eventData(event).localTrace)
        hold all
        plot(peakPos,peakAmp,'b*') 
        xlabel('frame')
        ylabel('local trace')
        title(['Event Nr. ' num2str(event)])
        subplot(2,1,2)
        plot(eventData(event).localDeriv)
        hold all
        plot(onsetTimeTotal,eventData(event).localDeriv(onsetTimeTotal),'r*')
        plot(offsetTimeTotal,eventData(event).localDeriv(offsetTimeTotal),'r*')
        xlabel('frame')
        ylabel('local derivative')
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDir, [AreaName ' Local Event Nr ' num2str(event)]))
    end
end
