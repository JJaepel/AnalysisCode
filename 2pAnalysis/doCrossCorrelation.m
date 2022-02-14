function doCrossCorrelation(events,AreaTrace,secAreaTrace,metadata,modus, saveDir, areaName)

%calculates how events in one are cross-correlated with the trace from the
%other area
%inputs
% - events: structure, containing timing information for events
% - AreaTrace: overall trace of the area from which the events are taken
% - secAreaTrace: overall trace of the other area that we want to match
% - metadata: contains timing information
% - modus: do you want to detect based on onset or peak of event?
% - saveDir: where do you want to save it?
% - areaName: events from which area are you correlating?

verbose = 0;

secWindow = ceil(1 / mean(diff(metadata.Imaging.time))); %take a 1 second time window
time = mean(diff(metadata.Imaging.time));

%pre-allocate variables
corrAll = zeros(length(events),4*secWindow+1);
firstTraceAll = zeros(length(events),2*secWindow+1);
secondTraceAll = zeros(length(events),2*secWindow+1);

lengthTrace = size(firstTraceAll,2);
for i = 1:length(events)
    switch modus
        case 'onset' %if you correlate based on the onset of the event, do the correlation starting 1 sec
            startFrame = round(max(1,events(i).onset-.5*secWindow));
            endFrame = startFrame+secWindow*2;
        case 'peak' %if you correlate based on the peak of the event, do the correlation starting 2 sec before peak 
            startFrame = max(1,events(i).peakTime-1*secWindow);
            endFrame = startFrame+secWindow*2;
    end
    firstTrace = AreaTrace(startFrame:endFrame);
    secondTrace = secAreaTrace(startFrame:endFrame);
    [corTrace,lag] = xcorr(firstTrace,secondTrace);%,'normalized');
    [~, maxCorrTrace] = max(corTrace);
    timeLagBetweenAreas(i) = lag(maxCorrTrace)*time;
   
    if length(firstTrace) < lengthTrace
        zerosToAdd = zeros(lengthTrace-length(firstTrace),1);
        firstTrace = [firstTrace zerosToAdd];
        secondTrace = [secondTrace zerosToAdd];
        corTrace = [zerosToAdd corTrace zerosToAdd];
    end
    
    %save all
    corrAll(i,:) =squeeze(corTrace);
    firstTraceAll(i,:) =squeeze(firstTrace);
    secondTraceAll(i,:) =squeeze(secondTrace);
end

lagInTime = lag*time;
traceTime = linspace(0,lengthTrace*time, lengthTrace);

%plot all correlations with the corresponding traces
if verbose
    for eventNr = 1:length(events)
        figure;
        subplot(2,1,1)
        plot(traceTime, firstTraceAll(eventNr,:), 'b')
        hold on
        plot(traceTime, secondTraceAll(eventNr,:), 'r')
        xlabel('time in s')
        ylabel('normalized amplitude')
        title([modus 'Correlation Events ' num2str(eventNr) ' in ' areaName])

        subplot(2,1,2)
        plot(lagInTime,corrAll(eventNr,:), 'g')
        xlabel('time in s')
        ylabel('correlation')
        set(gcf, 'color', 'w');
        saveas(gcf, fullfile(saveDir, [modus 'Correlation Events ' num2str(eventNr) ' in ' areaName]))
    end
end


%plot the mean correlation
figure;
errorbar(lagInTime,mean(corrAll,1),std(corrAll,0,1))
xlabel('time in s')
ylabel('correlation')
set(gcf, 'color', 'w');
title(['Mean ' modus ' Correlation: Events in ' areaName])
saveas(gcf, fullfile(saveDir, ['Mean ' modus 'Correlation_Events in ' areaName]))

%plot all correlations
figure
for event = 1:length(events)
    plot(lagInTime,corrAll(event,:));
    hold all
end
xlabel('time in s')
ylabel('correlation')
set(gcf, 'color', 'w');
title(['All ' modus ' Correlation: Events in ' areaName])
saveas(gcf, fullfile(saveDir, [modus 'Correlation_Events in ' areaName]))

%make histogram of max correlation timing 
figure
histogram(timeLagBetweenAreas, lengthTrace*2, 'normalization', 'probability')
xlabel('timelag in s')
ylabel('probability')
set(gcf, 'color', 'w');
title([modus ' Correlation: Timelag between events in ' areaName])
saveas(gcf, fullfile(saveDir, [modus 'Correlation_Timelag between events in ' areaName]))
