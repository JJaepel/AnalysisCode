function eventsArea = findEventSequence(eventsArea,eventsOtherArea, metadata, window)

% classify events as having simultaneous onset, preceding one event or
% following another event
% 
% inputs
% - eventsArea: structure containing all events in one area
% - eventsOtherArea: structure containing all events in the neighboring
% - metadata: containing information about timing
% area
% output:
% - eventsArea: structure containing all events in one area, adding one
% column in the structure that classifies event as 1) preceding other area,
% 2) simultaneous, 3) followig an event in the other area, 0) no event
% nearby

minFrameDifference = ceil(window / mean(diff(metadata.Imaging.time))); %event should start within half a second
otherAreaOnsets = [eventsOtherArea.onset];
if window < 1
    field = ['class0' num2str(window*10)];
else
    field = ['class' num2str(window*10)];
end

for i = 1:length(eventsArea)
    eventOnset = eventsArea(i).onset;
    [closesEventOtherArea, ind] = min(abs(otherAreaOnsets-eventOnset));
    if closesEventOtherArea == 0
        eventsArea(i).(field) = 1;
    elseif abs(closesEventOtherArea) < minFrameDifference
        eventDifference = eventsArea(i).onset - eventsOtherArea(ind).onset;        
        if eventDifference < 0
            eventsArea(i).(field) = 2;
        else
            eventsArea(i).(field) = 3;
        end
    else
        eventsArea(i).(field) = 0;
    end
end