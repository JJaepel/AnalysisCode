function [PrefAngle, curveSum, eventSum] = calcSumSpineDirResponses(Spines, IDs)


%create empty arrays
allResp = NaN(size(Spines(IDs(1)).funcData.meanResp,2), length(IDs));
allEvents = allResp;
allPrefDeg = NaN(length(IDs),1);

for i = 1:length(IDs)
    if Spines(IDs(i)).funcData.OSI > 0.1
        allPrefDeg(i) = Spines(IDs(i)).funcData.prefDir;

        meanResp = Spines(IDs(i)).funcData.meanResp;
        %meanRespNorm = (meanResp-min(meanResp))/(max(meanResp)-min(meanResp)); %normalize the resp
        meanRespNorm = (meanResp)/(max(meanResp)); %normalize the resp
        allResp(:,i) = meanRespNorm; %add it to the array

        eventsAllTrials = sum(Spines(IDs(i)).funcData.events,2);%sum over all trials
        eventsFolded = eventsAllTrials';
        %eventsNorm = (eventsFolded-min(eventsFolded))/(max(eventsFolded)-min(eventsFolded)); %normalize the events
        eventsNorm = (eventsFolded)/(max(eventsFolded)); %normalize the events
        allEvents(:,i) = eventsNorm; %add it to the array
    end
end
allPrefAng = deg2rad(allPrefDeg);
PrefAngle = rad2deg(circ_mean(allPrefAng));

if PrefAngle < 0
    PrefAngle = 90+PrefAngle;
end

meanInputResp = nanmean(allResp,2); %get the mean resp of all inputs
%curveSum = (meanInputResp-min(meanInputResp))/(max(meanInputResp)-min(meanInputResp)); %normalize the resp
curveSum = (meanInputResp)/(max(meanInputResp)); %normalize the resp

meanInputEvents = nanmean(allEvents,2); %get the mean events of all inputs
%eventSum = (meanInputEvents-min(meanInputEvents))/(max(meanInputEvents)-min(meanInputEvents)); %normalize the events
eventSum = (meanInputEvents)/(max(meanInputEvents)); %normalize the events