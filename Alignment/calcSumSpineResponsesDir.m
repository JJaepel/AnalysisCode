function [PrefAngleDir, curveSumDir eventSumDir] = calcSumSpineResponsesDir(Spines, IDs)


%create empty arrays
allResp = zeros(size(Spines(IDs(1)).funcData.meanResp,2), length(IDs));
allEvents = allResp;
allPrefDeg = zeros(length(IDs),1);

for i = 1:length(IDs)
    allPrefDeg(i) = Spines(IDs(i)).funcData.prefDir;
    
    meanResp = Spines(IDs(i)).funcData.meanResp;
    %meanRespNorm = (meanResp-min(meanResp))/(max(meanResp)-min(meanResp)); %normalize the resp
    meanRespNorm = meanResp/(max(meanResp));
    allResp(:,i) = meanRespNorm; %add it to the array
    
    eventsAllTrials = sum(Spines(IDs(i)).funcData.events,2)';%sum over all trials
    %eventsNorm = (eventsAllTrials-min(eventsAllTrials))/(max(eventsAllTrials)-min(eventsAllTrials)); %normalize the events
    eventsNorm = eventsAllTrials/(max(eventsAllTrials));
    allEvents(:,i) = eventsNorm; %add it to the array
end
allPrefAng = deg2rad(allPrefDeg);
PrefAngleDir = rad2deg(circ_mean(allPrefAng));

if PrefAngleDir < 0
    PrefAngleDir = 180+PrefAngleDir;
end

meanInputResp = mean(allResp,2); %get the mean resp of all inputs
curveSumDir = (meanInputResp-min(meanInputResp))/(max(meanInputResp)-min(meanInputResp)); %normalize the resp

meanInputEvents = mean(allEvents,2); %get the mean events of all inputs
eventSumDir = (meanInputEvents-min(meanInputEvents))/(max(meanInputEvents)-min(meanInputEvents)); %normalize the events