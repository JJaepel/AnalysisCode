function data = FitCellsToDendrite(analysisParams, metadata, data)

% cutting up data based on stimIDs and times
scanPeriod = 1/metadata.TwoPhoton.rate;
if analysisParams.level
    scanPeriod = scanPeriod*5;
    TwoPhotonTimes = metadata.TwoPhoton.time(1:5:end);
else 
    TwoPhotonTimes = metadata.TwoPhoton.time;
end
preStimPeriod = floor(analysisParams.pre /scanPeriod);
stimDur = floor(metadata.StimParams.stimDuration / scanPeriod);
analysisParams.post = 0.5;
postPeriod = floor(analysisParams.post / scanPeriod);
analysisLength = preStimPeriod + stimDur + postPeriod;

for cc = 1:length(data.roi)
    data.roi(cc).cyc = zeros(metadata.StimParams.uniqStims, metadata.StimParams.numTrials, analysisLength); 
    trialList = zeros(1, metadata.StimParams.uniqStims);

    for i = 1:metadata.StimParams.numberOfStims
        stimTimeStart = find(TwoPhotonTimes > metadata.StimParams.StimOnTimes(2,i),1);
        windowStart = stimTimeStart-preStimPeriod+1;
        windowStop = stimTimeStart + stimDur + postPeriod;
        stimTime = windowStart:windowStop;
        ind = find(metadata.StimParams.uniqStimIds == metadata.StimParams.StimOnTimes(1,i));
        trialList(ind) = trialList(ind)+1;
        f = data.roi(cc).dff(stimTime);
        data.roi(cc).cyc(ind, trialList(ind),:) = f;
    end
    
    data.roi(cc).rawRes = data.roi(cc).rawF;
    data.roi(cc).cycRes = data.roi(cc).cyc;
    data.roi(cc).slope = [];
    data.roi(cc).corr = -1;
end