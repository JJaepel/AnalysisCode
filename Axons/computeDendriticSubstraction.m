function data = computeDendriticSubstraction(analysisParams, metadata, data)

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
    
    
    if metadata.StimParams.StimOnTimes(2,end) > TwoPhotonTimes(end)
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
    else
       numTrial= floor(size(metadata.StimParams.StimOnTimes,2)/metadata.StimParams.uniqStims);
       if (metadata.StimParams.StimOnTimes(2,metadata.StimParams.uniqStims*numTrial)+ (stimDur + postPeriod)/metadata.TwoPhoton.rate) > TwoPhotonTimes(end)
           numTrial = numTrial-1;
       end
       for i = 1:metadata.StimParams.uniqStims*numTrial
            stimTimeStart = find(TwoPhotonTimes > metadata.StimParams.StimOnTimes(2,i),1);
            windowStart = stimTimeStart-preStimPeriod+1;
            windowStop = stimTimeStart + stimDur + postPeriod;
            stimTime = windowStart:windowStop;
            ind = find(metadata.StimParams.uniqStimIds == metadata.StimParams.StimOnTimes(1,i));
            trialList(ind) = trialList(ind)+1;
            f = data.roi(cc).dff(stimTime);
            data.roi(cc).cyc(ind, trialList(ind),:) = f;
        end
    end
end

for cc = 1:length(data.roi)
    if data.roi(cc).spine
        %find dendrite
        dendPoint = cc;
        while ~data.roi(dendPoint).dendrite
            dendPoint = dendPoint+1;
        end
        SpDff = data.roi(cc).dff;
        SpDff = SpDff(1:round(length(SpDff)*.9));
        SpDff(isinf(SpDff)) = 0;
        SpDff_sub = SpDff(SpDff < nanmedian(SpDff)+abs(min(SpDff)));
        %noiseM = nanmedian(SpDff_sub); %should be near 0
        noiseSD = nanstd(SpDff_sub);
        
        slope = robustfit(data.roi(dendPoint).cyc(:), data.roi(cc).cyc(:));
        
        % dendritic scalar applied (from robust fit) to cyc and
        % substraction
        data.roi(cc).slope = slope(2);
        data.roi(cc).cycRes = data.roi(cc).cyc - slope(2).*data.roi(dendPoint).cyc;
        data.roi(cc).cycRes(data.roi(cc).cycRes < 0) = 0;
        
        % calculate correlation
        rSp = data.roi(cc).dff;
        rDn = data.roi(dendPoint).dff;
        rSp = rSp - slope(2).*rDn;
        rSp(rSp < -noiseSD) = -noiseSD;
        rSp(isinf(rSp)) = 0;
        rDn(isinf(rDn)) = 0;
        data.roi(cc).rawRes = rSp;
        rSp(rSp <= 0) = nan;
        rSp(rSp <= 0) = nan;
        [r, ~] = corrcoef(rSp, rDn, 'rows', 'pairwise');
        data.roi(cc).corr = r(2);
    else
        data.roi(cc).rawRes = data.roi(cc).rawF;
        data.roi(cc).cycRes = data.roi(cc).cyc;
        data.roi(cc).slope = [];
        data.roi(cc).corr = -1;
    end
end