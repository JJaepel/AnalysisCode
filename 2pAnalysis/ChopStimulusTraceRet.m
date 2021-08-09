function [analysis, metadata, data] = ChopStimulusTraceRet(analysis,metadata,data,level, field, varargin)
    if mod(length(varargin),2) ~= 0
        error('You need to specify both the field and the argument')
    elseif strcmp(metadata.StimParams.type, 'Spontaneous')
        warning('Unable to chop up Spontaneous activity into meaningful chunks')
        return
    else
        % defaults for values
        windowStart= 0;  % choose the entire stim period
        windowStop = 0.5;    % choose the entire stim period
        tfield = field; % default to field name as analysis name
        pre=0;
        post=1;
        skipTrial=0;
        baseline=true;
        for i = 1:2:length(varargin)
            switch varargin{i}
                case 'analysisName'
                    tfield = varargin{i+1};
                case 'windowStart'
                    windowStart= varargin{i+1};
                case 'windowStop'
                    windowStop=varargin{i+1};
                case 'pre'
                    pre = varargin{i+1};
                case 'post'
                    post = varargin{i+1};
                case 'skipTrial'
                    skipTrial = varargin{i+1};
                case 'baseline'
                    baseline=varargin{i+1};
                otherwise
                    error(['You have specified an invalid field: ', varargin{i}])
            end
        end
    end
    if level == 1
        metadata.TwoPhoton.time = metadata.TwoPhoton.time(1:5:end);
        metadata.TwoPhoton.rate = 1/median(diff(metadata.TwoPhoton.time));
    end
    TwoPhotonRate= metadata.TwoPhoton.rate;
    analysis.(tfield)=struct;
    analysis.(tfield).preTrialTime = pre;
    analysis.(tfield).postTrialTime = post;
    offsetPre = round(metadata.TwoPhoton.rate * analysis.(tfield).preTrialTime);
    offsetPost = round(metadata.TwoPhoton.rate * analysis.(tfield).postTrialTime);
    analysis.(tfield).stimStart=offsetPre+1;
    analysis.(tfield).stimStop=round(analysis.(tfield).stimStart+metadata.TwoPhoton.rate * metadata.StimParams.stimDuration());
    analysis.(tfield).roi = struct;

    windowStartIdx= round(windowStart * metadata.TwoPhoton.rate +offsetPre+1);
    windowStopIdx = round(windowStop * TwoPhotonRate + offsetPost);
    if windowStartIdx == windowStopIdx
        windowStopIdx = analysis(tfield).stimStop;
    end
    analysis.(tfield).windowStart =windowStartIdx;
    analysis.(tfield).windowStop = windowStopIdx;

    %do a little clean up here to make stim numbers match
%     while (metadata.StimParams.StimOnTimes(2,metadata.StimParams.numberOfStims)+metadata.StimParams.stimDuration > metadata.TwoPhoton.time(end))
%         if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ==0
%             metadata.StimParams.numTrials = metadata.StimParams.numTrials -1;
%         else
%             metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
%         end
%         metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
%     end

    Stimtype = metadata.StimParams.type;
    switch Stimtype
        case 'Patch'
            numVisStimsPerTrial = metadata.StimParams.numAzimuth * metadata.StimParams.numElevation * 2;
            numBlanksPerTrial = floor(metadata.StimParams.numStimAzi * metadata.StimParams.numStimElev * metadata.StimParams.blankpercent *2);
        case 'RetWedge'
            numVisStimsPerTrial = metadata.StimParams.numWedges;
            numBlanksPerTrial = metadata.StimParams.nBlank;
        case 'Ret_Annulus'
            numVisStimsPerTrial = metadata.StimParams.numSizes;
            numBlanksPerTrial = metadata.StimParams.nBlank;
        case 'Retinotopy_2D'
            numVisStimsPerTrial = metadata.StimParams.numAzimuth * metadata.StimParams.numElevation;
            numBlanksPerTrial = 0;
            metadata.StimParams.uniqStimIds = metadata.StimParams.uniqStimIds(2:end);
        case 'PatchGrating'
            numVisStimsPerTrial = metadata.StimParams.numAzimuth * metadata.StimParams.numElevation;
            numBlanksPerTrial = floor(metadata.StimParams.numAzimuth * metadata.StimParams.numElevation * metadata.StimParams.blankpercent) * metadata.StimParams.doBlank;
        otherwise
            try
                numBlanksPerTrial = floor(metadata.StimParams.numAzimuth * metadata.StimParams.numElevation * metadata.StimParams.blankpercent);
            catch
                numBlanksPerTrial = metadata.StimParams.doBlank;
            end
            numVisStimsPerTrial = metadata.StimParams.numAzimuth * metadata.StimParams.numElevation;
    end
    
    numBlanks = numBlanksPerTrial * metadata.StimParams.numTrials;
    numVisStims = numVisStimsPerTrial * metadata.StimParams.numTrials;
    
    if (numBlanks+numVisStims) > metadata.StimParams.numberOfStims
        metadata.StimParams.numTrials = metadata.StimParams.numTrials-1;
        numBlanks = numBlanksPerTrial * metadata.StimParams.numTrials;
        numVisStims = numVisStimsPerTrial * metadata.StimParams.numTrials;
    elseif max(metadata.TwoPhoton.time) < max(metadata.StimParams.StimOnTimes(2,:))
        metadata.StimParams.numTrials = metadata.StimParams.numTrials-1;
        numBlanks = numBlanksPerTrial * metadata.StimParams.numTrials;
        numVisStims = numVisStimsPerTrial * metadata.StimParams.numTrials;
    end
    
    try 
        for i =1:size(metadata.StimParams.StimOnTimes,2)
            metadata.StimParams.stimStartIndex(i)= find(metadata.TwoPhoton.time > metadata.StimParams.StimOnTimes(2,i),1);
            metadata.StimParams.stimStopIndex(i)= floor(metadata.StimParams.stimStartIndex(i)+ metadata.StimParams.stimDuration* metadata.TwoPhoton.rate);
        end
    catch
        for i =1:(numBlanks+numVisStims)
            metadata.StimParams.stimStartIndex(i)= find(metadata.TwoPhoton.time > metadata.StimParams.StimOnTimes(2,i),1);
            metadata.StimParams.stimStopIndex(i)= floor(metadata.StimParams.stimStartIndex(i)+ metadata.StimParams.stimDuration* metadata.TwoPhoton.rate);
        end
    end

    stimPeriod=(analysis.(tfield).stimStart:analysis.(tfield).stimStop);
    analysisPeriod=(analysis.(tfield).windowStart:analysis.(tfield).windowStop);

    %correct for our specified window

    if windowStartIdx > max(stimPeriod) || windowStopIdx > max(stimPeriod)
        warning('You have choosen a period to average over that is longer than the stim period')
    end

    stimStarts=metadata.StimParams.stimStartIndex;
    stimStops=metadata.StimParams.stimStopIndex;

    for stimID= 1:size(metadata.StimParams.uniqStimIds,2)
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        StimonTimes= metadata.StimParams.StimOnTimes;
        StimonTimes=StimonTimes(:,1:metadata.StimParams.numberOfStims);
        stimIndices = find(StimonTimes(1,:)==stimulus);
        for trialNumber= 1:metadata.StimParams.numTrials
            stimIndex = stimIndices(trialNumber);
            for i=1:length(data.roi)
                selectedFramesTrace = (stimStarts(stimIndex)-offsetPre):(stimStops(stimIndex)+offsetPost);
                if max(selectedFramesTrace) <= length(data.roi(i).(field))
                    analysis.(tfield).roi(i).stimResponseTrace(stimID,trialNumber,:)= data.roi(i).(field)(selectedFramesTrace);
                else
                    analysis.(tfield).roi(i).stimResponseTrace(stimID,trialNumber,:) = zeros(1,1,length(selectedFramesTrace));
                end
            end
        end

        for i=1:length(data.roi)
            if skipTrial >0 && skipTrial <=size(analysis.(tfield).roi(i).stimResponseTrace,2)
                analysis.(tfield).roi(i).stimResponseTrace(:,skipTrial,:)=[];
            end
            traces=squeeze(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:));
            bltraces=repmat(mean(traces(:,1:(stimPeriod(1)-1)),2)', size(traces,2),1)';

            if baseline
                traces=traces- bltraces;
            end
            try
                analysis.(tfield).roi(i).stimResponseTraceRaw(stimID,:,:)=analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:);
            catch ME
                ME.message
            end
            analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:)= traces;
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = mean(traces, 1);
            n= size(analysis.(tfield).roi(i).stimResponseTrace,2);
            y=analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:);
            analysis.(tfield).roi(i).SEMResponseTrace(stimID,:) = std(y,[],2)/sqrt(n);
            analysis.(tfield).roi(i).avgStimResponse(stimID,:) = mean(analysis.(tfield).roi(i).avgResponseTrace(stimID,analysisPeriod),2);
            n= size(analysis.(tfield).roi(i).stimResponseTrace,2);
            y=analysis.(tfield).roi(i).avgStimResponse(stimID,:);
            analysis.(tfield).SEMCellResponses(i,stimID)= std(y,[],2)/sqrt(n);
            analysis.(tfield).avgCellResponses(i,stimID)= mean(analysis.(tfield).roi(i).avgStimResponse(stimID,:),2);
        end
    end

    for i = 1:length(data.roi)
        for stimID = 1:metadata.StimParams.uniqStims-1
            analysis.(tfield).roi(i).avgStimResponse(stimID) = analysis.(tfield).roi(i).avgStimResponse(stimID) -analysis.(tfield).roi(i).avgStimResponse(end) ;
            analysis.(tfield).roi(i).avgResponseTraceRaw(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:);
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:) - analysis.(tfield).roi(i).avgResponseTrace(end,:);
            analysis.(tfield).roi(i).Fb(stimID,:)= mean(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,1:analysisPeriod(1)-1),3);

        end
    end
    analysis.(tfield).preTrialTime=pre;
end    