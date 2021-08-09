function [analysis, metadata] = ChopStimulusTraceEpi(analysis,metadata,data,field)
    offsetPre  = round(metadata.Imaging.rate * metadata.StimParams.isi/2);
    offsetPost = round(metadata.Imaging.rate * metadata.StimParams.isi/2);
    stimDuration = ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration);
    blPeriod   = (1:offsetPre);

    % Collect the stimulus evoked traces
    imgStack = data.(field);
    [nRows,nCols,~] = size(imgStack);
    analysisPeriod = -offsetPre:(stimDuration+offsetPost);
    nFrames        = length(analysisPeriod);
    stimResponseTrace = zeros([metadata.StimParams.uniqStims,metadata.StimParams.numTrials-1,nRows,nCols,nFrames,],'single');
    for stimID= 1:metadata.StimParams.uniqStims
        % Get stim ID
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        disp(strcat('Processing StimID: ', num2str(stimulus)))
        stimIndices = find(metadata.StimParams.StimOnTimes(1,:)==stimulus);

        % Loop through each trial for the associated stim ID
        for trialNumber= 1:metadata.StimParams.numTrials-1
            startFrame = metadata.StimParams.stimStartIndex(stimIndices(trialNumber));
            stimResponseTrace(stimID,trialNumber,:,:,:) = single(data.(field)(:,:,startFrame+analysisPeriod));
        end
    end
    
    disp('Pre-Stimulus Blank')
    for stimID = 1:metadata.StimParams.uniqStims
        for trial = 1:metadata.StimParams.numTrials-1
            blmean = squeeze(nanmean(stimResponseTrace(stimID, trial, :,:, blPeriod), 5));
            for frameNum=1:nFrames
                stimResponseTrace(stimID, trial, :,:, frameNum)=(squeeze(stimResponseTrace(stimID, trial, :,:, frameNum))-blmean)./blmean;
            end
        end
    end
    
    analysis.(field).roi.stimResponseTrace = permute(stimResponseTrace,[5 2 1 3 4]);
    %defaultOrder = {'x','y','t','c','n'}; % DO NOT CHANGE this variable!!! This is the way that the code originally structures the stimResponseTrace
                
end 