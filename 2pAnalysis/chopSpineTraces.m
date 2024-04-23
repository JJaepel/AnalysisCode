function ce = chopSpineTraces(ce, metadata)

% Uses the metadata to chop the traces into chunks according to the
% stimulus that was show at that time
%
% Input: 
% - ce: structure containing all the ROIs with anatomical information, raw
% trace and dff trace
% - metadat: structure containing stimulus information 
%
% Output:
% - ce: now updated with cyc, which cotains the trace split up in uniqStims
% x Trials x timeWindow

%% Step 0: Get the stimIDs
stimID = metadata(1).copyStimID;

%% Step 1: Go through all ROIs
for cc = 1:length(ce)
    %prepare the structures
    ce(cc).cyc = zeros(length(metadata(1).uniqStims),metadata(1).ntrials,metadata(1).stimDur2+metadata(1).postPeriod2+metadata(1).prestimPeriod2);
    trialList = zeros(1,length(metadata(1).uniqStims));
    
    %for each stim
    for ii = 1:metadata(1).numStims
        %find the stim wnindow
        stimTime2 = metadata(1).stimOn2pFrame(ii)+1-metadata(1).prestimPeriod2:metadata(1).stimOn2pFrame(ii)+metadata(1).stimDur2+metadata(1).postPeriod2;
        %which of the stims was shown
        ind = find(metadata(1).uniqStims==stimID(ii));
        %which trials are wee looking at
        trialList(ind) = trialList(ind)+1;
        f = ce(cc).dff(stimTime2);
        %save it in there
        ce(cc).cyc(ind,trialList(ind),:) = f;% - median(f(1:2));
    end
end