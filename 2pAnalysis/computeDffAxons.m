function data = computeDffAxons(metadata, data)
%field source is usually 'rawF', fieldTarge = 'dff'


%% get all baseline data together

%find start and end frames of all baseline times
TimeBaselineStart = metadata.StimParams.StimOnTimes(2,:)-1-metadata.StimParams.isi;
TimeBaselineEnd = metadata.StimParams.StimOnTimes(2,:)-1;
TimeStimEnd = metadata.StimParams.StimOnTimes(2,:)+metadata.StimParams.stimDuration;

for base = 1:length(TimeBaselineEnd)
    FrameBaselineStart(base) = find(metadata.TwoPhoton.time > TimeBaselineStart(base),1);
    FrameBaselineEnd(base) = find(metadata.TwoPhoton.time > TimeBaselineEnd(base),1);
    FrameStimEnd(base) = find(metadata.TwoPhoton.time > TimeStimEnd(base),1);
end
recLength = length(data.roi(1).rawF);

for RoiNR = 1:length(data.roi)
    inputTrace = data.roi(RoiNR).rawF;

    for Stim = 1:length(TimeBaselineEnd)
        F0(Stim) = nanmean(inputTrace(FrameBaselineStart(Stim):FrameBaselineEnd(Stim)));
        x(Stim) = mean(FrameBaselineStart(Stim):FrameBaselineEnd(Stim));
    end

    %% find minima and interpolate F0 trace
    [~,~,~, minidx] = extrema(F0);
    minidx = unique([1 minidx length(F0)]); %just take the first and last F0 value to prevent stupid interpolation;

    mini = F0(sort(minidx))';

    F0_mini_interp  = interp1(x(sort(minidx)), smooth(mini,10,'rloess'), 1:recLength ,'linear', 'extrap');

    %% Using interpolated F0 to calculate dff
    DF = NaN(1,length(inputTrace));
    for stim = 1:length(FrameBaselineStart)
        F0New = nanmean(F0_mini_interp(FrameBaselineStart(Stim):FrameBaselineEnd(Stim)));
        DF(FrameBaselineStart(stim):FrameStimEnd(stim)) = inputTrace(FrameBaselineStart(stim):FrameStimEnd(stim))-F0New;
        DFF(FrameBaselineStart(stim):FrameStimEnd(stim)) = DF(FrameBaselineStart(stim):FrameStimEnd(stim))./F0New*100;

        %2nd zeroing
        %DFF(FrameBaselineStart(stim):FrameStimEnd) = DFF(FrameBaselineStart(stim):FrameStimEnd) - nanmean(DFF(FrameBaselineStart(stim):FrameBaselineEnd));
    end
    data.roi(RoiNR).dff = DFF;
    
    %% check trace
    %figure(2352)
    %plot(inputTrace);hold all
    %plot(F0_mini_interp); 
    %plot(DFF)
end