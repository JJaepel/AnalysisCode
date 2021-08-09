function [metadata, data] = getF0FromBaseline(metadata, data, fieldSource)
%field source is usually 'rawF'


%% get all baseline data together

%find start and end frames of all baseline times
TimeBaselineStart = metadata.StimParams.StimOnTimes(2,:)-1-metadata.StimParams.isi;
TimeBaselineEnd = metadata.StimParams.StimOnTimes(2,:)-1;
TimeStimEnd = metadata.StimParams.StimOnTimes(2,:)+metadata.StimParams.stimDuration;

for base = 1:length(TimeBaselineEnd)
    FrameBaselineStart(i) = find(metadata.TwoPhoton > TimeBaselineStart(i),1);
    FrameBaselineEnd(i) = find(metadata.TwoPhoton > TimeBaselineEnd(i),1);
    FrameStimEnd(i) = find(metadata.TwoPhoton > TimeStimEnd(i),1);
end

inputTrace = data.roi(i).(fieldSource)';

for i = 1:length(TimeBaselineEnd)
    F0(i) = nanmean(inputTrace(FrameBaselineStart(i):FrameBaselineEnd(i));
    x(i) = mean(FrameBaselineStart(i):FrameBaselineEnd(i));
end

%% find minima and interpolate F0 trace
[~,~,~, minidx] = extrema(F0);
minidx = unique([1 minidx length(F0)]); %just take the first and last F0 value to prevent stupid interpolation;

mini = F0(sort(minidx))';

F0_mini_interp  = interp1(x(sort(minidx)), smooth(mini,10,'rloess'), 1:reclength ,'linear', 'extrap');
F0_interp = 1;

%% reextract F0 values from interpolated F0 trace
DF = NaN(length(inputTrace));
for stim = 1:length(F0)
    DF = inputTrace(FrameBaselineStart(stim):FrameStimEnd)-F0;
end
