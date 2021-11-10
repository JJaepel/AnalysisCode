function metadata = createStimCodesEpi(metadata, analysisParams, stackSize)

numberOfConditions = metadata.StimParams.numberOfStims;
stimStartIndex = zeros(numberOfConditions,1,'double'); 
stimStopIndex  = zeros(numberOfConditions,1,'double');
metadata.StimParams.directions = linspace(0,360,metadata.StimParams.uniqStims);
metadata.StimParams.directions = metadata.StimParams.directions(1:end-1);

metadata.StimParams.stimDuration = 5;
metadata.StimParams.isi = 5;

for i=1:numberOfConditions
    stimStartIndex(i) = find(metadata.Imaging.time>=metadata.StimParams.StimOnTimes(2,i),1,'first');
    stimStopIndex(i) = find(metadata.Imaging.time>=metadata.StimParams.StimOnTimes(2,i)+metadata.StimParams.stimDuration,1,'first');
end
metadata.StimParams.stimStartIndex = stimStartIndex;
metadata.StimParams.stimStopIndex = stimStopIndex;

if analysisParams.downsample > 1
    metadata.Imaging.time = metadata.Imaging.time(1:analysisParams.downsample:stackSize);
    metadata.Imaging.rate = metadata.Imaging.rate/analysisParams.downsample;
    try
        metadata.StimParams.stimStartIndex = floor(metadata.StimParams.stimStartIndex/analysisParams.downsample);
        metadata.StimParams.stimStopIndex = floor(metadata.StimParams.stimStopIndex/analysisParams.downsample);
    catch
    end
end