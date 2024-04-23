function metadata = getMetadataAndTimes(analysisParams, info)

%This reads in the twophoton times and the stimulus information from the
%saved text files and adapts them for further use
%
% Input:
% - analysisParams:uses the analysis parameter settings such as prestim,
% stimduration and poststim period
% - info: contains all info about the experiment, such as animal name and
% folder directions
%
% Output:
% - metadata: structure containing stimulus information as well as
% information about when frames were acquired and when stims were on
% - Stiminfo.mat: file containing the metadata

%% Step 1: load frametimes from twophotontimes.txt -> twophotontimes
twoPhotonFile = 'twophotontimes.txt';
tpFullFile = fullfile(info.Sp2Dir, twoPhotonFile);
infotpFullFile = dir(tpFullFile);
%Failsave: if the file is nearly empty, read frametrigger instad
if infotpFullFile.bytes < 10
    twoPhotonFile = 'frametrigger.txt';
    tpFullFile = fullfile(path, twoPhotonFile);
end
twophotontimes = load(tpFullFile);

%% Step 2: load triggers and stimIDs from stimontimes.txt -> stimOn, stimID
S = load([info.Sp2Dir 'stimontimes.txt']);
stimOn = S(2:2:length(S));
stimID = S(1:2:length(S)-1);
if stimID(1)==0 % Check to see if the first StimID is 0.  if it is, then delete it 
    stimOn(1) = [];% (initialization error with serial port in psychopy)
    stimID(1) = [];
end
if sum(stimID==0)>1 %if you make a mistake and 0 is a stim code
    stimID = stimID+1;
end
uniqStims = unique(stimID);
disp(['Animal: ' char(info.animal) ', experiment: ' char(info.name)]);
disp(['Loaded ', num2str(length(uniqStims)), ' unique stimCodes'])
disp(['Loaded ', num2str(length(stimOn)), ' stim on times'])

%% Step 3: Save everything in the metadata structure
metadata = [];
metadata(1).twophotontimes = twophotontimes;
metadata(1).copyStimID = stimID;
metadata(1).copyStimOn = stimOn;
metadata(1).uniqStims = uniqStims;
metadata(1).scanPeriod = mean(diff(metadata(1).twophotontimes(1:10)));
metadata(1).rate = 1/metadata(1).scanPeriod;

metadata(1).prestimPeriod = ceil(analysisParams.prestimPeriod ./metadata(1).scanPeriod);
metadata(1).stimDur = ceil(analysisParams.stimDur ./metadata(1).scanPeriod);
metadata(1).postPeriod = ceil(analysisParams.postPeriod ./metadata(1).scanPeriod);

%combine ontimes and twophotontimes to find frame times
metadata(1).ntrials = floor(length(stimOn)/length(uniqStims)); %automatically removes unfinished trials
metadata(1).numStims = metadata(1).ntrials*length(uniqStims);
stimOn2pFrame = zeros(1,metadata(1).numStims);

%go through all stims
for ii = 1:metadata(1).numStims
    id1 = stimOn(ii)<twophotontimes;
    id2 = stimOn(ii)>twophotontimes;
    ind = id1.*id2;
    ind = find(ind==1);
    if ~isempty(ind)
        stimOn2pFrame(ii) = ind;
    else
        stimOn2pFrame(ii) = find(diff(id1)==1);
    end
end

%downsample 4x and get df/f cycles
metadata(1).stimOn2pFrame = floor(stimOn2pFrame./4);
metadata(1).prestimPeriod2 = round(metadata(1).prestimPeriod/4);
metadata(1).stimDur2 = round(metadata(1).stimDur/4);
metadata(1).postPeriod2 = round(metadata(1).postPeriod/4);

%% Step 5: Save it in a file
save([info.saveDir 'StimInfos.mat'], 'metadata', '-mat') 