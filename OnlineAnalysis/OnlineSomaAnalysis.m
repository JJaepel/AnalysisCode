% set up directory and files
animalnum = 'F2726';
file2p  =8;
fileSpk =8;
date = '2022-09-02'; 
useLessTrials = 0;    

prestimPeriod = 0;    
stimDur = 2; %stimdur
postPeriod = 0;  

% Spk2dir = ['C:\Spike2Data\',date2,'\'];


addpath(['Z:\Juliane\Data\2P_Data\',animalnum,'_',date,'\'])
a = xlsread(['Z:\Juliane\Data\2P_Data\',animalnum,'_',date,'\t',sprintf('%.5d',file2p),'_IntegrationRois_00001.csv']);
T = a(:,3:size(a,2));
timeinfo =  a(:,1:2);
%% get Spk2 data

% cd([Spk2dir,'t',sprintf('%.5d',fileSpk)])
cd(['Z:\Juliane\Data\Spike2Data\',animalnum,'_',date,'\t',sprintf('%.5d',fileSpk)])
spk2twophotontimes = load('frametrigger.txt');
% twophotontimes = load('frametrigger.txt');
% twophotontimes = timeinfo(2,1)+spk2twophotontimes;%(1);
twophotontimes = spk2twophotontimes;
S = load('stimontimes.txt');
stimOn = S(2:2:length(S));
stimID = S(1:2:length(S)-1);

% Check to see if the first StimID is 0.  if it is, then delete it
% (initialization error with serial port in psychopy)
if stimID(1)==0
    stimOn(1) = [];
    stimID(1) = [];
end

if sum(stimID==0)>1 %if you make a mistake and 0 is a stim code
    stimID = stimID+1;
end


% f1 = fopen('stimtimes.txt', 'r');
% ST = fscanf(f1, '%f');
uniqStims = unique(stimID);
disp(['Loaded ', num2str(length(uniqStims)), ' unique stimCodes.'])
disp(['Loaded ', num2str(length(stimOn)), ' stim on times'])
preVisStim = find( twophotontimes < stimOn(1));
% twophotontimes = twophotontimes(1:4:end); %comment this out for single
% plane imaging

% global ce
ce = [];
ce(1).twophotontimes = twophotontimes;
ce(1).copyStimID = stimID;
ce(1).copyStimOn = stimOn;
ce(1).uniqStims = uniqStims;
%% get excel data

% ind = find(isnan(T(1,:)),1); % first empyt column
% ind = 1;
% numROIs = ind-1;
numROIs = size(T,2);
for i = 1:numROIs
    raw = T(:,i);
    ce(i).raw = raw;
    
    if i==1
        % take the avg of first 10 frames to give actual img rate
        ce(1).scanPeriod = nanmean(diff(ce(1).twophotontimes(1:10)));
        ce(1).rate = 1/ce(1).scanPeriod;
    end
    
end   
disp('loaded excel integration data')    
%% align stimuli times and extract response

%ce(1).scanPeriod = nanmean(diff(ce(1).twophotontimes(1:10)));
scanPeriod = ce(1).scanPeriod;
if scanPeriod>0.037
    scanPeriod = 0.0341;
end
% resPeriod = resPeriod/2; for 60Hz imaging???

prestimPeriod2 = ceil(prestimPeriod./scanPeriod);
stimDur2 = ceil(stimDur./scanPeriod);
postPeriod2 = ceil(postPeriod./scanPeriod);

stimID = ce(1).copyStimID;
stimOn = ce(1).copyStimOn;
uniqStims = ce(1).uniqStims;
twophotontimes = ce(1).twophotontimes;

if sum(uniqStims==0)==1
    uniqStims = uniqStims(2:end);
end

ntrials = floor(length(stimOn)/length(uniqStims));
ntrials = ntrials-useLessTrials;

numStims = ntrials*length(uniqStims);
% numStims = numStims-2;
stimOn2pFrame = zeros(1,numStims);

%convert stim times into frame times
for ii = 2:numStims
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
    
xt = timeinfo(:,2);
ft = timeinfo(:,1);
xq = 1:timeinfo(end,2);
fq = timeinfo(1,1):mean(diff(timeinfo(1:100,1))):timeinfo(end,1);

stimOn2pFrame = floor(stimOn2pFrame./4);
stimDur2 = round(stimDur2/4);
postPeriod2 = round(postPeriod2/4);

disp 'Cutting up data based off stimIDs...'
for cc = 1:length(ce)
    
    ce(cc).cyc = zeros(length(uniqStims),ntrials,stimDur2+postPeriod2+prestimPeriod2);
    ce(cc).stimOn2pFrame = stimOn2pFrame;
    
        
    raw = ce(cc).raw;
    
    raw2 = interp1(xt,raw,xq);
    
    dff = filterBaseline_dFcomp(resample(raw2,1,4));
    ce(cc).dff = dff;
    
    trialList = zeros(1,length(uniqStims));
    for ii = 2:numStims-1
        stimTime2 = stimOn2pFrame(ii)+1-prestimPeriod2:stimOn2pFrame(ii)+stimDur2+postPeriod2;
        ind = find(uniqStims==stimID(ii));
        trialList(ind) = trialList(ind)+1;
        if stimTime2(1)>length(dff)
            f = zeros(length(stimTime2),1);
        else
            f = dff(stimTime2);
        end
        ce(cc).cyc(ind,trialList(ind),:) = f;% - median(f(1:2));
    end
    fprintf(num2str(cc))
end    
    
%%
for k = 1:numROIs
    figure;
    plotcycRes2(ce(k).cyc,1:size(ce(k).cyc,1))
    supertitle(['soma ',num2str(k)])    
end
    
    
    
    