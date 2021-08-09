function Patches_Jad(analysisParams)


close all

%% 0.) define folders and structures
if analysisParams.server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

base2pDirectory= [TwoPhontondir analysisParams.animal];
tifDirectory = [base2pDirectory filesep analysisParams.name];
Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep];
analysisParams.baseDirectory = base2pDirectory;


metadata = struct;
metadata.ROI = struct;
analysis = struct;

%% 1). load Data and metadata
disp('Loading data')
data = LoadRoisS2p(analysisParams);
number_of_ROIs = length(data.roi);

twoP_times = importdata ([Sp2dDirectory 'frametrigger.txt']);    % choose one and go to stim directory
twoP_times_new = twoP_times (1:end-2);
figure; plot(twoP_times,'r*');

stim_on_times_and_condition_nb = importdata ([Sp2dDirectory 'stimontimes.txt']); % go to corresponding directory

stim_condition = stim_on_times_and_condition_nb (1, 3:2:end)';  % change "3" to "1" for Kuo's setup
stim_on_times = stim_on_times_and_condition_nb (1, 4:2:end)';   % change "4" to "2" for Kuo's setup

%% 2.) Convert stim times into frame times
uniqStims = unique(stim_condition);
ntrials = floor(length(stim_on_times)/length(uniqStims));
numStims = ntrials*length(uniqStims);
stimOn2pFrame = zeros(1,numStims);
for ii = 1:numStims
    id1 = find(stim_on_times(ii)<twoP_times,1,'first');
    id2 = find(stim_on_times(ii)>twoP_times,1,'last');
    ind = id2;
    if ~isempty(ind)
        stimOn2pFrame(ii) = ind;
    elseif mean(diff(id1)~=0)
        stimOn2pFrame(ii) = find(diff(id1)==1);
    else
        stimOn2pFrame(ii) = 1;
    end
end

%% 3.) Cutting up data based off stim_conditions

prestimPeriod = 1;   % take last 1 second of 4 second ISI
stimDur = 2; %2 second stim
postPeriod = 1;
useLessTrials = 0; 


scanPeriod = mean(diff(twoP_times(1:10)));
rate = 1/scanPeriod;
prestimPeriod2 = ceil(prestimPeriod./scanPeriod);
stimDur2 = ceil(stimDur./scanPeriod);
postPeriod2 = ceil(postPeriod./scanPeriod);

ce = struct;
disp 'Cutting up data based off stimIDs...'

for cc = 1:number_of_ROIs
    
    ce(cc).cyc = zeros(length(uniqStims),ntrials,(stimDur2+postPeriod2));
    ce(cc).cyc2 = zeros(length(uniqStims),ntrials,(stimDur2+postPeriod2));
    ce(cc).stimOn2pFrame = stimOn2pFrame;
    ce(cc).raw_signal = data.roi(cc).rawF;
    
    % raw_signal = ce(cc).raw_signal;
    % dff = filterBaseline_dFcomp(raw);
    dff = ce(cc).raw_signal;
    ce(cc).dff = dff;
    
    trialList = zeros(1,length(uniqStims));
    dFtrace = [];
    for ii = 1:numStims
       
        prestimTime = stimOn2pFrame(ii)-prestimPeriod2:stimOn2pFrame(ii);
        stimTime = stimOn2pFrame(ii)+1:stimOn2pFrame(ii)+stimDur2+postPeriod2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%try manual dFF
        %check to see if there is prestimTime
     %   prestimTime = prestimTime(find(prestimTime>0));
%         if max(stimTime)>length(ce(cc).raw)
%             indt = find(stimTime>length(ce(cc).raw),1,'first')
%             stimTime = stimTime(1:indt-1);
%         end
        
        Ftrace = ce(cc).raw_signal(stimTime);
        Fo = nanmean(ce(cc).raw_signal(prestimTime));
        
        Fo_trace = ce(cc).raw_signal(prestimTime);    %% added by Jad on 8/6/2018
              
        
        if Fo==0 
            dFtrace = nan(size(Ftrace));
            dFo_trace = nan(size(Fo_trace));      %% added by Jad on 8/6/2018
        else
            dFtrace = (Ftrace - Fo)./Fo;
            dFo_trace = (Fo_trace - Fo)./Fo;      %% added by Jad on 8/6/2018
        end
        
        if length(uniqStims)==1
            ind=1;
        else
            ind = stim_condition(ii);
        end
        trialList(ind) = trialList(ind)+1;
        ce(cc).cyc(ind,trialList(ind),1:length(stimTime)) = dFtrace';
        
        ce(cc).Fotrace(ind,trialList(ind),1:length(prestimTime)) = dFo_trace';         %% added by Jad on 8/6/2018
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%try the slow filter dFF
        f = dff(stimTime);
        ce(cc).cyc2(ind,trialList(ind),1:length(stimTime)) = f - median(f(1:2));
    end
end

%% 4.) Average dF/F for each ROI for each trial for each condition

for cc = 1:number_of_ROIs
    ce(cc).trial_response = nanmean (ce(cc).cyc, 3);
    
    ce(cc).trial_pre_stim_response = nanmean (ce(cc).Fotrace, 3);
    ce(cc).trial_pre_stim_response_std = std (ce(cc).Fotrace, [], 3, 'omitnan');        
end


%% 5.) Trial average dF/F for each ROI for each condition

for cc = 1:number_of_ROIs
   
    ce(cc).condition_response = nanmean (ce(cc).trial_response, 2);
    ce(cc).condition_pre_stim_response = nanmean (ce(cc).trial_pre_stim_response, 2);
    
    ce(cc).condition_response_std = nanstd (ce(cc).trial_response')' ;
    ce(cc).condition_pre_stim_response_std = nanstd (ce(cc).trial_pre_stim_response')' ;
    
    ce(cc).condition_response_sem = (nanstd (ce(cc).trial_response') ./ sqrt (ntrials))' ;
    ce(cc).condition_pre_stim_response_sem = (nanstd (ce(cc).trial_pre_stim_response') ./ sqrt (ntrials))' ;
    
end

%% 6.) Response significance tests %%

% Response significance (2 SD from baseline) for each ROI for each condition (NOT for each trial)
for cc = 1:number_of_ROIs
    
    ce(cc).condition_baseline_plus_2SD = ce(cc).condition_pre_stim_response + 2* (ce(cc).condition_pre_stim_response_std);
    
    for i= 1: length (uniqStims (:,1))
     if ce(cc).condition_response (i,1) > ce(cc).condition_baseline_plus_2SD (i,1)
        ce(cc).condition_response_significance (i,1) = 1;
    else
        ce(cc).condition_response_significance (i,1) = 0;
     end
    end
end


% Response significance (2 SD from baseline) for each ROI  (ROI_responsiveness)
for cc = 1:number_of_ROIs
    ROI_responsiveness (cc,1) = any (ce(cc).condition_response_significance);
end


% Trial response significance (2 SD from baseline) for each ROI for each trial    
for cc = 1:number_of_ROIs
    
    ce(cc).trial_baseline_plus_2SD = ce(cc).trial_pre_stim_response + 2* (ce(cc).trial_pre_stim_response_std);
    
    for i= 1: length (ce(cc).trial_response (:,1))
        for ii= 1: length (ce(cc).trial_response (1,:))
     if ce(cc).trial_response (i,ii) > ce(cc).trial_baseline_plus_2SD (i,ii)
        ce(cc).trial_response_significance (i,ii) = 1;
    else
        ce(cc).trial_response_significance (i,ii) = 0;
     end
        end
    end
end
    

% paired t-test for response significance (comparing for each cell for each condition the trial responses pre and post stim onset)
for i = 1: number_of_ROIs
    for ii = 1: length (uniqStims)
        h = ttest2(ce(i).trial_pre_stim_response(ii,:), ce(i).trial_response(ii,:));
        paired_t_test_significance (i, ii) = h;
    end
end



%% 7.) Plot raw signal for 5x5 RF %%


for i = 54  % set ROI number
    figure; 
    for ii = 1: floor(length (uniqStims)/2)
        for iii = 1: ntrials 
            
            subplot (8,7,ii);
            
            temp_trace (1,:) = (ce(i).cyc(ii,iii,:));
            concatenated_temp_traces (iii,:) = temp_trace;
            plot (temp_trace, 'color', [0.5 0.5 0.5]);
            hold on;
            
            temp_F0_trace (1,:) = (ce(i).Fotrace(ii,iii,:));
            concatenated_temp_F0_traces (iii,:) = temp_F0_trace;
            plot ([-32:1:-1], temp_F0_trace, 'color', [0.5 0.5 0.5]);
            hold on;
            
            axis ('square');
            xlim ([-23 38]);
%           ylim ([-0.2 1]);
%           plot ([0 0],[-10 10], '--r');
%           plot ([30 30],[-10 10], '--r');
            
%             if ii == 26
%                 title ('Blank');
%             end
            
            if paired_t_test_significance (i,ii) == 1
                set (gca, 'color', [0.2, 0.75, 0.95]);
            end
        
        end
        
        Average_temp_trace = mean (concatenated_temp_traces,1);
        plot (Average_temp_trace, '-k', 'LineWidth', 2);           
            
        Average_temp_F0_trace = mean (concatenated_temp_F0_traces,1);
        plot ([-32:1:-1], Average_temp_F0_trace, '-k', 'LineWidth', 2);
        
    end
end

%% 8.) Calculate and display individual, or averaged RFs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:number_of_ROIs
    Resp_no_blank (i, 1:98) = ce(i).condition_response (1:end-1);
end

Resp_matrix = reshape(Resp_no_blank', [14,7,number_of_ROIs]);      
temp_Resp_matrix = rot90 (Resp_matrix, 3);
temp2_Resp_matrix = fliplr (temp_Resp_matrix);
Resp_matrix = temp2_Resp_matrix;

for ii = 1:number_of_ROIs   % finding peak response
    [Pref_location_resp(ii) Pref_location_index(ii)] = max (Resp_no_blank (ii,:));
end

% to look at single cell RF examples
figure;
for i =30:55     % number_of_ROIs
subplot (5,5,i);
imagesc (Resp_matrix(:,:,i));
colormap (gray);
axis ('square');
end

% to calculate and plot the average RF of all cells
Average_Resp_matrix = mean (Resp_matrix, 3);
imagesc (Average_Resp_matrix);
axis ('square');


% normalize reponse to peak
for i = 1:number_of_ROIs  
    Resp_no_blank_normalized (i, 1:98) = ce(i).condition_response (1:end-1) / ce(i).condition_response (Pref_location_index (1,i),1);
end

Resp_matrix_normalized = reshape(Resp_no_blank_normalized', [14,7,number_of_ROIs]);  
temp_Resp_matrix_normalized = rot90 (Resp_matrix_normalized, 3);
temp2_Resp_matrix_normalized = fliplr (temp_Resp_matrix_normalized);
Resp_matrix_normalized = temp2_Resp_matrix_normalized;

% to look at single cell RF examples normalized to each cells peak response
figure;
for i =30:55     % number_of_ROIs
subplot (5,5,i);
imagesc (Resp_matrix_normalized(:,:,i));
colormap (gray);
axis ('square');
end

% to calculate and plot the average normalized RF of all cells
Average_Resp_matrix_normalized = mean (Resp_matrix_normalized, 3);
imagesc (Average_Resp_matrix_normalized);
axis ('square');
