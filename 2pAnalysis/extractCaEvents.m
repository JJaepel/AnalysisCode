function ce = extractCaEvents(ce, info, verbose)

% Returns significant Ca events and the respective amplitude based on 
% ce.cyc(Res) size knows which stim conditions were used. events
% are only extracted for stimuli, not blank. so events will have size
% (cond x trials)
%
% Inputs:
% - ce: structure containing all information about the ROIs, each row for
% an individual ROI
% - info: all info about the specific experiment such as the save dir
% - verbose: do you want to plot?
%
% Outputs:
% - ce: updated with cycRes, containing the spine trace after dendritic
% substraction, & corr, containing the residual correlation with dendrite
% - ROIEvents.mat: .mat file with new ce structure

%% Step 0: Set parameters

if nargin<3
    verbose = 0;
end

SDCUTTOFF = 2;

for cc = 1:length(ce)
    if ce(cc).spine
        cyc = ce(cc).cycRes;
    else
        cyc = ce(cc).cyc;
    end
    
    %% Step 1: Calculate spine F noise before substraction for later use
    Spdff = filterBaseline_dFcomp(resample(ce(cc).raw,1,4));
    Spdff = Spdff(1:round(length(Spdff)*.9));
    Spdff(isinf(Spdff)) = 0;
    Spdff_sub = Spdff(Spdff < nanmedian(Spdff)+abs(min(Spdff)));
    noiseM = nanmedian(Spdff_sub); % should be near 0
    noiseSD = nanstd(Spdff_sub);
    
    %% Step 2:Identify spine 'events': decay is 400 ms time constant
    % (roughly following Konnerth)
    alpha = 0.4;
    ind = 1:size(cyc,1)-1; % dont go through blank
    
    events = zeros(length(ind),size(cyc,2));
    eventAmps = events;
    for jj = 1:size(cyc,2)
        for i = 1:length(ind)
            ii = ind(i);
            r = squeeze(cyc(ii,jj,:));
            %%%%%%%%%%%%%%%%%%%%%%%
            %exponential filter (r: dff trace)
            rf = filter(alpha, [1 alpha-1],r);

            [~,locs]= findpeaks(rf);
            %ignore trial beginning/end
            locs = locs(locs>1 & locs<length(rf)-3);
            amps = zeros(length(locs),1);
            %avg 2 pts (choose best from forward or backward)
            for z = 1:length(locs)
                amps(z) = max([mean(r(locs(z)-1:locs(z))) mean(r(locs(z):locs(z)+1))]);
            end

            locs = locs(amps>(noiseSD*SDCUTTOFF+noiseM));

            events(i,jj) = sum(locs)>0;

            if isempty(locs)
                eventAmps(i,jj) = 0;
            else
                eventAmps(i,jj) = max(amps);
            end

            if verbose

                figure(jj)
                subplot(1,size(cyc,1),ii)
                hold off
                plot(1:length(r),r,'b',1:length(rf),rf,'r')
                hold on
                plot(1:length(r),ones(1,length(r)).*(noiseSD*SDCUTTOFF +noiseM),'--k')
                ylim([-noiseSD 1])


            end
        end
    end
    ce(cc).events = events;
    ce(cc).eventAmps = eventAmps;
    
end

%% Step 3: Save in mat file
save([info.saveDir 'ROIsEvents.mat'], 'ce', '-mat') 