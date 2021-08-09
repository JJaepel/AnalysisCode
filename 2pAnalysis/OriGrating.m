function OriGrating(analysisParams)

close all
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
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
if ~exist(ROIsaveDirectory, 'dir')
    mkdir(ROIsaveDirectory);  
end

disp('Loading data')
field = analysisParams.field;

%% load Data and metadata
if analysisParams.reloadData
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    data = LoadRoisS2p(analysisParams);
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');

    save(fullfile(saveDirectory, 'oriGrating.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
else
    load(fullfile(saveDirectory, 'oriGrating.mat'), 'data', 'metadata', 'analysis');
end

%% chop traces
disp('Chopping Traces')
[analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,analysisParams.level, field, 'pre', analysisParams.pre, 'post',metadata.StimParams.isi,'windowStart',analysisParams.windowStart, 'windowStop',analysisParams.windowStop);

%% find maxResponses and significant responses
disp('Calculating significant responses')
for i = 1:length(analysis.(field).roi)
    pretrialTime= analysis.(field).preTrialTime;
    preTrialIndex= (1:floor(pretrialTime * metadata.TwoPhoton.rate));
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    %collect our pretrial interval
    analysis.(field).roi(i).isRespSignificant = false;
    analysis.(field).roi(i).respThreshold = [];
    baselines=analysis.(field).roi(i).stimResponseTrace(:,:,preTrialIndex);
    analysis.(field).roi(i).baselineSD = std(baselines,[],3);
    analysis.(field).roi(i).baselineMean = mean(baselines,3);
    analysis.(field).roi(i).baselineMean(analysis.(field).roi(i).baselineMean < 0) = mean(mean(analysis.(field).roi(i).baselineMean,2),1);
    analysisPeriod=(analysis.(field).windowStart:analysis.(field).windowStop);
    stimResp = analysis.(field).roi(i).stimResponseTrace(:,:,analysisPeriod);
    analysis.(field).roi(i).peaks = max(stimResp,[],3);
    analysis.(field).roi(i).zscore = ([analysis.(field).roi(i).peaks]-[analysis.(field).roi(i).baselineMean])./[analysis.(field).roi(i).baselineSD];
    analysis.(field).roi(i).crosser = sum(analysis.(field).roi(i).zscore > analysisParams.zThresh,2);
    analysis.(field).roi(i).respStim = analysis.(field).roi(i).crosser >= ((metadata.StimParams.numTrials)*analysisParams.fraction);
    if sum(analysis.(field).roi(i).respStim) > 0
        analysis.(field).roi(i).isResponseSignificant = 1;
    else 
        analysis.(field).roi(i).isResponseSignificant = 0;
    end
end

disp('Calculating ROI properties')
numberOfStims = metadata.StimParams.uniqStims;
numStims=numberOfStims-1;
metadata.StimParams.numSf = 1;
metadata.StimParams.numTf = 1;
metadata.StimParams.theta= [0:2*pi/(numStims) : 2*pi-2*pi/(numStims)];
theta = metadata.StimParams.theta;
theta=mod(theta, pi);
orientations = linspace(0,180,numStims/2+1);
directions = linspace(0,360,numStims+1);
metadata.StimParams.directions = directions(1:end-1);
orientations = orientations(1:end-1);
metadata.StimParams.numOrientations = length(orientations);
lastori = length(orientations);
for i = 1:length(data.roi)
    %%find preferred orientation for all cells
    medResponse = analysis.(field).roi(i).stimResponseTrace(1:end-1, :, :);
    medResponse = mean(medResponse(:,:,stimWindow),3);
    Respones = medResponse;
    medResponseA(1:size(medResponse,1)/2, 1:size(medResponse,2)) = medResponse(1:size(medResponse,1)/2,:);
    medResponseA(1:size(medResponse,1)/2, 1+size(medResponse,2):2*size(medResponse,2))= medResponse(1+size(medResponse,1)/2:size(medResponse,1),:);
    medResponseA= median(medResponseA,2);
    medResponse = medResponseA(:)';
    theta=theta(1:lastori);
    theta= repmat(theta, 1, length(medResponse)/length(theta));
    
    [analysis.(field).roi(i).preferredOrientation,analysis.(field).roi(i).coeffOr,analysis.(field).roi(i).rsqOr] = ComputePreferredOrientations(medResponse, theta);
    analysis.(field).roi(i).preferredOrientation =mod(analysis.(field).roi(i).preferredOrientation, 180);
    analysis.(field).roi(i).rsqOr= Rsquared(medResponse, vonMisesFit(analysis.(field).roi(i).coeffOr, theta*2), true);
    analysis.(field).roi(i).OSIFit = computeOSIFit(analysis.(field).roi(i).coeffOr);
    [analysis.(field).roi(i).OSI, ~] = computeOSI(lastori,medResponse);
    [analysis.(field).roi(i).cohensD] = computeCohensD(lastori,Respones);
    [~, ind] = min(abs(orientations-analysis.(field).roi(i).preferredOrientation));
    analysis.(field).roi(i).prefOriStim = orientations(ind);
    analysis.(field).roi(i).prefOriStimInd = ind;
    medResponse = mean(analysis.(field).roi(i).avgStimResponse, 2)';
    analysis.(field).vsOrientation(i)= wrapTo2Pi(angle(vectorSum(abs(medResponse(1:lastori)),2)))/2 * 180/pi;
    analysis.(field).orientationMag(i) = abs(vectorSum(abs(medResponse(1:lastori))/max(medResponse(1:lastori)),2));
    %find preferred direction for all cells
    dirResponse = analysis.(field).roi(i).stimResponseTrace(1:end-1, :, :);
    dirResponse = mean(dirResponse(:,:,stimWindow),3);
    dirResponse = median(dirResponse,2)';
    bestCoeff=[NaN, NaN, NaN, NaN];
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    bestCoeff=fitVonMisesLinSumFunction(dirResponse, metadata.StimParams.theta, bestCoeff, options);
    xo=[0:pi/40:2*pi];
    fit=vonMisesLinSum(bestCoeff, xo);
    if bestCoeff(1) > bestCoeff(2)
        analysis.(field).roi(i).preferredDirection = rad2deg(bestCoeff(4));
    else
        analysis.(field).roi(i).preferredDirection = mod(rad2deg(bestCoeff(4))+180,360);
    end
    [analysis.(field).roi(i).DSI, ~] = computeDSI(lastori*2, dirResponse, 'preferred');
    [analysis.(field).roi(i).CircVar] = TT_CircularVariance(dirResponse);
end

[zoom, setup] = getzoom(tifDirectory);
metadata.zoom = zoom;
if setup == 1
    fieldofview = 1000/zoom;
    umperpixel = fieldofview/512;
    disp('Setup Ben')
elseif setup == 2
    umperpixel = 2.73/zoom;
end

%% plot results
coc_prop = cbrewer('qual', 'Paired', 12);

figure
imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
axis image
hold on
for i = 1:length(data.roi)
    xpos= data.roi(i).xPos;
    ypos= data.roi(i).yPos;
    plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor', 'blue');
    text(xpos, ypos, num2str(i),'HorizontalAlignment','center', 'VerticalAlignment','middle','color', 'white') 
end
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'ROI_positions.png'))

%% plot prefOri on top of template
for types = 1:3
    if types == 1
        rois_ori = linspace(1,length(analysis.dff.roi),length(analysis.dff.roi));
        alloriprefs = [analysis.dff.roi.preferredOrientation];
        alldirprefs = [analysis.dff.roi.preferredDirection];
        rois_dir = rois_ori;
    elseif types == 2
        rois_ori = find([analysis.dff.roi.isResponseSignificant] == 1);
        alloriprefs = [analysis.dff.roi(rois_ori).preferredOrientation];
        alldirprefs = [analysis.dff.roi(rois_ori).preferredDirection];
        rois_dir = rois_ori;
    elseif types == 3
        rois_ori = find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alloriprefs = [analysis.dff.roi(rois_ori).preferredOrientation];
        rois_dir = find([analysis.dff.roi.DSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alldirprefs = [analysis.dff.roi(rois_dir).preferredDirection];
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    PlotPrefOnTemplateOri(analysis, data, metadata,1, field,data.template, rois_ori)

    subplot(2,2,2)
    histogram(alloriprefs,linspace(0,180,lastori+1), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel(sprintf('Orientation preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)]/2)
    axis square;
    set(gca,'Box','off');

    subplot(2,2,3)
    PlotPrefOnTemplateOri(analysis, data, metadata,2, field,data.template, rois_dir)

    subplot(2,2,4)
    histogram(alldirprefs,linspace(0,360,2*lastori+1), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel(sprintf('Direction preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)])
    axis square;
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    if types == 1
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_all_cells.png'))
    elseif types == 2
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_resp_cells.png'))
    elseif types == 3
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_selective_cells.png'))
    end
    %close all
end
%% 
figure
rois_ori = find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
colormap(hsv)
LUT = hsv(180);
hold on
for i = 1:length(rois_ori)
    l = rois_ori(i);
    xpos= data.roi(l).xPos;
    ypos= data.roi(l).yPos;
    try
        plot(xpos,ypos,'ok','MarkerSize', 20,'MarkerFaceColor',LUT(1+floor(analysis.dff.roi(l).preferredOrientation),:));
    catch
    end
end
saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_forAI.eps'))

figure
subplot(1,3,1)
all = length(analysis.dff.roi);
non_resp = length(find([analysis.dff.roi.isResponseSignificant] == 0)) ./all;
resp = length(find([analysis.dff.roi.isResponseSignificant] == 1)) ./all;
h = pie([non_resp resp]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', coc_prop(7,:));
try set(hp(2), 'FaceColor', coc_prop(8,:)); end
title('Responsive')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,2)
ori = length(find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1)) ./all;
non_ori = 1- ori - non_resp;
h = pie([non_resp non_ori ori]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', coc_prop(7,:));
try set(hp(3), 'FaceColor', coc_prop(2,:)); end
try set(hp(2), 'FaceColor', coc_prop(1,:)); end
title('Orientation-selective')
legend({'Non-resp', 'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,3)
dir = length(find([analysis.dff.roi.DSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1)) ./all;
non_dir = 1-dir-non_resp;
h = pie([non_resp non_dir dir]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', coc_prop(7,:));
try set(hp(3), 'FaceColor', coc_prop(4,:)); end
try set(hp(2), 'FaceColor', coc_prop(3,:)); end
title('Direction-selective')
set(gca, 'box', 'off')
legend({'Non-resp', 'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'resp_cells.png'))
%close gcf

figure 
allOSI = [analysis.dff.roi.OSIFit];
resp_cells = find([analysis.dff.roi.isResponseSignificant] == 1);
respOSI = [analysis.dff.roi(resp_cells).OSIFit];
allDSI = [analysis.dff.roi.DSI];
respDSI = [analysis.dff.roi(resp_cells).DSI];
allCohensD = [analysis.dff.roi.cohensD];
respCohensD = [analysis.dff.roi(resp_cells).cohensD];

subplot(2,3,1)
histogram(allOSI, 10, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlabel('OSI');
title('All cells')
set(gca,'Box','off');

subplot(2,3,2)
histogram(allDSI, 10, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlabel('DSI');
title('All cells')
set(gca,'Box','off');

subplot(2,3,3)
hdi(1) = cdfplot(allCohensD);  hold all;
set(hdi(1), 'Color', coc_prop(9,:), 'LineWidth', 2);
grid off;
ylabel('Cumulative probablity');
xlabel('CohensD');
xlim([0 5]);
title('All cells')
set(gca,'Box','off');

subplot(2,3,4)
histogram(respOSI,10,'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlabel('OSI');
title('All responsive cells')
set(gca,'Box','off');

subplot(2,3,5)
histogram(respDSI,10, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlabel('DSI');
title('All responsive cells')
set(gca,'Box','off');

subplot(2,3,6)
try
    hdi(1) = cdfplot(respCohensD);  hold all;
    set(hdi(1), 'Color', coc_prop(9,:), 'LineWidth', 2);
    grid off;
    ylabel('Cumulative probablity');
    xlabel('CohensD');
    xlim([0 5]);
    title('All responsive cells')
    set(gca,'Box','off');

    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, 'OSI_DSI_CohensD_distribtion.png'))
    %close gcf
catch
end

disp('Calculating map parameters')
%% plot delta Ori vs. distance of ROIs
ori_sel = find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
dist_ori = zeros(length(ori_sel), length(ori_sel));
delta_ori = zeros(length(ori_sel), length(ori_sel));
for A = 1:length(ori_sel)
    for B = 1:length(ori_sel)
        if analysisParams.level == 1
            if data.roi(A).plane == data.roi(B).plane
                dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
                dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
                dist_ori(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
                delta_ori(A,B) = abs(analysis.dff.roi(A).preferredOrientation - analysis.dff.roi(B).preferredOrientation);
                if delta_ori(A,B) > 90
                    delta_ori(A,B) = 180 - delta_ori(A,B);
                end
            else
                dist_ori(A,B) = NaN;
                delta_ori(A,B) = NaN; 
            end
        else
            dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
            dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
            dist_ori(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
            delta_ori(A,B) = abs(analysis.dff.roi(A).preferredOrientation - analysis.dff.roi(B).preferredOrientation);
            if delta_ori(A,B) > 90
                delta_ori(A,B) = 180 - delta_ori(A,B);
            end
        end
    end
end
L = triu(dist_ori);
L = reshape(L,1,length(ori_sel)^2);
L(L == 0) = NaN;
dist_ori_clean = L((~isnan(L)));
delta_ori = reshape(delta_ori,1,length(ori_sel)^2);
delta_ori = delta_ori((~isnan(L)));
%plot(dist_ori_clean,delta_ori, '*')
edges = linspace(0, 700, 15);
[n, edges, bin_dist] = histcounts(dist_ori_clean,edges);
edge = edges(1:end-1)+25;
mean_deltaOri = zeros(length(edge),1);
SEM_deltaOri = zeros(length(edge),1);
for bin = 1:length(edge)
    mean_deltaOri(bin) = nanmean(delta_ori(bin_dist == bin));
    SEM_deltaOri(bin) = nanstd(delta_ori(bin_dist == bin))/sqrt(n(bin));
end

bin_deltaOri_shuffle = zeros(analysisParams.shufflenum,bin);
sem_deltaOri_shuffle = zeros(analysisParams.shufflenum,bin);
for rep = 1:analysisParams.shufflenum
    delta_ori_shuffle = delta_ori(randperm(length(delta_ori)));
    for bin = 1:length(edge)
        bin_deltaOri_shuffle(rep,bin) = nanmean(delta_ori_shuffle(bin_dist == bin));
        sem_deltaOri_shuffle(rep,bin) = nanstd(delta_ori_shuffle(bin_dist == bin))/sqrt(n(bin));
    end
end
mean_deltaOri_shuffle = nanmean(bin_deltaOri_shuffle);
SEM_deltaOri_shuffle = nanmean(sem_deltaOri_shuffle);

for i = 1:length(delta_ori)
    analysis.dff.ori_pair(i).distance = dist_ori_clean(i);
    analysis.dff.ori_pair(i).deltaOri = delta_ori(i);
end
    
figure
errorbar(edge,mean_deltaOri,SEM_deltaOri, 'o-', 'Color', coc_prop(2,:), 'MarkerFaceColor', coc_prop(2,:))
hold all
errorbar(edge,mean_deltaOri_shuffle,SEM_deltaOri_shuffle, 'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel('Distance in \mum')
ylabel('\DeltaOrientation preference (\circ)')
ylim([0 90])
xlim([0 700])
legend('Data', 'Shuffled')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'OSI_distance.png'))

% calculate HetereogeneityIndex
radius = linspace(0,250,6);
if analysisParams.level == 1 
    plane_ori = [data.roi(ori_sel).plane];
end
HI = zeros(length(ori_sel),length(radius)-1);
preferences = [analysis.dff.roi(ori_sel).preferredOrientation];
for dis = 1:length(radius)-1
    for cellID = 1:length(ori_sel)
        distances = dist_ori(cellID,:);
        withinRadius1 = distances >= radius(dis);
        withinRadius2 = distances < radius (dis+1);
        withinRadius = withinRadius1&withinRadius2;
        if analysisParams.level == 1
            sameplane = [plane_ori == data.roi(cellID).plane];
            withinRadius = withinRadius&sameplane;
        end
        withinRadius(cellID)=0;
        thresholdOrientations=exp(2*1i*preferences(withinRadius));
        HI(cellID, dis) = 1-(abs(1/sum(withinRadius) * sum(thresholdOrientations)));
        if HI(cellID,dis) == 0
            HI(cellID,dis) = NaN;
        end
    end
end
analysis.dff.ori_cells.HomeogeneityIndex = HI;

figure
subplot(1,3,1)
plot(nanmedian(HI(:,1)), max(histcounts(HI(:,1))),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(5,:),'MarkerFaceColor',coc_prop(5,:)); hold on
text(nanmedian(HI(:,1)), 1.05 * max(histcounts(HI(:,1))),num2str(round(100*nanmedian(HI(:,1)))/100),'HorizontalAlignment','center', 'Color', coc_prop(6,:), 'FontSize', 12')
histogram(HI(:,1),10,'FaceColor', coc_prop(5,:), 'EdgeColor', coc_prop(6,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (50 \mum)');
set(gca,'Box','off');
try
    subplot(1,3,2)
    hdl(1) = cdfplot(HI(:,1)); hold all
    hdl(2) = cdfplot(HI(:,2)); hold all
    hdl(3) = cdfplot(HI(:,3)); hold all
    hdl(4) = cdfplot(HI(:,4)); hold all
    hdl(5) = cdfplot(HI(:,5)); hold all
    set(hdl(1), 'Color', coc_prop(1,:), 'LineWidth', 2)
    set(hdl(2), 'Color', coc_prop(2,:), 'LineWidth', 2)
    set(hdl(3), 'Color', coc_prop(3,:), 'LineWidth', 2)
    set(hdl(4), 'Color', coc_prop(4,:), 'LineWidth', 2)
    set(hdl(5), 'Color', coc_prop(5,:), 'LineWidth', 2)
    grid off
    xlabel('Homogeneity Index')
    xlim([0 1])
    ylabel('Cumulative fraction of cells')
    legend('0-50', '50-100', '100-150', '150-200', '200-250', 'Location', 'NorthWest'); legend('boxoff')
    set(gca,'Box','off');
    title('')
end

subplot(1,3,3)
avg_HI = nanmean(HI);
sem_HI = nanstd(HI)/sqrt(length(HI));
try
    errorbar(radius(2:end), avg_HI, sem_HI,'o-', 'Color', coc_prop(6,:), 'MarkerFaceColor', coc_prop(6,:))
end
xlabel('Maximal radius in \mum')
ylabel('Homeogeneity Index')
xlim([0 250])
ylim([0 1])
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'HomeogeneityIndex.png'))

%plot all ROIs
if analysisParams.plotROIs
    for i = 1:length(data.roi)
        PlotAvgStimResponseOri(metadata, analysis, field, i)
        saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
        close gcf
    end
    for i = 1:length(data.roi)
        PlotTrialStimResponse(metadata, analysis, field, i)
        saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
        close gcf
    end
end

%%
if analysisParams.predictor
    disp('Predicting stimulus from cell responses')
    % Get response data for cell using the stimResponseTrace field    
    interval = analysisPeriod;
    nCells = length(analysis.dff.roi);
    [nConds, nTrials, ~] = size(analysis.dff.roi(1).stimResponseTrace);
    isSig     = false(nCells,1);
    cellResps = zeros(nCells,(nConds-1)/2,2*nTrials);
    for n = 1:nCells
        dff = data.roi(n).dff;
        trialResps = analysis.dff.roi(n).stimResponseTrace;
        trialResps = nanmean(trialResps(1:(end-1),:,interval),3); %this is for decoding the direction, size numberDir x numTrials
        OriTrialResps = cat(2,trialResps(1:(nConds-1)/2,:),trialResps(((nConds-1)/2+1):end,:)); %this is for decoding the orientation, size numberOri x numTrials
        zscoredResps = (trialResps-nanmean(dff))./nanstd(dff);
        zscoredOriResps = (OriTrialResps-nanmean(dff))./nanstd(dff);
        cellResps(n,:,:) = zscoredOriResps; 
        isSig(n) = analysis.dff.roi(n).isResponseSignificant;
    end
    % Setup MATLAB anomynous functions to evalute decoding on real/shuffled data
    decodingStyle = 'similarity'; % 1-using similarity metric from Vince et al.,2-using correlation similarity
    switch decodingStyle
        case 'similarity'
            decoder = @(x,y) nansum(x(:).*y(:))./(sqrt(nansum(x(:).^2))*sqrt(nansum(y(:).^2)));
        case 'correlation'
            decoder = @(x,y) nansum(zscore(x(:)).*zscore(y(:)));
    end
    collapse3dMatrix = @(matrix,matrixSize) reshape(matrix,[matrixSize(1) matrixSize(2)*matrixSize(3)]);
    shuffle2dMatrix  = @(matrix) matrix(:,randperm(size(matrix,2)));
    restore3dMatrix  = @(matrix,matrixSize) reshape(matrix,[matrixSize(1) matrixSize(2) matrixSize(3)]);
    shuffle3dMatrix  = @(matrix) restore3dMatrix(shuffle2dMatrix(collapse3dMatrix(matrix,size(matrix))),size(matrix));

    % Process data with loops
    nBootstraps = 1000;
    nCellsUsed  = 50; % we used 50 in our paper, but could be nCells
    rng(1);
    medianTemplateResps = nanmedian(cellResps,3);
    
    decodingData  = {};
    for isShuffle = 0:1
        decoderPerformance  = zeros(nBootstraps,1);
        parfor(loop = 1:nBootstraps) % Can make a par for loopp
            % Pick a random subset of cells
            cellSubset = ceil(nCells*rand(nCells,1));
            cellSubset = cellSubset(1:nCellsUsed); 

            % Apply decoder on the selected population responses to each
            % stimulus trial
            if isShuffle
                x = shuffle3dMatrix(cellResps);
                y = nanmedian(shuffle3dMatrix(cellResps),3);
                %disp('Shuffled data')
            else
                x = cellResps;
                y = medianTemplateResps;
                %disp('Actual Data')
            end
            decodingMatrix = zeros((nConds-1)/2,nTrials,(nConds-1)/2);
            for cond = 1:(nConds-1)/2
                for trial = 1:nTrials
                    for refCond = 1:(nConds-1)/2
                        decodingMatrix(cond,trial,refCond) = decoder(x(cellSubset,cond,trial),y(cellSubset,refCond));
                    end
                end
            end

            % Evaluate decoding performance across all stimulus trials
            [~,predictedTrialLabel] = max(decodingMatrix,[],3);
            actualTrialLabel = repmat(1:(nConds-1)/2,[nTrials 1])';
            correctTrials = predictedTrialLabel==actualTrialLabel;  
            decoderPerformance(loop) = mean(correctTrials(:));
        end
        decodingData{isShuffle+1}=decoderPerformance; % Save the decoder for actual data and shuffle
    end
    figure    
    DecoderPerformanceAll = [decodingData{1}, decodingData{2}];
    groups = {'Data', 'Shuffle'};
    boxplot(DecoderPerformanceAll, 'labels', groups);
    ylim([0 1])
    ylabel('Fraction of correctly decoded trials');
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, 'Decoder.png'))
end

save(fullfile(saveDirectory, 'oriGratingAna.mat'), 'data', 'metadata', 'analysisParams', 'analysis');

end

