function Ori_Sf_S2p_Axon(animal, expt_id, sp2id, name, reloadData, plotROIs, server)
% 
% animal = 'F2363_2019-09-13'; 
% expt_id = 't00040';
% sp2id = expt_id;
% name = 'Fov_3';
% reloadData = 0;
% plotROIs = 0;
close all
if server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end


TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

windowStop=2;
windowStart=0;
pre=1;
field = 'dff';

z_thresh = 5;
fraction = 0.5;
shufflenum = 100;
predictor = 1;
plotROIsResps = plotROIs;

base2pDirectory= [TwoPhontondir animal];
tifDirectory = [base2pDirectory filesep name];
Sp2dDirectory = [Sp2dir animal filesep sp2id filesep];
saveDirectory = [savedir animal filesep expt_id filesep];
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
if ~exist(ROIsaveDirectory, 'dir')
    mkdir(ROIsaveDirectory);  
end
if ~exist(ROIRespsaveDirectory, 'dir')
    mkdir(ROIRespsaveDirectory);  
end
if ~exist(ROINonRespsaveDirectory, 'dir')
    mkdir(ROINonRespsaveDirectory);  
end

disp('Loading data')

%% load Data and metadata
if reloadData
    sliceparams = struct;
    sliceparams.expt_id = expt_id;
    sliceparams.baseDirectory = base2pDirectory;
    metadata.StimParams=Load_stimparams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    metadata.StimParams.series=expt_id;
    data = Suite2pAxonExtractorFct(animal,name, server);
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');
    metadata.ROI = struct;
    analysis = struct;
    save(fullfile(saveDirectory, 's1_ori_sf.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
else
    load(fullfile(saveDirectory, 's1_ori_sf.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
end

%% chop traces
disp('Chopping Traces')
[analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, 'pre', pre, 'post',metadata.StimParams.isi,'windowStart',windowStart, 'windowStop',windowStop);

%% find maxResponses and significant responses
disp('Calculating significant responses')
sigVector=zeros(1,length(analysis.(field).roi));

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
    analysisPeriod=(analysis.(field).windowStart:analysis.(field).windowStop);
    stimResp = analysis.(field).roi(i).stimResponseTrace(:,:,analysisPeriod);
    analysis.(field).roi(i).peaks = max(stimResp,[],3);
    analysis.(field).roi(i).zscore = ([analysis.(field).roi(i).peaks]-[analysis.(field).roi(i).baselineMean])./[analysis.(field).roi(i).baselineSD];
    analysis.(field).roi(i).crosser = sum(analysis.(field).roi(i).zscore > z_thresh,2);
    analysis.(field).roi(i).respStim = analysis.(field).roi(i).crosser > ((metadata.StimParams.numTrials)*fraction);
    if sum(analysis.(field).roi(i).respStim) > 0
        analysis.(field).roi(i).isResponseSignificant = 1;
    else 
        analysis.(field).roi(i).isResponseSignificant = 0;
    end
end

disp('Calculating ROI properties over all stimuli')
numStims= metadata.StimParams.numOrientations;
metadata.StimParams.theta= [0:2*pi/(numStims) : 2*pi-2*pi/(numStims)];

%% find preferred orientation for all cells
theta = metadata.StimParams.theta;
theta=mod(theta, pi);
lastOrientation = metadata.StimParams.numOrientations;
orientations = linspace(0,180,numStims/2+1);
orientations = orientations(1:end-1);
%find out if it is across several directions or if it is a sweep with only
%few orientations
if max(metadata.StimParams.directions > 180)
    sweep = 0;
else
    sweep = 1;
end
metadata.StimParams.sweep = sweep;
for i = 1:length(data.roi)
    %collapse all sf for the same orientation 
    spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
    metadata.StimParams.numSf = size(spatFreq_cell,2);
    medResponse = zeros(metadata.StimParams.numOrientations,metadata.StimParams.numTrials,size(analysis.(field).roi(i).stimResponseTrace,3));
    endnum = metadata.StimParams.numSf * metadata.StimParams.numOrientations+1;
    for ori = 1:metadata.StimParams.numOrientations
        medResponse(ori,:,:) = mean(analysis.(field).roi(i).stimResponseTrace(ori:metadata.StimParams.numOrientations:endnum, 1:metadata.StimParams.numTrials, :));
    end
    %calculate parameters
    medResponse = mean(medResponse(:,:,stimWindow),3);
    medResponse_dir = median(medResponse,2);
    medResponseA(1:size(medResponse,1)/2, 1:size(medResponse,2)) = medResponse(1:size(medResponse,1)/2,:);
    medResponseA(1:size(medResponse,1)/2, 1+size(medResponse,2):2*size(medResponse,2))= medResponse(1+size(medResponse,1)/2:size(medResponse,1),:);
    medResponseA= median(medResponseA,2);
    medResponse = medResponseA(:)';
    theta=theta(1:lastOrientation/2);
    theta= repmat(theta, 1, length(medResponse)/length(theta));
    
    %find preferred direction for all cells
    bestCoeff=[NaN, NaN, NaN, NaN];
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    bestCoeff=fitVonMisesLinSumFunction(medResponse_dir', metadata.StimParams.theta, bestCoeff, options);

    if sweep
        xo=[0:pi/40:pi];
        fit=vonMisesLinSum(bestCoeff, xo);
        if bestCoeff(1) > bestCoeff(2)
            analysis.(field).roi(i).preferredDirection = rad2deg(bestCoeff(4));
        else
            analysis.(field).roi(i).preferredDirection = mod(rad2deg(bestCoeff(4))+90,180);
        end
        analysis.(field).roi(i).preferredOrientation = analysis.(field).roi(i).preferredDirection;
        analysis.(field).roi(i).OSIFit = 1;
        interval = metadata.StimParams.directions(2) - metadata.StimParams.directions(1);
        ind = round(analysis.(field).roi(i).preferredOrientation / interval);
        ind = ind + 1;
        if ind > length(metadata.StimParams.directions) || isnan(ind)
            ind = 1;
        end
        analysis.(field).roi(i).prefOriStimInd = ind;
    else
        xo=[0:pi/40:2*pi];
        fit=vonMisesLinSum(bestCoeff, xo);
        if bestCoeff(1) > bestCoeff(2)
            analysis.(field).roi(i).preferredDirection = rad2deg(bestCoeff(4));
        else
            analysis.(field).roi(i).preferredDirection = mod(rad2deg(bestCoeff(4))+180,360);
        end
        [analysis.(field).roi(i).preferredOrientation,analysis.(field).roi(i).coeffOr,analysis.(field).roi(i).rsqOr] = ComputePreferredOrientations(medResponse, theta);
        analysis.(field).roi(i).preferredOrientation =mod(analysis.(field).roi(i).preferredOrientation, 180);
        analysis.(field).roi(i).rsqOr= Rsquared(medResponse, vonMisesFit(analysis.(field).roi(i).coeffOr, theta*2), true);
        analysis.(field).roi(i).OSIFit = computeOSIFit(analysis.(field).roi(i).coeffOr);
        [~, ind] = min(abs(orientations-analysis.(field).roi(i).preferredOrientation));
        analysis.(field).roi(i).prefOriStim = orientations(ind);
        analysis.(field).roi(i).prefOriStimInd = ind;
        medResponse = mean(analysis.(field).roi(i).avgStimResponse, 2)';
        analysis.(field).vsOrientation(i)= wrapTo2Pi(angle(vectorSum(abs(medResponse(1:lastOrientation)),2)))/2 * 180/pi;
        analysis.(field).orientationMag(i) = abs(vectorSum(abs(medResponse(1:lastOrientation))/max(medResponse(1:lastOrientation)),2));                

    end
end

%% find preferred sF for each cell by looking at preferred Orientation
for i = 1:length(data.roi)
    Response_sf_prefOri = zeros(metadata.StimParams.numSf,metadata.StimParams.numTrials,size(analysis.(field).roi(i).stimResponseTrace,3));
    for sf = 1:metadata.StimParams.numSf
        temp = analysis.(field).roi(i).stimResponseTrace(1+(sf-1)*ori:sf*ori, 1:metadata.StimParams.numTrials, :);
        if sweep
            Response_sf_prefOri(sf, :, :) = temp(analysis.(field).roi(i).prefOriStimInd,:,:);
        else
            Response_sf_prefOri(sf, :, :) = (temp(analysis.(field).roi(i).prefOriStimInd,:,:)+temp(analysis.(field).roi(i).prefOriStimInd+numStims/2,:,:)/2);
        end
    end
    avgResponse_sf_prefOri = mean(Response_sf_prefOri(:,:,stimWindow),3);
    medResponse_sf_prefOri = mean(avgResponse_sf_prefOri,2)';
    PrefSFInd = find(medResponse_sf_prefOri == max(medResponse_sf_prefOri));
    maxResp = max(medResponse_sf_prefOri);
    minResp = min(medResponse_sf_prefOri);
    SFSI = (maxResp - minResp) ./ (maxResp + minResp);
    analysis.(field).roi(i).prefSfStim = PrefSFInd;
    analysis.(field).roi(i).prefSf = spatFreq_cell{PrefSFInd};
    analysis.(field).roi(i).SFSI = SFSI;
end

%% compute DSI and OSI at prefered SF
disp('Calculating ROI properties for each sf')
if sweep
    for i= 1:length(data.roi)
        prefSfResponse = analysis.(field).roi(i).stimResponseTrace((analysis.(field).roi(i).prefSfStim-1)*lastOrientation+1:analysis.(field).roi(i).prefSfStim*lastOrientation, 1:metadata.StimParams.numTrials, :);
        prefSfResponse = mean(prefSfResponse(:,:,stimWindow),3);
        prefSfResponse_dir = median(prefSfResponse,2)';
        [analysis.(field).roi(i).DSI, ~] = computeDSI(lastOrientation, prefSfResponse_dir, 'preferred');
        analysis.(field).roi(i).OSI = analysis.(field).roi(i).DSI;
        analysis.(field).roi(i).OSIFit = analysis.(field).roi(i).DSI;
    end
else
    for i=1:length(data.roi)
        prefSfResponse = analysis.(field).roi(i).stimResponseTrace((analysis.(field).roi(i).prefSfStim-1)*lastOrientation+1:analysis.(field).roi(i).prefSfStim*lastOrientation, 1:metadata.StimParams.numTrials, :);
        prefSfResponse = mean(prefSfResponse(:,:,stimWindow),3);
        prefSfResponse_dir = median(prefSfResponse,2)';
        prefSfResponse(1:size(prefSfResponse,1)/2, 1:size(prefSfResponse,2)) = prefSfResponse(1:size(prefSfResponse,1)/2,:);
        prefSfResponse(1:size(prefSfResponse,1)/2, 1+size(prefSfResponse,2):2*size(prefSfResponse,2))= prefSfResponse(1+size(prefSfResponse,1)/2:size(prefSfResponse,1),:);
        prefSfResponse= median(prefSfResponse,2);
        prefSfResponse_ori = prefSfResponse(:)';
        [analysis.(field).roi(i).OSI, ~] = computeOSI(lastOrientation/2,prefSfResponse_ori(1:lastOrientation/2));
        [analysis.(field).roi(i).DSI, ~] = computeDSI(lastOrientation, prefSfResponse_dir, 'preferred');
    end
end

%% compute ori preference at each sf

prefDir_allSf = zeros(metadata.StimParams.numSf,1);
prefOri_allSf = zeros(metadata.StimParams.numSf,1);
OSIFit_allSf = zeros(metadata.StimParams.numSf,1);
OSI_allSf = zeros(metadata.StimParams.numSf,1);
DSI_allSf = zeros(metadata.StimParams.numSf,1);
isResp_allSf = zeros(metadata.StimParams.numSf,1);

for i = 1:length(data.roi)
    for sf = 1:metadata.StimParams.numSf
        medResponse = analysis.(field).roi(i).stimResponseTrace((sf-1)*lastOrientation+1:sf*lastOrientation, :,:);
        medResponse = mean(medResponse(:,:,stimWindow),3);
        dirResponse = median(medResponse,2)';
        medResponseA(1:size(medResponse,1)/2, 1:size(medResponse,2)) = medResponse(1:size(medResponse,1)/2,:);
        medResponseA(1:size(medResponse,1)/2, 1+size(medResponse,2):2*size(medResponse,2))= medResponse(1+size(medResponse,1)/2:size(medResponse,1),:);
        medResponseA= median(medResponseA,2);
        medResponse = medResponseA(:)';
        
        bestCoeff=[NaN, NaN, NaN, NaN];
        options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
        bestCoeff=fitVonMisesLinSumFunction(dirResponse, metadata.StimParams.theta, bestCoeff, options);
        if sweep
           xo=[0:pi/40:pi];
           fit=vonMisesLinSum(bestCoeff, xo);
           if bestCoeff(1) > bestCoeff(2)
               prefDir_allSf(sf) = rad2deg(bestCoeff(4));
           else
               prefDir_allSf(sf) = mod(rad2deg(bestCoeff(4))+90,180);
           end
           DSI_allSf(sf) = computeDSI(lastOrientation, dirResponse, 'preferred');
           prefOri_allSf(sf) = prefDir_allSf(sf);
           OSIFit_allSf(sf) = DSI_allSf(sf);
           OSI_allSf(sf) = DSI_allSf(sf);
        else
            xo=[0:pi/40:2*pi];
            fit=vonMisesLinSum(bestCoeff, xo);
            if bestCoeff(1) > bestCoeff(2)
                prefDir_allSf(sf) = rad2deg(bestCoeff(4));
            else
                prefDir_allSf(sf) = mod(rad2deg(bestCoeff(4))+180,360);
            end
            DSI_allSf(sf) = computeDSI(lastOrientation, dirResponse, 'preferred');
            [prefOri_allSf(sf), coeffOr, ~] = ComputePreferredOrientations(medResponse, theta);
            prefOri_allSf(sf) =mod(prefOri_allSf(sf), 180);
            OSIFit_allSf(sf) = computeOSIFit(coeffOr);
            OSI_allSf(sf) = computeOSI(lastOrientation/2,medResponse(1:lastOrientation/2));
        end
        tempCrosser = analysis.(field).roi(i).crosser((sf-1)*lastOrientation+1:sf*lastOrientation);
        tempRespStim = tempCrosser > ((metadata.StimParams.numTrials)*fraction);
        if sum(tempRespStim) > 0
            isResp_allSf(sf) = 1;
        else
            isResp_allSf(sf) = 0;
        end
    end
    analysis.(field).roi(i).preferredOrientation_allSf = prefOri_allSf;
    analysis.(field).roi(i).preferredDirection_allSf = prefDir_allSf;
    analysis.(field).roi(i).OSIFit_allSf = OSIFit_allSf;
    analysis.(field).roi(i).OSI_allSf = OSI_allSf;
    analysis.(field).roi(i).DSI_allSf = DSI_allSf;
    analysis.(field).roi(i).isResp_allSf = isResp_allSf;
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

%plot prefOri on top of template
for types = 1:3
    if types == 1
        rois_ori = linspace(1,length(analysis.dff.roi),length(analysis.dff.roi));
        alloriprefs = [analysis.dff.roi.preferredOrientation];
        alldirprefs = [analysis.dff.roi.preferredDirection];
        allsfpref = [analysis.dff.roi.prefSfStim];
        rois_dir = rois_ori;
        rois_sf = rois_ori;
    elseif types == 2
        rois_ori = find([analysis.dff.roi.isResponseSignificant] == 1); 
        alloriprefs = [analysis.dff.roi(rois_ori).preferredOrientation];
        alldirprefs = [analysis.dff.roi(rois_ori).preferredDirection];
        allsfpref = [analysis.dff.roi(rois_ori).prefSfStim];
        rois_dir = rois_ori;
        rois_sf = rois_ori;
    elseif types == 3
        rois_ori = find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alloriprefs = [analysis.dff.roi(rois_ori).preferredOrientation];
        rois_dir = find([analysis.dff.roi.DSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alldirprefs = [analysis.dff.roi(rois_dir).preferredDirection];
        rois_sf = find([analysis.dff.roi.SFSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        allsfpref = [analysis.dff.roi(rois_sf).prefSfStim];
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    %plot prefOri on top of template
    subplot(2,3,1)
    PlotPrefOnTemplate(analysis, data, metadata,1, field,data.template, rois_ori)
    
    subplot(2,3,4)
    histogram(alloriprefs,linspace(0,180,5), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel(sprintf('Orientation preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)]/2)
    axis square;
    set(gca,'Box','off');
    
    subplot(2,3,2)
    PlotPrefOnTemplate(analysis, data, metadata,3, field,data.template, rois_dir)

    subplot(2,3,5)
    histogram(alldirprefs,linspace(0,360,lastOrientation+1), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel(sprintf('Direction preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)])
    axis square;
    set(gca,'Box','off');
    
    %plot prefSF on top of template
    subplot(2,3,3)
    PlotPrefOnTemplate(analysis, data, metadata,2, field,data.template, rois_sf)

    subplot(2,3,6)
    numSf = size(spatFreq_cell,2);
    sf_cat = cellfun(@str2double,spatFreq_cell);
    sf_counts = histcounts(allsfpref, numSf);
    bar(sf_cat, sf_counts, 'FaceColor', coc_prop(5,:), 'EdgeColor', coc_prop(6,:));
    title('Histogram')
    ylabel('Cells');
    xlabel(sprintf('Spatial frequency preference'));
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

for types = 1:2
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    rows = ceil(sqrt(metadata.StimParams.numSf));
    for sf = 1:metadata.StimParams.numSf
        if metadata.StimParams.numSf > 6
            subplot(rows,rows,sf)
        else
            subplot(2,rows,sf)
        end
        tempResp = [analysis.dff.roi.isResp_allSf];
        isResponsive = find(tempResp(sf,:) == 1);
        if sum(isResponsive) > 0
            if types == 1
                tempOSI= [analysis.dff.roi.OSIFit_allSf];
                tempOriSelect = find(tempOSI(sf,:) > 0.2);
                rois = intersect(tempOriSelect, isResponsive);
                if ~isempty(rois)
                    tempPrefOri = [analysis.dff.roi(rois).preferredOrientation_allSf];
                    allprefs = tempPrefOri(sf,:);
                    tempFit = [analysis.dff.roi(rois).OSIFit_allSf];
                    allFit = tempFit(sf,:);
                    PlotPrefOnTemplateSF(allprefs, allFit,data, metadata,1, data.template,rois)
                    title(['Ori pref map, ' num2str(spatFreq_cell{sf}) ' cpd'])
                else
                    imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                    colormap(hsv)
                    title(['Ori pref map, ' num2str(spatFreq_cell{sf}) ' cpd'])
                    caxis([0 180])
                    colorbar('Location', 'southoutside');
                end
            elseif types == 2
                tempDSI= [analysis.dff.roi.DSI_allSf];
                tempDirSelect = find(tempDSI(sf,:) > 0.2);
                rois = intersect(tempDirSelect, isResponsive);
                if ~isempty(rois)
                    tempPrefDir = [analysis.dff.roi(rois).preferredDirection_allSf];
                    allprefs = tempPrefDir(sf,:);
                    tempFit = [analysis.dff.roi(rois).DSI_allSf];
                    allFit = tempFit(sf,:);
                    PlotPrefOnTemplateSF(allprefs, allFit,data, metadata,3, data.template,rois)
                    title(['Dir pref map, ' num2str(spatFreq_cell{sf}) ' cpd'])
                else
                    imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                    colormap(hsv)
                    title(['Dir pref map, ' num2str(spatFreq_cell{sf}) ' cpd'])
                    caxis([0 360])
                    colorbar('Location', 'southoutside');
                end
            end
        else
            imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
            colormap(hsv)
            if types == 1
                title(['Ori pref map, ' num2str(spatFreq_cell{sf}) ' cpd'])
                caxis([0 180])
            elseif types == 2
                title(['Dir pref map, ' num2str(spatFreq_cell{sf}) ' cpd'])
                caxis([0 360])
            end
             colorbar('Location', 'southoutside');
        end
    end
    if types == 1
        saveas(gcf, fullfile(saveDirectory, 'OrientationPreference_allSF.png'))
    elseif types == 2
        saveas(gcf, fullfile(saveDirectory, 'DirectionPreference_allSF.png'))
    end
end

%plot all ROIs
if plotROIs
    for i = 1:length(data.roi)
        PlotAvgStimResponse(metadata, analysis, field, i)
        saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
        close gcf
    end
    for i = 1:length(data.roi)
        PlotTrialStimResponse(metadata, analysis, field, i)
        saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
        close gcf
    end
end
if plotROIsResps
    for i = 1:length(data.roi)
        if analysis.(field).roi(i).isResponseSignificant == 1
            PlotTrialStimResponse(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        else
            PlotTrialStimResponse(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        end
    end
end

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
ori = length(find([analysis.dff.roi.OSI] > 0.33 & [analysis.dff.roi.isResponseSignificant] == 1)) ./all;
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
dir = length(find([analysis.dff.roi.DSI] > 0.33 & [analysis.dff.roi.isResponseSignificant] == 1)) ./all;
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

%% plot delta Ori vs. distance of ROIs
[zoom, setup] = getzoom(tifDirectory);
metadata.zoom = zoom;
if setup == 1
    fieldofview = 1000/zoom;
    umperpixel = fieldofview/512;
    disp('Setup Ben')
elseif setup == 2
    umperpixel = 2.73/zoom;
end


ori_sel = find([analysis.dff.roi.OSI] > 0.33 & [analysis.dff.roi.isResponseSignificant] == 1);
dist_ori = zeros(length(ori_sel), length(ori_sel));
delta_ori = zeros(length(ori_sel), length(ori_sel));
for A = 1:length(ori_sel)
    for B = 1:length(ori_sel)
        dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
        dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
        dist_ori(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
        delta_ori(A,B) = abs(analysis.dff.roi(A).preferredOrientation - analysis.dff.roi(B).preferredOrientation);
        if delta_ori(A,B) > 90
            delta_ori(A,B) = 180 - delta_ori(A,B);
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
bin_deltaOri_shuffle = zeros(shufflenum,bin);
sem_deltaOri_shuffle = zeros(shufflenum,bin);
for rep = 1:shufflenum
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
%calculate distance of ROIs to each other and check whether their preferred
%sf is the same
%%

% 
% sf_sel = find([analysis.dff.roi.SFSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
% dist_sf = zeros(length(sf_sel), length(sf_sel));
% sfPref_sf_sel = zeros(length(sf_sel), length(sf_sel));
% for A = 1:length(sf_sel)
%     for B = 1:length(sf_sel)
%         dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
%         dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
%         dist_sf(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
%         if analysis.(field).roi(A).prefSf == analysis.(field).roi(B).prefSf
%             sfPref_sf_sel(A,B) = 1;
%         else
%             sfPref_sf_sel(A,B) = 0;
%         end
%     end
% end
% 
% %calculate HeterogeneityIndex for sf as follow: 
% %HI = 1 - cells with same sf/all cells for different distances
% %if HI = 1 -> all cell with same sf, if HI = 0, all cells with different sf
% max_dist = max(max(dist_sf));
% num_dist = ceil(max_dist/50);
% radius = linspace(0,num_dist*50,num_dist+1);
% HI_sf = zeros(length(sf_sel),length(radius)-1);
% HI_sf_shuffle = zeros(length(sf_sel),length(radius)-1);
% sameSf_shuffle_all = zeros(1, shufflenum);
% for dis = 1:length(radius)-1
%     for cellID = 1:length(sf_sel)
%         distances = dist_sf(cellID,:);
%         withinRadius1 = distances >= radius(dis);
%         withinRadius2 = distances < radius (dis+1);
%         withinRadius = withinRadius1&withinRadius2;
%         withinRadius(cellID)=0;
%         sameSf = sfPref_sf_sel(cellID,:) & withinRadius;
%         for shuffle = 1:shufflenum
%             sfPref_sf_sel_shuffle = sfPref_sf_sel(randperm(length(sfPref_sf_sel(cellID,:))));
%             sameSf_shuffle_all(shuffle) = sum(sfPref_sf_sel_shuffle & withinRadius);
%         end
%         sameSf_shuffle = mean(sameSf_shuffle_all);
%         HI_sf(cellID, dis) = 1 - sum(sameSf)/sum(withinRadius);
%         HI_sf_shuffle(cellID, dis) = 1 - sum(sameSf_shuffle)/sum(withinRadius);
%     end
% end
% analysis.dff.ori_cells.HomeogeneityIndexSf = HI_sf;
% 
% figure
% avg_HI_sf = nanmean(HI_sf);
% sem_HI_sf = nanstd(HI_sf)/sqrt(length(HI_sf));
% avg_HI_sf_shuffle = nanmean(HI_sf_shuffle);
% sem_HI_sf_shuffle = nanstd(HI_sf_shuffle)/sqrt(length(HI_sf_shuffle));
% errorbar(radius(2:end), avg_HI_sf, sem_HI_sf,'o-', 'Color', coc_prop(6,:), 'MarkerFaceColor', coc_prop(6,:)); hold all
% errorbar(radius(2:end), avg_HI_sf_shuffle, sem_HI_sf_shuffle,'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
% xlabel('Maximal radius in \mum')
% ylabel('Homeogeneity Index Spatial Frequency')
% xlim([0 num_dist*50])
% ylim([0 1])
% legend('Data', 'Shuffled')
% legend('boxoff')
% set(gca,'Box','off');
% set(gcf, 'color', 'w');
% saveas(gcf, fullfile(saveDirectory, 'HeterogeneityIndexSF.png'))

save(fullfile(saveDirectory, 's1_ori_sf_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');

end

%% additional functions

function [twophoton] = LoadFrameTimes(path)
    %ExtractFrameTimes Reads in Frame Triggers from frametrigger.txt
    %  	Args:
    %     path: path of the metadata file to extract
    %   Returns:
    %     twophoton (struct):   .time :trigger times (s)
    %                           .rate :frame rate (Hz)

    twoPhotonFile = 'twophotontimes.txt';
    counter = 1;

    if ~exist(path, 'dir')
        error('%s does not exist', path);
    end
    
    twophoton = struct;
    twophoton.time = [];
    tpFullFile = fullfile(path, twoPhotonFile);
    infotpFullFile = dir(tpFullFile);
    if infotpFullFile.bytes < 10
        twoPhotonFile = 'frametrigger.txt';
        tpFullFile = fullfile(path, twoPhotonFile);
    end
    while ~exist(tpFullFile, 'file') || isempty(twophoton.time)
        if exist(tpFullFile, 'file')
            disp(['Loading... ', tpFullFile])
            twophoton.time = load(tpFullFile);
            twophoton.rate = 1/median(diff(twophoton.time));

        else
            if counter > size(defaultFileNames,1)
               warning ('Could not get frame acquisition triggers')
               return;
            end
            tpFullFile = fullfile(twoPhotonPath, char(defaultFileNames{counter}));
            counter = counter + 1;
        end
    end
end
function [stim_params]=Load_stimparams(exptpath)     
    stim_params=struct;
    stim_matrix=[];
    disp(strcat('Loading....',exptpath, '\stimontimes.txt'))
    stimtimes=load(strcat(exptpath, '\stimontimes.txt'));
    i=1;
    while i <length(stimtimes)
        idx= stimtimes(i);
        stim_matrix(1, floor(i/2)+1)= stimtimes(i);
        stim_matrix(2, floor(i/2)+1)= stimtimes(i+1);
        i=i+2;
    end
    stim_params.StimOnTimes=stim_matrix;
    stim_params.StimOnTimes=stim_params.StimOnTimes(:, 2:end); %discard the very first stim because its garbage
    files=dir(fullfile(exptpath, 'frametrigger.txt'));
    files={files.name};
    if ~isempty(files)
        frametimes=load(fullfile(exptpath, 'frametrigger.txt'));
    else
        frametimes=[];
    end
    stim_params.frameTimes=frametimes;
    disp(strcat('Loading....',exptpath, '\stimorientations.txt'))
    stimorientations=load(strcat(exptpath, '\stimorientations.txt'));
    stim_params.directions = stimorientations;
    if isempty(stim_params.StimOnTimes)
        stim_params.type = 'Spontaneous';
        return;
    else
        stim_params.uniqStimIds= unique(stim_params.StimOnTimes(1,1:end)); %always ignore the first stim code
        stim_params.uniqStims= length(stim_params.uniqStimIds);
        stim_params.numberOfStims=length(stim_params.StimOnTimes);
        files=dir(strcat(exptpath, '\*.py'));
        files={files.name};

        for i=1:length(files) % find the type of stimulus
            if isempty(strfind(char(files{i}), 'serialTriggerDaqOut'))
                stim_params.type= strrep(char(files{i}), '.py', '');
                stim_params.file= strcat(exptpath,'\',char(files{i}));
            end
        end
        if isempty(stim_params.type) || strcmp(stim_params.type, 'blackScreen') % check if there wasn't a
            stim_params.type= 'Spontaneous';
        end

        %extract all the relevant field for the stimulation
        stimFields= stimField(stim_params.type);
        for i=1:length(stimFields) %
            stim_params.(stimFields{i})= FindStimulusParam(stim_params.file, char(stimFields{i}));
        end
        %correct the stim ids for driftingGrating to ignore the
        if isfield(stim_params, 'flashInterval')
            stim_params.stimDuration = stim_params.flashInterval;
        elseif isfield(stim_params, 'stimDuration')
            stim_params.stimDuration = stim_params.stimDuration;
        else 
            stim_params.stimDuration = NaN;
        end
    end

    if strcmp(stim_params.type, 'driftingGrating')
        stims = stim_params.uniqStims;
        if stim_params.doBlank
            stims = stims-1;

        end
          stim_params.theta = (0:2*pi/stims:2*pi - 2*pi/stims);
    end
end
function y = percentileFilt1(x,percentile,n,blksz,DIM)
    narginchk(1,5);
    if nargin < 4, blksz = []; end
    if nargin < 5, DIM = []; end

    % Check the input data type. Single precision is not supported.
    % try
    %     chkinputdatatype(x,n,blksz,DIM);
    % catch ME
    %     throwAsCaller(ME);
    % end

    % Check if the input arguments are valid
    if isempty(n)
      n = 3;
    end

    if ~isempty(DIM) && DIM > ndims(x)
        error(message('signal:medfilt1:InvalidDimensions'))
    end

    % Reshape x into the right dimension.
    if isempty(DIM)
        % Work along the first non-singleton dimension
        [x, nshifts] = shiftdim(x);
    else
        % Put DIM in the first (row) dimension (this matches the order 
        % that the built-in filter function uses)
        perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
        x = permute(x,perm);
    end

    % Verify that the block size is valid.
    siz = size(x);
    if isempty(blksz)
        blksz = siz(1); % siz(1) is the number of rows of x (default)
    else
        blksz = blksz(:);
    end

    % Initialize y with the correct dimension
    y = zeros(siz); 

    % Call medfilt1D (vector)
    for i = 1:prod(siz(2:end))
        y(:,i) = prctilefilt1d(x(:,i),n,blksz,percentile);
    end

    % Convert y to the original shape of x
    if isempty(DIM)
        y = shiftdim(y, -nshifts);
    else
        y = ipermute(y,perm);
    end
end
function y = percentile(x, k)
    x = sort(x);
    n = size(x,1);

    p = 1 + (n-1) * k / 100;

    if p == fix(p)
        y = x(p);
    else
        r1 = floor(p); r2 = r1+1;
        y = x(r1) + (x(r2)-x(r1)) * k / 100;
    end
end
function y = prctilefilt1d(x,n,blksz,percentile)
    %PRCTILEFILT1D  One dimensional median filter.
    %
    % Inputs:
    %   x     - vector
    %   n     - order of the filter
    %   blksz - block size

    nx = length(x);
    if rem(n,2)~=1    % n even
        m = n/2;
    else
        m = (n-1)/2;
    end
    X = [zeros(m,1); x; zeros(m,1)];
    y = zeros(nx,1);

    % Work in chunks to save memory
    indr = (0:n-1)';
    indc = 1:nx;
    for i=1:blksz:nx
        ind = indc(ones(1,n),i:min(i+blksz-1,nx)) + ...
              indr(:,ones(1,min(i+blksz-1,nx)-i+1));
        xx = reshape(X(ind),n,min(i+blksz-1,nx)-i+1);
        y(i:min(i+blksz-1,nx)) = prctile(xx,percentile,1);
    end
end
function [zoom, setup] = getzoom(tifDirectory)
    filename = [tifDirectory filesep '*.tif'];
    files = dir(filename);
    filepath = [tifDirectory filesep files(1).name];
    InfoImage = imfinfo(filepath);
    try 
        a = InfoImage(1).Software;
        setup = 1;
    catch 
        a = InfoImage(1).ImageDescription;
        setup = 2;
    end
    zoom = regexp(a,'(?<=scanZoomFactor = )\d+\.?\d*', 'match');
    zoom = str2num(zoom{1});
end
function [fields] = stimField(Stimtype)
    switch Stimtype
        case 'driftingGrating'
            fields = {'numTrials',...
                'doBlank',...
                'changeDirectionAt',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'temporalFreq',...
                'spatialFreq',...
                'contrast',...
                'textureType',...
                'maxv',...
                'minv',...
                'dutyCycle',...
                'foregroundColor',...
                'startingPhase',...
                'centerPoint',...
                'stimSize',...
                'animalOrientation'
                };
        case 'driftingGrating_ori_sf'
            fields = {'numTrials',...
                    'doBlank',...
                    'changeDirectionAt',...
                    'stimDuration',...
                    'isi',...
                    'isRandom',...
                    'initialDelay',...
                    'temporalFreq',...
                    'spatialFreq',...
                    'contrast',...
                    'textureType',...
                    'numOrientations',...
                    'dutyCycle',...
                    'foregroundColor',...
                    'startingPhase',...
                    'centerPoint',...
                    'stimSize',...
                    'orientations',...
                    };
        otherwise
            fields ={};
    end
end
function [val]= FindStimulusParam(fname,paramKey)
    fid=fopen(fname);
    try
        % get the rows for the paramKey
        fileText= textscan(fid, '%s', 'delimiter', '\n');
        fileText= fileText{1};
        fileText= strrep(fileText, ' ', ''); % delete all whitespace
        keys = strfind(fileText, strcat(paramKey, '='));
        keys= find(~cellfun(@isempty, keys));
        line = fileText(keys(1));
        line = strsplit(char(line), '#');
        line = strsplit(char(line(1)), '=');
        val= str2num(char(line(2)));
        if isempty(val)
            val = char(line(2));
        end
    catch ME
        val= '';
    end
    fclose(fid);
end
function [metadata, data]= baselinePercentileFilter(metadata,data, fieldSource,fieldTarget,filteredCutoff, desiredPercentileRank)            
    numSeries = 1;
    fps= metadata.TwoPhoton.rate; %Hz
    if ~isfield(metadata, fieldTarget)
        metadata.baseline= struct;
    end
    if ~isfield(metadata.baseline, fieldSource)
        metadata.baseline.(fieldSource)= struct;
    end
    metadata.baseline.(fieldSource).filteredCutoff=filteredCutoff;
    metadata.baseline.(fieldSource).desiredPercentileRank= desiredPercentileRank;
    for i=1:length(data.roi)
        for j = 1:numSeries
            disp({'baselining cell #',  num2str(i)});
            if isfield(metadata.StimParams, 'seriesBreaks')
                startIdx = 1+ metadata.StimParams.seriesBreaks(j);
                stopIdx = metadata.StimParams.seriesBreaks(j+1);
                inputTrace = data.roi(i).(fieldSource)(startIdx:stopIdx)';
            else
                inputTrace = data.roi(i).(fieldSource)';
            end
            paddingLength = ceil(length(inputTrace)/1);

            paddedTrace   = [inputTrace(paddingLength:-1:1); inputTrace; inputTrace(paddingLength:-1:1)];
            filteredTrace = percentileFilt1(paddedTrace, desiredPercentileRank, round(filteredCutoff*fps));
            filteredTrace = filteredTrace(paddingLength+(1:length(inputTrace)));

            butterWorthOrder = 1;
            Wn = (1/filteredCutoff) / (fps/2);
            [b,a] = butter(butterWorthOrder, Wn, 'low');
            highpassFilteredTrace = filtfilt(b,a,[filteredTrace(paddingLength:-1:1); filteredTrace; filteredTrace(1:paddingLength)]);
            highpassFilteredTrace = highpassFilteredTrace(paddingLength+[1:length(inputTrace)]);
            if j==1
                data.roi(i).(fieldTarget)=highpassFilteredTrace';
            else
                data.roi(i).(fieldTarget)= horzcat(obj.data.roi(i).(fieldTarget), highpassFilteredTrace');
            end
        end
    end
end
function [data] = computeDff(data, sourceSignal, sourceBaseline, fieldTarget)
    for i=1:length(data.roi)
        data.roi(i).(fieldTarget)= (data.roi(i).(sourceSignal)-data.roi(i).(sourceBaseline))./abs(data.roi(i).(sourceBaseline));
    end
end
function [analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, varargin)
    if mod(length(varargin),2) ~= 0
        error('You need to specify both the field and the argument')
    elseif strcmp(metadata.StimParams.type, 'Spontaneous')
        warning('Unable to chop up Spontaneous activity into meaningful chunks')
        return
    else
        % defaults for values
        windowStart= 0;  % choose the entire stim period
        windowStop = 0.5;    % choose the entire stim period
        tfield = field; % default to field name as analysis name
        pre=0;
        post=1;
        skipTrial=0;
        baseline=false;
        for i = 1:2:length(varargin)
            switch varargin{i}
                case 'analysisName'
                    tfield = varargin{i+1};
                case 'windowStart'
                    windowStart= varargin{i+1};
                case 'windowStop'
                    windowStop=varargin{i+1};
                case 'pre'
                    pre = varargin{i+1};
                case 'post'
                    post = varargin{i+1};
                case 'skipTrial'
                    skipTrial = varargin{i+1};
                case 'baseline'
                    baseline=varargin{i+1};
                otherwise
                    error(['You have specified an invalid field: ', varargin{i}])
            end
        end
    end

    TwoPhotonRate= metadata.TwoPhoton.rate;
    analysis.(tfield)=struct;
    analysis.(tfield).preTrialTime = pre;
    analysis.(tfield).postTrialTime = post;
    offsetPre = round(metadata.TwoPhoton.rate * analysis.(tfield).preTrialTime);
    offsetPost = round(metadata.TwoPhoton.rate * analysis.(tfield).postTrialTime);
    analysis.(tfield).stimStart=offsetPre+1;
    analysis.(tfield).stimStop=round(analysis.(tfield).stimStart+metadata.TwoPhoton.rate * metadata.StimParams.stimDuration());
    analysis.(tfield).roi = struct;

    windowStartIdx= round(windowStart * metadata.TwoPhoton.rate +offsetPre+1);
    windowStopIdx = round(windowStop * TwoPhotonRate + offsetPre);
    if windowStartIdx == windowStopIdx
        windowStopIdx = analysis(tfield).stimStop;
    end
    analysis.(tfield).windowStart =windowStartIdx;
    analysis.(tfield).windowStop = windowStopIdx;

    %do a little clean up here to make stim numbers match
%     while (metadata.StimParams.StimOnTimes(2,metadata.StimParams.numberOfStims)+metadata.StimParams.stimDuration > metadata.TwoPhoton.time(end))
%         if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ==0
%             metadata.StimParams.numTrials = metadata.StimParams.numTrials;
%         else
%             metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
%         end
%         metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
%     end
% 
    if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ~=0
        metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
    end

    for i =1:metadata.StimParams.numberOfStims
        metadata.StimParams.stimStartIndex(i)= find(metadata.TwoPhoton.time > metadata.StimParams.StimOnTimes(2,i),1);
        metadata.StimParams.stimStopIndex(i)= floor(metadata.StimParams.stimStartIndex(i)+ metadata.StimParams.stimDuration* metadata.TwoPhoton.rate);
    end

    stimPeriod=(analysis.(tfield).stimStart:analysis.(tfield).stimStop);
    analysisPeriod=(analysis.(tfield).windowStart:analysis.(tfield).windowStop);

    %correct for our specified window

    if windowStartIdx > max(stimPeriod) || windowStopIdx > max(stimPeriod)
        warning('You have choosen a period to average over that is longer than the stim period')
    end

    stimStarts=metadata.StimParams.stimStartIndex;
    stimStops=metadata.StimParams.stimStopIndex;
    
    StimonTimes= metadata.StimParams.StimOnTimes;
    %StimonTimes=StimonTimes(:,1:metadata.StimParams.numberOfStims);
    
    for stimID = 1:metadata.StimParams.uniqStims
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        stimIndexTrials(stimID) = length(find(StimonTimes(1,:)==stimulus));
    end
    if metadata.StimParams.numTrials > min(stimIndexTrials)
       metadata.StimParams.numTrials =  min(stimIndexTrials);
    end
    
    for stimID= 1:metadata.StimParams.uniqStims
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        stimIndices = find(StimonTimes(1,:)==stimulus);
        stimIndices = stimIndices(1:metadata.StimParams.numTrials);
        for trialNumber= 1:metadata.StimParams.numTrials 
            stimIndex = stimIndices(trialNumber);
            for i=1:length(data.roi)
                selectedFramesTrace = (stimStarts(stimIndex)-offsetPre):(stimStops(stimIndex)+offsetPost);
                if max(selectedFramesTrace) <= length(data.roi(i).(field))
                    analysis.(tfield).roi(i).stimResponseTrace(stimID,trialNumber,:)= data.roi(i).(field)(selectedFramesTrace);
                else
                    analysis.(tfield).roi(i).stimResponseTrace(stimID,trialNumber,:) = zeros(1,1,length(selectedFramesTrace));
                end
            end
        end

        for i=1:length(data.roi)
            if skipTrial >0 && skipTrial <=size(analysis.(tfield).roi(i).stimResponseTrace,2)
                analysis.(tfield).roi(i).stimResponseTrace(:,skipTrial,:)=[];
            end
            traces=squeeze(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:));
            bltraces=repmat(mean(traces(:,1:(stimPeriod(1)-1)),2)', size(traces,2),1)';

            if baseline
                traces=traces- bltraces;
            end
            try
                analysis.(tfield).roi(i).stimResponseTraceRaw(stimID,:,:)=analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:);
            catch ME
                ME.message
            end
            analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:)= traces;
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = mean(traces, 1);
            n= size(analysis.(tfield).roi(i).stimResponseTrace,2);
            y=analysis.(tfield).roi(i).stimResponseTrace(stimID,:,:);
            analysis.(tfield).roi(i).SEMResponseTrace(stimID,:) = std(y,[],2)/sqrt(n);
            analysis.(tfield).roi(i).avgStimResponse(stimID,:) = mean(analysis.(tfield).roi(i).avgResponseTrace(stimID,analysisPeriod),2);
            n= size(analysis.(tfield).roi(i).stimResponseTrace,2);
            y=analysis.(tfield).roi(i).avgStimResponse(stimID,:);
            analysis.(tfield).SEMCellResponses(i,stimID)= std(y,[],2)/sqrt(n);
            analysis.(tfield).avgCellResponses(i,stimID)= mean(analysis.(tfield).roi(i).avgStimResponse(stimID,:),2);
        end
    end

    for i = 1:length(data.roi)
        for stimID = 1:length(metadata.StimParams.uniqStimIds)
            analysis.(tfield).roi(i).avgStimResponse(stimID) = analysis.(tfield).roi(i).avgStimResponse(stimID) -analysis.(tfield).roi(i).avgStimResponse(end) ;
            analysis.(tfield).roi(i).avgResponseTraceRaw(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:);
            analysis.(tfield).roi(i).avgResponseTrace(stimID,:) = analysis.(tfield).roi(i).avgResponseTrace(stimID,:) - analysis.(tfield).roi(i).avgResponseTrace(end,:);
            analysis.(tfield).roi(i).Fb(stimID,:)= mean(analysis.(tfield).roi(i).stimResponseTrace(stimID,:,1:analysisPeriod(1)-1),3);

        end
    end
    analysis.(tfield).preTrialTime=pre;
end        
function [prefOri, coeffOr, rsqPOr] = ComputePreferredOrientations(medResponse, theta)
    theta = theta *2;
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    coeffOr=[NaN, NaN, NaN, NaN];
    [coeffOr, rsqPOr]=fitVonMisesFunction(medResponse, theta, coeffOr, options);
    prefOri = rad2deg(coeffOr(3)/2);
    hold off;
    plot(rad2deg(theta/2), medResponse, 'bo')
    hold on;
    plot(rad2deg([0:pi/32:pi]), vonMisesFit(coeffOr, [0:pi/16:2*pi]), 'r');
    title([num2str(prefOri),':',num2str(rsqPOr)]);
    drawnow()
end
function [bestCoeff, bestrsq]=fitVonMisesFunction(ydata, tdata, bestCoeff, options)
    stims=length(ydata);
    fun=@(x, tdata) vonMisesFit(x,tdata);

    prefdir=mod(tdata(ydata == max(ydata(:))), pi);
    if length(prefdir) >1
        prefdir=prefdir(1);
    end
    %coeff0=[abs(max(ydata(:))), abs(max(ydata(:))), pi/2, pi/stims];
    
    lb=[0,0,0,0];
    ub=[abs(1.5*max(ydata(:))),abs(max(ydata(:))*2),2*pi,2.5*pi];
    coeff0=(ub-lb)./2 + lb;
    rsq=0;
    bestrsq=0;
    
    iterations=0;
    try
    if sum(isnan(lb))==0 && sum(isnan(ub)) ==0
        while rsq < 0.7 && iterations <5
            

            [currentCoeff, ~, ~,exitflag]=lsqcurvefit(fun, coeff0, tdata, ydata, lb, ub, options);

            [rsq, ~]=Rsquared(ydata, vonMisesFit(currentCoeff, tdata), true);

             
            if rsq > bestrsq
                bestCoeff= currentCoeff;
                bestrsq=rsq;
            end
            iterations=iterations+1;
            coeff0= ub-lb .* rand(1) + lb;
        end
        %disp(['Ran ', num2str(iterations),  ' iterations']);
    else
        bestCoeff=[NaN, NaN, NaN, NaN];
    end 
    catch
        bestCoeff=[NaN,NaN,NaN,NaN];
        rsq=0;
    end
        
end
function bestCoeff=fitVonMisesLinSumFunction (ydata, tdata, bestCoeff,  options)
    stims=length(ydata);
    fun=@(x, tdata)vonMisesLinSum(x,tdata);

    prefdir=mod(tdata(ydata == max(ydata(:))), pi);
    if length(prefdir) >1
        prefdir=prefdir(1);
    end
    coeff0=[abs(max(ydata(:))), abs(max(ydata(:))),abs(max(ydata(:))), prefdir, 2*pi/stims,2*pi/stims];

    lb=[min(abs(min(ydata(:))),0),min(abs(min(ydata(:))),0),0,0,0,0];
    ub=[abs(1.5*max(ydata(:))),abs(1.5*max(ydata(:))),abs(max(ydata(:))*2),pi,10,10];
    rsq=0;
    
    % if we've undersampled- double the length of ydata and tdata
    if length(ydata) < length(lb)
        tdata_old = tdata;
        ydata_old = ydata;
        tdata = linspace(0, tdata_old(end), length(lb)*length(ydata));
        ydata = interp1(tdata_old, ydata_old, tdata);
    end
    

    iterations=0;
    if sum(isnan(lb))==0 && sum(isnan(ub)) ==0
        while rsq < 0.7 && iterations <10
            lastrsq=rsq;

            [currentCoeff, ~, ~,exitflag]=lsqcurvefit(fun, coeff0, tdata, ydata, lb, ub, options);

            [rsq, ~]=Rsquared(ydata, vonMisesLinSum(currentCoeff, tdata), true);
            if rsq > lastrsq
                bestCoeff= currentCoeff;
            end
            iterations=iterations+1;
            coeff0= ub-lb .* rand(1) + lb;
        end
        %disp(['Ran ', num2str(iterations),  ' iterations']);
    else
        bestCoeff=[NaN, NaN, NaN, NaN, NaN,NaN];
    end
end
function [r2 rmse] = Rsquared(y,f,varargin)
    if isempty(varargin); c = true; 
    elseif length(varargin)>1; error 'Too many input arguments';
    elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
    else c = varargin{1}; 
    end

    % Compare inputs
    if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

    % Check for NaN
    tmp = ~or(isnan(y),isnan(f));
    y = y(tmp);
    f = f(tmp);

    if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
    else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
        if r2<0
        % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
            warning('Consider adding a constant term to your model') %#ok<WNTAG>
            r2 = 0;
        end
    end

    rmse = sqrt(mean((y(:) - f(:)).^2));
end
function [OSIFit] = computeOSIFit(coeffOr)
    Rpref= vonMisesFit(coeffOr, coeffOr(3));
    Rnull = vonMisesFit(coeffOr, coeffOr(3) + pi);
    OSIFit = (Rpref - Rnull) / (Rpref + Rnull);
end
function [ out ] = vonMisesFit(x, tdata)
    %VONMISESFIT Linear Sum of two vonMises functions that's constrained to be
    %pi radians apart
    %   x is a vector of the parameters of the linear sum and vonMises
    %   parameters
    %   tdata is the x range over which to compute the vonMises
    %   out = [Alpha Beta mu kappa]
    if sum(isnan(x)) > 0 || length(x) <4
        out = NaN * ones(size(tdata));
        return
    end
    A= x(1);
    B= x(2);  % a dc component
    mu1=x(3); % center
    kappa1=x(4); % width
    out = A*vonMisesFunction(kappa1, mu1, tdata ) + B;
end
function [ out ] = vonMisesLinSum(x, tdata)
%VONMISESFIT Linear Sum of two vonMises functions that's constrained to be
%pi radians apart
%   x is a vector of the parameters of the linear sum and vonMises
%   parameters
%   tdata is the x range over which to compute the vonMises

if sum(isnan(x)) > 0 || length(x) <6
    out = NaN * ones(size(tdata));
    return
end

A= x(1); % first direction
B=x(2);  % 2nd direction
D=x(3);  % a dc component
mu1=x(4);
kappa1=x(5);
kappa2=x(6);



mu2=mod(mu1 + pi, 2*pi);

out = A*vonMisesFunction(kappa1, mu1, tdata ) + B*vonMisesFunction(kappa2, mu2, tdata )+D;

end
function [ out ] = vonMisesFunction(kappa, mu, x )
%VONMISESFUNCTION Summary of this function goes here
%   Detailed explanation goes here

    out = exp(kappa * cos (x-mu))/ (2 * pi * besseli(0,kappa));
end
function summedVector = vectorSum(array,harmonic,dim)
    % summedVector = vectorSum(stack,harmonic,dim)
    % vector sum of input array with 2nd harmonic response 
    % along the last dimension of the array (unless 
    % otherwise specified by user input). 

    if(nargin<3), dim = length(size(array)); end
    if(nargin<2), harmonic = 2;              end

    stackSize = size(array);
    phaseArraySize = 1+0*stackSize; phaseArraySize(dim) = stackSize(dim);

    vectorArraySize = stackSize;    vectorArraySize(dim) = 1;
    phaseValues = linspace(0,360,stackSize(dim)+1); 
    summedVector = sum(array.*repmat(reshape(exp(2*pi*1i*harmonic*mod(phaseValues(1:(end-1)),360/harmonic)/360),phaseArraySize),vectorArraySize),dim);
end
function [dsi, di]= computeDSI(numStims, responses, type)
    %numStims = number of stimulations
    %responses is a 1xnumStims array containing the average response
    if nargin < 3
        type = 'preferred';
    end
    
    
    [~, prefIdx]= max(responses);
    nullIdx = mod(numStims/2 + prefIdx -1, numStims)+1; %this wraps our indices from 1-N to find the null directio
    
    switch type
        case 'preferred'
            Rnull = responses(nullIdx);
            Rpref = responses(prefIdx);

        case 'hemifield'
            
            hemifield = -numStims/4 : numStims/4;
            
            nullHemi = mod(nullIdx + hemifield -1, numStims)+1;
            prefHemi = mod(prefIdx + hemifield -1, numStims)+1;
            
                       
            Rpref = sum(responses(prefHemi));
            Rnull = sum(responses(nullHemi));
    end
    
    dsi= (Rpref-Rnull)/(Rpref+Rnull);
    di = (Rpref-Rnull)/Rpref;
    
    %restrict dsi [-1,1]
    dsi(dsi>1)= 1;
    dsi(dsi<-1) = -1;
    
end
function [osi, oi] = computeOSI(numStims, responses)
    %numStims = number of stimulations
    %responses is a 1xnumStims array containing the average response

    [Rpref, preferredDirInd]= max(responses);
    orthIdx=[ (mod(numStims/4 + preferredDirInd -1, numStims) +1), (mod( preferredDirInd -1- numStims/4 , numStims) +1)];

    Rorth= mean(responses(orthIdx));
    osi= (Rpref-Rorth)/(Rpref+Rorth);
    oi = (Rpref-Rorth)/Rpref;
    
    
    %restrict osi [-1,1]
    osi(osi>1)= 1;
    osi(osi<-1) = -1;
    
end
function PlotPrefOnTemplate(analysis, data, metadata,type, field,template,rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,template,template,template)/prctile(template(:),99));
    colormap(hsv) 
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        title('Spatial frequency preference map')
        numsf = (metadata.StimParams.uniqStims-1)/metadata.StimParams.numOrientations;
        coc_sf = cbrewer('qual', 'Set1', numsf);
    elseif type == 3
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        l = rois(i);
        xpos= data.roi(l).xPos;
        ypos= data.roi(l).yPos;
        try
            if type == 1
                sizeMarker = ceil(20 * analysis.(field).roi(l).OSIFit);
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(1+floor(analysis.(field).roi(l).preferredOrientation),:));
            elseif type == 2
                sizeMarker = ceil(20 * analysis.(field).roi(l).SFSI);
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',coc_sf(analysis.(field).roi(l).prefSfStim,:));
            elseif type == 3
                sizeMarker = ceil(20 * analysis.(field).roi(l).DSI);
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(1+floor(analysis.(field).roi(l).preferredDirection),:));
            end
        catch
        end
    end
end
function PlotPrefOnTemplateSF(pref, SI, data, metadata,type, template,rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,template,template,template)/prctile(template(:),99));
    colormap(hsv) 
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        title('Spatial frequency preference map')
        numsf = (metadata.StimParams.uniqStims-1)/metadata.StimParams.numOrientations;
        coc_sf = cbrewer('qual', 'Set1', numsf);
    elseif type == 3
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        xpos= data.roi(rois(i)).xPos;
        ypos= data.roi(rois(i)).yPos;
        sizeMarker = ceil(20 * SI(i));
        try
            if type == 1
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor((pref(i))),:));
            elseif type == 2
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',coc_sf((pref(i))));
            elseif type == 3
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor((pref(i))),:));
            end
        catch
        end
    end
end
function PlotAvgStimResponse(metadata, analysis, field, roi)
    %PlotAvgStimResponse Plots Avg + SEM response for stimuli
    %  	Args:
    %     obj (Population Imaging): class
    %     field (char): analysis field key
    %     roi (int): cell index number
    %   Returns:
    %     NONE
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    colorlevels = metadata.StimParams.numSf;
    coc_sf = cbrewer('qual', 'Set1', colorlevels);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    for i=1:metadata.StimParams.numOrientations
        ax=subplot(1,metadata.StimParams.numOrientations, i);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
            bl=x(stimWindow);
            patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        for sf = 1:metadata.StimParams.numSf
            id = metadata.StimParams.numOrientations*(sf-1)+i;
            err=analysis.(field).roi(roi).SEMResponseTrace(id,:);
            hold all
            patch([x fliplr(x)],[y+err fliplr(y-err)], coc_sf(sf,:), 'LineStyle', 'none');
            %patch([bl fliplr(bl)],[zeros(1,length(bl)), ones(1,length(bl))],[0 0 0]);
            plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', coc_sf(sf,:), 'LineWidth', 3);
            hold all
        end
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        title(metadata.StimParams.direction(i))
        axis off
    end
    hold off
    set(gcf, 'Color', 'w')
end
function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% 
% [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging), 'qual' (qualitative)
%   - cname: name of colortable. It changes depending on ctype.
%   - ncol:  number of color in the table. It changes according to ctype and
%            cname
%   - interp_method: interpolation method (see interp1.m). Default is "cubic" )
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of 
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types and
%           names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
% ------------------------------
% 18.09.2015  Minor fixes, fixed a bug where the 'spectral' color table did not appear in the preview


    % load colorbrewer data
    load('C:\Users\jaepelj\Dropbox\Work\Code\colorbrewer.mat')
    % initialise the colormap is there are any problems
    colormap=[];
    if (~exist('interp_method', 'var'))
        interp_method='cubic';
    end

    % If no arguments
    if (~exist('ctype', 'var') | ~exist('cname', 'var') | ~exist('ncol', 'var'))
        disp(' ')
        disp('[colormap] = cbrewer(ctype, cname, ncol [, interp_method])')
        disp(' ')
        disp('INPUT:')
        disp('  - ctype: type of color table *seq* (sequential), *div* (divergent), *qual* (qualitative)')
        disp('  - cname: name of colortable. It changes depending on ctype.')
        disp('  - ncol:  number of color in the table. It changes according to ctype and cname')
        disp('  - interp_method:  interpolation method  (see interp1.m). Default is "cubic" )')

        disp(' ')
        disp('Sequential tables:')
        z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
                 'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
        disp(z')     

        disp('Divergent tables:')
        z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
        disp(z')

        disp(' ')
        disp('Qualitative tables:')
        %getfield(colorbrewer, 'qual')
        z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
        disp(z')

        plot_brewer_cmap
        return
    end

    % Verify that the input is appropriate
    ctype_names={'div', 'seq', 'qual'};
    if (~ismember(ctype,ctype_names))
        disp('ctype must be either: *div*, *seq* or *qual*')
        colormap=[];
        return
    end

    if (~isfield(colorbrewer.(ctype),cname))
        disp(['The name of the colortable of type *' ctype '* must be one of the following:'])
        getfield(colorbrewer, ctype)
        colormap=[];
        return
    end

    if (ncol>length(colorbrewer.(ctype).(cname)))
    %     disp(' ')
    %     disp('----------------------------------------------------------------------')
    %     disp(['The maximum number of colors for table *' cname '* is ' num2str(length(colorbrewer.(ctype).(cname)))])
    %     disp(['The new colormap will be extrapolated from these ' num2str(length(colorbrewer.(ctype).(cname))) ' values'])
    %     disp('----------------------------------------------------------------------')
    %     disp(' ')
        cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
        colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
        colormap=colormap./255;
        return
    end

    if (isempty(colorbrewer.(ctype).(cname){ncol}))

        while(isempty(colorbrewer.(ctype).(cname){ncol}))
            ncol=ncol+1;
        end        
        disp(' ')
        disp('----------------------------------------------------------------------')
        disp(['The minimum number of colors for table *' cname '* is ' num2str(ncol)])
        disp('This minimum value shall be defined as ncol instead')
        disp('----------------------------------------------------------------------')
        disp(' ')
    end

    colormap=(colorbrewer.(ctype).(cname){ncol})./255;

 end
function PlotTrialStimResponse(metadata, analysis, field, roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    colorlevels = metadata.StimParams.numSf;
    colorlevels_tr = metadata.StimParams.numTrials+1;
    coc_sf = cbrewer('qual', 'Set1', colorlevels);
    coc_tr = [linspace(0,1,colorlevels_tr);linspace(0,1,colorlevels_tr);linspace(0,1,colorlevels_tr)]';
    coc_tr = coc_tr(2:end,:);

    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    if metadata.StimParams.sweep
        for i=1:metadata.StimParams.numOrientations
            for sf = 1:metadata.StimParams.numSf
                id = metadata.StimParams.numOrientations*(sf-1)+i;
                plotnr = metadata.StimParams.numSf*(i-1)+sf;
                ax=subplot(metadata.StimParams.numOrientations, metadata.StimParams.numSf,plotnr);
                cla(ax)
                y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(id,:),3);
                x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
                if stimWindow(1)~=0
                        bl=x(stimWindow);
                        patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
                end
                hold all
                for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                    plot(ax,x, squeeze(analysis.(field).roi(roi).stimResponseTrace(id,trial,:)), 'Color', coc_tr(trial,:));
                    hold all
                end
                plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', coc_sf(sf,:), 'LineWidth', 3);
                hold all
                ylim([ymin, ymax])
                xlim([min(x) max(x)])
                title(metadata.StimParams.directions(i))
                axis off
            end
        end
    else
        for sf = 1:metadata.StimParams.numSf
            for i=1:metadata.StimParams.numOrientations
                id = metadata.StimParams.numOrientations*(sf-1)+i;
                ax=subplot(metadata.StimParams.numSf,metadata.StimParams.numOrientations, id);
                cla(ax)
                y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(id,:),3);
                x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
                if stimWindow(1)~=0
                        bl=x(stimWindow);
                        patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
                end
                hold all
                for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                    plot(ax,x, squeeze(analysis.(field).roi(roi).stimResponseTrace(id,trial,:)), 'Color', coc_tr(trial,:));
                    hold all
                end
                plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', coc_sf(sf,:), 'LineWidth', 3);
                hold all
                ylim([ymin, ymax])
                xlim([min(x) max(x)])
                title(metadata.StimParams.directions(i))
                axis off
            end
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end