function OriTf(analysisParams)
% 
close all
if server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

field = analysisParams.field;
plotROIsResps = analysisParams.plotROIs;

base2pDirectory= [TwoPhontondir analysisParams.animal];
tifDirectory = [base2pDirectory filesep analysisParams.name];
Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep];
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
if analysisParams.reloadData
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    metadata.StimParams.series=expt_id;
    data = LoadRoisS2p(analysisParams);
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');
    metadata.ROI = struct;
    analysis = struct;
    save(fullfile(saveDirectory, 'oriTf.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
else
    load(fullfile(saveDirectory, 'oriTf.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
end

disp('Loading data')

%% chop traces
disp('Chopping Traces')
[analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, 'pre', analysisParams.pre, 'post',metadata.StimParams.isi,'windowStart',analysisParams.windowStart, 'windowStop',analysisParams.windowStop);

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
    analysisPeriod=(analysis.(field).windowStart:analysis.(field).windowStop);
    stimResp = analysis.(field).roi(i).stimResponseTrace(:,:,analysisPeriod);
    analysis.(field).roi(i).peaks = max(stimResp,[],3);
    analysis.(field).roi(i).zscore = ([analysis.(field).roi(i).peaks]-[analysis.(field).roi(i).baselineMean])./[analysis.(field).roi(i).baselineSD];
    analysis.(field).roi(i).crosser = sum(analysis.(field).roi(i).zscore > analysisParams.zThresh,2);
    analysis.(field).roi(i).respStim = analysis.(field).roi(i).crosser > ((metadata.StimParams.numTrials)*analysisParams.fraction);
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
    tempFreq_cell = strsplit(metadata.StimParams.temporalFreq(2:end-1), ',');
    metadata.StimParams.numTf = size(tempFreq_cell,2);
    metadata.StimParams.numSf = 1;
    medResponse = zeros(metadata.StimParams.numOrientations,metadata.StimParams.numTrials,size(analysis.(field).roi(i).stimResponseTrace,3));
    endnum = metadata.StimParams.numTf * metadata.StimParams.numOrientations+1;
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

%% find preferred tF for each cell by looking at preferred Orientation
for i = 1:length(data.roi)
    Response_tf_prefOri = zeros(metadata.StimParams.numTf,metadata.StimParams.numTrials,size(analysis.(field).roi(i).stimResponseTrace,3));
    for tf = 1:metadata.StimParams.numTf
        temp = analysis.(field).roi(i).stimResponseTrace(1+(tf-1)*ori:tf*ori, 1:metadata.StimParams.numTrials, :);
        if sweep
            Response_tf_prefOri(tf, :, :) = temp(analysis.(field).roi(i).prefOriStimInd,:,:);
        else
            Response_tf_prefOri(tf, :, :) = (temp(analysis.(field).roi(i).prefOriStimInd,:,:)+temp(analysis.(field).roi(i).prefOriStimInd+numStims/2,:,:)/2);
        end
    end
    avgResponse_tf_prefOri = mean(Response_tf_prefOri(:,:,stimWindow),3);
    medResponse_tf_prefOri = mean(avgResponse_tf_prefOri,2)';
    PrefTFInd = find(medResponse_tf_prefOri == max(medResponse_tf_prefOri));
    maxResp = max(medResponse_tf_prefOri);
    minResp = min(medResponse_tf_prefOri);
    TFSI = (maxResp - minResp) ./ (maxResp + minResp);
    analysis.(field).roi(i).prefTfStim = PrefTFInd;
    analysis.(field).roi(i).prefTf = tempFreq_cell{PrefTFInd};
    analysis.(field).roi(i).TFSI = TFSI;
end

%% compute DSI and OSI at prefered TF
disp('Calculating ROI properties for each tf')
if sweep
    for i= 1:length(data.roi)
        prefTfResponse = analysis.(field).roi(i).stimResponseTrace((analysis.(field).roi(i).prefTfStim-1)*lastOrientation+1:analysis.(field).roi(i).prefTfStim*lastOrientation, 1:metadata.StimParams.numTrials, :);
        prefTfResponse = mean(prefTfResponse(:,:,stimWindow),3);
        prefTfResponse_dir = median(prefTfResponse,2)';
        [analysis.(field).roi(i).DSI, ~] = computeDSI(lastOrientation, prefTfResponse_dir, 'preferred');
        analysis.(field).roi(i).OSI = analysis.(field).roi(i).DSI;
        analysis.(field).roi(i).OSIFit = analysis.(field).roi(i).DSI;
    end
else
    for i=1:length(data.roi)
        prefTfResponse = analysis.(field).roi(i).stimResponseTrace((analysis.(field).roi(i).prefTfStim-1)*lastOrientation+1:analysis.(field).roi(i).prefTfStim*lastOrientation, 1:metadata.StimParams.numTrials, :);
        prefTfResponse = mean(prefTfResponse(:,:,stimWindow),3);
        prefTfResponse_dir = median(prefTfResponse,2)';
        prefTfResponse(1:size(prefTfResponse,1)/2, 1:size(prefTfResponse,2)) = prefTfResponse(1:size(prefTfResponse,1)/2,:);
        prefTfResponse(1:size(prefTfResponse,1)/2, 1+size(prefTfResponse,2):2*size(prefTfResponse,2))= prefTfResponse(1+size(prefTfResponse,1)/2:size(prefTfResponse,1),:);
        prefTfResponse= median(prefTfResponse,2);
        prefTfResponse_ori = prefTfResponse(:)';
        [analysis.(field).roi(i).OSI, ~] = computeOSI(lastOrientation/2,prefTfResponse_ori);
        [analysis.(field).roi(i).DSI, ~] = computeDSI(lastOrientation, prefTfResponse_dir, 'preferred');
    end
end

%% compute ori preference at each sf

prefDir_allTf = zeros(metadata.StimParams.numTf,1);
prefOri_allTf = zeros(metadata.StimParams.numTf,1);
OSIFit_allTf = zeros(metadata.StimParams.numTf,1);
OSI_allTf = zeros(metadata.StimParams.numTf,1);
DSI_allTf = zeros(metadata.StimParams.numTf,1);
isResp_allTf = zeros(metadata.StimParams.numTf,1);

for i = 1:length(data.roi)
    for tf = 1:metadata.StimParams.numTf
        medResponse = analysis.(field).roi(i).stimResponseTrace((tf-1)*lastOrientation+1:tf*lastOrientation, :,:);
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
               prefDir_allTf(tf) = rad2deg(bestCoeff(4));
           else
               prefDir_allTf(tf) = mod(rad2deg(bestCoeff(4))+90,180);
           end
           DSI_allTf(tf) = computeDSI(lastOrientation, dirResponse, 'preferred');
           prefOri_allTf(tf) = prefDir_allTf(tf);
           OSIFit_allTf(tf) = DSI_allTf(tf);
           OSI_allTf(tf) = DSI_allTf(tf);
        else
            xo=[0:pi/40:2*pi];
            fit=vonMisesLinSum(bestCoeff, xo);
            if bestCoeff(1) > bestCoeff(2)
                prefDir_allTf(tf) = rad2deg(bestCoeff(4));
            else
                prefDir_allTf(tf) = mod(rad2deg(bestCoeff(4))+180,360);
            end
            DSI_allTf(tf) = computeDSI(lastOrientation, dirResponse, 'preferred');
            [prefOri_allTf(tf), coeffOr, ~] = ComputePreferredOrientations(medResponse, theta);
            prefOri_allTf(tf) =mod(prefOri_allTf(tf), 180);
            OSIFit_allTf(tf) = computeOSIFit(coeffOr);
            OSI_allTf(tf) = computeOSI(lastOrientation/2,medResponse(1:lastOrientation/2));
        end
        tempCrosser = analysis.(field).roi(i).crosser((tf-1)*lastOrientation+1:tf*lastOrientation);
        tempRespStim = tempCrosser > ((metadata.StimParams.numTrials)*fraction);
        if sum(tempRespStim) > 0
            isResp_allTf(tf) = 1;
        else
            isResp_allTf(tf) = 0;
        end
    end
    analysis.(field).roi(i).preferredOrientation_allTf = prefOri_allTf;
    analysis.(field).roi(i).preferredDirection_allTf = prefDir_allTf;
    analysis.(field).roi(i).OSIFit_allTf = OSIFit_allTf;
    analysis.(field).roi(i).OSI_allTf = OSI_allTf;
    analysis.(field).roi(i).DSI_allTf = DSI_allTf;
    analysis.(field).roi(i).isResp_allTf = isResp_allTf;
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
        alltfpref = [analysis.dff.roi.prefTfStim];
        rois_dir = rois_ori;
        rois_tf = rois_ori;
    elseif types == 2
        rois_ori = find([analysis.dff.roi.isResponseSignificant] == 1); 
        alloriprefs = [analysis.dff.roi(rois_ori).preferredOrientation];
        alldirprefs = [analysis.dff.roi(rois_ori).preferredDirection];
        alltfpref = [analysis.dff.roi(rois_ori).prefTfStim];
        rois_dir = rois_ori;
        rois_tf = rois_ori;
    elseif types == 3
        rois_ori = find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alloriprefs = [analysis.dff.roi(rois_ori).preferredOrientation];
        rois_dir = find([analysis.dff.roi.DSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alldirprefs = [analysis.dff.roi(rois_dir).preferredDirection];
        rois_tf = find([analysis.dff.roi.TFSI] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
        alltfpref = [analysis.dff.roi(rois_tf).prefTfStim];
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    %plot prefOri on top of template
    subplot(2,3,1)
    PlotPrefOnTemplateOri(analysis, data, metadata,1, field,data.template, rois_ori)
    
    subplot(2,3,4)
    histogram(alloriprefs,linspace(0,180,5), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel(sprintf('Orientation preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)]/2)
    axis square;
    set(gca,'Box','off');
    
    subplot(2,3,2)
    PlotPrefOnTemplateOri(analysis, data, metadata,2, field,data.template, rois_dir)

    subplot(2,3,5)
    histogram(alldirprefs,linspace(0,360,lastOrientation+1), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel(sprintf('Direction preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)])
    axis square;
    set(gca,'Box','off');

    %plot prefTF on top of template
    subplot(2,3,3)
    PlotPrefOnTemplateOri(analysis, data, metadata,4, field,data.template, rois_tf)

    subplot(2,3,6)
    numTf = size(tempFreq_cell,2);
    tf_cat = cellfun(@str2double,tempFreq_cell);
    tf_counts = histcounts(alltfpref, numTf);
    bar(tf_cat, tf_counts);
    title('Histogram')
    ylabel('Cells');
    xlabel(sprintf('Temporal frequency preference'));
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
    rows = ceil(sqrt(metadata.StimParams.numTf));
    for tf = 1:metadata.StimParams.numTf
        if metadata.StimParams.numTf > 6
            subplot(rows,rows,tf)
        else
            subplot(2,rows,tf)
        end
        tempResp = [analysis.dff.roi.isResp_allTf];
        isResponsive = find(tempResp(tf,:) == 1);
        if sum(isResponsive) > 0
            if types == 1
                tempOSI= [analysis.dff.roi.OSIFit_allTf];
                tempOriSelect = find(tempOSI(tf,:) > 0.2);
                rois = intersect(tempOriSelect, isResponsive);
                if ~isempty(rois)
                    tempPrefOri = [analysis.dff.roi(rois).preferredOrientation_allTf];
                    allprefs = tempPrefOri(tf,:);
                    tempFit = [analysis.dff.roi(rois).OSIFit_allTf];
                    allFit = tempFit(tf,:);
                    PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
                    title(['Ori pref map, ' num2str(tempFreq_cell{tf}) ' Hz'])
                else
                    imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                    colormap(hsv)
                    title(['Ori pref map, ' num2str(tempFreq_cell{tf}) ' Hz'])
                    caxis([0 180])
                    colorbar('Location', 'southoutside');
                end
            elseif types == 2
                tempDSI= [analysis.dff.roi.DSI_allTf];
                tempDirSelect = find(tempDSI(tf,:) > 0.2);
                rois = intersect(tempDirSelect, isResponsive);
                if ~isempty(rois)
                    tempPrefDir = [analysis.dff.roi(rois).preferredDirection_allTf];
                    allprefs = tempPrefDir(tf,:);
                    tempFit = [analysis.dff.roi(rois).DSI_allTf];
                    allFit = tempFit(tf,:);
                    PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
                    title(['Dir pref map, ' num2str(tempFreq_cell{tf}) ' Hz'])
                else
                    imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                    colormap(hsv)
                    title(['Dir pref map, ' num2str(tempFreq_cell{tf}) ' Hz'])
                    caxis([0 360])
                    colorbar('Location', 'southoutside');
                end
            end
        else
            imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
            colormap(hsv)
            if types == 1
                title(['Ori pref map, ' num2str(tempFreq_cell{tf}) ' Hz'])
                caxis([0 180])
            elseif types == 2
                title(['Dir pref map, ' num2str(tempFreq_cell{tf}) ' Hz'])
                caxis([0 180])
            end
             colorbar('Location', 'southoutside');
        end
    end
    if types == 1
        saveas(gcf, fullfile(saveDirectory, 'OrientationPreference_allTF.png'))
    elseif types == 2
        saveas(gcf, fullfile(saveDirectory, 'DirectionPreference_allTF.png'))
    end
end
%plot all ROIs
if plotROIs
    for i = 1:length(data.roi)
        PlotAvgStimResponseOri(metadata, analysis, field, i)
        saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
        close gcf
    end
    for i = 1:length(data.roi)
        PlotTrialStimResponseOri(metadata, analysis, field, i)
        saveas(gcf, fullfile(ROIsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
        close gcf
    end
end
if plotROIsResps
    for i = 1:length(data.roi)
        if analysis.(field).roi(i).isResponseSignificant == 1
            PlotTrialStimResponseOri(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        else
            PlotTrialStimResponseOri(metadata, analysis, field, i)
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

save(fullfile(saveDirectory, 'oriTfAna.mat'), 'data', 'metadata', 'sliceparams', 'analysis');

end
