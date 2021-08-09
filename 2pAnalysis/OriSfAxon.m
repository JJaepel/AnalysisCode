function OriSfAxon(analysisParams)

close all
if analysisParams.server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
Sp2dir = [drive '\Data\Spike2Data\'];
savedir = [drive '\Data\ImageAnalysis\'];

field = 'dff';
plotROIsResps = analysisParams.plotROIs;

base2pDirectory= [TwoPhontondir analysisParams.animal];
tifDirectory = [base2pDirectory filesep analysisParams.name];
Sp2dDirectory = [Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [savedir analysisParams.animal filesep analysisParams.expID filesep];
%ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
% if ~exist(ROIsaveDirectory, 'dir')
%     mkdir(ROIsaveDirectory);  

disp('Loading data')

%% load Data and metadata
if analysisParams.reloadData
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    if analysisParams.makeROIs
        if analysisParams.level == 1
            Suite2pAxonTifCombiner(analysisParams)
        end
        data = Suite2pAxonExtractorFct(analysisParams);
    else
        if analysisParams.level == 1
            
            dataFile = ([base2pDirectory filesep analysisParams.expID filesep 'suite2p\combined\data.mat']);
        else
            dataFile = ([base2pDirectory filesep analysisParams.expID filesep 'suite2p\plane0\data.mat']);
        end
        File = load(dataFile, 'data');
        data = File.data;
        clear File
    end
    [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
    data = computeDff(data, 'rawF', 'baseline', 'dff');
    if analysisParams.checkDFF
        figure
        subplot(3,1,1)
        plot(data.roi(1).rawF)
        axis off
        subplot(3,1,2)
        plot(data.roi(1).baseline)
        axis off
        subplot(3,1,3)
        plot(data.roi(1).dff)
        axis off
    end
    metadata.ROI = struct;
    analysis = struct;
    save(fullfile(saveDirectory, 'oriSf.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
else
    load(fullfile(saveDirectory, 'oriSf.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
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
metadata.StimParams.orientations = orientations;
directions = linspace(0,360,numStims+1);
metadata.StimParams.directions = directions(1:end-1);
%find out if it is across several directions or if it is a sweep with only
%few orientations
try
    if max(metadata.StimParams.directions > 180)
        sweep = 0;
    else
        sweep = 1;
    end
catch
    sweep = 0;
end
metadata.StimParams.sweep = sweep;
for i = 1:length(data.roi)
    %collapse all sf for the same orientation 
    spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
    metadata.StimParams.numSf = size(spatFreq_cell,2);
    metadata.StimParams.numTf = 1;
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
        tempRespStim = tempCrosser > ((metadata.StimParams.numTrials)*analysisParams.fraction);
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
    
    %plot prefSF on top of template
    subplot(2,3,3)
    PlotPrefOnTemplateOri(analysis, data, metadata,3, field,data.template, rois_sf)

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
                    PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
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
                    PlotPrefOnTemplateCon(allprefs, allFit,data,types, data.template,rois)
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

%% plot ROIs
if analysisParams.plotROIs
    if ~exist(ROIRespsaveDirectory, 'dir')
        % make new file directory
        mkdir(ROIRespsaveDirectory); 
    else
        filePattern = fullfile(ROIRespsaveDirectory, '*.png'); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        for k = 1 : length(theFiles)
          baseFileName = theFiles(k).name;
          fullFileName = fullfile(ROIRespsaveDirectory, baseFileName);
          delete(fullFileName);
        end
        %remove old files
    end
    if ~analysisParams.plotRespROIsOnly
        if ~exist(ROINonRespsaveDirectory, 'dir')
            mkdir(ROINonRespsaveDirectory);
        else
            filePattern = fullfile(ROIRespsaveDirectory, '*.png'); % Change to whatever pattern you need.
            theFiles = dir(filePattern);
            for k = 1 : length(theFiles)
              baseFileName = theFiles(k).name;
              fullFileName = fullfile(ROIRespsaveDirectory, baseFileName);
              delete(fullFileName);
            end
        end
    end

    for i = 1:length(data.roi)
        if analysis.(field).roi(i).isResponseSignificant == 1
            PlotAvgStimResponseOri(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
            close gcf
            PlotTrialStimResponseOri(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        else
            if ~analysisParams.plotRespROIsOnly
                PlotAvgStimResponseOri(metadata, analysis, field, i)
                saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
                close gcf
                PlotTrialStimResponseOri(metadata, analysis, field, i)
                saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                close gcf
            end
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
direction = length(find([analysis.dff.roi.DSI] > 0.33 & [analysis.dff.roi.isResponseSignificant] == 1)) ./all;
non_dir = 1-direction-non_resp;
h = pie([non_resp non_dir direction]);
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

%% save analysis file
save(fullfile(saveDirectory, 'oriSfAna.mat'), 'data', 'metadata', 'analysisParams', 'analysis');

end