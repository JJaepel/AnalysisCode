function GratingAnalysis(analysisParams)

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
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];

if ~exist(saveDirectory, 'dir')
    % make new file directory
    mkdir(saveDirectory); 
end

metadata = struct;
metadata.ROI = struct;
analysis = struct;

%% 1.) load data and metadata
disp('Loading data')
if analysisParams.reloadData
    % make sure that you also do analysis afterwards
    analysisParams.reanalyse = 1;
    
    % load stimulus parameters and spike2data
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParams(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    
    % load tif data
    switch analysisParams.dataType
        case 1
            data = LoadRoisS2p(analysisParams);
            
            %do baseline filtering and compute dff
            if strcmp(analysisParams.field, 'dff')
                [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                data = computeDff(data, 'rawF', 'baseline', 'dff');
            end
            
        case 2
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
            %do baseline filtering and compute dff
            if strcmp(analysisParams.field, 'dff')
                %etadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                %data = computeDff(data, 'rawF', 'baseline', 'dff');[m
                data = computeDffAxons(metadata, data);
            end
            
        case 3
            %do spine data reloading here
            data = SpineROIExtractFct(analysisParams);
            data = computeDffSpines(data);
            structureType = analysisParams.region;
            switch structureType
                case 'cells'
                    data = FitCellsToDendrite(analysisParams,metadata, data);
                case 'dendrite'
                    data = computeDendriticSubstraction(analysisParams, metadata, data);    
            end
            analysisParams.field = 'rawRes';
    end
    
    
    %save file
    switch analysisParams.stimType
        case 1
            save(fullfile(saveDirectory, 'oriGrating.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
        case 2
            save(fullfile(saveDirectory, 'oriSf.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
        case 5
            save(fullfile(saveDirectory, 'oriBino.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
    end

else
    switch analysisParams.stimType
        case 1
            load(fullfile(saveDirectory, 'oriGrating.mat'), 'data', 'metadata', 'analysis');
        case 2
            load(fullfile(saveDirectory, 'oriSf.mat'), 'data', 'metadata', 'analysis');
        case 5
            save(fullfile(saveDirectory, 'oriBino.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
    end
end

%% 2.) Do analysis
if analysisParams.reanalyse
    % chop traces
    disp('Chopping Traces')
    [analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,analysisParams.level, analysisParams.field, 'pre', analysisParams.pre, 'post',metadata.StimParams.isi,'windowStart',analysisParams.windowStart, 'windowStop',analysisParams.windowStop);
    
    %clean up properties
    numStims=metadata.StimParams.uniqStims-1; %we have -1 number of stims (=Blank)
    try
        spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
        metadata.StimParams.numSf = size(spatFreq_cell,2);
    catch
        metadata.StimParams.numSf = length(metadata.StimParams.spatialFreq);
    end
    try
        tempFreq_cell = strsplit(metadata.StimParams.tempFreq(2:end-1), ',');
        metadata.StimParams.numTf = size(tempFreq_cell,2);
    catch
        metadata.StimParams.numTf = length(metadata.StimParams.temporalFreq);
    end
    if analysisParams.stimType == 5
        metadata.StimParams.numDirections = metadata.StimParams.numOrientations;
    else
        metadata.StimParams.numDirections = numStims/(metadata.StimParams.numSf*metadata.StimParams.numTf);
    end
    metadata.StimParams.orientations = linspace(0,180,metadata.StimParams.numDirections/2+1); %this would give us a all orientations from 0 to 180
    metadata.StimParams.orientations = metadata.StimParams.orientations(1:end-1); %180 == 0, so we need to substract htat one
    metadata.StimParams.directions = linspace(0,360,metadata.StimParams.numDirections+1); %this would give us a all orientations from 0 to 180
    metadata.StimParams.directions = metadata.StimParams.directions(1:end-1); %180 == 0, so we need to substract htat one
    metadata.StimParams.numOrientations = length(metadata.StimParams.orientations); %how many oris do we have?
    
    theta= [0:2*pi/(metadata.StimParams.numDirections) : 2*pi-2*pi/(metadata.StimParams.numDirections)];%theta is needed for the fitting
    theta=mod(theta, pi); 
    thetaOri = theta(1:metadata.StimParams.numOrientations);
    stimWindow=(analysis.(analysisParams.field).windowStart: analysis.(analysisParams.field).windowStop);
    
    %calculate indices
    disp('Calculating indices')
    for i = 1:length(data.roi) %for each ROI
        %for multiple sf or tfs, collapse them 
        if metadata.StimParams.numTf * metadata.StimParams.numSf > 1
            dirResponse = zeros(metadata.StimParams.numOrientations,metadata.StimParams.numTrials,size(analysis.(analysisParams.field).roi(i).stimResponseTrace,3));
            endnum = metadata.StimParams.numSf * metadata.StimParams.numTf * metadata.StimParams.numDirections+1;
            for ori = 1:metadata.StimParams.numDirections
                dirResponse(ori,:,:) = mean(analysis.(analysisParams.field).roi(i).stimResponseTrace(ori:metadata.StimParams.numOrientations:endnum, 1:metadata.StimParams.numTrials, :));
            end
        %alternatively, collapse for ipsi, contra and binocular response
        elseif analysisParams.stimType == 5
            dirResponse = zeros(metadata.StimParams.numOrientations,metadata.StimParams.numTrials,size(analysis.(analysisParams.field).roi(i).stimResponseTrace,3));
            endnum = 3 * metadata.StimParams.numDirections+1;
            for ori = 1:metadata.StimParams.numDirections
                dirResponse(ori,:,:) = mean(analysis.(analysisParams.field).roi(i).stimResponseTrace(ori:metadata.StimParams.numOrientations:endnum, 1:metadata.StimParams.numTrials, :));
            end
        else
            dirResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace(1:end-1, :, :);
        end
        
        %calculate the median response curve for each direction and orientation
        dirResponse = mean(dirResponse(:,:,stimWindow),3); %average over the stimWindow
        Responses = dirResponse; %this contains the variability over the trials and is needed for Cohen's D

        medResponse(1:size(dirResponse,1)/2, 1:size(dirResponse,2)) = dirResponse(1:size(dirResponse,1)/2,:); %rearrange so that the two corresponding directions can be summed
        medResponse(1:size(dirResponse,1)/2, 1+size(dirResponse,2):2*size(dirResponse,2))= dirResponse(1+size(dirResponse,1)/2:size(dirResponse,1),:);
        allOriResponses = medResponse;
        medResponse= median(medResponse,2);
        oriResponse = medResponse(:)';  %this is the median response for each orientation
        dirResponse = median(dirResponse,2)'; %this is the median response for each direction

        %OSI/DSI: (Pref Ori - Orth Ori)/ (Pref Ori + Orth Ori)
        [analysis.(analysisParams.field).roi(i).OSI] = computeOSI(metadata.StimParams.numOrientations ,oriResponse);
        [analysis.(analysisParams.field).roi(i).DSI] = computeDSI(metadata.StimParams.numOrientations *2, dirResponse, 'preferred');

        %CircVar, DirCircVar
        [analysis.(analysisParams.field).roi(i).OriCircVar] = 1-OriCircularVariance(dirResponse);
        [analysis.(analysisParams.field).roi(i).DirCircVar] = 1-DirCircularVariance(dirResponse);

        %Cohen's D
        [analysis.(analysisParams.field).roi(i).cohensD] = computeCohensD(metadata.StimParams.numOrientations ,Responses);

        %Find preferred orientation based on gaussian fit
        [analysis.(analysisParams.field).roi(i).preferredOrientation,analysis.(analysisParams.field).roi(i).coeffOr,analysis.(analysisParams.field).roi(i).rsqOr] = ComputePreferredOrientations(oriResponse, thetaOri);
        analysis.(analysisParams.field).roi(i).OSIFit = computeOSIFit(analysis.(analysisParams.field).roi(i).coeffOr); %compute OSI from fit

        %Find preferred direction based on double gaussian fit
        try
            [analysis.(analysisParams.field).roi(i).preferredDirection] = ComputePreferredDirection(dirResponse, theta);
        catch
            [analysis.(analysisParams.field).roi(i).preferredDirection] = NaN;
        end
        
        interval = metadata.StimParams.directions(2) - metadata.StimParams.directions(1);
        ind = round(analysis.(analysisParams.field).roi(i).preferredOrientation / interval);
        ind = ind + 1;
        if ind > metadata.StimParams.numOrientations || isnan(ind)
            ind = 1;
        end
        analysis.(analysisParams.field).roi(i).prefOriStimInd = ind;
        
        %check if significantly ori/dir-selective
        vectorRespOri = (allOriResponses'*transpose(exp(sqrt(-1)*2*mod(metadata.StimParams.orientations*pi/180,pi))));
        try
            [analysis.(analysisParams.field).roi(i).OriSignH, analysis.(analysisParams.field).roi(i).pOriSign] =  hotellingt2test(vectorRespOri, [0 0]);
            [analysis.(analysisParams.field).roi(i).DirSignH, analysis.(analysisParams.field).roi(i).pDirSign] =  computeDirSignificanceDotProduct(metadata.StimParams.directions, Responses');
        catch
            analysis.(analysisParams.field).roi(i).OriSignH = NaN; analysis.(analysisParams.field).roi(i).pOriSign = NaN;
            analysis.(analysisParams.field).roi(i).DirSignH = NaN; analysis.(analysisParams.field).roi(i).pDirSign = NaN;
        end
    end
    
    if metadata.StimParams.numTf * metadata.StimParams.numSf > 1 || analysisParams.stimType == 5
        for i= 1:length(data.roi)
            % find preferred condition for each cell by looking at preferred Orientation
            if analysisParams.stimType == 5
                metadata.StimParams.numCon = 3;
            else
                metadata.StimParams.numCon = max([metadata.StimParams.numSf metadata.StimParams.numTf]);
            end
            ResponseConPrefOri = zeros(metadata.StimParams.numCon,metadata.StimParams.numTrials,size(analysis.(analysisParams.field).roi(i).stimResponseTrace,3));
            for con = 1:metadata.StimParams.numCon
                temp = analysis.(analysisParams.field).roi(i).stimResponseTrace(1+(con-1)*ori:con*ori, 1:metadata.StimParams.numTrials, :);
                ResponseConPrefOri(con, :, :) = (temp(analysis.(analysisParams.field).roi(i).prefOriStimInd,:,:)+temp(analysis.(analysisParams.field).roi(i).prefOriStimInd+metadata.StimParams.numOrientations,:,:)/2);
            end
            avgResponseConPrefOri = mean(ResponseConPrefOri(:,:,stimWindow),3);
            medResponseConprefOri = mean(avgResponseConPrefOri,2)';
            PrefConInd = find(medResponseConprefOri == max(medResponseConprefOri));
            maxResp = max(medResponseConprefOri);
            minResp = min(medResponseConprefOri);
            ConSI = (maxResp - minResp) ./ (maxResp + minResp);
            switch analysisParams.stimType
                case 2
                    analysis.(analysisParams.field).roi(i).prefSfStim = PrefConInd;
                    analysis.(analysisParams.field).roi(i).prefSf = spatFreq_cell{PrefConInd};
                    analysis.(analysisParams.field).roi(i).SFSI = ConSI;
                case 4
                    analysis.(analysisParams.field).roi(i).prefTfStim = PrefConInd;
                    analysis.(analysisParams.field).roi(i).prefTf = tempFreq_cell{PrefConInd};
                    analysis.(analysisParams.field).roi(i).TFSI = ConSI;
                case 5
                    analysis.(analysisParams.field).roi(i).ODI = ConSI;
                    analysis.(analysisParams.field).roi(i).prefCon = PrefConInd;
            end
        
            % compute DSI and OSI at preferred condition
            prefConResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace((PrefConInd-1)* metadata.StimParams.numDirections +1:PrefConInd*metadata.StimParams.numDirections, 1:metadata.StimParams.numTrials, :);
            prefConResponse = mean(prefConResponse(:,:,stimWindow),3);
            prefConResponse_dir = median(prefConResponse,2)';
            prefConResponse(1:size(prefConResponse,1)/2, 1:size(prefConResponse,2)) = prefConResponse(1:size(prefConResponse,1)/2,:);
            prefConResponse(1:size(prefConResponse,1)/2, 1+size(prefConResponse,2):2*size(prefConResponse,2))= prefConResponse(1+size(prefConResponse,1)/2:size(prefConResponse,1),:);
            prefConResponse= median(prefConResponse,2);
            prefConResponse_ori = prefConResponse(:)';
            
            [analysis.(analysisParams.field).roi(i).OSI] = computeOSI(metadata.StimParams.numOrientations,prefConResponse_ori(1:metadata.StimParams.numOrientations));
            [analysis.(analysisParams.field).roi(i).DSI] = computeDSI(metadata.StimParams.numDirections, prefConResponse_dir, 'preferred');
            
            %CircVar, DirCircVar
            [analysis.(analysisParams.field).roi(i).OriCircVar] = 1-OriCircularVariance(prefConResponse);
            [analysis.(analysisParams.field).roi(i).DirCircVar] = 1-DirCircularVariance(prefConResponse);
            
            %Cohen's D
            [analysis.(analysisParams.field).roi(i).cohensD] = computeCohensD(metadata.StimParams.numOrientations,Responses);
        
            % compute ori preference at each condition
            for con = 1:metadata.StimParams.numCon
                dirResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace((con-1)*metadata.StimParams.numDirections+1:con*metadata.StimParams.numDirections, :,:);
                dirResponse = mean(dirResponse(:,:,stimWindow),3); %average over the stimWindow
                medResponse(1:size(dirResponse,1)/2, 1:size(dirResponse,2)) = dirResponse(1:size(dirResponse,1)/2,:); %rearrange so that the two corresponding directions can be summed
                medResponse(1:size(dirResponse,1)/2, 1+size(dirResponse,2):2*size(dirResponse,2))= dirResponse(1+size(dirResponse,1)/2:size(dirResponse,1),:);
                medResponse= median(medResponse,2);
                oriResponse = medResponse(:)';  %this is the median response for each orientation
                dirResponse = median(dirResponse,2)';
                
                %Find preferred orientation based on gaussian fit
                [prefOri_allCon(con),~,~] = ComputePreferredOrientations(oriResponse, thetaOri);
                OSIFit_allCon(con) = computeOSIFit(analysis.(analysisParams.field).roi(i).coeffOr); %compute OSI from fit
                
                %Find preferred direction based on double gaussian fit
                [prefDir_allCon(con)] = ComputePreferredDirection(dirResponse, theta);
                [DSI_allCon(con)] = computeDSI(metadata.StimParams.numDirections, dirResponse, 'preferred');
            end
            switch analysisParams.stimType
                case 2
                    analysis.(analysisParams.field).roi(i).preferredOrientation_allSf = prefOri_allCon;
                    analysis.(analysisParams.field).roi(i).OSIFit_allSf = OSIFit_allCon;
                    analysis.(analysisParams.field).roi(i).preferredDirection_allSf = prefDir_allCon;
                    analysis.(analysisParams.field).roi(i).DSI_allSf = DSI_allCon;
                case 3
                    analysis.(analysisParams.field).roi(i).preferredOrientation_allTf = prefOri_allCon;
                    analysis.(analysisParams.field).roi(i).OSIFit_allTf = OSIFit_allCon;
                    analysis.(analysisParams.field).roi(i).preferredDirection_allTf = prefDir_allCon;
                    analysis.(analysisParams.field).roi(i).DSI_allTf = DSI_allCon;
                case 5
                    analysis.(analysisParams.field).roi(i).preferredOrientation_allEyes = prefOri_allCon;
                    analysis.(analysisParams.field).roi(i).OSIFit_allEyes = OSIFit_allCon;
                    analysis.(analysisParams.field).roi(i).preferredDirection_allEyes = prefDir_allCon;
                    analysis.(analysisParams.field).roi(i).DSI_allEyes = DSI_allCon;
            end
            
        end
    else
        metadata.StimParams.numCon = 1;
    end
    
    
    % find maxResponses and significant responses
    disp('Calculating significant responses')
    disp(['Criteria: threshold (z  = ' num2str(analysisParams.zThresh) ') needs to be crossed ' num2str(analysisParams.fraction *100) ' % of all trials'])
    
    for i = 1:length(data.roi)
        pretrialTime= analysis.(analysisParams.field).preTrialTime;
        preTrialIndex= (1:floor(pretrialTime * metadata.TwoPhoton.rate));
        baselines=analysis.(analysisParams.field).roi(i).stimResponseTrace(:,:,preTrialIndex);
        analysis.(analysisParams.field).roi(i).baselineSD = std(baselines,[],3);
        analysis.(analysisParams.field).roi(i).baselineMean = mean(baselines,3);
        analysis.(analysisParams.field).roi(i).baselineMean(analysis.(analysisParams.field).roi(i).baselineMean < 0) = mean(mean(analysis.(analysisParams.field).roi(i).baselineMean,2),1);
        analysisPeriod=(analysis.(analysisParams.field).windowStart:analysis.(analysisParams.field).windowStop);
        stimResp = analysis.(analysisParams.field).roi(i).stimResponseTrace(:,:,analysisPeriod);
        analysis.(analysisParams.field).roi(i).peaks = max(stimResp,[],3);
        analysis.(analysisParams.field).roi(i).zscore = ([analysis.(analysisParams.field).roi(i).peaks]-[analysis.(analysisParams.field).roi(i).baselineMean])./[analysis.(analysisParams.field).roi(i).baselineSD];
        analysis.(analysisParams.field).roi(i).crosser = sum(analysis.(analysisParams.field).roi(i).zscore > analysisParams.zThresh,2);
        analysis.(analysisParams.field).roi(i).respStim = analysis.(analysisParams.field).roi(i).crosser >= ((metadata.StimParams.numTrials)*analysisParams.fraction);
        if sum(analysis.(analysisParams.field).roi(i).respStim) > 0
            analysis.(analysisParams.field).roi(i).isResponseSignificant = 1;
        else 
            analysis.(analysisParams.field).roi(i).isResponseSignificant = 0;
        end
        
        if analysisParams.stimType == 2 || analysisParams.stimType == 3 || analysisParams.stimType == 5
            for con = 1:metadata.StimParams.numCon
                tempCrosser = analysis.(analysisParams.field).roi(i).crosser((con-1)*metadata.StimParams.numOrientations+1:con*metadata.StimParams.numOrientations);
                tempRespStim = tempCrosser > ((metadata.StimParams.numTrials)*analysisParams.fraction);
                if sum(tempRespStim) > 0
                    isResp_allCon(con) = 1;
                else
                    isResp_allCon(con) = 0;
                end
            end
            analysis.(analysisParams.field).roi(i).isResp_allCon = isResp_allCon;

        end
    end
    
    % calculate deltaOri vs. distance of ROIs
    disp('Calculating distance between ROIs, deltaOri and HI')
    
    [zoom, setup] = getzoom(tifDirectory); %first get um per pixel for the experiment
    metadata.zoom = zoom;
    if setup == 1
        fieldofview = 1000/zoom;
        umperpixel = fieldofview/512;
        disp('Setup Ben')
    elseif setup == 2
        umperpixel = 2.73/zoom;
    end
    
    ori_sel = find([analysis.(analysisParams.field).roi.OSIFit] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    [distOri, deltaOri, distROIs] = calcDistDeltaOri(data, analysis, ori_sel,umperpixel, analysisParams.level, analysisParams.field);
    for i = 1:length(deltaOri)
        analysis.(analysisParams.field).ori_pair(i).distance = distOri(i);
        analysis.(analysisParams.field).ori_pair(i).deltaOri = deltaOri(i);
    end
    analysis.(analysisParams.field).distROIs = distROIs;
    
    % calculate HeterogeneityIndex
    analysis.(analysisParams.field).ori_cells.HomeogeneityIndex = calculateHI(analysis, ori_sel, analysisParams.level, data, analysisParams.field);
    
    % try to predict stimulus from the data
    if analysisParams.predictor
        disp('Predicting stimulus from cell responses')
        [analysis.(analysisParams.field).decoderData] = decoderOri(analysis, data, analysisParams.field);
    end
    
    % save file
    disp('Saving analyzed data')
    save(fullfile(saveDirectory, 'AnaData.mat'), 'data', 'metadata', 'analysisParams', 'analysis');

else
    disp('Loading analyzed data')
    load(fullfile(saveDirectory, 'AnaData.mat'), 'data', 'metadata', 'analysis');
end
metadata.StimParams.sweep = 0;
%% 3.) Plot data
coc_prop = cbrewer('qual', 'Paired', 12);

% 3.1 Plot ROI positions
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

% 3.2 Plot prefOri on top of template
for types = 1:3
    if types == 1
        rois_ori = linspace(1,length(analysis.(analysisParams.field).roi),length(analysis.(analysisParams.field).roi));
        rois_dir = rois_ori;
        rois_con = rois_ori;
    elseif types == 2
        rois_ori = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        rois_dir = rois_ori;
        rois_con = rois_ori;
    elseif types == 3
        rois_ori = find([analysis.(analysisParams.field).roi.OSIFit] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        rois_dir = find([analysis.(analysisParams.field).roi.DSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        if analysisParams.stimType == 2
            rois_con = find([analysis.(analysisParams.field).roi.SFSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        elseif analysisParams.stimType == 3
            rois_con = find([analysis.(analysisParams.field).roi.TFSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        end
    end
    
    alloriprefs = [analysis.(analysisParams.field).roi(rois_ori).preferredOrientation];
    alldirprefs = [analysis.(analysisParams.field).roi(rois_ori).preferredDirection];
    if analysisParams.stimType == 2
        allConpref = [analysis.(analysisParams.field).roi(rois_con).prefSfStim];
    elseif analysisParams.stimType == 3
        allConpref = [analysis.(analysisParams.field).roi(rois_con).prefTfStim];
    elseif analysisParams.stimType == 5
        allConpref = [analysis.(analysisParams.field).roi(rois_con).prefCon];
    end

    h=figure('units','normalized','outerposition',[0 0 1 1]);
    if analysisParams.stimType == 1
        subplot(2,2,1)
    else 
        subplot(2,3,1)
    end
    PlotPrefOnTemplateOri(analysis, data, 1, analysisParams.field,data.template, rois_ori)
    
    if analysisParams.stimType == 1
        subplot(2,2,2)
    else 
        subplot(2,3,4)
    end
    histogram(alloriprefs,linspace(0,180,metadata.StimParams.numOrientations +1), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel(sprintf('Orientation preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)]/2)
    axis square;
    set(gca,'Box','off');
    
    if analysisParams.stimType == 1
        subplot(2,2,3)
    else 
        subplot(2,3,2)
    end
    PlotPrefOnTemplateOri(analysis, data, 2, analysisParams.field,data.template, rois_dir)

    if analysisParams.stimType == 1
        subplot(2,2,4)
    else 
        subplot(2,3,5)
    end
    histogram(alldirprefs,linspace(0,360,2*metadata.StimParams.numOrientations +1), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel(sprintf('Direction preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)])
    axis square;
    set(gca,'Box','off');
    
    if analysisParams.stimType == 2
        subplot(2,3,3)
        PlotPrefOnTemplateOri(analysis, data, 3, analysisParams.field,data.template, rois_con)
        subplot(2,3,6)
        spatFreq_cell = strsplit(metadata.StimParams.spatialFreq(2:end-1), ',');
        numSf = size(spatFreq_cell,2);
        sf_cat = cellfun(@str2double,spatFreq_cell);
        sf_counts = histcounts(allConpref, numSf);
        bar(sf_cat, sf_counts, 'FaceColor', coc_prop(5,:), 'EdgeColor', coc_prop(6,:));
        title('Histogram')
        ylabel('Cells');
        xlabel(sprintf('Spatial frequency preference'));
        axis square;
        set(gca,'Box','off');
    elseif analysisParams.stimType == 3
        subplot(2,3,3)
        PlotPrefOnTemplateOri(analysis, data, 4, analysisParams.field,data.template, rois_con)
        subplot(2,3,6)
        tempFreq_cell = strsplit(metadata.StimParams.temporalFreq(2:end-1), ',');
        numTf = size(tempFreq_cell,2);
        tf_cat = cellfun(@str2double,tempFreq_cell);
        tf_counts = histcounts(allConpref, numTf);
        bar(tf_cat, tf_counts, 'FaceColor', coc_prop(5,:), 'EdgeColor', coc_prop(6,:));
        title('Histogram')
        ylabel('Cells');
        xlabel(sprintf('Temporal frequency preference'));
        axis square;
        set(gca,'Box','off');
    end    
    
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

if analysisParams.stimType == 2 || analysisParams.stimType == 3 || analysisParams.stimType == 5
    %plot preferred ori/dir for all conditions
    for types = 1:2
        h=figure('units','normalized','outerposition',[0 0 1 1]);
        rows = ceil(sqrt(metadata.StimParams.numCon));
        for con = 1:metadata.StimParams.numCon
            if metadata.StimParams.numCon > 6
                subplot(rows,rows,con)
            else
                subplot(2,rows,con)
            end
            tempResp = [analysis.(analysisParams.field).roi.isResp_allCon];
            tempResp = tempResp(con:metadata.StimParams.numCon:end);
            isResponsive = find(tempResp == 1);
            if sum(isResponsive) > 0
                if types == 1
                    if analysisParams.stimType == 2  
                        tempOSI= [analysis.(analysisParams.field).roi.OSIFit_allSf];
                    elseif analysisParams.stimType == 3
                        tempOSI= [analysis.(analysisParams.field).roi.OSIFit_allTf];
                    elseif analysisParams.stimType == 5
                            tempOSI = [analysis.(analysisParams.field).roi.OSIFit_allEyes];
                    end
                    tempOSI = tempOSI(con:metadata.StimParams.numCon:end);
                    tempOriSelect = find(tempOSI > 0.2);
                    rois = intersect(tempOriSelect, isResponsive);
                    if ~isempty(rois)
                        if analysisParams.stimType == 2
                            tempPrefOri = [analysis.(analysisParams.field).roi(rois).preferredOrientation_allSf];
                            tempFit = [analysis.(analysisParams.field).roi(rois).OSIFit_allSf];
                        elseif analysisParams.stimType == 3
                                tempPrefOri = [analysis.(analysisParams.field).roi(rois).preferredOrientation_allTf];
                                tempFit = [analysis.(analysisParams.field).roi(rois).OSIFit_allTf];
                        elseif analysisParams.stimType == 5
                                tempPrefOri = [analysis.(analysisParams.field).roi(rois).preferredOrientation_allEyes];
                                tempFit = [analysis.(analysisParams.field).roi(rois).OSIFit_allEyes];
                        end
                        allprefs = tempPrefOri(con:metadata.StimParams.numCon:end);
                        allFit = tempFit(con:metadata.StimParams.numCon:end);
                        PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
                        if analysisParams.stimType == 2
                            title(['Ori pref map, ' num2str(spatFreq_cell{con}) ' cpd'])
                        elseif analysisParams.stimType == 3
                            title(['Ori pref map, ' num2str(tempFreq_cell{con}) ' Hz'])
                        elseif analysisParams.stimType == 5
                            if con == 1
                                title('Ori pref map, contra')
                            elseif con == 2
                                title('Ori pref map, ipsi')
                            elseif con == 3
                                title('Ori pref map, binocular')
                            end
                        end
                    else
                        imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                        colormap(hsv)
                        if analysisParams.stimType == 2
                            title(['Ori pref map, ' num2str(spatFreq_cell{con}) ' cpd'])
                        elseif analysisParams.stimType == 3
                            title(['Ori pref map, ' num2str(tempFreq_cell{con}) ' Hz'])
                        elseif analysisParams.stimType == 5
                            if con == 1
                                title('Ori pref map, contra')
                            elseif con == 2
                                title('Ori pref map, ipsi')
                            elseif con == 3
                                title('Ori pref map, binocular')
                            end
                        end
                        caxis([0 180])
                        colorbar('Location', 'southoutside');
                    end
                elseif types == 2
                    if analysisParams.stimType == 2 
                        tempDSI= [analysis.(analysisParams.field).roi.DSI_allSf];
                    elseif analysisParams.stimType == 3
                        tempDSI= [analysis.(analysisParams.field).roi.DSI_allTf];
                    elseif analysisParams.stimType == 5
                        tempDSI= [analysis.(analysisParams.field).roi.DSI_allEyes];
                    end
                    tempDSI = tempDSI(con:metadata.StimParams.numCon:end);
                    tempDirSelect = find(tempDSI > 0.2);
                    rois = intersect(tempDirSelect, isResponsive);
                    if ~isempty(rois)
                        if analysisParams.stimType == 2
                            tempPrefDir = [analysis.(analysisParams.field).roi(rois).preferredDirection_allSf];
                            tempFit = [analysis.(analysisParams.field).roi(rois).DSI_allSf];
                        elseif analysisParams.stimType == 3
                            tempPrefDir = [analysis.(analysisParams.field).roi(rois).preferredDirection_allTf];
                            tempFit = [analysis.(analysisParams.field).roi(rois).DSI_allTf];
                        elseif analysisParams.stimType == 5
                            tempPrefDir = [analysis.(analysisParams.field).roi(rois).preferredDirection_allEyes];
                            tempFit = [analysis.(analysisParams.field).roi(rois).DSI_allEyes];
                        end
                        allprefs = tempPrefDir(con:metadata.StimParams.numCon:end);
                        allFit = tempFit(con:metadata.StimParams.numCon:end);
                        PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
                        if analysisParams.stimType == 2
                            title(['Dir pref map, ' num2str(spatFreq_cell{con}) ' cpd'])
                        elseif analysisParams.stimType == 3
                            title(['Dir pref map, ' num2str(tempFreq_cell{con}) ' Hz'])
                        elseif analysisParams.stimType == 5
                            if con == 1
                                title('Dir pref map, contra')
                            elseif con == 2
                                title('Dir pref map, ipsi')
                            elseif con == 3
                                title('Dir pref map, binocular')
                            end
                        end
                    else
                        imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                        colormap(hsv)
                        if analysisParams.stimType == 2
                            title(['Dir pref map, ' num2str(spatFreq_cell{con}) ' cpd'])
                        elseif analysisParams.stimType == 3
                            title(['Dir pref map, ' num2str(tempFreq_cell{con}) ' Hz'])
                        elseif analysisParams.stimType == 5
                            if con == 1
                                title('Dir pref map, contra')
                            elseif con == 2
                                title('Dir pref map, ipsi')
                            elseif con == 3
                                title('Dir pref map, binocular')
                            end
                        end
                        caxis([0 360])
                        colorbar('Location', 'southoutside');
                    end
                end
            else
                imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                colormap(hsv)
                if types == 1
                    if analysisParams.stimType == 2
                        title(['Ori pref map, ' num2str(spatFreq_cell{con}) ' cpd'])
                    elseif analysisParams.stimType == 3
                        title(['Ori pref map, ' num2str(tempFreq_cell{con}) ' Hz'])
                    elseif analysisParams.stimType == 5
                            if con == 1
                                title('Ori pref map, contra')
                            elseif con == 2
                                title('Ori pref map, ipsi')
                            elseif con == 3
                                title('Ori pref map, binocular')
                            end
                    end
                    caxis([0 180])
                elseif types == 2
                    if analysisParams.stimType == 2
                        title(['Dir pref map, ' num2str(spatFreq_cell{con}) ' cpd'])
                    elseif analysisParams.stimType == 3
                        title(['Dir pref map, ' num2str(tempFreq_cell{con}) ' Hz'])
                    elseif analysisParams.stimType == 5
                            if con == 1
                                title('Dir pref map, contra')
                            elseif con == 2
                                title('Dir pref map, ipsi')
                            elseif con == 3
                                title('Dir pref map, binocular')
                            end
                    end
                    caxis([0 360])
                end
                 colorbar('Location', 'southoutside');
            end
        end
        if types == 1
            saveas(gcf, fullfile(saveDirectory, 'OrientationPreference_allCon.png'))
        elseif types == 2
            saveas(gcf, fullfile(saveDirectory, 'DirectionPreference_allCon.png'))
        end
    end 
end



% 3.3 Plot percentages of all ROIs
figure
subplot(1,3,1)
all = length(analysis.(analysisParams.field).roi);
non_resp = length(find([analysis.(analysisParams.field).roi.isResponseSignificant] == 0)) ./all;
resp = length(find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1)) ./all;
h = pie([non_resp resp]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', coc_prop(7,:));
try set(hp(2), 'FaceColor', coc_prop(8,:)); end
title('Responsive')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,2)
ori = length(find([analysis.(analysisParams.field).roi.OSIFit] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1)) ./all;
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
dirResp = length(find([analysis.(analysisParams.field).roi.DSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1)) ./all;
non_dir = 1-dirResp-non_resp;
h = pie([non_resp non_dir dirResp]);
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

% 3.4 Plot OSI/DSI metrices
% for all cells
%allOSIFit = [analysis.(analysisParams.field).roi.OSIFit];

figure
subplot(1,6,1)
distributionPlot([analysis.(analysisParams.field).roi.OSI]','color', coc_prop(2,:)); hold all
boxplot([analysis.(analysisParams.field).roi.OSI])
title('OSI')
subplot(1,6,2)
distributionPlot([analysis.(analysisParams.field).roi.OSIFit]','color', coc_prop(2,:)); hold all
boxplot([analysis.(analysisParams.field).roi.OSI])
title('OSI after fitting')
subplot(1,6,3)
distributionPlot([analysis.(analysisParams.field).roi.OriCircVar]','color', coc_prop(2,:)); hold all
boxplot([analysis.(analysisParams.field).roi.OriCircVar])
title('CircVar')
subplot(1,6,4)
distributionPlot([analysis.(analysisParams.field).roi.cohensD]','color', coc_prop(2,:)); hold on
boxplot([analysis.(analysisParams.field).roi.cohensD])
title('cohensD')
subplot(1,6,5)
distributionPlot([analysis.(analysisParams.field).roi.DSI]','color', coc_prop(4,:)); hold on
boxplot([analysis.(analysisParams.field).roi.cohensD])
title('DSI')
subplot(1,6,6)
distributionPlot([analysis.(analysisParams.field).roi.DSI]','color', coc_prop(4,:)); hold on
boxplot([analysis.(analysisParams.field).roi.DirCircVar])
title('DirCircVar')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'OSI_DSI_distribution_all.png'))

% for resp cells only
resp_cells = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
%respOSI = [analysis.(analysisParams.field).roi(resp_cells).OSIFit];
if ~isempty(resp_cells)
    figure
    subplot(1,6,1)
    distributionPlot([analysis.(analysisParams.field).roi(resp_cells).OSI],'color', coc_prop(2,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(resp_cells).OSI])
    title('OSI')
    subplot(1,6,2)
    distributionPlot([analysis.(analysisParams.field).roi(resp_cells).OSIFit],'color', coc_prop(2,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(resp_cells).OSI])
    title('OSI after fitting')
    subplot(1,6,3)
    distributionPlot([analysis.(analysisParams.field).roi(resp_cells).OriCircVar],'color', coc_prop(2,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(resp_cells).OriCircVar])
    title('CircVar')
    subplot(1,6,4)
    distributionPlot([analysis.(analysisParams.field).roi(resp_cells).cohensD],'color', coc_prop(2,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(resp_cells).cohensD])
    title('cohensD')
    subplot(1,6,5)
    distributionPlot([analysis.(analysisParams.field).roi(resp_cells).DSI],'color', coc_prop(4,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(resp_cells).DSI])
    title('DSI')
    subplot(1,6,6)
    distributionPlot([analysis.(analysisParams.field).roi(resp_cells).DirCircVar],'color', coc_prop(4,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(resp_cells).DirCircVar])
    title('DirCircVar')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, 'OSI_DSI_distribution_resp.png'))
end

% 3.5 Plot deltaOri vs. distance of ROIs
try 
    delta_ori = [analysis.(analysisParams.field).ori_pair.deltaOri]; 
    maxDis = floor(max([analysis.(analysisParams.field).ori_pair.distance]));
    edges = linspace(0, maxDis, 15);
    [n, edges, bin_dist] = histcounts([analysis.(analysisParams.field).ori_pair.distance],edges);
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

    figure
    errorbar(edge,mean_deltaOri,SEM_deltaOri, 'o-', 'Color', coc_prop(2,:), 'MarkerFaceColor', coc_prop(2,:))
    hold all
    errorbar(edge,mean_deltaOri_shuffle,SEM_deltaOri_shuffle, 'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
    xlabel('Distance in \mum')
    ylabel('\DeltaOrientation preference (\circ)')
    ylim([0 90])
    xlim([0 maxDis])
    legend('Data', 'Shuffled')
    legend('boxoff')
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, 'OSI_distance.png'))
end

% 3.6 Plot homogeneity index 
try 
    HI = analysis.(analysisParams.field).ori_cells.HomeogeneityIndex;
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
end

% 3.7 Plot decoder performance
if analysisParams.predictor
    figure    
    DecoderPerformanceAll = [analysis.(analysisParams.field).decoderData{1}, analysis.(analysisParams.field).decoderData{2}];
    groups = {'Data', 'Shuffle'};
    boxplot(DecoderPerformanceAll, 'labels', groups);
    ylim([0 1])
    ylabel('Fraction of correctly decoded trials');
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, 'Decoder.png'))
end

% 3.8 Plot ROIs
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
    metadata.StimParams.numCon = metadata.StimParams.numSf*metadata.StimParams.numTf;
    for i = 1:length(data.roi)
        if analysis.(analysisParams.field).roi(i).isResponseSignificant == 1
            try
                PlotAvgStimResponseOri(metadata, analysis, analysisParams.field, i)
                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
                close gcf
                PlotTrialStimResponseOri(metadata, analysis, analysisParams.field, i)
                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                close gcf
            end
        else
            try
                if ~analysisParams.plotRespROIsOnly
                    PlotAvgStimResponseOri(metadata, analysis, analysisParams.field, i)
                    saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp.png']))
                    close gcf
                    PlotTrialStimResponseOri(metadata, analysis, analysisParams.field, i)
                    saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                    close gcf
                end
            end
        end
    end

end

end