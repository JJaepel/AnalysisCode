close all
clear all
animal = 'F2342_2019-05-28\';
TwoPhontondir = 'E:\Data\2P_Data\';
Sp2dir = 'E:\Data\Spike2Data\';
savedir = 'E:\Data\ImageAnalysis\';
expt_id='t00006';

reloadData = 1;
plotROIs = 1;
z_thresh = 2;
fraction = 0.25;
selective = 1;
shufflenum = 100;

totalSlices=1;
channelLoad= 1;
roiFiles='';
windowStop=2;
windowStart=0;
pre=1;
field = 'dff';

base2pDirectory= [TwoPhontondir animal];
Sp2dDirectory = [Sp2dir animal filesep expt_id filesep];
saveDirectory = [savedir animal filesep expt_id filesep];
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
if ~exist(ROIsaveDirectory, 'dir')
    mkdir(ROIsaveDirectory);  
end

for slice = 1:4
    %%load Data and metadata
    if reloadData
        sliceparams = struct;
        sliceparams.expt_id = expt_id;
        sliceparams.channel = channelLoad ;
        sliceparams.slicenum = slice;
        sliceparams.totalSlices = totalSlices;
        sliceparams.baseDirectory = base2pDirectory;
        metadata.StimParams=Load_stimparams(Sp2dDirectory);
        metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
        metadata.StimParams.path=fullfile(Sp2dDirectory);
        metadata.StimParams.series=expt_id;
        data = Load_rois(sliceparams, roiFiles);
        [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
        data = computeDff(data, 'rawF', 'baseline', 'dff');
        metadata.ROI = struct;
        analysis = struct;
        save(fullfile(saveDirectory, ['s' num2str(slice) '_ori_Grating.mat']), 'data', 'metadata', 'sliceparams', 'analysis');
    else
        load(fullfile(saveDirectory, ['s' num2str(slice) '_ori_Grating.mat']), 'data', 'metadata', 'sliceparams', 'analysis');
    end

    %%chop traces
    disp('Chopping Traces')
    [analysis, metadata, data] = ChopStimulusTrace(analysis,metadata,data,field, 'pre', pre, 'post',metadata.StimParams.isi,'windowStart',windowStart, 'windowStop',windowStop);

    %%find maxResponses and significant responses
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
        analysis.(field).roi(i).baselineMean(analysis.(field).roi(i).baselineMean < 0) = mean(mean(analysis.(field).roi(i).baselineMean,2),1);
        analysisPeriod=(analysis.(field).windowStart:analysis.(field).windowStop);
        stimResp = analysis.(field).roi(i).stimResponseTrace(:,:,analysisPeriod);
        analysis.(field).roi(i).peaks = max(stimResp,[],3);
        analysis.(field).roi(i).zscore = ([analysis.(field).roi(i).peaks]-[analysis.(field).roi(i).baselineMean])./[analysis.(field).roi(i).baselineSD];
        analysis.(field).roi(i).crosser = sum(analysis.(field).roi(i).zscore > z_thresh,2);
        analysis.(field).roi(i).respStim = analysis.(field).roi(i).crosser >= ((metadata.StimParams.numTrials)*fraction);
        if sum(analysis.(field).roi(i).respStim) > 0
            analysis.(field).roi(i).isResponseSignificant = 1;
        else 
            analysis.(field).roi(i).isResponseSignificant = 0;
        end
    end

    numberOfStims = metadata.StimParams.uniqStims;
    numStims=numberOfStims-1;
    metadata.StimParams.theta= [0:2*pi/(numStims) : 2*pi-2*pi/(numStims)];
    theta = metadata.StimParams.theta;
    theta=mod(theta, pi);
    orientations = linspace(0,180,numStims/2+1);
    directions = linspace(0,360,numStims+1);
    metadata.StimParams.directions = directions(1:end-1);
    orientations = orientations(1:end-1);
    lastori = length(orientations);
    for i = 1:length(data.roi)
        %%find preferred orientation for all cells
        medResponse = analysis.(field).roi(i).stimResponseTrace(1:end-1, :, :);
        medResponse = mean(medResponse(:,:,stimWindow),3);
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
        [analysis.(field).roi(i).CircVar] = 1-TT_CircularVariance(dirResponse);
    end

    %%load template
    imgFileName = [expt_id filesep 'Registered' filesep 'slice' num2str(slice) filesep 'slice' num2str(slice) '_c1_MasterTemplate.tif'];
    imgFile=fullfile(base2pDirectory,imgFileName);
    template=read_Tiffs(imgFile,1);
    template = double(template);
    data.template = template./prctile(template(:),99.9);

    [zoom, setup] = getzoom(base2pDirectory, expt_id);
    metadata.zoom = zoom;
    if setup == 1
        fieldofview = 1000/zoom;
        umperpixel = fieldofview/512;
        disp('Setup Ben')
    elseif setup == 2
        umperpixel = 2.73/zoom;
    end

    coc_prop = cbrewer('qual', 'Paired', 12);

    %plot prefOri on top of template
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
        PlotPrefOnTemplate(analysis, data, metadata,1, field,template, rois_ori)

        subplot(2,2,2)
        histogram(alloriprefs,linspace(0,180,lastori+1), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
        ylabel('Cells');
        xlabel(sprintf('Orientation preference (%s)',char(145)));
        xlim([-22.5 (360+22.5)]/2)
        axis square;
        set(gca,'Box','off');

        subplot(2,2,3)
        PlotPrefOnTemplate(analysis, data, metadata,3, field,template, rois_dir)

        subplot(2,2,4)
        histogram(alldirprefs,linspace(0,360,2*lastori+1), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
        ylabel('Cells');
        xlabel(sprintf('Direction preference (%s)',char(145)));
        xlim([-22.5 (360+22.5)])
        axis square;
        set(gca,'Box','off');
        set(gcf, 'color', 'w');
        if types == 1
            saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_Overlaymaps_all_cells.png']))
        elseif types == 2
            saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_Overlaymaps_resp_cells.png']))
        elseif types == 3
            saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_Overlaymaps_selective_cells.png']))
        end
        %close all
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
    saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_resp_cells.png']))
    %close gcf

    figure 
    allOSI = [analysis.dff.roi.OSIFit];
    resp_cells = find([analysis.dff.roi.isResponseSignificant] == 1);
    respOSI = [analysis.dff.roi(resp_cells).OSIFit];
    allDSI = [analysis.dff.roi.DSI];
    respDSI = [analysis.dff.roi(resp_cells).DSI];

    subplot(2,2,1)
    histogram(allOSI, 10, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel('OSI');
    title('All cells')
    set(gca,'Box','off');

    subplot(2,2,2)
    histogram(allDSI, 10, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel('DSI');
    title('All cells')
    set(gca,'Box','off');

    subplot(2,2,3)
    histogram(respOSI,10,'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel('OSI');
    title('All responsive cells')
    set(gca,'Box','off');

    subplot(2,2,4)
    histogram(respDSI,10, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel('DSI');
    title('All responsive cells')
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_OSI_DSI_distribtion.png']))
    %close gcf


    %plot delta Ori vs. distance of ROIs
    ori_sel = find([analysis.dff.roi.OSIFit] > 0.2 & [analysis.dff.roi.isResponseSignificant] == 1);
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
    saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_OSI_distance.png']))

    % calculate HetereogeneityIndex
    radius = linspace(0,250,6);
    HI = zeros(length(ori_sel),length(radius)-1);
    preferences = [analysis.dff.roi(ori_sel).preferredOrientation];
    for dis = 1:length(radius)-1
        for cellID = 1:length(ori_sel)
            distances = dist_ori(cellID,:);
            withinRadius1 = distances >= radius(dis);
            withinRadius2 = distances < radius (dis+1);
            withinRadius = withinRadius1&withinRadius2;
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

    subplot(1,3,3)
    avg_HI = nanmean(HI);
    sem_HI = nanstd(HI)/sqrt(length(HI));
    errorbar(radius(2:end), avg_HI, sem_HI,'o-', 'Color', coc_prop(6,:), 'MarkerFaceColor', coc_prop(6,:))
    xlabel('Maximal radius in \mum')
    ylabel('Homeogeneity Index')
    xlim([0 250])
    ylim([0 1])
    set(gca,'Box','off');
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, ['s' num2str(slice) '_HomeogeneityIndex.png']))

    %plot all ROIs
    if plotROIs
        for i = 1:length(data.roi)
            PlotAvgStimResponse(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIsaveDirectory, ['s' num2str(slice) '_ROI_Nr_' num2str(i) '_AvgStimResp.png']))
            close gcf
        end
        for i = 1:length(data.roi)
            PlotTrialStimResponse(metadata, analysis, field, i)
            saveas(gcf, fullfile(ROIsaveDirectory, ['s' numst2(slice) '_ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        end
    end
end

save(fullfile(saveDirectory, ['s' num2str(slice) '_ori_Grating_ana.mat']), 'data', 'metadata', 'sliceparams', 'analysis');

function [data]= Load_rois(sliceparams, roiFiles)
    %IRConfigs = ConfigurationParser.ReadImageRegistrationConfig;     
    slicenum = sliceparams.slicenum;
    baseDirectory = sliceparams.baseDirectory;
    channelLoad = sliceparams.channel;
    expt_id = sliceparams.expt_id;

    RM= ij.plugin.frame.RoiManager;
    RC= RM.getInstance();
    count= RC.getCount();
    if count > 0
        RC.setSelectedIndexes([0:count-1]);
        RC.runCommand('Delete');
    end
    roiFileName = [expt_id filesep 'Registered' filesep 'slice' num2str(slicenum) filesep 'RoiSet.zip'];
    
    if isempty(roiFiles)
        roiFiles= dir(fullfile(baseDirectory, roiFileName));
        roiFiles= {roiFiles.name};
        if ~isempty(roiFiles)
            roiFile=fullfile(baseDirectory,roiFileName);
        end
    end
    disp(['Loading ROI from: ', roiFile]);

    if RC.getCount() == 0 && isempty(roiFile)
        % check if there is a zip file in the expected directory
        fprintf('Could not find ROI file in %s', roiPath, '\n');
        fprintf('Load files for %s', expt_id, '\n');
        RC.runCommand('Open','');
    elseif RC.getCount() == 0
        RC.runCommand('Open', roiFile)
    end

    if RC.getCount() == 0
        error('There are still no ROIs defined')
    end
    
    metadata.File=struct;
    data.roi=[];
            
    for chIdx=1:length(channelLoad)
        channel=channelLoad(chIdx);
        sliceparams.channel= channel;
        metadata.File.baseDirectory = baseDirectory;
        metadata.File.series = expt_id;
        metadata.File.slicenum = slicenum;
        registeredImageDirectory = [baseDirectory filesep expt_id filesep 'Registered' filesep 'slice' num2str(slicenum) filesep];
        metadata.File.ImgData = registeredImageDirectory;
        registeredFileName = 'stack*.tif';

        disp(['Loading from ...' registeredImageDirectory])
        disp(['Looking for...' fullfile(registeredImageDirectory, registeredFileName)])
        files= dir(fullfile(registeredImageDirectory, registeredFileName));
        files=sort({files.name}'); 

        for fnum = 1:size(files,1)
           files{fnum} = fullfile(registeredImageDirectory, char(files{fnum})); 
        end

        NumberofTifs= size(files, 1);
        disp(['Found ' num2str(NumberofTifs) ' files'])
        % this will save the ROIs to a specified folder
        
        count = RC.getCount();
        for i =1:count
            tempROI= RC.getRoi(i-1);
            tempPos= tempROI.getContourCentroid();
            data.roi(i).xPos= tempPos(1);
            data.roi(i).yPos= tempPos(2);
            data.roi(i).mask= tempROI.getContainedPoints();
            data.roi(i).name= tempROI.getName;
        end

        %% Export ROI fluorescence traces
        for i=1:NumberofTifs
            if i==1
                for j=1:RC.getCount()
                    if channel ==1
                        data.roi(j).rawF=[];
                    end
                end
            end
            if i < 10
                registeredFileName = ['stack_c1_0' num2str(i) '.tif'];
            else
                registeredFileName = ['stack_c1_' num2str(i) '.tif'];
            end

            filename = fullfile(registeredImageDirectory, registeredFileName);
            if(~isempty(find(ismember(files, filename), 1)))
                data.roi=ExtractRawFMIJ(RC, data, filename, channel);
            end
        end
    end
end
function [cells]=ExtractRawFMIJ(RC, data,fullfilename, channel)
    %ExtractRawF generates an array of structs that contains the centroid of
    %ROI and the Raw fluorscence
    %   Args:
    %           RC: handler for ROI Manager
    %           cells(array): array of strucs
    %   Returns:
    %           cells(array): reassignments as below
    %               cells(i).name: ROI name
    %               cells(i).x: x location
    %               cells(i).y: y location
    %               cells(i).name: ROI name
    %               cells(i).name: ROI name
    %ij.IJ.run('Image Sequence...', strrep(sprintf('open=%s sort',filename),'\','\\'));
    cells=data.roi;
    disp(strcat('Loading...', fullfilename));
    MIJ.run('Open...', strrep(sprintf('path=%s', fullfilename), '\', '\\'));
    count = RC.getCount();
    targetField='rawF';
    if channel > 1
        targetField=strcat(targetField, '_c', num2str(channel));
    end
    for i = 1:count
        RC.select(i-1)
        cells(i).name = char(RC.getName(i-1));
        % This is the way the position was extracted by using ROI names
        %this method figures out x,y position by computing the centroid
        % Pull the fluorescence plot
        if ~isfield(cells(i), targetField)
            cells(i).(targetField)=[];
        end

            a=ij.plugin.ZAxisProfiler();
            a.run('ZPlot');
            zplot=a.getPlot();

            tmp = double(zplot.getYValues());
            MIJ.run('Close','');
            cells(i).(targetField) = [cells(i).(targetField), tmp(:,1)'];
    end
    MIJ.run('Close', '');
end
function [twophoton] = LoadFrameTimes(path)
    %ExtractFrameTimes Reads in Frame Triggers from frametrigger.txt
    %  	Args:
    %     path: path of the metadata file to extract
    %   Returns:
    %     twophoton (struct):   .volume :trigger times of volume(s)
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
            twophoton.volume = twophoton.time(1:5:end);

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
function [zoom, setup] = getzoom(base2pDirectory, expt_id)
    filename = [base2pDirectory filesep expt_id '*.tif'];    
    files = dir(filename);
    filepath = [base2pDirectory files(1).name];
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
        case 'Spontaneous'
            fields ={};
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
        case 'fullScreenFlash'
            fields = {'initialDelay',...
                'lum0',...
                'lum1',...
                'flashInterval',...
                'isi',...
                'numTrials'...
                };
        case 'waveletStim_forJoe'
            fields = {'Iname',...
                'repeats',...
                'interRepeatInt',...
                'interMovieInt',...
                'stimFramesPerTrial',...
                'movies'

                };
        case 'continuousEdge_withOrientationsAndTriggers'
            fields = {'numOrientations', ...
                'animalOrientation',...
                'barColor',...
                'backGroundColor',...
                'backGroundBarColor',...
                'isRandom',...
                'edgeWidth',...
                'centerPoint',...
                'numberOfTrials',...
                'initialDelay',...
                'movementPeriod',...
                'animalName'
                };
        case 'MovingDots'
            fields = {'numTrials',...
                'doBlank',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'nDots',...
                'dotSize',...
                'contrast',...
                'fieldShape',...
                'dotLife',...
                'minv',...
                'fieldSize',...
                'speed',...
                'coherence',...
                'theta',...
                };
        case 'Patch'
            fields = {'numTrials',...
                'doBlank',...
                'stimDuration',...
                'isi',...
                'isRandom',...
                'initialDelay',...
                'centerPoint',...
                'stimSize',...
                'numStimElev',...
                'numStimAzi',...
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
                    };
        case 'driftingGrating_ori_tf'
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

    TwoPhotonRate= metadata.TwoPhoton.rate/5;
    metdata.TwoPhoton.rate = metadata.TwoPhoton.rate/5;
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
    while (metadata.StimParams.StimOnTimes(2,metadata.StimParams.numberOfStims)+metadata.StimParams.stimDuration > metadata.TwoPhoton.volume(end))
        if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ==0
            metadata.StimParams.numTrials = metadata.StimParams.numTrials -1;
        else
            metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
        end
        metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
    end

    if mod(metadata.StimParams.numberOfStims, metadata.StimParams.uniqStims) ~=0
        metadata.StimParams.numTrials = floor(metadata.StimParams.numberOfStims/metadata.StimParams.uniqStims);
        metadata.StimParams.numberOfStims = metadata.StimParams.uniqStims * metadata.StimParams.numTrials;
    end

    for i =1:metadata.StimParams.numberOfStims
        metadata.StimParams.stimStartIndex(i)= find(metadata.TwoPhoton.volume > metadata.StimParams.StimOnTimes(2,i),1);
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

    for stimID= 1:metadata.StimParams.uniqStims
        stimulus = metadata.StimParams.uniqStimIds(stimID);
        StimonTimes= metadata.StimParams.StimOnTimes;
        StimonTimes=StimonTimes(:,1:metadata.StimParams.numberOfStims);
        stimIndices = find(StimonTimes(1,:)==stimulus);
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
        for stimID = 1:metadata.StimParams.uniqStims-1
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
function [CircularVariance, ResultantLength, ResultantAngle] = TT_CircularVariance( TuningCurve )
% function [CircularVariance, ResultantLength, ResultantAngle] = ...
%                         tpd_circularVariance( TuningCurve )
%
% Calculates circular variance, resultant length and angle of the 
% resultant based on the tuning curve. Assumes a full 360 degree measured
% tuningcurve.
%
% Input:
% - TuningCurve: Array with mean response per direction stimulus
%
% Ouput:
% - CircularVariance: Circular variance of the tuning curve
% - ResultantLength:  Length of the tuning curve resultant
% - ResultantAngle:   Angle of the tuning curve resultant
%
% Written by Pieter Goltstein
% Version 1.0: July 25th, 2011
%

    % get number of datapoints on full circle
    [a, nDataPoints] = size(TuningCurve);
    
    % remove negatives from tuning curve
    TuningCurve( TuningCurve < 0 ) = 0;

    % convert theta to complex plane
    Directions = (1:nDataPoints) * (360/nDataPoints);
    Theta = (Directions/360) * 2 * pi;
    Theta = exp(1i*Theta);  

    % calculate complex resultant
    Resultant = sum(TuningCurve.*Theta)/sum(TuningCurve);

    % calculate resultant length
    ResultantLength = abs(Resultant);

    % calculate resultant angle
    ResultantAngle = (angle(Resultant)/(2*pi))*360;

    %imaginary numbers give negative angles, so convert to positive ones
    while ResultantAngle < 0
        ResultantAngle = ResultantAngle + 360;
    end

    % calculates circular variance
    CircularVariance = 1 - abs(Resultant);

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
        title('Temporal frequency preference map')
        numtf = (metadata.StimParams.uniqStims-1)/metadata.StimParams.numOrientations;
        LUT = hsv(numtf);
        caxis([0 numtf]);colorbar('Location', 'southoutside');
    elseif type == 3
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        xpos= data.roi(i).xPos;
        ypos= data.roi(i).yPos;
        try
           if type == 1
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(1+floor(analysis.(field).roi(i).preferredOrientation),:));
            elseif type == 2
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(analysis.(field).roi(i).prefTfStim,:));
            elseif type == 3
                plot(xpos,ypos,'ok','MarkerSize',12,'MarkerFaceColor',LUT(1+floor(analysis.(field).roi(i).preferredDirection),:));
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
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    stimNum = metadata.StimParams.uniqStims-1;
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    for i=1:stimNum
        ax=subplot(1,stimNum, i);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
            bl=x(stimWindow);
            patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        err=analysis.(field).roi(roi).SEMResponseTrace(i,:);
        hold on
        patch([x fliplr(x)],[y+err fliplr(y-err)],[.5 .5 .5], 'LineStyle', 'none');
        %patch([bl fliplr(bl)],[zeros(1,length(bl)), ones(1,length(bl))],[0 0 0]);
        %plot(ax,x,y, 'k-')
        plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTraceRaw(i,:)), 'k');
        hold off
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        title(metadata.StimParams.directions(i))
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
    load('C:\Users\jaepelj\Dropbox\Work\colorbrewer.mat')
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
    colorlevels_tr = metadata.StimParams.numTrials;
    coc_tr = cbrewer('div', 'Spectral', colorlevels_tr);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).stimResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).stimResponseTrace(:));

    for i=1:metadata.StimParams.uniqStims-1
        ax=subplot(1,metadata.StimParams.uniqStims-1, i);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
                bl=x(stimWindow);
                patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        hold all
        for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
            plot(ax,x, smooth(squeeze(analysis.(field).roi(roi).stimResponseTrace(i,trial,:))), 'Color', coc_tr(trial,:));
            hold all
        end
        plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(i,:)), 'k');
        hold all
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        title(metadata.StimParams.directions(i))
        axis off
    end
    hold off
    set(gcf, 'Color', 'w')
end