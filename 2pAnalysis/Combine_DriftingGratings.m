%% set parameters
clear all
close all
TwoPhontondir = 'E:\Data\2P_Data\';
Sp2dir = 'E:\Data\Spike2Data\';
savedir = 'E:\Data\ImageAnalysis\';

animal = 'F2363_2019-09-13';
exp_id{1}  = 't00021'; exp_id{2}  = 't00023'; %exp_id{3}  = 't00041';
name = 'FoV1';

windowStop=2;
windowStart=0;
pre=1;
field = 'dff';

z_thresh = 5;
fraction = 0.5;
plotROIsResps = 1;

saveDirectory = [savedir animal filesep name filesep];
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
%% load Data

for exp = 1:length(exp_id)
    datapath = [savedir animal filesep exp_id{exp} filesep];
    filename = [datapath filesep '*ana.mat'];
    files = dir(filename);
    master{exp} = load(fullfile(datapath, files(1).name), 'data', 'metadata', 'analysis');
    disp(['loaded Data from exp ' num2str(exp_id{exp})])
end

combined = struct;

%% asumption: all exp have the same ROIs and the same number of directions that are shown
%create stimcodes -> how many different conditions do we have per exp?
disp('Calculating significant responses')

combined.metadata.numDirections = length(master{1}.metadata.StimParams.directions);
combined.metadata.theta = master{1}.metadata.StimParams.theta;
combined.metadata.directions = master{1}.metadata.StimParams.directions;
combined.metadata.template = master{1}.data.template;
combined.metadata.preTrialTime = master{1}.analysis.dff.preTrialTime;
combined.metadata.postTrialTime = master{1}.analysis.dff.postTrialTime;
combined.metadata.windowStart = master{1}.analysis.dff.windowStart;
combined.metadata.windowStop = master{1}.analysis.dff.windowStop;
combined.metadata.TwoPhoton.rate = master{1}.metadata.TwoPhoton.rate;
combined.analysis.roi = struct;

for i = 1:length(master{1}.data.roi)
    combined.analysis.roi(i).xPos = master{1}.data.roi(i).xPos;
    combined.analysis.roi(i).yPos = master{1}.data.roi(i).yPos;
    combined.analysis.roi(i).name = master{1}.data.roi(i).name;
end

counter = 1;
sumTrials = 1;
for exp = 1:length(exp_id)
    num_cond = (master{exp}.metadata.StimParams.uniqStims-1)/combined.metadata.numDirections;
    for cond = 1:num_cond
        try 
            tempFreq_cell = strsplit(master{exp}.metadata.StimParams.temporalFreq(2:end-1), ',');
            combined.metadata.tempFreq(counter) = str2num(tempFreq_cell{cond});
        catch combined.metadata.tempFreq(counter) = master{exp}.metadata.StimParams.temporalFreq; end
        try 
            spatFreq_cell = strsplit(master{exp}.metadata.StimParams.spatialFreq(2:end-1), ',');
            combined.metadata.spatialFreq(counter) = str2num(spatFreq_cell{cond});
        catch combined.metadata.spatialFreq(counter) = master{exp}.metadata.StimParams.spatialFreq; end
        combined.metadata.stimDuration(counter) = master{exp}.metadata.StimParams.stimDuration;
        combined.metadata.isi(counter) = master{exp}.metadata.StimParams.isi;
        combined.metadata.numTrials(counter) = master{exp}.metadata.StimParams.numTrials;
        counter = counter +1;
    end
end
numConditions = counter-1;
combined.metadata.numConditions = numConditions;

for i = 1:length(combined.analysis.roi)
    combined.analysis.roi(i).peaks = master{1}.analysis.dff.roi(i).peaks(1:end-1,1:combined.metadata.numTrials(1));
    combined.analysis.roi(i).zscore = master{1}.analysis.dff.roi(i).zscore(1:end-1,1:combined.metadata.numTrials(1));
    combined.analysis.roi(i).crosser = master{1}.analysis.dff.roi(i).crosser(1:end-1);
    combined.analysis.roi(i).respStim = master{1}.analysis.dff.roi(i).respStim(1:end-1);
    combined.analysis.roi(i).avgStimResponse = master{1}.analysis.dff.roi(i).avgStimResponse(1:end-1);
    combined.analysis.roi(i).stimResponseTrace = master{1}.analysis.dff.roi(i).stimResponseTrace(1:end-1,1:combined.metadata.numTrials(1),:);
    combined.analysis.roi(i).avgResponseTrace = master{1}.analysis.dff.roi(i).avgResponseTrace(1:end-1,:);
    combined.analysis.roi(i).SEMResponseTrace = master{1}.analysis.dff.roi(i).SEMResponseTrace(1:end-1,:);
    for exp = 2:length(exp_id)
        combined.analysis.roi(i).stimResponseTrace =  cat(1,combined.analysis.roi(i).stimResponseTrace,master{exp}.analysis.dff.roi(i).stimResponseTrace(1:end-1,1:combined.metadata.numTrials(1),:));
        combined.analysis.roi(i).avgResponseTrace =  cat(1,combined.analysis.roi(i).avgResponseTrace,master{exp}.analysis.dff.roi(i).avgResponseTrace(1:end-1,:));
        combined.analysis.roi(i).SEMResponseTrace =  cat(1,combined.analysis.roi(i).SEMResponseTrace,master{exp}.analysis.dff.roi(i).SEMResponseTrace(1:end-1,:));
        combined.analysis.roi(i).peaks =  cat(1,combined.analysis.roi(i).peaks,master{exp}.analysis.dff.roi(i).peaks(1:end-1,1:combined.metadata.numTrials(1)));
        combined.analysis.roi(i).zscore =  cat(1,combined.analysis.roi(i).zscore,master{exp}.analysis.dff.roi(i).zscore(1:end-1,1:combined.metadata.numTrials(1)));
        combined.analysis.roi(i).crosser =  [combined.analysis.roi(i).crosser; master{exp}.analysis.dff.roi(i).crosser(1:end-1)];
        combined.analysis.roi(i).respStim =  [combined.analysis.roi(i).respStim; master{exp}.analysis.dff.roi(i).respStim(1:end-1)];
        combined.analysis.roi(i).avgStimResponse =  [combined.analysis.roi(i).avgStimResponse; master{exp}.analysis.dff.roi(i).avgStimResponse(1:end-1)];
    end
    if sum(combined.analysis.roi(i).respStim) > 0
        combined.analysis.roi(i).isResponseSignificant = 1;
    else 
        combined.analysis.roi(i).isResponseSignificant = 0;
    end
end
clear master

%% find preferred orientation for all cells
disp('Calculating ROI properties over all stimuli')

orientations = linspace(0, 180, combined.metadata.numDirections/2+1);
lastOrientation = combined.metadata.numDirections/2;
theta = combined.metadata.theta;
theta=mod(theta, pi);
stimWindow = (combined.metadata.windowStart:combined.metadata.windowStop);

for i = 1:length(combined.analysis.roi)
    medResponse_ori = zeros(combined.metadata.numDirections,combined.metadata.numTrials(1),size(combined.analysis.roi(i).stimResponseTrace,3));
    endnum = numConditions * combined.metadata.numDirections;
    %collapse all conditions for the same direction
    for dir = 1:combined.metadata.numDirections
        medResponse_ori(dir,:,:) = mean(combined.analysis.roi(i).stimResponseTrace(dir:combined.metadata.numDirections:endnum, 1:combined.metadata.numTrials(1), :));
    end
    %calculate parameters
    medResponse_ori = mean(medResponse_ori(:,:,stimWindow),3);
    medResponse_dir = median(medResponse_ori,2);
    medResponseA(1:size(medResponse_ori,1)/2, 1:size(medResponse_ori,2)) = medResponse_ori(1:size(medResponse_ori,1)/2,:);
    medResponseA(1:size(medResponse_ori,1)/2, 1+size(medResponse_ori,2):2*size(medResponse_ori,2))= medResponse_ori(1+size(medResponse_ori,1)/2:size(medResponse_ori,1),:);
    medResponseA= median(medResponseA,2);
    medResponse_ori = medResponseA(:)';
    theta=theta(1:lastOrientation);
    theta=repmat(theta, 1, length(medResponse_ori)/length(theta));
    
    [combined.analysis.roi(i).preferredOrientation,combined.analysis.roi(i).coeffOr,combined.analysis.roi(i).rsqOr] = ComputePreferredOrientations(medResponse_ori, theta);
    combined.analysis.roi(i).preferredOrientation =mod(combined.analysis.roi(i).preferredOrientation, 180);
    combined.analysis.roi(i).rsqOr= Rsquared(medResponse_ori, vonMisesFit(combined.analysis.roi(i).coeffOr, theta*2), true);
    combined.analysis.roi(i).OSIFit = computeOSIFit(combined.analysis.roi(i).coeffOr);
    [~, ind] = min(abs(orientations-combined.analysis.roi(i).preferredOrientation));
    combined.analysis.roi(i).prefOriStim = orientations(ind);
    combined.analysis.roi(i).prefOriStimInd = ind;
    medResponse_ori = mean(combined.analysis.roi(i).avgStimResponse, 2)';
    combined.analysis.vsOrientation(i)= wrapTo2Pi(angle(vectorSum(abs(medResponse_ori(1:lastOrientation)),2)))/2 * 180/pi;
    combined.analysis.orientationMag(i) = abs(vectorSum(abs(medResponse_ori(1:lastOrientation))/max(medResponse_ori(1:lastOrientation)),2));
    
    %find preferred direction for all cells
    bestCoeff=[NaN, NaN, NaN, NaN];
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    bestCoeff=fitVonMisesLinSumFunction(medResponse_dir', combined.metadata.theta, bestCoeff, options);
    xo=[0:pi/40:2*pi];
    fit=vonMisesLinSum(bestCoeff, xo);
    if bestCoeff(1) > bestCoeff(2)
        combined.analysis.roi(i).preferredDirection = rad2deg(bestCoeff(4));
    else
        combined.analysis.roi(i).preferredDirection = mod(rad2deg(bestCoeff(4))+180,360);
    end
end

%% find preferred condition for each cell by looking at preferred Orientation
disp('Calculating ROI properties for each condition')
for i = 1:length(combined.analysis.roi)
    ResponseCondPrefOri = zeros(numConditions,combined.metadata.numTrials(1),size(combined.analysis.roi(i).stimResponseTrace,3));
    for con = 1:numConditions
        temp = combined.analysis.roi(i).stimResponseTrace(1+(con-1)*lastOrientation*2:con*lastOrientation*2, 1:combined.metadata.numTrials(1), :);
        ResponseCondPrefOri(con, :, :) = (temp(combined.analysis.roi(i).prefOriStimInd,:,:)+temp(combined.analysis.roi(i).prefOriStimInd+lastOrientation/2,:,:)/2);
    end
    avgResponseConprefOri = mean(ResponseCondPrefOri(:,:,stimWindow),3);
    medResponseConprefOri = mean(avgResponseConprefOri,2)';
    semResponseConprefOri = (std(avgResponseConprefOri,[],2)/sqrt(size(medResponseConprefOri,2)))';
    PrefConInd = find(medResponseConprefOri == max(medResponseConprefOri));
    maxResp = max(medResponseConprefOri);
    minResp = min(medResponseConprefOri);
    CSI = (maxResp - minResp) ./ (maxResp + minResp);
    combined.analysis.roi(i).prefConStim = PrefConInd;
    combined.analysis.roi(i).prefSf = combined.metadata.spatialFreq(PrefConInd);
    combined.analysis.roi(i).prefTf = combined.metadata.tempFreq(PrefConInd);
    combined.analysis.roi(i).CSI = CSI;
end

%% compute OSI and DSI at prefered condition
for i=1:length(combined.analysis.roi)
    prefConResponse = combined.analysis.roi(i).stimResponseTrace((combined.analysis.roi(i).prefConStim-1)*lastOrientation*2+1:combined.analysis.roi(i).prefConStim*lastOrientation*2, 1:combined.metadata.numTrials(1), :);
    prefConResponse = mean(prefConResponse(:,:,stimWindow),3);
    prefConResponse_dir = median(prefConResponse,2)';
    prefConResponse(1:size(prefConResponse,1)/2, 1:size(prefConResponse,2)) = prefConResponse(1:size(prefConResponse,1)/2,:);
    prefConResponse(1:size(prefConResponse,1)/2, 1+size(prefConResponse,2):2*size(prefConResponse,2))= prefConResponse(1+size(prefConResponse,1)/2:size(prefConResponse,1),:);
    prefConResponse= median(prefConResponse,2);
    prefConResponse_ori = prefConResponse(:)';
    [combined.analysis.roi(i).OSI, ~] = computeOSI(lastOrientation,prefConResponse_ori(1:lastOrientation));
    [combined.analysis.roi(i).DSI, ~] = computeDSI(lastOrientation*2, prefConResponse_dir, 'preferred');
end

%% compute ori preference at each condition
prefDir_allCon = zeros(numConditions,1);
prefOri_allCon = zeros(numConditions,1);
OSIFit_allCon = zeros(numConditions,1);
OSI_allCon = zeros(numConditions,1);
DSI_allCon = zeros(numConditions,1);
isResp_allCon = zeros(numConditions,1);

for i = 1:length(combined.analysis.roi)
    for con = 1:numConditions
        medResponse = combined.analysis.roi(i).stimResponseTrace((con-1)*lastOrientation*2+1:con*lastOrientation*2, :,:);
        medResponse = mean(medResponse(:,:,stimWindow),3);
        dirResponse = median(medResponse,2)';
        medResponseA(1:size(medResponse,1)/2, 1:size(medResponse,2)) = medResponse(1:size(medResponse,1)/2,:);
        medResponseA(1:size(medResponse,1)/2, 1+size(medResponse,2):2*size(medResponse,2))= medResponse(1+size(medResponse,1)/2:size(medResponse,1),:);
        medResponseA= median(medResponseA,2);
        medResponse = medResponseA(:)';
        
        [prefOri_allCon(con), coeffOr, ~] = ComputePreferredOrientations(medResponse, theta);
        prefOri_allCon(con) =mod(prefOri_allCon(con), 180);
        OSIFit_allCon(con) = computeOSIFit(coeffOr);
        
        bestCoeff=fitVonMisesLinSumFunction(dirResponse, combined.metadata.theta, bestCoeff, options);
        xo=[0:pi/40:2*pi];
        fit=vonMisesLinSum(bestCoeff, xo);
        if bestCoeff(1) > bestCoeff(2)
            prefDir_allCon(con) = rad2deg(bestCoeff(4));
        else
            prefDir_allCon(con) = mod(rad2deg(bestCoeff(4))+180,360);
        end
        
        OSI_allCon(con) = computeOSI(lastOrientation,medResponse(1:lastOrientation));
        DSI_allCon(con) = computeDSI(lastOrientation*2, dirResponse, 'preferred');
        
        tempCrosser = combined.analysis.roi(i).crosser((con-1)*lastOrientation*2+1:con*lastOrientation*2);
        tempRespStim = tempCrosser > ((combined.metadata.numTrials(1))*fraction);
        if sum(tempRespStim) > 0
            isResp_allCon(con) = 1;
        else
            isResp_allCon(con) = 0;
        end
    end
    combined.analysis.roi(i).preferredOrientation_allCon = prefOri_allCon;
    combined.analysis.roi(i).preferredDirection_allCon = prefDir_allCon;
    combined.analysis.roi(i).OSIFit_allCon = OSIFit_allCon;
    combined.analysis.roi(i).OSI_allCon = OSI_allCon;
    combined.analysis.roi(i).DSI_allCon = DSI_allCon;
    combined.analysis.roi(i).isResp_allCon = isResp_allCon;
end

%% plot results
coc_prop = cbrewer('qual', 'Paired', 12);

%plot prefOri on top of template
for types = 1:3
    if types == 1
        rois_ori = linspace(1,length(combined.analysis.roi),length(combined.analysis.roi));
        alloriprefs = [combined.analysis.roi.preferredOrientation];
        alldirprefs = [combined.analysis.roi.preferredDirection];
        allconpref = [combined.analysis.roi.prefConStim];
        rois_dir = rois_ori;
        rois_con = rois_ori;
    elseif types == 2
        rois_ori = find([combined.analysis.roi.isResponseSignificant] == 1); 
        alloriprefs = [combined.analysis.roi(rois_ori).preferredOrientation];
        alldirprefs = [combined.analysis.roi(rois_ori).preferredDirection];
        allconpref = [combined.analysis.roi(rois_ori).prefConStim];
        rois_dir = rois_ori;
        rois_con = rois_ori;
    elseif types == 3
        rois_ori = find([combined.analysis.roi.OSIFit] > 0.2 & [combined.analysis.roi.isResponseSignificant] == 1);
        alloriprefs = [combined.analysis.roi(rois_ori).preferredOrientation];
        rois_dir = find([combined.analysis.roi.DSI] > 0.2 & [combined.analysis.roi.isResponseSignificant] == 1);
        alldirprefs = [combined.analysis.roi(rois_dir).preferredDirection];
        rois_con = find([combined.analysis.roi.CSI] > 0.2 & [combined.analysis.roi.isResponseSignificant] == 1);
        allconpref = [combined.analysis.roi(rois_con).prefConStim];
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    %plot prefOri on top of template
    subplot(2,3,1)
    PlotPrefOnTemplate(combined,1, rois_ori)
    
    subplot(2,3,4)
    histogram(alloriprefs,linspace(0,180,5), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
    ylabel('Cells');
    xlabel(sprintf('Orientation preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)]/2)
    axis square;
    set(gca,'Box','off');
    
    subplot(2,3,2)
    PlotPrefOnTemplate(combined, 3, rois_dir)

    subplot(2,3,5)
    histogram(alldirprefs,linspace(0,360,lastOrientation+1), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
    ylabel('Cells');
    xlabel(sprintf('Direction preference (%s)',char(145)));
    xlim([-22.5 (360+22.5)])
    axis square;
    set(gca,'Box','off');
    
    %plot prefCon on top of template
    subplot(2,3,3)
    PlotPrefOnTemplate(combined,2, rois_con)

    subplot(2,3,6)
    conCat = linspace(1, numConditions,numConditions);
    numCounts = histcounts(allconpref, numConditions);
    bar(conCat, numCounts, 'FaceColor', coc_prop(5,:), 'EdgeColor', coc_prop(6,:));
    title('Histogram')
    ylabel('Cells');
    xlabel(sprintf('Condition preference'));
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
    rows = ceil(sqrt(numConditions));
    for con = 1:numConditions
        subplot(rows,rows,con)
        tempResp = [combined.analysis.roi.isResp_allCon];
        isResponsive = find(tempResp(con,:) == 1);
        if sum(isResponsive) > 0
            if types == 1
                tempOSI= [combined.analysis.roi.OSIFit_allCon];
                tempOriSelect = find(tempOSI(con,:) > 0.2);
                rois = intersect(tempOriSelect, isResponsive);
                tempPrefOri = [combined.analysis.roi(rois).preferredOrientation_allCon];
                allprefs = tempPrefOri(con,:);
                PlotPrefOnTemplateCon(allprefs, combined,1,rois)
                title(['Ori pref, ' num2str(combined.metadata.spatialFreq(con)) ' cpd, ' num2str(combined.metadata.tempFreq(con)) ' Hz'])
            elseif types == 2
                tempDSI= [combined.analysis.roi.DSI_allCon];
                tempDirSelect = find(tempDSI(con,:) > 0.2);
                rois = intersect(tempDirSelect, isResponsive);
                tempPrefDir = [combined.analysis.roi(rois).preferredDirection_allCon];
                allprefs = tempPrefDir(con,:);
                PlotPrefOnTemplateCon(allprefs, combined,3,rois)
                title(['Dir pref, ' num2str(combined.metadata.spatialFreq(con)) ' cpd, ' num2str(combined.metadata.tempFreq(con)) ' Hz'])
            end
        else
            imshow(cat(3,combined.metadata.template,combined.metadata.template,combined.metadata.template)/prctile(combined.metadata.template(:),99));
            colormap(hsv)
            if types == 1
                title(['Ori pref, ' num2str(combined.metadata.spatialFreq(con)) ' cpd, ' num2str(combined.metadata.tempFreq(con)) ' Hz'])
                caxis([0 180])
            elseif types == 2
                title(['Dir pref, ' num2str(combined.metadata.spatialFreq(con)) ' cpd, ' num2str(combined.metadata.tempFreq(con)) ' Hz'])
                caxis([0 180])
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

%% plot all ROIs
if plotROIsResps
    for i = 1:length(combined.analysis.roi)
        if combined.analysis.roi(i).isResponseSignificant == 1
            PlotTrialStimResponse(combined,i)
            saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        else
            PlotTrialStimResponse(combined, i)
            saveas(gcf, fullfile(ROINonRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
            close gcf
        end
    end
end

save(fullfile(saveDirectory, 'OriCombined_ana.mat'), 'combined');

%% additional functions

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
function PlotPrefOnTemplate(combined, type, rois)
    imshow(cat(3,combined.metadata.template,combined.metadata.template,combined.metadata.template)/prctile(combined.metadata.template(:),99));
    colormap(hsv) 
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        title('Condition preference map')
        numCon = combined.metadata.numConditions;
        cocCon = cbrewer('qual', 'Set1', numCon);
    elseif type == 3
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        l = rois(i);
        xpos= combined.analysis.roi(l).xPos;
        ypos= combined.analysis.roi(l).yPos;
        try
            if type == 1
                sizeMarker = ceil(15 * combined.analysis.roi(l).OSIFit);
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor(combined.analysis.roi(l).preferredOrientation),:));
            elseif type == 2
                sizeMarker = ceil(15 * combined.analysis.roi(l).CSI);
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',cocCon(combined.analysis.roi(l).prefConStim,:));
            elseif type == 3
                sizeMarker = ceil(15 * combined.analysis.roi(l).DSI);
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor(combined.analysis.roi(l).preferredDirection),:));
            end
        catch
        end
    end
end
function PlotPrefOnTemplateCon(pref, combined, type, rois)
    %h=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(cat(3,combined.metadata.template,combined.metadata.template,combined.metadata.template)/prctile(combined.metadata.template(:),99));
    colormap(hsv) 
    if type == 1
        LUT = hsv(180);
        title('Orientation preference map')
        caxis([0 180]); colorbar('Location', 'southoutside');
    elseif type == 2
        title('Condition preference map')
        numCon = combined.metadata.numConditions;
        cocCon = cbrewer('qual', 'Set1', numCon);
    elseif type == 3
        LUT = hsv(360);
        title('Direction preference map')
        caxis([0 360]); colorbar('Location', 'southoutside');
    end
    axis image
    hold on
    for i = 1:length(rois)
        xpos= combined.analysis.roi(rois(i)).xPos;
        ypos= combined.analysis.roi(rois(i)).yPos;
        try
            if type == 1
                sizeMarker = ceil(20 * combined.analysis.roi(rois(i)).OSIFit);
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor((pref(i))),:));
            elseif type == 2
                sizeMarker = ceil(20 * combined.analysis.roi(rois(i)).CSI);
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',cocCon((pref(i))));
            elseif type == 3
                sizeMarker = ceil(20 * combined.analysis.roi(rois(i)).DSI);
                plot(xpos,ypos,'ok','MarkerSize',sizeMarker,'MarkerFaceColor',LUT(1+floor((pref(i))),:));
            end
        catch
        end
    end
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
function PlotTrialStimResponse(combined,roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    colorlevels = combined.metadata.numConditions;
    colorlevels_tr = combined.metadata.numTrials(1)+1;
    coc_sf = cbrewer('qual', 'Set1', colorlevels);
    coc_tr = [linspace(0,1,colorlevels_tr);linspace(0,1,colorlevels_tr);linspace(0,1,colorlevels_tr)]';
    coc_tr = coc_tr(2:end,:);

    
    pretrialTime= combined.metadata.preTrialTime;
    stimWindow=(combined.metadata.windowStart:combined.metadata.windowStop);
    ymax=max(combined.analysis.roi(roi).avgResponseTrace(:)+combined.analysis.roi(roi).SEMResponseTrace(:));
    ymin=min(combined.analysis.roi(roi).avgResponseTrace(:)-combined.analysis.roi(roi).SEMResponseTrace(:));

    for con = 1:combined.metadata.numConditions
        for i=1:combined.metadata.numDirections
            id = combined.metadata.numDirections*(con-1)+i;
            if combined.metadata.numDirections > 2
                ax=subplot(combined.metadata.numConditions,combined.metadata.numDirections, id);
            else
               ax=subplot(combined.metadata.numDirections,combined.metadata.numConditions, id);
            end
            cla(ax)
            y=medfilt1(combined.analysis.roi(roi).avgResponseTrace(id,:),3);
            x=(1:length(y))./combined.metadata.TwoPhoton.rate-pretrialTime;
            if stimWindow(1)~=0
                    bl=x(stimWindow);
                    patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
            end
            hold all
            for trial =1:size(combined.analysis.roi(roi).stimResponseTrace,2)
                plot(ax,x, squeeze(combined.analysis.roi(roi).stimResponseTrace(id,trial,:)), 'Color', coc_tr(trial,:));
                hold all
            end
            plot(ax,x, squeeze(combined.analysis.roi(roi).avgResponseTrace(id,:)), 'Color', coc_sf(con,:), 'LineWidth', 3);
            hold all
            ylim([ymin, ymax])
            xlim([min(x) max(x)])
            title(combined.metadata.directions(i))
            axis off
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end