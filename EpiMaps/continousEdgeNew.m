close all
clear all

animal = 'F2528_2021-05-28';
expt_id = 3;
sp2id = 3;
quantification = 0;

EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
SaveDir ='Z:\Juliane\Data\ImageAnalysis\';

cmap = jet;
OnSet=0;
OffSet = 9;
field = 'rawF';
intrinsic = 0;
responseSign = 1; % -1 for epi imaging, 1 for intrinsic imaging


EpiDirectory = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id) filesep];
if sp2id > 9
    Sp2dDirectory = [Sp2Dir animal filesep 't000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't000' num2str(expt_id) filesep];
else
    Sp2dDirectory = [Sp2Dir animal filesep 't0000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't0000' num2str(expt_id) filesep];
end

if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end

expParam = struct;
expParam.animal = animal;
expParam.expt_id = expt_id;
expParam.sp2id = sp2id;
expParam.baseDirectory = EpiDirectory;
expParam.Sp2dDirectory = Sp2dDirectory;
expParam.saveDirectory = saveDirectory;

%% load metadata
metadata.StimParams=LoadStimParams(Sp2dDirectory);
metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
metadata.StimParams.path=fullfile(Sp2dDirectory);
metadata.StimParams.series=expt_id;
switch metadata.StimParams.type
    case 'continuousEdge_withOrientationsAndTriggers'
        if isa(metadata.StimParams.startPoint, 'double')
            if metadata.StimParams.orientations == 90 || metadata.StimParams.orientations == 270
                metadata.StimParams.visualSpace = 'azimuth';
            elseif metadata.StimParams.orientations == 0 || metadata.StimParams.orientations == 180
                metadata.StimParams.visualSpace = 'elevation';
            end
        else
            if isempty(metadata.StimParams.right_screen)
                metadata.StimParams.startPoint = str2double(regexp({metadata.StimParams.startPoint},'-?[0-9]*', 'match', 'once'));
                metadata.StimParams.endPoint = str2double(regexp({metadata.StimParams.endPoint},'-?[0-9]*', 'match', 'once'));
            else
                if metadata.StimParams.orientations == 90
                    metadata.StimParams.startPoint = -metadata.StimParams.left_screen;
                    metadata.StimParams.endPoint = metadata.StimParams.right_screen;
                    metadata.StimParams.visualSpace = 'azimuth';
                elseif metadata.StimParams.orientations == 270
                    metadata.StimParams.startPoint = metadata.StimParams.right_screen;
                    metadata.StimParams.endPoint = -metadata.StimParams.left_screen;
                    metadata.StimParams.visualSpace = 'azimuth';
                elseif metadata.StimParams.orientations == 0
                    metadata.StimParams.startPoint = -metadata.StimParams.lower_screen;
                    metadata.StimParams.endPoint = metadata.StimParams.upper_screen; 
                    metadata.StimParams.visualSpace = 'elevation';
                elseif metadata.StimParams.orientations == 180
                    metadata.StimParams.startPoint = metadata.StimParams.upper_screen;
                    metadata.StimParams.endPoint =-metadata.StimParams.lower_screen; 
                    metadata.StimParams.visualSpace = 'elevation';
                end
            end
        end
    case 'continuousBar_withTriggers'
        if metadata.StimParams.stimID == 1
            metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
            metadata.StimParams.endPoint = - metadata.StimParams.endPoint; 
            metadata.StimParams.visualSpace = 'elevation';
        elseif metadata.StimParams.stimID == 3
            metadata.StimParams.visualSpace = 'azimuth';
        elseif    metadata.StimParams.stimID == 5
            metadata.StimParams.visualSpace = 'elevation';
        elseif   metadata.StimParams.stimID == 7
            metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
            metadata.StimParams.endPoint = - metadata.StimParams.endPoint; 
            metadata.StimParams.visualSpace = 'azimuth';
        end        
    case 'continuousBar_withTriggers_down'
        metadata.StimParams.orientations = 180;
        metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
        metadata.StimParams.endPoint = - metadata.StimParams.endPoint;
        metadata.StimParams.visualSpace = 'elevation';
    case 'continuousBar_withTriggers_up'
        metadata.StimParams.orientations = 0;
        metadata.StimParams.visualSpace = 'elevation';
    case 'continuousBar_withTriggers_left'
        metadata.StimParams.orientations = 90;
        metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
        metadata.StimParams.endPoint = - metadata.StimParams.endPoint;
        metadata.StimParams.visualSpace = 'azimuth';
    case 'continuousBar_withTriggers_right'
        metadata.StimParams.orientations = 270;
        metadata.StimParams.visualSpace = 'azimtuh';
    case 'continuousEdge_azimuth'
        if metadata.StimParams.orientations == 90
           metadata.StimParams.startPoint = -metadata.StimParams.left_screen;
           metadata.StimParams.endPoint = metadata.StimParams.right_screen;
        elseif metadata.StimParams.orientations == 90
            metadata.StimParams.startPoint = - metadata.StimParams.right_screen;
            metadata.StimParams.endPoint = metadata.StimParams.left_screen;
        end
        metadata.StimParams.visualSpace = 'azimuth';
    case 'continousEdge_elevation'
        if metadata.StimParams.orientations == 0
            metadata.StimParams.startPoint = -metadata.StimParams.lower_screen;
            metadata.StimParams.endPoint = metadata.StimParams.upper_screen;
        elseif metadata.StimParams.orientations == 180
            metadata.StimParams.startPoint = metadata.StimParams.upper_screen;
            metadata.StimParams.endPoint =-metadata.StimParams.lower_screen;
        end
        metadata.StimParams.visualSpace = 'elevation';
    case 'Retinotopy_LR'
        Direction = regexp(metadata.StimParams.Direction,'\w*','match', 'once');
        try 
            if Direction == 'left'
                metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
                metadata.StimParams.endPoint = -metadata.StimParams.endPoint;
            end
        catch            
        end
        metadata.StimParams.visualSpace = 'azimuth';
    case 'Retinotopy_UD'
        Direction = regexp(metadata.StimParams.Direction,'\w*','match', 'once');
        try 
            if Direction == 'down'
                metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
                metadata.StimParams.endPoint = -metadata.StimParams.endPoint;
            end
        catch
        end
        metadata.StimParams.visualSpace = 'elevation';
    case 'continuousEdge_left_RHField'
        metadata.StimParams.visualSpace = 'azimuth';
end
%% load tiffs
t0=tic;
tifStack = readingImagingData(EpiDirectory);
toc(t0)
expParam.ROI =true( [size(tifStack,1),size(tifStack,2)]); 
expParam.rawFMeanImg = mean(tifStack,3);
expParam.baseImg = mean(tifStack(:,:,1:50),3);
expParam.gaussMeanImg = imgaussfilt(mean(tifStack, 3), 4);

%% create stimCodes
metadata.Imaging.startTime = metadata.Imaging.time(1);
metadata.Imaging.offsetFrameTimes = metadata.Imaging.time-metadata.Imaging.startTime;
metadata.StimParams.stimConditionOnsetTime = metadata.StimParams.StimOnTimes(2:2:end);
metadata.Imaging.meanCCDTime = median(diff(metadata.Imaging.offsetFrameTimes));
stimConditionOnsetTimeCorr = metadata.StimParams.stimConditionOnsetTime-metadata.Imaging.startTime;

%% make fourierMaps
metadata.Imaging.frameRate = 1/metadata.Imaging.meanCCDTime;
metadata.Imaging.stimConditionOnsetTimeInFrames = ceil(stimConditionOnsetTimeCorr*metadata.Imaging.frameRate);

%Read Imaging Data For Each Cycle
time=tic;
displayOnlineMaps = true;
harmonicsToTest = [1 2];
numberOfTrials    = length(metadata.Imaging.stimConditionOnsetTimeInFrames);
freqData   = zeros(size(tifStack,1),size(tifStack,2),length(harmonicsToTest),numberOfTrials);
Win_startFrames = floor(OnSet*metadata.Imaging.frameRate);
Win_endFrames = floor(OffSet*metadata.Imaging.frameRate);
framesPerCycle    = mode(diff(metadata.Imaging.stimConditionOnsetTimeInFrames));

for currentCycle = 1:numberOfTrials
    selectedFrames = metadata.Imaging.stimConditionOnsetTimeInFrames(currentCycle)+Win_startFrames:...
        (metadata.Imaging.stimConditionOnsetTimeInFrames(currentCycle)+framesPerCycle+Win_startFrames+Win_endFrames-1);
    temp = zeros(size(tifStack,1),size(tifStack,2),length(selectedFrames));
    if(selectedFrames(end)>size(tifStack,3)),selectedFrames = selectedFrames(1):size(tifStack,3); end
    temp(:,:,1:length(selectedFrames)) = tifStack(:,:,selectedFrames);
    temp = double(temp);
    
    % Create Temporally Encoded Maps at the Desired Frequencies using DFT
    meanTrialFrame = mean(temp,3); % Used as a Cocktail Blank
    currentRelativeHarmonic = 0;
    for currentStimFreq = harmonicsToTest
        currentRelativeHarmonic = currentRelativeHarmonic+1;
        freqData(:,:,currentRelativeHarmonic,currentCycle) = freqData(:,:,currentRelativeHarmonic,currentCycle)...
            + sum((temp-repmat(meanTrialFrame,[1 1 size(temp,3)]))...
            .* repmat(reshape(exp(2*pi*1i*currentStimFreq*([1:size(temp,3)]/framesPerCycle)),[1 1 size(temp,3)]),[size(temp,1) size(temp,2) 1]),3);
    end
    
    % Display Online Maps
    if(displayOnlineMaps && currentCycle == 1)
        h=figure; drawnow;
    end
    if(displayOnlineMaps)
        figure(h);
        imagesc(((polarMapNew(mean(freqData(:,:,1,1:currentCycle),4))))); title('1st Order Harmonic'); axis image; colormap(hsv); colorbar;
        axis equal
        set(gca,'xtick',[],'ytick',[]);
        drawnow; pause(0.25);
    end
    
    disp(['  *Finished Reading and Constructing Imaging Data for Trial ' num2str(currentCycle) ' - Time Elapsed: ' num2str(toc(time)) ' seconds']);
end
analysis.maps.polar = polarMapNew(mean(squeeze(freqData(:,:,1,:)),3),cmap,2,[0 1]);
z = mean(squeeze(freqData(:,:,1,:)),3);
analysis.maps.phaseMap = angle(z)*180/pi;
analysis.maps.angleMap = wrapTo360(analysis.maps.phaseMap)./360;
clear temp
clear tifStack
clear freqData
saveas(gcf, fullfile(saveDirectory, 'PolarMap.png'))
%% make mask
figure
subplot(2,2,1)
imagesc(analysis.maps.polar)
title('Polar map'); colormap(jet); cbh=colorbar;
axis off

%make bv mask
Mask = expParam.baseImg;
Mean = nanmean(nanmean(Mask,2));
Std = std(std(Mask));
threshold = Mean + .5 * Std;
Mask(Mask<threshold) = nan;
Mask(~isnan(Mask)) = 1;
Mask(isnan(Mask)) = 0;


%make threshold mask
Mask2 = expParam.gaussMeanImg;
threshold = Mean - .9 * Std;
Mask2(Mask2<threshold) = nan;
Mask2(~isnan(Mask2)) = 1;
Mask2(isnan(Mask2)) = 0;

[labeledObject,nPolygons] = bwlabel(Mask2);
MaskSize = zeros(nPolygons,1);
for p = 1:nPolygons
    MaskSize(p) = sum(labeledObject(:)==p);
end
[~,LargestObjNr] = max(MaskSize);
Mask3 = zeros(size(Mask2,1), size(Mask2,2));
Mask3(labeledObject == LargestObjNr) = 1;
analysis.maps.mask = times(Mask3, Mask);

subplot(2,2,2)
analysis.maps.polarMask = times(analysis.maps.mask, analysis.maps.polar);
imagesc(analysis.maps.polarMask)
title('Polar map masked'); colormap(jet); cbh=colorbar;
axis off

%make phaseMap
subplot(2,2,3)
imagesc(analysis.maps.angleMap)
title('Angle'); colormap(jet); cbh=colorbar;
axis off

subplot(2,2,4)
analysis.maps.angleMapMask = times(analysis.maps.mask, analysis.maps.angleMap);
analysis.maps.angleMapMask(analysis.maps.angleMapMask == 0) = nan;
imagesc(analysis.maps.angleMapMask)
title('Angle masked'); colormap(jet); cbh=colorbar;
axis off
set(gcf,'color','w');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2)-4*100 6*100, 4*100]);
saveas(gcf, fullfile(saveDirectory, 'Mask&Angle.png'))

%% get retinotopic profile
if quantification
    rotate = 0;
    fig = figure;
    imagesc(analysis.maps.angleMapMask)
    hold all
    colormap(jet); cbh=colorbar; 

    %draw line over V3 & get profile
    [xV3, yV3]  = getline(fig);
    analysis.quant.slope = (yV3(2) - yV3(1))/(xV3(2)-xV3(1));
    
    if analysis.quant.slope > 1.2 || analysis.quant.slope < -1.2
        temp = imrotate(analysis.maps.angleMapMask,90);
        fig = figure;
        imagesc(temp)
        hold all
        colormap(jet); cbh=colorbar; axis off
        [xV3, yV3]  = getline(fig);
        analysis.quant.slope = (yV3(2) - yV3(1))/(xV3(2)-xV3(1));
        rotate = 1;
    end

    analysis.quant.interceptV3 = yV3(1) - analysis.quant.slope * xV3(1);
    analysis.quant.xV3= [1:1:size(analysis.maps.angleMapMask,2)];
    analysis.quant.yV3= analysis.quant.slope * analysis.quant.xV3 + analysis.quant.interceptV3;
    plot(analysis.quant.xV3, analysis.quant.yV3, 'w', 'LineWidth', 5)

    %make parallel line in V1
    [xP, yP]  = getpts(fig);
    plot(xP, yP, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    analysis.quant.interceptV1 = yP-analysis.quant.slope*xP;
    analysis.quant.xV1= [1:1:size(analysis.maps.angleMapMask,2)];
    analysis.quant.yV1= analysis.quant.slope * analysis.quant.xV1 + analysis.quant.interceptV1;
    plot(analysis.quant.xV1, analysis.quant.yV1, 'k', 'LineWidth', 5)
    saveas(gcf, fullfile(saveDirectory, 'Profile.png'))
    
    %rotate by 90 deg and get profile along that line
    if rotate == 0
        rotImg = imrotate(analysis.maps.angleMapMask,90);
    else
        rotImg = imrotate(analysis.maps.angleMapMask,180);
    end
    fig = figure;
    imagesc(rotImg)
    hold all
    xPNew = yP;
    yPNew = size(analysis.maps.angleMapMask,2) - xP;
    plot(xPNew, yPNew, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    
    interceptNew = yPNew-analysis.quant.slope*xPNew;
    xV1New= [1:1:size(analysis.maps.angleMapMask,1)];
    yV1New = analysis.quant.slope * xV1New + interceptNew;
    plot(xV1New, yV1New, 'k', 'LineWidth', 5)    

    %get profiles
    analysis.quant.V3Profile = improfile(analysis.maps.angleMapMask,[analysis.quant.xV3(1) analysis.quant.xV3(end)],[analysis.quant.yV3(1) analysis.quant.yV3(end)]);
    analysis.quant.V1Profile = improfile(analysis.maps.angleMapMask,[analysis.quant.xV1(1) analysis.quant.xV1(end)],[analysis.quant.yV1(1) analysis.quant.yV1(end)]);
    analysis.quant.AcrossAreasProfile = improfile(rotImg,[xV1New(1) xV1New(end)],[yV1New(1) yV1New(end)]);

    %change to deg and um
    metadata.StimParams.coverageDegree = abs(metadata.StimParams.startPoint-metadata.StimParams.endPoint);
    analysis.quant.V1ProfileDeg = analysis.quant.V1Profile * metadata.StimParams.coverageDegree + metadata.StimParams.startPoint;
    analysis.quant.V3ProfileDeg = analysis.quant.V3Profile * metadata.StimParams.coverageDegree + metadata.StimParams.startPoint;
    analysis.quant.AcrossAreasProfileDeg = analysis.quant.AcrossAreasProfile * metadata.StimParams.coverageDegree + metadata.StimParams.startPoint;
    umperpixel = 8300 / size(analysis.maps.angleMapMask,2);
    lengthProfileTotal = sqrt(size(analysis.maps.angleMapMask,2)^2 + (analysis.quant.yV3(1)-analysis.quant.yV3(end))^2);
    lengthUM = lengthProfileTotal*umperpixel;
    analysis.quant.lenghtProfile = linspace(0, lengthUM, size(analysis.maps.angleMapMask,2));
    analysis.quant.lengthProfileOtherDir = linspace(0, lengthUM, size(analysis.maps.angleMapMask,1));
    
    figure
    plot(analysis.quant.lenghtProfile, analysis.quant.V3ProfileDeg','b')
    hold on
    plot(analysis.quant.lenghtProfile, analysis.quant.V1ProfileDeg', 'g')
    try
        if metadata.StimParams.visualSpace == 'azimuth'
            title('Retinotopic position - Azimuth')
        elseif metadata.StimParams.visualSpace == 'elevation'
            title('Retinotopic position - Elevation') 
        end
    catch
        title('Retinotopic position')
    end
    ylabel('Position in Deg')
    xlabel('Position in \mum')
    legend({'V3', 'V1'}, 'Location', 'best')
    legend('boxoff')
    set(gcf,'color','w');
    saveas(gcf, fullfile(saveDirectory, 'Quantification.png'))
    
    %plot change across areas
    figure
    smoothAcross = smoothdata(analysis.quant.AcrossAreasProfileDeg', 'movmedian',15)';
    plot(analysis.quant.lengthProfileOtherDir,smoothAcross ,'r')
    hold on
    ylabel('Position in Deg')
    xlabel('Position in \mum')
    set(gcf,'color','w');
    %see if you can find reversal & correct for that
    if metadata.StimParams.visualSpace == 'azimuth'
        prompt = 'Do you want to look for 1) minima or 2) maxima?';
        choice = input(prompt);
        if choice == 1
            [~, LocalMinProm]  = islocalmin(smoothAcross);
        elseif choice == 2
            [~, LocalMinProm]  = islocalmax(smoothAcross);
        else
            error('invalid input');
        end
            
        [~, posMin] = max(LocalMinProm);
        plot(analysis.quant.lengthProfileOtherDir(posMin),smoothAcross(posMin), 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
        saveas(gcf, fullfile(saveDirectory, 'Quantification_acrossAreas.png'))
        
        zeroPointDeg = smoothAcross(posMin);
        disp(['shift by ' num2str(zeroPointDeg)])
        if choice == 1
            analysis.quant.V1ProfileDegCorrected = analysis.quant.V1ProfileDeg - zeroPointDeg;
            analysis.quant.V3ProfileDegCorrected = analysis.quant.V3ProfileDeg - zeroPointDeg;
            analysis.quant.AcrossAreasProfileDegCorrected = analysis.quant.AcrossAreasProfileDeg - zeroPointDeg;
            
        else
            analysis.quant.V1ProfileDegCorrected = analysis.quant.V1ProfileDeg + zeroPointDeg;
            analysis.quant.V3ProfileDegCorrected = analysis.quant.V3ProfileDeg + zeroPointDeg;
            analysis.quant.AcrossAreasProfileDegCorrected = analysis.quant.AcrossAreasProfileDeg + zeroPointDeg;
        end
               
        figure
        plot(analysis.quant.lenghtProfile, smoothdata(analysis.quant.V3ProfileDegCorrected', 'movmedian',15),'b')
        hold on
        plot(analysis.quant.lenghtProfile, smoothdata(analysis.quant.V1ProfileDegCorrected', 'movmedian',15), 'g')
        title('Retinotopic position - Azimuth')
        ylabel('Position in Deg')
        xlabel('Position in \mum')
        legend({'V3', 'V1'}, 'Location', 'best')
        legend('boxoff')
        set(gcf,'color','w');
        saveas(gcf, fullfile(saveDirectory, 'Quantification_Corrected.png'))
        
        figure
        smoothAcrossCorr = smoothdata(analysis.quant.AcrossAreasProfileDegCorrected', 'movmedian',15)';
        plot(analysis.quant.lengthProfileOtherDir,smoothAcrossCorr ,'r')
        hold on
        ylabel('Position in Deg')
        xlabel('Position in \mum')
        set(gcf,'color','w');
        saveas(gcf, fullfile(saveDirectory, 'Quantification_acrossAreasCorrected.png'))
    else
        saveas(gcf, fullfile(saveDirectory, 'Quantification_acrossAreas.png'))
    end
end

%% save file
save(fullfile(saveDirectory, 'continousEdge_ana.mat'), 'expParam', 'metadata','analysis');