function analysis = makeContourMaps(metadata, analysis, analysisParams)

saveDirectory = analysisParams.saveDirectory;
clippingPercentile = 0.2;
mapsWithContours = 0;

includedFrames = [round(metadata.Imaging.rate * metadata.StimParams.isi/2)+1:round(metadata.Imaging.rate * metadata.StimParams.isi/2)+ceil(metadata.Imaging.rate * metadata.StimParams.stimDuration)];
stimResponseTrace = permute(analysis.(analysisParams.field).roi.stimResponseTrace, [4 5 3 2 1]);
indRespMaps = mean(stimResponseTrace(:,:,:,:,includedFrames),5);
for stim=1:size(indRespMaps,3)
    for trial = 1:size(indRespMaps,4)
        temp = indRespMaps(:,:,stim, trial);
        temp(~analysis.maskBV(:)) = NaN;
        indRespMaps(:,:,stim,trial)=temp;
    end
end

clipValue = prctile(indRespMaps(:),[clippingPercentile 100-clippingPercentile]); 


%z scoring data
disp('Z scoring data')
mu = nanmean(indRespMaps(:));
sd = nanstd(indRespMaps(:));
analysis.(analysisParams.field).zScore = (indRespMaps - mu)/sd;

% figure('units','normalized','outerposition',[0 0 1 1])
% i = 1;
% for stim=1:size(indRespMaps,3)-1
%     for trial = 1:size(indRespMaps,4)
%         subplot(size(indRespMaps,3)-1, size(indRespMaps,4),i)
%         imagesc(analysis.(analysisParams.field).zScore(:,:,stim,trial))
%         if stim == 1
%            title(['Trial ' num2str(trial)]);
%         end
%         axis off
%         axis equal
%         if trial == 1
%             h = text(-20, size(indRespMaps,1),[num2str(metadata.StimParams.directions(stim)) ' deg']);
%             set(h, 'rotation', 90)
%         end
%         i = i+1;
%         colormap('gray')
%     end
% end
% set(gcf, 'color', 'w');
% saveas(gcf, fullfile(saveDirectory, 'ZscoredTrials.png'))

%making  contours on trial masks
contourMaps = zeros(size(indRespMaps,1), size(indRespMaps,2), size(indRespMaps,3),size(indRespMaps,4));
for stim = 1:size(indRespMaps,3)-1
   for trial=1:size(indRespMaps,4)
       [c, ~] = contour(analysis.(analysisParams.field).zScore(:,:,stim,trial));
       cdata = contourdata(c);
       levelData = find([cdata.level] ==0);
       analysis.(analysisParams.field).zeroContourXData{stim,trial} = {cdata(levelData).xdata};
       analysis.(analysisParams.field).zeroContourYData{stim,trial} = {cdata(levelData).ydata};
       counter = 1;
       for l = 1:size(analysis.(analysisParams.field).zeroContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).zeroContourXData{stim,trial}{l},1) > 5
                mask(:,:,counter) = poly2mask(analysis.(analysisParams.field).zeroContourXData{stim,trial}{l}, analysis.(analysisParams.field).zeroContourYData{stim,trial}{l}, size(indRespMaps,1), size(indRespMaps,2));
                counter = counter + 1;          
            end
       end
       temp = sum(mask,3);
       temp(temp > 1) = 1;
       contourMaps(:,:,stim, trial) = temp;
   end
end

%making contours of average masks
AvgRespMaps = zeros(size(indRespMaps,1), size(indRespMaps,2), size(indRespMaps,3));
contourAvgMaps = zeros(size(indRespMaps,1), size(indRespMaps,2), size(indRespMaps,3));
for stim = 1:size(analysis.(analysisParams.field).trialAveragedMaps,3)
    temp = analysis.(analysisParams.field).trialAveragedMaps(:,:,stim);
    temp(~analysis.maskBV(:)) = NaN;
    AvgRespMaps(:,:,stim) = temp;
    [c, ~] = contour(temp);
    cdata = contourdata(c);
    levelData = find([cdata.level] ==0);
    analysis.(analysisParams.field).zeroContourAvgXData{stim} = {cdata(levelData).xdata};
    analysis.(analysisParams.field).zeroContourAvgYData{stim} = {cdata(levelData).ydata};
    counter = 1;
    for l = 1:size(analysis.(analysisParams.field).zeroContourAvgXData{stim},2)
        if size(analysis.(analysisParams.field).zeroContourAvgXData{stim}{l},1) > 5
            mask(:,:,counter) = poly2mask(analysis.(analysisParams.field).zeroContourAvgXData{stim}{l}, analysis.(analysisParams.field).zeroContourAvgYData{stim}{l}, size(indRespMaps,1), size(indRespMaps,2));
            counter = counter + 1;          
        end
    end
    temp = sum(mask,3);
    temp(temp > 1) = 1;
    contourAvgMaps(:,:,stim) = temp;
    close gcf
end

%calculating overlap
analysis.overlap = zeros(size(indRespMaps,3), size(indRespMaps, 4)-1);
analysis.overlapA19 = zeros(size(indRespMaps,3), size(indRespMaps, 4)-1);
analysis.overlapV1 = zeros(size(indRespMaps,3), size(indRespMaps, 4)-1);
for stim = 1:size(indRespMaps,3)-1
   for trial=1:size(indRespMaps,4)
       overlapMask = contourMaps(:,:,stim, trial) & contourAvgMaps(:,:,stim);
       analysis.overlap(stim, trial) = nnz(overlapMask)/nnz(contourAvgMaps(:,:,stim));
       analysis.overlapA19(stim,trial) = nnz(overlapMask & analysis.maskA19)/(nnz(contourAvgMaps(:,:,stim)&analysis.maskA19));
       analysis.overlapV1(stim,trial) = nnz(overlapMask & analysis.maskV1)/(nnz(contourAvgMaps(:,:,stim)&analysis.maskV1));
   end
end

%put a mask over Ori and dir mask
OriMap = analysis.(analysisParams.field).orientationMap;
OriMap(~analysis.maskBV(:)) = NaN;
DirMap = analysis.(analysisParams.field).directionMap;
DirMap(~analysis.maskBV(:)) = NaN;

conC = cbrewer('qual', 'Paired',metadata.StimParams.numTrials);
cocV1 = cbrewer('seq', 'RdPu',5);
cocA19 = cbrewer('seq', 'PuBuGn',5);
if metadata.StimParams.numTrials > 6
    if metadata.StimParams.numTrials > 10
        nRows = 5; nCol = 5;
        plotNrTrials = [1 2 3 6 7 8 11 12 13 16 17 18 21 22 23];
        plotNrAll= [4:5 9:10];
        plotNrAverage = [14 19]; 
        plotNrQuant = [15 20]; 
    else
        nRows = 4; nCol = 5;
        plotNrTrials = [1 2 3 6 7 8 11 12 13 16 17 18];
        plotNrAll = [4:5 9:10];
        plotNrAverage = [14 19]; 
        plotNrQuant = [15 20]; 
    end
else
    nRows = 4; nCol = 4;
    plotNrTrials = [1 2 5 6 9 10];
    plotNrAll = [3:4 7:8];
    plotNrAverage = [11 15];
    plotNrQuant = [12 16];
end

disp('Drawing countours for individual trials')
if mapsWithContours
    for type = 1:2
        for stim = 1:size(indRespMaps,3)-1
            cRGBmap = zeros(size(indRespMaps,1),size(indRespMaps,2),3);
            figure('units','normalized','outerposition',[0 0 1 1])
            for trial=1:size(indRespMaps,4)
                if type == 1
                    %plot trial z-scored map
                    subplot(nRows,nCol,plotNrTrials(trial))
                    zScoreMapTrial = analysis.(analysisParams.field).zScore(:,:,stim,trial);
                    zScoreMapTrial(isnan(zScoreMapTrial)) = 1;
                    imagesc(zScoreMapTrial);
                    axis off
                    axis equal
                    colormap('gray')
                else
                    %plot trial z-scored map with contour on top
                    subplot(nRows,nCol,plotNrTrials(trial))
                    zScoreMapTrial = analysis.(analysisParams.field).zScore(:,:,stim,trial);
                    zScoreMapTrial(isnan(zScoreMapTrial)) = 1;
                    for c=1:3
                        cRGBmap(:,:,c) = zScoreMapTrial;
                    end
                    imagesc(cRGBmap)
                    hold on
                    for i = 1:size(analysis.(analysisParams.field).zeroContourXData{stim,trial},2)
                        if size(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i},1) > 5
                            plot(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i}, analysis.(analysisParams.field).zeroContourYData{stim,trial}{i}, 'color', conC(trial,:))
                            hold all
                        end
                    end
                    axis off
                    axis equal
                end
            end

            %plot average map for this stim
            subplot(nRows, nCol, plotNrAverage)
            imagesc(AvgRespMaps(:,:,stim))
            hold on
            colormap('gray')
            title('trial-averaged map')

            if type == 2
                for i = 1:size(analysis.(analysisParams.field).zeroContourAvgXData{stim},2)
                    if size(analysis.(analysisParams.field).zeroContourAvgXData{stim}{i},1) > 5
                        plot(analysis.(analysisParams.field).zeroContourAvgXData{stim}{i}, analysis.(analysisParams.field).zeroContourAvgYData{stim}{i}, 'color', [1,1,1])
                        hold all
                    end
                    axis off
                    axis equal
                end

                subplot(nRows, nCol, plotNrAll)
                for trial=1:size(indRespMaps,4)
                    for i = 1:size(analysis.(analysisParams.field).zeroContourXData{stim,trial},2)
                        if size(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i},1) > 5
                            plot(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i}, analysis.(analysisParams.field).zeroContourYData{stim,trial}{i}, 'color', conC(trial,:))
                            hold all
                        end
                    end
                end
                title('overlayed contours')
                set(gca,'ydir','reverse')
            end
            axis off
            axis equal

            set(gcf, 'color', 'w');
            supertitle([num2str(metadata.StimParams.directions(stim)) ' deg'])
            if type == 1
                saveas(gcf, fullfile(saveDirectory, ['SingleTrialZScore_Stim' num2str(stim) '.png']))
            else
                saveas(gcf, fullfile(saveDirectory, ['SingleTrialContour_Stim' num2str(stim) '.png']))
            end
            close gcf
        end
    end
end

disp('Drawing zero contour lines for individual trials')

for stim = 1:size(indRespMaps,3)-1
    figure('units','normalized','outerposition',[0 0 1 1])
    for trial=1:size(indRespMaps,4)
        h = subplot(nRows,nCol,plotNrTrials(trial));
        for i = 1:size(analysis.(analysisParams.field).zeroContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i}, analysis.(analysisParams.field).zeroContourYData{stim,trial}{i}, 'color', conC(trial,:))
                hold all
            end
        end
        set(gca,'ydir','reverse')
        axis off
        axis equal
    
        subplot(nRows,nCol,plotNrAll)
        for i = 1:size(analysis.(analysisParams.field).zeroContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i}, analysis.(analysisParams.field).zeroContourYData{stim,trial}{i}, 'color', conC(trial,:))
                hold all
            end
        end
        set(gca,'ydir','reverse')
        axis off
        axis equal
    end
    subplot(nRows, nCol, plotNrAverage)
    imagesc(AvgRespMaps(:,:,stim))
    hold on
    colormap('gray')
    axis off
    axis equal
    title('trial-averaged map')
    
    subplot(nRows,nCol, plotNrQuant)
    boxplot([analysis.overlapA19(stim,:); analysis.overlapV1(stim,:)]', 'Labels',{'A19','V1'})
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),cocA19(4,:),'FaceAlpha',.5);
    patch(get(h(2),'XData'),get(h(2),'YData'),cocV1(4,:),'FaceAlpha',.5);
    box off
    ylabel('percentage overlap')
    
    set(gcf, 'color', 'w');
    supertitle([num2str(metadata.StimParams.directions(stim)) ' deg'])
    saveas(gcf, fullfile(saveDirectory, ['SingleTrialContourOnly_Stim' num2str(stim) '.png']))
    close gcf
end

if length(metadata.StimParams.directions) > 8
    nRows = 4;
    nCol = 5;
    plotNrTrials = [1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19];
    plotNrOriMap = 5;
    plotNrDirMap = 10;
    plotNrQuant = 15;
else
    nRows = 3;
    nCol = 4;
    plotNrTrials = [1 2 3 5 6 7 9 10 11];
    plotNrOriMap = 4;
    plotNrDirMap = 8;
    plotNrQuant = 12;
end

figure('units','normalized','outerposition',[0 0 1 1])
for stim = 1:length(metadata.StimParams.directions)
    subplot(nRows, nCol,plotNrTrials(stim))
    for trial=1:size(indRespMaps,4)
        for i = 1:size(analysis.(analysisParams.field).zeroContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).zeroContourXData{stim,trial}{i}, analysis.(analysisParams.field).zeroContourYData{stim,trial}{i}, 'color', conC(trial,:))
                hold all
            end
        end
    end
    set(gca,'ydir','reverse')
    axis off
    axis equal
    title([num2str(metadata.StimParams.directions(stim)) ' deg'])
end
subplot(nRows,nCol, plotNrOriMap)
cmap=hsv;
imagesc(polarMapEpi(OriMap, [0 clipValue(2)]));
colormap(cmap), cbh=colorbar; axis image ;
axis off
axis equal
title('Orientation map')

subplot(nRows,nCol, plotNrDirMap)
imagesc(polarMapEpi(DirMap, [0 clipValue(2)]));
colormap(cmap), cbh=colorbar; axis image ;
axis off
axis equal
title('Direction map')

subplot(nRows,nCol, plotNrQuant)
boxplot([analysis.overlapA19(:), analysis.overlapV1(:)], 'Labels',{'A19','V1'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),cocA19(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocV1(4,:),'FaceAlpha',.5)
box off
ylabel('percentage overlap')
set(gcf, 'color', 'w');

analysis.contourMaps = contourMaps;
analysis.contourAvgMaps = contourAvgMaps;
saveas(gcf, fullfile(saveDirectory, ['SingleTrialCorrelationAllStims' num2str(stim) '.png']))