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
respData = indRespMaps(:,:,1:end-1);
mu = nanmean(respData(:));
sd = nanstd(respData(:));
analysis.(analysisParams.field).zScore = (indRespMaps - mu)/sd;

%making  contours on trial masks
contourMaps = zeros(size(indRespMaps,1), size(indRespMaps,2), size(indRespMaps,3),size(indRespMaps,4));
for stim = 1:size(indRespMaps,3)-1
   for trial=1:size(indRespMaps,4)
        [c, ~] = contour(analysis.(analysisParams.field).zScore(:,:,stim,trial));
        cdata = contourdata(c);
        levelData = find([cdata.level] ==0);
        analysis.(analysisParams.field).zeroContourXData{stim,trial} = {cdata(levelData).xdata};
        analysis.(analysisParams.field).zeroContourYData{stim,trial} = {cdata(levelData).ydata};
        temp = analysis.(analysisParams.field).zScore(:,:,stim,trial);
        temp(temp>0) = 1;
        temp(temp<0) = 0;
        temp(isnan(temp)) = 0;
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
    temp = zeros(size(indRespMaps,1), size(indRespMaps,2));
    for i = 1:size(analysis.(analysisParams.field).zeroContourAvgXData{stim},2)
        if size(analysis.(analysisParams.field).zeroContourAvgXData{stim}{i},1) > 5
            xValues = round(analysis.(analysisParams.field).zeroContourAvgXData{stim}{i});
            yValues = round(analysis.(analysisParams.field).zeroContourAvgYData{stim}{i});
            for dots = 1:length(xValues)
                temp(yValues(dots), xValues(dots)) = 1;
                temp(yValues(dots)-1, xValues(dots)) = 1;
                temp(yValues(dots)+1, xValues(dots)) = 1;
                temp(yValues(dots), xValues(dots)-1) = 1;
                temp(yValues(dots), xValues(dots)+1) = 1;
                temp(yValues(dots)-1, xValues(dots)-1) = 1;
                temp(yValues(dots)+1, xValues(dots)-1) = 1;
                temp(yValues(dots)-1, xValues(dots)-1) = 1;
                temp(yValues(dots)-1, xValues(dots)+1) = 1;
                temp(yValues(dots)-1, xValues(dots)+1) = 1;
                temp(yValues(dots)+1, xValues(dots)+1) = 1;
                temp(yValues(dots)+1, xValues(dots)-1) = 1;
                temp(yValues(dots)+1, xValues(dots)+1) = 1;
            end
        end
    end
    contourAvgMaps(:,:,stim) = temp;
end

%calculating overlap
analysis.overlap = zeros(size(indRespMaps,3), size(indRespMaps, 4)-1);
analysis.overlapA19 = zeros(size(indRespMaps,3), size(indRespMaps, 4)-1);
analysis.overlapV1 = zeros(size(indRespMaps,3), size(indRespMaps, 4)-1);
for stim = 1:size(indRespMaps,3)-1
   for trial=1:size(indRespMaps,4)
       overlapMask = contourMaps(:,:,stim, trial) & contourMaps(:,:,stim, 1);
       analysis.overlap(stim, trial) = nnz(overlapMask)/nnz(contourMaps(:,:,stim, 1));
       analysis.overlapA19(stim,trial) = nnz(overlapMask & analysis.maskA19)/(nnz(contourMaps(:,:,stim, 1)& analysis.maskA19));
       analysis.overlapV1(stim,trial) = nnz(overlapMask & analysis.maskV1)/(nnz(contourMaps(:,:,stim, 1)& analysis.maskV1));
       analysis.corr(stim, trial) = corr2(overlapMask,contourMaps(:,:,stim, 1));
       analysis.corrA19(stim,trial)= corr2(contourMaps(:,:,stim, trial) & analysis.maskA19,contourMaps(:,:,stim, 1) & analysis.maskA19);
       analysis.corrV1(stim,trial)= corr2(contourMaps(:,:,stim, trial) & analysis.maskV1,contourMaps(:,:,stim, 1) & analysis.maskV1); 
   end
end

%calculating chance overlap
analysis.overlapRandom = zeros(size(indRespMaps,3),1000);
analysis.overlapRandom = zeros(size(indRespMaps,3),1000);
analysis.overlapRandom = zeros(size(indRespMaps,3),1000);
for stim = 1:size(indRespMaps,3)-1
   for drawing = 1:1000
       randomStim = randi([1 size(indRespMaps,3)-1]);
       randomTrial= randi([1 size(indRespMaps,4)]);
       overlapMask = contourMaps(:,:,randomStim, randomTrial) & contourMaps(:,:,stim, 1);
       OverlapRandom(drawing) = nnz(overlapMask)/nnz(contourMaps(:,:,stim, 1));
       OverlapA19Random(drawing)= nnz(overlapMask & analysis.maskA19)/(nnz(contourMaps(:,:,stim, 1)&analysis.maskA19));
       OverlapV1Random(drawing)= nnz(overlapMask & analysis.maskV1)/(nnz(contourMaps(:,:,stim, 1)&analysis.maskV1));
       CorrA19Random(drawing)= corr2(contourMaps(:,:,randomStim, randomTrial) & analysis.maskA19,contourMaps(:,:,stim, 1)&analysis.maskA19);
       CorrV1Random(drawing)= corr2(contourMaps(:,:,randomStim, randomTrial) & analysis.maskV1,contourMaps(:,:,stim, 1)&analysis.maskV1);

   end
   analysis.overlapRandom(stim,:) =  OverlapRandom;
   analysis.overlapA19Random(stim,:) = OverlapA19Random;
   analysis.overlapV1Random(stim,:) = OverlapV1Random;
   analysis.corrA19Random(stim,:) = CorrA19Random;
   analysis.corrlapV1Random(stim,:) = CorrV1Random;
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
                        plot(analysis.(analysisParams.field).zeroContourAvgXData{stim}{i}, analysis.(analysisParams.field).zeroContourAvgYData{stim}{i}, 'color', [0,0,0])
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
    allOverlap = [analysis.overlapA19(stim,:)'; analysis.overlapA19Random(stim,:)'; analysis.overlapV1(stim,:)'; analysis.overlapV1Random(stim,:)'];
    boxHelp = [zeros(length(analysis.overlapA19(stim,:)), 1); ones(length(analysis.overlapA19Random(stim,:)), 1); 2*ones(length(analysis.overlapV1(stim,:)), 1); 3*ones(length(analysis.overlapV1Random(stim,:)), 1)];
    boxplot(allOverlap, boxHelp, 'Labels',{'A19','A19 Random' 'V1', 'V1 Random'})
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),cocA19(3,:),'FaceAlpha',.75)
    patch(get(h(2),'XData'),get(h(2),'YData'),cocA19(5,:),'FaceAlpha',.75)
    patch(get(h(3),'XData'),get(h(3),'YData'),cocV1(3,:),'FaceAlpha',.75)
    patch(get(h(4),'XData'),get(h(4),'YData'),cocV1(5,:),'FaceAlpha',.75)
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

allOverlap = [analysis.overlapA19(:); analysis.overlapA19Random(:); analysis.overlapV1(:); analysis.overlapV1Random(:)];
boxHelp = [zeros(length(analysis.overlapA19(:)), 1); ones(length(analysis.overlapA19Random(:)), 1); 2*ones(length(analysis.overlapV1(:)), 1); 3*ones(length(analysis.overlapV1Random(:)), 1)];
boxplot(allOverlap, boxHelp, 'Labels',{'A19','A19 Random' 'V1', 'V1 Random'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),cocA19(3,:),'FaceAlpha',.75)
patch(get(h(2),'XData'),get(h(2),'YData'),cocA19(5,:),'FaceAlpha',.75)
patch(get(h(3),'XData'),get(h(3),'YData'),cocV1(3,:),'FaceAlpha',.75)
patch(get(h(4),'XData'),get(h(4),'YData'),cocV1(5,:),'FaceAlpha',.75)
box off
ylabel('percentage overlap')
set(gcf, 'color', 'w');

analysis.contourMaps = contourMaps;
analysis.contourAvgMaps = contourAvgMaps;
saveas(gcf, fullfile(saveDirectory, ['SingleTrialCorrelationAllStims' num2str(stim) '.png']))