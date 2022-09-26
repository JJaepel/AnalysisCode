function plotSingleTrialContours(analysisParams, analysis, metadata, type)

cRGBmap = zeros(size(analysis.(analysisParams.field).zScore,1),size(analysis.(analysisParams.field).zScore,2),3);
figure('units','normalized','outerposition',[0 0 1 1])
for trial=1:size(analysis.(analysisParams.field).zScore,4)
    if type == 1
        %plot trial z-scored map
        subplot(analysisParams.nRows,analysisParams.nCol,analysisParams.plotNrTrials(trial))
        zScoreMapTrial = analysis.(analysisParams.field).zScore(:,:,stim,trial);
        zScoreMapTrial(isnan(zScoreMapTrial)) = 1;
        imagesc(zScoreMapTrial);
        axis off
        axis equal
        colormap('gray')
    else
        %plot trial z-scored map with contour on top
        subplot(analysisParams.nRows,analysisParams.nCol,analysisParams.plotNrTrials(trial))
        zScoreMapTrial = analysis.(analysisParams.field).zScore(:,:,stim,trial);
        zScoreMapTrial(isnan(zScoreMapTrial)) = 1;
        for c=1:3
            cRGBmap(:,:,c) = zScoreMapTrial;
        end
        imagesc(cRGBmap)
        hold on
        for i = 1:size(analysis.(analysisParams.field).ContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).ContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).ContourXData{stim,trial}{i}, analysis.(analysisParams.field).ContourYData{stim,trial}{i}, 'color', conC(trial,:))
                hold all
            end
        end
        axis off
        axis equal
    end
end

%plot average map for this stim
subplot(analysisParams.nRows, analysisParams.nCol, analysisParams.plotNrAverage)
imagesc(analysis.(analysisParams.field).AvgRespMaps(:,:,stim))
hold on
colormap('gray')
title('trial-averaged map')

if type == 2
    for i = 1:size(analysis.(analysisParams.field).ContourXData{stim},2)
        if size(analysis.(analysisParams.field).ContourXData{stim}{i},1) > 5
            plot(analysis.(analysisParams.field).ContourXData{stim}{i}, analysis.(analysisParams.field).ContourYData{stim}{i}, 'color', [0,0,0])
            hold all
        end
        axis off
        axis equal
    end

    subplot(analysisParams.nRows, analysisParams.nCol, analysisParams.plotNrAll)
    for trial=1:size(analysis.(analysisParams.field).zScore,4)
        for i = 1:size(analysis.(analysisParams.field).ContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).ContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).ContourXData{stim,trial}{i}, analysis.(analysisParams.field).ContourYData{stim,trial}{i}, 'color', conC(trial,:))
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
    saveas(gcf, fullfile(analysisParams.saveDirectory, ['SingleTrialZScore_Stim' num2str(stim) '.png']))
else
    saveas(gcf, fullfile(analysisParams.saveDirectory, ['SingleTrialContour_Stim' num2str(stim) '.png']))
end
close gcf