function plotAllContoursSingleTrial(analysisParams, analysis, metadata)

for stim = 1:size(analysis.(analysisParams.field).zScore,3)-1
    figure('units','normalized','outerposition',[0 0 1 1])
    for trial=1:size(analysis.(analysisParams.field).zScore,4)
        subplot(analysisParams.nRows,analysisParams.nCol,analysisParams.plotNrTrials(trial));
        for i = 1:size(analysis.(analysisParams.field).ContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).ContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).ContourXData{stim,trial}{i}, analysis.(analysisParams.field).ContourYData{stim,trial}{i}, 'color', analysisParams.conC(trial,:))
                hold all
            end
        end
        set(gca,'ydir','reverse')
        axis off
        axis equal
    
        subplot(analysisParams.nRows,analysisParams.nCol,analysisParams.plotNrAll)
        for i = 1:size(analysis.(analysisParams.field).ContourXData{stim,trial},2)
            if size(analysis.(analysisParams.field).ContourXData{stim,trial}{i},1) > 5
                plot(analysis.(analysisParams.field).ContourXData{stim,trial}{i}, analysis.(analysisParams.field).ContourYData{stim,trial}{i}, 'color', analysisParams.conC(trial,:))
                hold all
            end
        end
        set(gca,'ydir','reverse')
        axis off
        axis equal
    end
    subplot(analysisParams.nRows, analysisParams.nCol, analysisParams.plotNrAverage)
    imagesc(analysis.(analysisParams.field).trialAveragedMaps(:,:,stim))
    hold on
    colormap('gray')
    axis off
    axis equal
    title('trial-averaged map')
    
    subplot(analysisParams.nRows,analysisParams.nCol, analysisParams.plotNrQuant(1))
    allOverlap = [analysis.(analysisParams.field).perOverlapA19(stim,:)'; analysis.(analysisParams.field).perOverlapA19Rnd(stim,:)'; analysis.(analysisParams.field).perOverlapV1(stim,:)'; analysis.(analysisParams.field).perOverlapV1Rnd(stim,:)'];
    boxHelp = [zeros(length(analysis.(analysisParams.field).perOverlapA19(stim,:)), 1); ones(length(analysis.(analysisParams.field).perOverlapA19Rnd(stim,:)), 1); 2*ones(length(analysis.(analysisParams.field).perOverlapV1(stim,:)), 1); 3*ones(length(analysis.(analysisParams.field).perOverlapV1Rnd(stim,:)), 1)];
    boxplot(allOverlap, boxHelp, 'Labels',{'A19','A19 Random', 'V1', 'V1 Random'})
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),analysisParams.cocA19(3,:),'FaceAlpha',.75)
    patch(get(h(2),'XData'),get(h(2),'YData'),analysisParams.cocA19(5,:),'FaceAlpha',.75)
    patch(get(h(3),'XData'),get(h(3),'YData'),analysisParams.cocV1(3,:),'FaceAlpha',.75)
    patch(get(h(4),'XData'),get(h(4),'YData'),analysisParams.cocV1(5,:),'FaceAlpha',.75)
    box off
    ylabel('percentage overlap')
    
    subplot(analysisParams.nRows,analysisParams.nCol, analysisParams.plotNrQuant(2))
    allCorr = [analysis.(analysisParams.field).corrMapsA19(stim,:)'; analysis.(analysisParams.field).corrMapsA19Rnd(stim,:)'; analysis.(analysisParams.field).corrMapsV1(stim,:)'; analysis.(analysisParams.field).corrMapsV1Rnd(stim,:)'];
    boxHelp = [zeros(length(analysis.(analysisParams.field).corrMapsA19(stim,:)), 1); ones(length(analysis.(analysisParams.field).corrMapsA19Rnd(stim,:)), 1); 2*ones(length(analysis.(analysisParams.field).corrMapsV1(stim,:)), 1); 3*ones(length(analysis.(analysisParams.field).corrMapsV1Rnd(stim,:)), 1)];
    boxplot(allCorr, boxHelp, 'Labels',{'A19','A19 Random', 'V1', 'V1 Random'})
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),analysisParams.cocA19(3,:),'FaceAlpha',.75)
    patch(get(h(2),'XData'),get(h(2),'YData'),analysisParams.cocA19(5,:),'FaceAlpha',.75)
    patch(get(h(3),'XData'),get(h(3),'YData'),analysisParams.cocV1(3,:),'FaceAlpha',.75)
    patch(get(h(4),'XData'),get(h(4),'YData'),analysisParams.cocV1(5,:),'FaceAlpha',.75)
    box off
    ylabel('correlation')
    
    subplot(analysisParams.nRows,analysisParams.nCol, analysisParams.plotNrQuant(3))
    allConF = [analysis.(analysisParams.field).convFactA19(stim,:)'; analysis.(analysisParams.field).convFactV1(stim,:)'];
    boxHelp = [zeros(length(analysis.(analysisParams.field).convFactA19(stim,:)), 1); ones(length(analysis.(analysisParams.field).convFactV1(stim,:)), 1)];
    boxplot(allConF, boxHelp, 'Labels',{'A19','V1'})
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),analysisParams.cocA19(3,:),'FaceAlpha',.75)
    patch(get(h(2),'XData'),get(h(2),'YData'),analysisParams.cocV1(3,:),'FaceAlpha',.75)
    box off
    ylabel('Convergence Factor')
    
    set(gcf, 'color', 'w');
    supertitle([num2str(metadata.StimParams.directions(stim)) ' deg'])
    saveas(gcf, fullfile(analysisParams.saveDirectory, ['SingleTrialContourOnly_Stim' num2str(stim) '.png']))
    close gcf
end
