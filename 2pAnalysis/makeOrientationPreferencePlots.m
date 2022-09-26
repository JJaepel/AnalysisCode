function makeOrientationPreferencePlots(analysisParams, analysis, data, metadata, types, saveDirectory)

if types == 1
    rois_ori = linspace(1,length(analysis.(analysisParams.field).roi),length(analysis.(analysisParams.field).roi));
    rois_dir = rois_ori;
    rois_con = rois_ori;
    number = 1;
elseif types == 2
    rois_ori = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    rois_dir = rois_ori;
    rois_con = rois_ori;
    number = 2;
elseif types == 3
    rois_ori = find([analysis.(analysisParams.field).roi.OSIFit] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    rois_dir = find([analysis.(analysisParams.field).roi.DSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    if analysisParams.stimType == 2
        rois_con = find([analysis.(analysisParams.field).roi.SFSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    elseif analysisParams.stimType == 3
        rois_con = find([analysis.(analysisParams.field).roi.TFSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    end
    number = 3;
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

h=figure(101+number);
set(h, 'units','normalized','outerposition',[0 0 1 1]);
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
histogram(alloriprefs,linspace(0,180,metadata.StimParams.numOrientations +1), 'FaceColor', analysisParams.coc_prop(1,:), 'EdgeColor', analysisParams.coc_prop(2,:));
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
histogram(alldirprefs,linspace(0,360,2*metadata.StimParams.numOrientations +1), 'FaceColor', analysisParams.coc_prop(3,:), 'EdgeColor', analysisParams.coc_prop(4,:));
ylabel('Cells');
xlabel(sprintf('Direction preference (%s)',char(145)));
xlim([-22.5 (360+22.5)])
axis square;
set(gca,'Box','off');

if analysisParams.stimType == 2
    subplot(2,3,3)
    PlotPrefOnTemplateOri(analysis, data, 3, analysisParams.field,data.template, rois_con)
    subplot(2,3,6)
    sf_counts = histcounts(allConpref, metadata.StimParams.numSf);
    bar(metadata.StimParams.spatialFreq, sf_counts, 'FaceColor', analysisParams.coc_prop(5,:), 'EdgeColor', analysisParams.coc_prop(6,:));
    title('Histogram')
    ylabel('Cells');
    xlabel(sprintf('Spatial frequency preference'));
    axis square;
    set(gca,'Box','off');
elseif analysisParams.stimType == 3
    subplot(2,3,3)
    PlotPrefOnTemplateOri(analysis, data, 4, analysisParams.field,data.template, rois_con)
    subplot(2,3,6)
    tf_counts = histcounts(allConpref, metadata.StimParams.numTf);
    bar(metadata.StimParams.temporalFreq, tf_counts, 'FaceColor', analysisParams.coc_prop(5,:), 'EdgeColor', analysisParams.coc_prop(6,:));
    title('Histogram')
    ylabel('Cells');
    xlabel(sprintf('Temporal frequency preference'));
    axis square;
    set(gca,'Box','off');
end    

set(gcf, 'color', 'w');
if types == 1
    saveas(gcf, fullfile(saveDirectory, '102_Overlaymaps_all_cells.png'))
elseif types == 2
    saveas(gcf, fullfile(saveDirectory, '103_Overlaymaps_resp_cells.png'))
elseif types == 3
    saveas(gcf, fullfile(saveDirectory, '104_Overlaymaps_selective_cells.png'))
end
%close all