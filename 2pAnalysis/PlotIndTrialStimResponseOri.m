function PlotIndTrialStimResponseOri(metadata, analysis, field, roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    cocTr= cbrewer('qual', 'Set3', metadata.StimParams.numTrials);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    if (ymax-ymin) > 3
        refBarY = 1;
    else
        refBarY = 0.1;
    end
    for i=1:metadata.StimParams.numDirections
        for trial = 1:metadata.StimParams.numTrials
            id = metadata.StimParams.numDirections*(trial-1)+i;
            ax=subplot(metadata.StimParams.numTrials, metadata.StimParams.numDirections,id);
            cla(ax)
            y=medfilt1(analysis.(field).roi(roi).stimResponseTrace(i,trial,:),3);
            x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
            if stimWindow(1)~=0
                bl=x(stimWindow);
                patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
            end
            plot(ax,x, squeeze(y), 'Color', cocTr(trial,:), 'LineWidth', 3);
            hold all
            ylim([ymin, ymax])
            xlim([min(x) max(x)])
            title(metadata.StimParams.directions(i))
            axis off
            if id == 1
                %add 0.1 F/F0 as reference
                patch([min(x) min(x)], [ymin ymin+refBarY],[1 1 1], 'LineWidth', 3)
                %add 1 s as reference
                patch([min(x) min(x)+1], [ymin ymin], [1 1 1], 'LineWidth',3);
            end
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end