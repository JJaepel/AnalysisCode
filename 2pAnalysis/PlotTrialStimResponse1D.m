function PlotTrialStimResponse1D(metadata, analysis, field, roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);    
    colorlevelsTr = metadata.StimParams.numTrials;
    cocTr = cbrewer('div', 'Spectral', colorlevelsTr);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).stimResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).stimResponseTrace(:));

    if (ymax-ymin) > 3
        refBarY = 1;
    else
        refBarY = 0.1;
    end
    
    for p = 1:metadata.StimParams.numPatches
        ax=subplot(1,metadata.StimParams.numPatches, p);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(p,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
                bl=x(stimWindow);
                patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        hold all
        for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
            plot(ax,x, smooth(squeeze(analysis.(field).roi(roi).stimResponseTrace(p,trial,:))), 'Color', cocTr(trial,:));
            hold all
        end
        plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(p,:)), 'k');
        hold all
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        axis off
        if p == 1
            %add 0.1 F/F0 as reference
            patch([min(x) min(x)], [ymin ymin+refBarY],[1 1 1], 'LineWidth', 3)
            %add 1 s as reference
            patch([min(x) min(x)+1], [ymin ymin], [1 1 1], 'LineWidth',3);
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end