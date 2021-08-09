function PlotAvgStimResponse1D(metadata, analysis, field, roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    cocPatch = cbrewer('qual', 'Set1', metadata.StimParams.numPatches);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));
    
    if (ymax-ymin) > 3
        refBarY = 1;
    else
        refBarY = 0.1;
    end

    for p=1:metadata.StimParams.numPatches
        ax=subplot(1,metadata.StimParams.numPatches, p);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(p,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
            bl=x(stimWindow);
            %add stimWindow in pink
            patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        err=analysis.(field).roi(roi).SEMResponseTrace(p,:);
        hold on
        %add stimTrace
        patch([x fliplr(x)],[y+err fliplr(y-err)],[.5 .5 .5], 'LineStyle', 'none');
        %add avgResponseTrace 
        plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(p,:)), 'Color', cocPatch(p,:), 'LineWidth', 3);
        hold off
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