function PlotAvgStimResponseOri(metadata, analysis, field, roi)
    %PlotAvgStimResponse Plots Avg + SEM response for stimuli
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    numCon = metadata.StimParams.numCon;
    cocCon = cbrewer('qual', 'Set1', numCon);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=nanmax(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=nanmin(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));
    
    if (ymax-ymin) > 7
        refBarY = 5;
    elseif (ymax-ymin) > 3
        refBarY = 1;
    elseif (ymax-ymin) > 0.5
        refBarY = 0.1;
    end
    
    for i=1:metadata.StimParams.numDirections 
        ax=subplot(1,metadata.StimParams.numDirections , i);
        cla(ax)
        y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(i,:),3);
        x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
        if stimWindow(1)~=0
            bl=x(stimWindow);
            patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
        end
        for tf = 1:numCon
            id = metadata.StimParams.numDirections*(tf-1)+i;
            err=analysis.(field).roi(roi).SEMResponseTrace(id,:);
            hold all
            patch([x fliplr(x)],[y+err fliplr(y-err)],[.5 .5 .5], 'LineStyle', 'none');
            plot(ax,x, smooth(squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:))), 'Color', cocCon(tf,:));
            hold all
            if id == 1
                %add 0.1 F/F0 as reference
                patch([min(x) min(x)], [ymin ymin+refBarY],[1 1 1], 'LineWidth', 3)
                %add 1 s as reference
                patch([min(x) min(x)+1], [ymin ymin], [1 1 1], 'LineWidth',3);
                text(min(x)-1,ymin, num2str(refBarY))            
            end
        end
        ylim([ymin, ymax])
        xlim([min(x) max(x)])
        title(metadata.StimParams.directions(i))
        axis off
    end
    hold off
    set(gcf, 'Color', 'w')
end