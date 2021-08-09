function PlotStimResponse2D(metadata, analysis, field, roi)
    %set variables 
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).stimResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).stimResponseTrace(:));


    h=figure('units','normalized','outerposition',[0 0 1 1]);   
    Stim = 1;
    for elev=1:metadata.StimParams.numElevation
        for azi = 1:metadata.StimParams.numAzimuth
            subplot(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, Stim)
            avgTempTrace = squeeze(analysis.(field).roi(roi).avgResponseTrace(Stim,:));
            
            %shade stimulus on times
            x=(1:length(avgTempTrace ))./metadata.TwoPhoton.rate-pretrialTime;
            bl=x(stimWindow);
            if analysis.(field).roi(roi).respStim(Stim,1) == 1
                patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [.9 .9 1], 'LineStyle', 'none');
            else
                patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
            end
            hold all
            
            %plot trial traces
            for trial = 1:metadata.StimParams.numTrials
                tempTrace = squeeze(analysis.(field).roi(roi).stimResponseTrace(Stim,trial,:));
                plot(x, tempTrace, 'color', [0.5 0.5 0.5]);
                hold all;
            end
            
            %plot averageTrace
            plot (x,avgTempTrace, '-k', 'LineWidth', 2);
            hold all
            
            axis ('square');
            ylim([ymin ymax])
            xlim([min(x) max(x)])

            ax = gca;
            ax.FontSize = 6; 
            title([num2str(metadata.StimParams.stimPosY(elev)) ' deg Elev, ' num2str(metadata.StimParams.stimPosX(azi)) ' deg Azi'],'FontSize', 7')
            Stim = Stim + 1;
        end
    end
    hold off
    set(gcf, 'Color', 'w')
end