function PlotTrialStimResponseOri(metadata, analysis, field, roi)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    numCon = metadata.StimParams.numCon;
    colorlevels = numCon;
    colorlevelsTr = metadata.StimParams.numTrials+1;
    cocCon = cbrewer('qual', 'Set1', colorlevels);
    cocTr= [linspace(0,1,colorlevelsTr);linspace(0,1,colorlevelsTr);linspace(0,1,colorlevelsTr)]';
    cocTr = cocTr(2:end,:);
    
    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);
    ymax=max(analysis.(field).roi(roi).avgResponseTrace(:)+analysis.(field).roi(roi).SEMResponseTrace(:));
    ymin=min(analysis.(field).roi(roi).avgResponseTrace(:)-analysis.(field).roi(roi).SEMResponseTrace(:));

    if (ymax-ymin) > 3
        refBarY = 1;
    else
        refBarY = 0.1;
    end
  
    if metadata.StimParams.sweep
        for i=1:metadata.StimParams.numDirections
            for con = numCon
                id = metadata.StimParams.numDirections*(con-1)+i;
                plotnr = metadata.StimParams.numTf*(i-1)+con;
                ax=subplot(metadata.StimParams.numOrientations, numCon,plotnr);
                cla(ax)
                y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(id,:),3);
                x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
                if stimWindow(1)~=0
                        bl=x(stimWindow);
                        patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
                end
                hold all
                for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                    plot(ax,x, squeeze(analysis.(field).roi(roi).stimResponseTrace(id,trial,:)), 'Color', cocTr(trial,:));
                    hold all
                end
                plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', cocCon(con,:), 'LineWidth', 3);
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
    else
        for con = 1:numCon
            for i=1:metadata.StimParams.numDirections
                id = metadata.StimParams.numDirections*(con-1)+i;
                ax=subplot(numCon,metadata.StimParams.numDirections, id);
                cla(ax)
                y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(id,:),3);
                x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
                if stimWindow(1)~=0
                        bl=x(stimWindow);
                        patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [1 .9 .9], 'LineStyle', 'none');
                end
                hold all
                for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                    plot(ax,x, squeeze(analysis.(field).roi(roi).stimResponseTrace(id,trial,:)), 'Color', cocTr(trial,:));
                    hold all
                end
                plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(id,:)), 'Color', cocCon(con,:), 'LineWidth', 3);
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
    end
    hold off
    set(gcf, 'Color', 'w')
end