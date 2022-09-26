function PlotTrialStimResponseOri(metadata, analysis, field, roi, saveDir)
    
    numCon = metadata.StimParams.numCon;
    
    conTr = zeros(metadata.StimParams.numTrials+5,3,6);
    conTr(:,:,1)= cbrewer('seq', 'Reds', metadata.StimParams.numTrials+5);
    conTr(:,:,2)= cbrewer('seq', 'Blues', metadata.StimParams.numTrials+5);
    conTr(:,:,3)= cbrewer('seq', 'Greens', metadata.StimParams.numTrials+5);
    conTr(:,:,4)= cbrewer('seq', 'Purples', metadata.StimParams.numTrials+5);
    conTr(:,:,5)= cbrewer('seq', 'Oranges', metadata.StimParams.numTrials+5);
    conTr(:,:,6)= cbrewer('seq', 'Greys', metadata.StimParams.numTrials+5);
    conTr(conTr>1)= 1; conTr(conTr<0) = 0;

    pretrialTime= analysis.(field).preTrialTime;
    stimWindow=(analysis.(field).stimStart: analysis.(field).stimStop);
    ymaxTrial=max(analysis.(field).roi(roi).stimResponseTrace(:));
    ymaxAvg = max(analysis.(field).roi(roi).avgResponseTrace(:));
    ymax = max(ymaxTrial, ymaxAvg);
    yminTrial=min(analysis.(field).roi(roi).stimResponseTrace(:));
    yminAvg = min(analysis.(field).roi(roi).avgResponseTrace(:));
    ymin = max(yminTrial, yminAvg);

    if (ymax-ymin) > 7
        refBarY = 5;
    elseif (ymax-ymin) > 3
        refBarY = 1;
    elseif (ymax-ymin) > 0.5
        refBarY = 0.1;
    end
  
    if metadata.StimParams.sweep
        for i=1:metadata.StimParams.numDirections
            for con = numCon
                stimID = metadata.StimParams.numDirections*(con-1)+i;
                plotnr = metadata.StimParams.numTf*(i-1)+con;
                ax=subplot(metadata.StimParams.numOrientations, numCon,plotnr);
                cla(ax)
                y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(stimID,:),3);
                x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
                if stimWindow(1)~=0
                        bl=x(stimWindow);
                        patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [.8 .8 .8], 'LineStyle', 'none');
                end
                hold all
                for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                    plot(ax,x, squeeze(analysis.(field).roi(roi).stimResponseTrace(stimID,trial,:)), 'Color', conTr(trial+1,:,con));
                    hold all
                end
                plot(ax,x, squeeze(analysis.(field).roi(roi).avgResponseTrace(stimID,:)), 'Color', [0 0 0], 'LineWidth', 3);
                hold all
                ylim([ymin, ymax])
                xlim([min(x) max(x)])
                title(metadata.StimParams.directions(i))
                axis off
                if stimID == 1
                    %add 0.1 F/F0 as reference
                    patch([min(x) min(x)], [ymin ymin+refBarY],[1 1 1], 'LineWidth', 3)
                    text(min(x),ymin-0.025*(ymax-ymin),'1 s')  
                    %add 1 s as reference
                    patch([min(x) min(x)+1], [ymin ymin], [1 1 1], 'LineWidth',3);
                    text(min(x)-1,ymin, num2str(refBarY))            
                end
            end

        end
    else
        for con = 1:numCon
            h=figure('units','normalized','outerposition',[0 0 1 1]);
            for i=1:metadata.StimParams.numDirections
                stimID = metadata.StimParams.numDirections*(con-1)+i;
                for trial =1:size(analysis.(field).roi(roi).stimResponseTrace,2)
                    plotID = i+metadata.StimParams.numDirections*(trial-1);
                    ax=subplot(metadata.StimParams.numTrials,metadata.StimParams.numDirections, plotID);
                    cla(ax)
                    y=squeeze(medfilt1(analysis.(field).roi(roi).stimResponseTrace(stimID,trial,:),3));
                    x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
                    
                    if stimWindow(1)~=0
                        bl=x(stimWindow);
                        patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [0.8 0.8 0.8], 'LineStyle', 'none');
                    end
                    hold all
                    plot(ax,x,y, 'Color', conTr(trial+5,:,con));
                    hold all
                    ylim([ymin, ymax])
                    xlim([min(x) max(x)])
                    if trial == 1
                        title(metadata.StimParams.directions(i))
                    end
                    axis off
                end
                hold off
            end
            set(gcf, 'Color', 'w')
            saveas(gcf, fullfile(saveDir, ['TrialStimResp_ROI_Nr_' num2str(roi) '_Con_' num2str(con) '.png']))
            %close gcf
        end
    end

end