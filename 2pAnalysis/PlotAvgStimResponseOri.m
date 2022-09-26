function PlotAvgStimResponseOri(metadata, analysis, field, roi,saveDir)
    %PlotAvgStimResponse Plots Avg + SEM response for stimuli
    h=figure(1000+roi);
    set(h,'units','normalized','outerposition',[0 0 1 1]);
    numCon = metadata.StimParams.numCon;
    if numCon > 1
        if metadata.StimParams.numTf > 1
            conLabels = metadata.StimParams.temporalFreq;
        elseif metadata.StimParams.numSf > 1
            conLabels = metadata.StimParams.spatialFreq;
        else
            conLabels = ['ipsi','contra', 'binocular']; %check label order
        end
    end
    cocCon = cbrewer('qual', 'Set1', max(numCon,3));
    
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
    else 
        refBarY = 0.01;
    end
    id = 1;
    for tf = 1:numCon
        for i=1:metadata.StimParams.numDirections 
            ax=subplot(numCon,metadata.StimParams.numDirections, id);
            cla(ax)
            y=medfilt1(analysis.(field).roi(roi).avgResponseTrace(id,:),3);
            x=(1:length(y))./metadata.TwoPhoton.rate-pretrialTime;
            if stimWindow(1)~=0
                bl=x(stimWindow); %make a patch of the stimulus window
                patch([bl fliplr(bl)], [ymin*ones(1,length(bl)) ymax*ones(1,length(bl))], [0.8 0.8 0.8], 'LineStyle', 'none');
            end
            err=analysis.(field).roi(roi).SEMResponseTrace(id,:);
            hold all
            shadedErrorBar(x,y,err, 'lineProps',{'-','Color',cocCon(tf,:), 'lineWidth', 3})
            hold all
            if id == 1
                %add 0.1 F/F0 as reference
                patch([min(x) min(x)], [ymin ymin+refBarY],[1 1 1], 'LineWidth', 3)
                text(min(x),ymin-0.025*(ymax-ymin),'1 s')  
                %add 1 s as reference
                patch([min(x) min(x)+1], [ymin ymin], [1 1 1], 'LineWidth',3);
                text(min(x)-1,ymin, num2str(refBarY))            
            end
            
            %add condition label at start of plotting for each condition
            try
                ylim([ymin, ymax])
            end
            xlim([min(x) max(x)])
            title(metadata.StimParams.directions(i))
            axis off
            if mod(id,metadata.StimParams.numDirections) == 1
               if metadata.StimParams.numTf > 1
                   t = text(min(x)-1,(ymax-ymin)*0.2, [num2str(conLabels(tf)) ' Hz'], 'FontSize', 12);
               elseif metadata.StimParams.numSf > 1
                   t = text(min(x)-1,(ymax-ymin)*0.2, [num2str(conLabels(tf)) ' cpd'], 'FontSize', 12);
               elseif metadata.StimParams.numCon > 3
                   t = text(min(x)-1,(ymax-ymin)*0.2, conLabels(tf), 'FontSize', 12);
               else
                   t = text(min(x)-1,(ymax-ymin)*0.2, [num2str(metadata.StimParams.temporalFreq) ' Hz, ' num2str(metadata.StimParams.spatialFreq) ' cpd'], 'FontSize', 12);
               end
               set(t,'Rotation',90);
            end
            id = id+1;
        end
    end
    hold off
    set(gcf, 'Color', 'w')
    saveas(gcf, fullfile(saveDir, ['AvgStimResp_ROI_Nr_' num2str(roi) '.png']))
    close gcf
end