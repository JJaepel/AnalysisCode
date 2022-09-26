function makeGratingPreferencePlots(analysisParams, analysis, data, metadata,types, saveDirectory)
if types == 1
    number = 5;
else
    number =6;
end
h=figure(100+number);
set(h,'units','normalized','outerposition',[0 0 1 1]);
rows = ceil(sqrt(metadata.StimParams.numCon));
for con = 1:metadata.StimParams.numCon
    if metadata.StimParams.numCon > 6
        subplot(rows,rows,con)
    else
        subplot(2,rows,con)
    end
    tempResp = [analysis.(analysisParams.field).roi.isResp_allCon];
    tempResp = tempResp(con:metadata.StimParams.numCon:end);
    isResponsive = find(tempResp == 1);
    if sum(isResponsive) > 0
        if types == 1
            if analysisParams.stimType == 2  
                tempOSI= [analysis.(analysisParams.field).roi.OSIFit_allSf];
            elseif analysisParams.stimType == 3
                tempOSI= [analysis.(analysisParams.field).roi.OSIFit_allTf];
            elseif analysisParams.stimType == 5
                    tempOSI = [analysis.(analysisParams.field).roi.OSIFit_allEyes];
            end
            tempOSI = tempOSI(con:metadata.StimParams.numCon:end);
            tempOriSelect = find(tempOSI > 0.2);
            rois = intersect(tempOriSelect, isResponsive);
            if ~isempty(rois)
                if analysisParams.stimType == 2
                    tempPrefOri = [analysis.(analysisParams.field).roi(rois).preferredOrientation_allSf];
                    tempFit = [analysis.(analysisParams.field).roi(rois).OSIFit_allSf];
                elseif analysisParams.stimType == 3
                        tempPrefOri = [analysis.(analysisParams.field).roi(rois).preferredOrientation_allTf];
                        tempFit = [analysis.(analysisParams.field).roi(rois).OSIFit_allTf];
                elseif analysisParams.stimType == 5
                        tempPrefOri = [analysis.(analysisParams.field).roi(rois).preferredOrientation_allEyes];
                        tempFit = [analysis.(analysisParams.field).roi(rois).OSIFit_allEyes];
                end
                allprefs = tempPrefOri(con:metadata.StimParams.numCon:end);
                allFit = tempFit(con:metadata.StimParams.numCon:end);
                PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
                if analysisParams.stimType == 2
                    title(['Ori pref map, ' num2str(metadata.StimParams.spatialFreq(con)) ' cpd'])
                elseif analysisParams.stimType == 3
                    title(['Ori pref map, ' num2str(metadata.StimParams.temporalFreq(con)) ' Hz'])
                elseif analysisParams.stimType == 5
                    if con == 1
                        title('Ori pref map, contra')
                    elseif con == 2
                        title('Ori pref map, ipsi')
                    elseif con == 3
                        title('Ori pref map, binocular')
                    end
                end
            else
                imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                colormap(hsv)
                if analysisParams.stimType == 2
                    title(['Ori pref map, ' num2str(metadata.StimParams.spatialFreq(con)) ' cpd'])
                elseif analysisParams.stimType == 3
                    title(['Ori pref map, ' num2str(metadata.StimParams.temporalFreq(con))  ' Hz'])
                elseif analysisParams.stimType == 5
                    if con == 1
                        title('Ori pref map, contra')
                    elseif con == 2
                        title('Ori pref map, ipsi')
                    elseif con == 3
                        title('Ori pref map, binocular')
                    end
                end
                caxis([0 180])
                colorbar('Location', 'southoutside');
            end
        elseif types == 2
            if analysisParams.stimType == 2 
                tempDSI= [analysis.(analysisParams.field).roi.DSI_allSf];
            elseif analysisParams.stimType == 3
                tempDSI= [analysis.(analysisParams.field).roi.DSI_allTf];
            elseif analysisParams.stimType == 5
                tempDSI= [analysis.(analysisParams.field).roi.DSI_allEyes];
            end
            tempDSI = tempDSI(con:metadata.StimParams.numCon:end);
            tempDirSelect = find(tempDSI > 0.2);
            rois = intersect(tempDirSelect, isResponsive);
            if ~isempty(rois)
                if analysisParams.stimType == 2
                    tempPrefDir = [analysis.(analysisParams.field).roi(rois).preferredDirection_allSf];
                    tempFit = [analysis.(analysisParams.field).roi(rois).DSI_allSf];
                elseif analysisParams.stimType == 3
                    tempPrefDir = [analysis.(analysisParams.field).roi(rois).preferredDirection_allTf];
                    tempFit = [analysis.(analysisParams.field).roi(rois).DSI_allTf];
                elseif analysisParams.stimType == 5
                    tempPrefDir = [analysis.(analysisParams.field).roi(rois).preferredDirection_allEyes];
                    tempFit = [analysis.(analysisParams.field).roi(rois).DSI_allEyes];
                end
                allprefs = tempPrefDir(con:metadata.StimParams.numCon:end);
                allFit = tempFit(con:metadata.StimParams.numCon:end);
                PlotPrefOnTemplateCon(allprefs, allFit,data, types, data.template,rois)
                if analysisParams.stimType == 2
                    title(['Dir pref map, ' num2str(metadata.StimParams.spatialFreq(con))  ' cpd'])
                elseif analysisParams.stimType == 3
                    title(['Dir pref map, ' num2str(metadata.StimParams.temporalFreq(con))  ' Hz'])
                elseif analysisParams.stimType == 5
                    if con == 1
                        title('Dir pref map, contra')
                    elseif con == 2
                        title('Dir pref map, ipsi')
                    elseif con == 3
                        title('Dir pref map, binocular')
                    end
                end
            else
                imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
                colormap(hsv)
                if analysisParams.stimType == 2
                    title(['Dir pref map, ' num2str(metadata.StimParams.spatialFreq(con))  ' cpd'])
                elseif analysisParams.stimType == 3
                    title(['Dir pref map, ' num2str(metadata.StimParams.temporalFreq(con))  ' Hz'])
                elseif analysisParams.stimType == 5
                    if con == 1
                        title('Dir pref map, contra')
                    elseif con == 2
                        title('Dir pref map, ipsi')
                    elseif con == 3
                        title('Dir pref map, binocular')
                    end
                end
                caxis([0 360])
                colorbar('Location', 'southoutside');
            end
        end
    else
        imshow(cat(3,data.template,data.template,data.template)/prctile(data.template(:),99));
        colormap(hsv)
        if types == 1
            if analysisParams.stimType == 2
                title(['Ori pref map, ' num2str(metadata.StimParams.spatialFreq(con))  ' cpd'])
            elseif analysisParams.stimType == 3
                title(['Ori pref map, ' num2str(metadata.StimParams.temporalFreq(con))  ' Hz'])
            elseif analysisParams.stimType == 5
                    if con == 1
                        title('Ori pref map, contra')
                    elseif con == 2
                        title('Ori pref map, ipsi')
                    elseif con == 3
                        title('Ori pref map, binocular')
                    end
            end
            caxis([0 180])
        elseif types == 2
            if analysisParams.stimType == 2
                title(['Dir pref map, ' num2str(metadata.StimParams.spatialFreq(con))  ' cpd'])
            elseif analysisParams.stimType == 3
                title(['Dir pref map, ' num2str(metadata.StimParams.temporalFreq(con))  ' Hz'])
            elseif analysisParams.stimType == 5
                    if con == 1
                        title('Dir pref map, contra')
                    elseif con == 2
                        title('Dir pref map, ipsi')
                    elseif con == 3
                        title('Dir pref map, binocular')
                    end
            end
            caxis([0 360])
        end
         colorbar('Location', 'southoutside');
    end
end
if types == 1
    saveas(gcf, fullfile(saveDirectory, '105_OrientationPreference_allCon.png'))
elseif types == 2
    saveas(gcf, fullfile(saveDirectory, '106_DirectionPreference_allCon.png'))
end