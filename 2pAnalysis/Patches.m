function Patches(analysisParams)

close all

%% 0.) define folders and structures
if analysisParams.server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
analysisParams.Sp2dir = [drive '\Data\Spike2Data\'];
analysisParams.savedir = [drive '\Data\ImageAnalysis\'];

base2pDirectory= [TwoPhontondir analysisParams.animal];
tifDirectory = [base2pDirectory filesep analysisParams.name];
Sp2dDirectory = [analysisParams.Sp2dir analysisParams.animal filesep analysisParams.sp2ID filesep];
saveDirectory = [analysisParams.savedir analysisParams.animal filesep analysisParams.expID filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];

if ~exist(saveDirectory, 'dir')
    % make new file directory
    mkdir(saveDirectory);  
end

metadata = struct;
metadata.ROI = struct;
analysis = struct;

%% 1). load Data and metadata
disp('Loading data')
if analysisParams.reloadData
    % make sure that you also do analysis afterwards
    analysisParams.reanalyse = 1;
    
    %load stimulus parameters and spike2data
    analysisParams.baseDirectory = base2pDirectory;
    metadata.StimParams=LoadStimParamsRet(Sp2dDirectory);
    metadata.TwoPhoton=LoadFrameTimes(Sp2dDirectory);
    metadata.StimParams.path=fullfile(Sp2dDirectory);
    
    Stimtype = metadata.StimParams.type;
    switch Stimtype
        case 'SparseNoise'
            metadata.StimParams.numTrials = metadata.StimParams.all;
        case 'rotatingGratingPatch'
            metadata.StimParams.stimDuration=metadata.StimParams.spinPeriodInSeconds;
    end
    
    % load tif data
    switch analysisParams.dataType
        case 1
            data = LoadRoisS2p(analysisParams);
            
            %do baseline filtering and compute dff
            if strcmp(analysisParams.field, 'dff')
                [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                data = computeDff(data, 'rawF', 'baseline', 'dff');
            end
        case 2
            if analysisParams.makeROIs
                if analysisParams.level == 1
                    Suite2pAxonTifCombiner(analysisParams)
                end
                data = Suite2pAxonExtractorFct(analysisParams);
            else
                if analysisParams.level == 1

                    dataFile = ([base2pDirectory filesep analysisParams.expID filesep 'suite2p\combined\data.mat']);
                else
                    dataFile = ([base2pDirectory filesep analysisParams.expID filesep 'suite2p\plane0\data.mat']);
                end
                File = load(dataFile, 'data');
                data = File.data;
                clear File
            end
            %do baseline filtering and compute dff
            if strcmp(analysisParams.field, 'dff')
                [metadata, data] = baselinePercentileFilter(metadata, data,'rawF', 'baseline', 60, 30);
                data = computeDff(data, 'rawF', 'baseline', 'dff');
            end
            
            
        case 3 
            %do spine data reloading here
            data = SpineROIExtractFct(analysisParams);
            data = computeDffSpines(data);
            structureType = analysisParams.region;
            switch structureType
                case 'cells'
                    data = FitCellsToDendrite(analysisParams,metadata, data);
                case 'dendrite'
                    data = computeDendriticSubstraction(analysisParams, metadata, data);    
            end
            analysisParams.field = 'rawRes';     
    end
    
    %save file
    save([saveDirectory, 'Patches.mat'], 'data', 'metadata', 'analysisParams', 'analysis');
%     if analysisParams.checkDFF
%         figure
%         subplot(3,1,1)
%         plot(data.roi(1).rawF)
%         axis off
%         subplot(3,1,2)
%         plot(data.roi(1).baseline)
%         axis off
%         subplot(3,1,3)
%         plot(data.roi(1).dff)
%         axis off
%     end
else
    load([saveDirectory, 'Patches.mat'], 'data', 'metadata', 'analysisParams', 'analysis');
end

%% 2.) Specify stimulus settings
disp('Specifiying stimulus settings')
Stimtype = metadata.StimParams.type;
[analysisParams, metadata] = getStimParamsPatches(analysisParams, metadata);

%% 3.) Run analysis
switch Stimtype
    case 'SparseNoise'
        [analysisParams, metadata, analysis] = calculateSignificantPatches(analysisParams, metadata, data, analysis);
    case 'Patch'
        [analysisParams, metadata, analysis] = calculateSignificantPatches(analysisParams, metadata, data, analysis);
        [analysisParams, metadata, data, analysis] = calcROIPropertiesPatches(analysisParams, metadata, data, analysis);
    case 'rotatingGratingPatch'
        [analysisParams, metadata, analysis] = calculateSignificantPatches(analysisParams, metadata, data, analysis);
        [analysisParams, metadata, data, analysis] = calcROIPropertiesPatches(analysisParams, metadata, data, analysis);
    case 'Retinotopy_2D'
        [analysisParams, metadata, analysis] = calculateSignificantPatches(analysisParams, metadata, data, analysis);
        [analysisParams, metadata, data, analysis] = calcROIPropertiesPatches(analysisParams, metadata, data, analysis);
    otherwise
        [analysisParams, metadata, data, analysis] = calcROIPropertiesPatches(analysisParams, metadata, data, analysis);
end



%% 4.) Plot results
coc_prop = cbrewer('qual', 'Paired', 12);

%% 4.a) Plot Preference on top of template
switch Stimtype
    case 'SparseNoise'
    otherwise
        for types = 1:2
    if types == 1
        rois = linspace(1,length(analysis.(analysisParams.field).roi),length(analysis.(analysisParams.field).roi));
        allPatchPrefs = [analysis.(analysisParams.field).roi.prefPatch];
        if metadata.StimParams.TwoDStim == 1
            allElevPrefs = [analysis.(analysisParams.field).roi.prefElevDeg];
            allAziPrefs = [analysis.(analysisParams.field).roi.prefAziDeg];
        end
    elseif types == 2
        rois = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        allPatchPrefs = [analysis.(analysisParams.field).roi(rois).prefPatch];
        if metadata.StimParams.TwoDStim == 1
            allElevPrefs = [analysis.(analysisParams.field).roi(rois).prefElevDeg];
            allAziPrefs = [analysis.(analysisParams.field).roi(rois).prefAziDeg];
        end
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    
    if metadata.StimParams.TwoDStim == 1
        subplot(2,3,1)
        PlotPrefOnTemplateRet(analysis, data, metadata,2, analysisParams.field,data.template, rois)

        subplot(2,3,4)
        histogram(allElevPrefs,linspace(metadata.StimParams.stimPosY(1), metadata.StimParams.stimPosY(end),metadata.StimParams.numElevation), 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
        ylabel('Cells');
        xlabel(sprintf('Elevation preference (%s)',char(145)));
        axis square;
        set(gca,'Box','off');

        subplot(2,3,2)
        PlotPrefOnTemplateRet(analysis, data, metadata,3, analysisParams.field,data.template, rois)

        subplot(2,3,5)
        histogram(allAziPrefs,linspace(metadata.StimParams.stimPosX(1), metadata.StimParams.stimPosX(end),metadata.StimParams.numAzimuth), 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
        ylabel('Cells');
        xlabel(sprintf('Azimuth preference (%s)',char(145)));
        axis square;
        set(gca,'Box','off');

        subplot(2,3,3)
        PlotPrefOnTemplateRet(analysis, data, metadata,1, analysisParams.field,data.template, rois)

        subplot(2,3,6)
        [NprefPatch, ~] = histcounts(allPatchPrefs, metadata.StimParams.numPatches);
        try
            NprefPatch = reshape(NprefPatch, metadata.StimParams.numElevation, metadata.StimParams.numAzimuth);
        catch
            NprefPatch = reshape(NprefPatch, metadata.StimParams.numElevation*2, metadata.StimParams.numAzimuth);

        end
        imagesc(NprefPatch)
        colorbar
        colormap('jet')
        axis off
        axis square
    
        ylabel('Cells');
        xlabel(sprintf('Patch preference'));
        axis square;
        set(gca,'Box','off');
    else
        subplot(2,2,[1 3])
        PlotPrefOnTemplateRet(analysis, data, metadata,1, analysisParams.field,data.template, rois)
        
        subplot(2,2,2)
        CM = colormap('jet');
        IX = round(linspace(1,64,metadata.StimParams.numPatches));
        for cc = 1:length(IX)
            C{cc} = CM(IX(cc),:);
        end
        switch Stimtype
            case 'RetWedge'
                stimPatches = [6 1; 5 2; 4 3];
                imagesc(stimPatches)
            case 'Ret_Annulus'
                for sizes = 1:metadata.StimParams.numSizes
                    plot(0,0,'o','MarkerSize',30*sizes, 'Color', C{sizes})
                    hold on
                end
        end
        axis off
        axis square
        title('Patch positions')
        
        subplot(2,2,4)
        [NprefPatch, ~] = histcounts(allPatchPrefs, metadata.StimParams.numPatches);
        for sizes = 1:metadata.StimParams.numPatches
            h = bar(sizes, NprefPatch(sizes));
            set(h,'FaceColor', C{sizes})
            hold all
        end
        
%         switch Stimtype
%             case 'RetWedge'
%                 NprefPatch = reshape(NprefPatch([6 5 4 1 2 3]),3,2);
%                 imagesc(NprefPatch)
%                 colorbar('Location', 'SouthOutside')
%             case 'Ret_Annulus'
%                 for sizes = 1:metadata.StimParams.numSizes
%                     h = bar(sizes, NprefPatch(sizes))
%                     set(h,'FaceColor', C{sizes})
%                     hold all
%                 end             
%         end

        axis off
        axis square
        title('Patch preference')
        
    end
    
    set(gcf, 'color', 'w');
    if types == 1
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_all_cells.png'))
    elseif types == 2
        saveas(gcf, fullfile(saveDirectory, 'Overlaymaps_resp_cells.png'))
    end
        end
end

%% 4. b) Plot overview of significant ROIs
switch Stimtype
    case 'SparseNoise'
    case 'Patch'
        signROIs = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        figure('units','normalized','outerposition',[0 0 1 1]); 
        i = 1; j = 1;
        for roi = 1:length(signROIs)
            ax(1)= subplot(8,10,i);
            imagesc(reshape(analysis.(analysisParams.field).roi(roi).avgRespNormalized(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth))
            colormap(ax(1),parula)
            axis ('square')
            axis off
            title(['ROI ' num2str(signROIs(roi)) ', Black'],'FontSize', 7')
            
            ax(2) = subplot(8,10,i+1);
            imagesc(reshape(analysis.(analysisParams.field).roi(roi).avgRespNormalized(metadata.StimParams.numElevation * metadata.StimParams.numAzimuth+1:end-1), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth))
            colormap(ax(2),hot)
            axis ('square')
            axis off
            title(['ROI ' num2str(signROIs(roi)) ', White'],'FontSize', 7')
            i = i+2;
            if i > 80
                set(gcf, 'Color', 'w')
                saveas(gcf, fullfile(saveDirectory, ['SignROIs_Fields_' num2str(j) '.png']))
                close gcf 
                j = j+1; i= 1;
                figure('units','normalized','outerposition',[0 0 1 1]); 
            end
        end
        set(gcf, 'Color', 'w')
        saveas(gcf, fullfile(saveDirectory, ['SignROIs_Fields_' num2str(j) '.png']))
        close gcf
    otherwise
        signROIs = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
        figure
        i = 1; j = 1;
        for roi = 1:length(signROIs)
            subplot(8,10,i)
            imagesc(reshape(analysis.(analysisParams.field).roi(roi).avgRespNormalized(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth))
            axis ('square')
            axis off
            colormap(hot)
            title(['ROI ' num2str(signROIs(roi))],'FontSize', 7')
            i = i+1;
            if i > 80
                set(gcf, 'Color', 'w')
                saveas(gcf, fullfile(saveDirectory, ['SignROIs_Fields_' num2str(j) '.png']))
                close gcf 
                j = j+1; i= 1;
                figure('units','normalized','outerposition',[0 0 1 1]); 
            end
        end
        set(gcf, 'Color', 'w')
        saveas(gcf, fullfile(saveDirectory, ['SignROIs_Fields_' num2str(j) '.png']))
        close gcf
end

%% 4. c) Calculate and Plot the average normalized RF of all significant ROIs and compare across the two analysis types
switch Stimtype
    case 'SparseNoise'
    case 'Patch'
        RespMatrixW = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
        RespMatrixB = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
        for roi = 1:length(signROIs)
            RespMatrixW(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgStimResponse(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
            RespMatrixB(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgStimResponse(metadata.StimParams.numElevation * metadata.StimParams.numAzimuth+1:end-1), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
        end
        AvgRespMatrixW = mean(RespMatrixW, 3);
        AvgRespMatrixB = mean(RespMatrixB, 3);
        
        
        RespMatrixNormW = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
        RespMatrixNormB = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
        for roi = 1:length(signROIs)
            RespMatrixNormW(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgRespNormalized(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
            RespMatrixNormB(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgRespNormalized(metadata.StimParams.numElevation * metadata.StimParams.numAzimuth+1:end-1), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
        end
        AvgRespMatrixNormW = mean(RespMatrixNormW, 3);
        AvgRespMatrixNormB = mean(RespMatrixNormB, 3);
        
        figure
        ax(1) = subplot(2,2,1);
        imagesc(AvgRespMatrixW)
        colormap(ax(1),parula)
        axis ('square')
        axis off
        title('Average Response RF, White','FontSize', 10')
        
        ax(2) = subplot(2,2,2);
        imagesc(AvgRespMatrixB)
        colormap(ax(2),hot)
        axis ('square')
        axis off
        title('Average Response RF, Black','FontSize', 10')
        
        ax(1) = subplot(2,2,3);
        imagesc(AvgRespMatrixNormW)
        colormap(ax(1),parula)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF, White','FontSize', 10')
        
        ax(2) = subplot(2,2,4);
        imagesc(AvgRespMatrixNormB)
        colormap(ax(2),hot)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF, Black','FontSize', 10')
        
        set(gcf, 'Color', 'w')
        saveas(gcf, fullfile(saveDirectory, 'SignROIs_AvgRF.png'))
        close gcf
    otherwise
        RespMatrix = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
        for roi = 1:length(signROIs)
            RespMatrix(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgStimResponse(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
        end
        AvgRespMatrix = mean(RespMatrix, 3);
        
        RespMatrixNorm = zeros(metadata.StimParams.numElevation,metadata.StimParams.numAzimuth, length(signROIs));
        for roi = 1:length(signROIs)
            RespMatrixNorm(:,:,roi) = reshape(analysis.(analysisParams.field).roi(roi).avgRespNormalized(1:metadata.StimParams.numElevation * metadata.StimParams.numAzimuth), metadata.StimParams.numElevation,metadata.StimParams.numAzimuth);
        end
        AvgRespMatrixNorm = mean(RespMatrixNorm, 3);
        
        figure
        subplot(2,1,1)
        imagesc(AvgRespMatrix)
        colormap(hot)
        axis ('square')
        axis off
        title('Average Response RF','FontSize', 10')
        
        subplot(2,1,2)
        imagesc(AvgRespMatrixNorm)
        colormap(hot)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF','FontSize', 10')
        set(gcf, 'Color', 'w')
        saveas(gcf, fullfile(saveDirectory, 'SignROIs_AvgRF.png'))
        close gcf
end

switch Stimtype
    case 'Patches'
        figure
        ax(1) = subplot(2,2,1);
        imagesc(AvgRespMatrixNormW)
        colormap(ax(1),parula)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF, ON','FontSize', 10')
        
        figure
        ax(1) = subplot(2,2,2);
        imagesc(AvgRespMatrixNormB)
        colormap(ax(1),parula)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF, OFF','FontSize', 10')
    case 'rotatingGratingPatch'
        figure
        ax(1) = subplot(2,2,1);
        imagesc(AvgRespMatrixNorm)
        colormap(ax(1),hot)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF','FontSize', 10')
    case 'Retinotopy_2D'
        figure
        ax(1) = subplot(2,2,1);
        imagesc(AvgRespMatrixNorm)
        colormap(ax(1),hot)
        axis ('square')
        axis off
        title('Norm. Avg. Response RF','FontSize', 10')
    otherwise
end

%% 4. d) Plot individual ROIs
switch Stimtype
    case 'SparseNoise'
    otherwise
        if analysisParams.plotROIs
            if ~exist(ROIRespsaveDirectory, 'dir')
                mkdir(ROIRespsaveDirectory);
            else
                %remove old files
                filePattern = fullfile(ROIRespsaveDirectory, '*.png'); % Change to whatever pattern you need.
                theFiles = dir(filePattern);
                for k = 1 : length(theFiles)
                  baseFileName = theFiles(k).name;
                  fullFileName = fullfile(ROIRespsaveDirectory, baseFileName);
                  delete(fullFileName);
                end
            end
            if ~analysisParams.plotRespROIsOnly
                if ~exist(ROINonRespsaveDirectory, 'dir')
                    mkdir(ROINonRespsaveDirectory);
                else
                    %remove old files
                    filePattern = fullfile(ROIRespsaveDirectory, '*.png'); % Change to whatever pattern you need.
                    theFiles = dir(filePattern);
                    for k = 1 : length(theFiles)
                      baseFileName = theFiles(k).name;
                      fullFileName = fullfile(ROIRespsaveDirectory, baseFileName);
                      delete(fullFileName);
                    end
                end
            end
            for i = 1:length(data.roi)
                if analysis.(analysisParams.field).roi(i).isResponseSignificant == 1
                    if metadata.StimParams.TwoDStim == 1
                        switch Stimtype
                            case 'Patch'
                                PlotStimResponsePatch(metadata, analysis, analysisParams.field, i)
                                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_StimResp.png']))
                                close gcf
                            otherwise
                                PlotStimResponse2D(metadata, analysis, analysisParams.field, i)
                                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_StimResp.png']))
                                close gcf
                        end                  
                    else
                        PlotTrialStimResponse1D(metadata, analysis, analysisParams.field, i)
                        saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                        close gcf
                        PlotAvgStimResponse1D(metadata, analysis, analysisParams.field, i)
                        saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                        close gcf
                    end
                else
                    if ~analysisParams.plotRespROIsOnly
                        if metadata.StimParams.TwoDStim == 1
                            try
                                PlotTrialStimResponse2D(metadata, analysis, analysisParams.field, i)
                                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                                close gcf
                                PlotAvgStimResponse2D(metadata, analysis, analysisParams.field, i)
                                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                                close gcf
                            catch
                                close gcf
                            end
                        else
                            try
                                PlotTrialStimResponse1D(metadata, analysis, analysisParams.field, i)
                                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_TrialStimResp_.png']))
                                close gcf
                                PlotAvgStimResponse1D(metadata, analysis, analysisParams.field, i)
                                saveas(gcf, fullfile(ROIRespsaveDirectory, ['ROI_Nr_' num2str(i) '_AvgStimResp_.png']))
                                close gcf
                            catch
                                close gcf
                            end
                        end
                    end
                end
            end
        end
end

save(fullfile(saveDirectory, 'PatchesAna.mat'), 'data', 'metadata', 'analysisParams', 'analysis');
end
%%

    
