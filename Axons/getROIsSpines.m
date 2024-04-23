function getROIsSpines(analysisParams, exp_info, ind)

Miji
for i = ind
    %set folder
    if analysisParams.server
        baseDir = ['Z:\Juliane\Data\2P_data\' char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    else 
        baseDir = ['F:\Data\2P_data\' char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    end

    % if it is a volume, combine the planes (if not already done)
    if exp_info.vol{i}
        RegTifsDir = [baseDir '\combined\'];
        if ~exist(RegTifsDir, 'dir')
            Suite2pSpineTifCombiner(analysisParams.server, char(exp_info.animal{i}), char(exp_info.name{i}))
        end
        ROIDir = [baseDir '\combined\Projection\'];
    else
        ROIDir = [baseDir '\slice1\Projection\'];
    end
    
    %% make ROIs if there is no ROI file yet
    ROIFiles = dir([ROIDir 'ROIs.mat']);
    if isempty(ROIFiles)
        %read in the tiff files
        tifStack = [];
        tifFiles = dir([baseDir '\slice1\*.tif']);
        readFiles = min([3, length(tifFiles)]);
        for currentFile = 1:readFiles  
            % Specify stack name
            filePath = [dirName tifFiles(currentFile).name];

            % Read images into tifStack
            tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
        end
        
        % make projection for spine ROIing
        cd(saveDir)
        meanImg = uint16(mean(tifStack(:,:,1:1000),3));
        imwrite(meanImg, saveFile, 'tiff', 'writemode', 'overwrite', 'compression', 'none')
        avg = mijread([saveDir filename]);
        
        %run Fiji plugin
        MIJ.run('SpineROIs') %opens the ROI manager
        f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
        h = uicontrol('Position',[10 10 190 30],'String','Save and Next Experiment?','Callback','uiresume(gcbf)');
        uiwait(gcf);
        MIJ.run('AxonROIs') %saves the ROIs
        MIJ.run('Close All');
        close gcf

    end
end