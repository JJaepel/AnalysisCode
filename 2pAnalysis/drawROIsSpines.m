function drawROIsSpines(analysisParams, expInfo, ind)

%This automatically opens the projection in Fiji, so that the user can mark
%all ROIs that are then automatically saved and converted to MATLAB
%readadable file. 
%
%Important: Requires the following Miji Plugins to be installed -  needs to
%be in the folder Fiji.app -> Plugins ->AxonROIs.ijm
%- SpineROIs: Opens the ROI manager
% run("ROI Manager...");
%- AxonROIs: Saves the ROIs
% roiManger("Save", "RoiSet.zip");
%
%Input: 
% - analysisParams: contains all info about the parameters of the analysis,
% here important for knowing which server to use
% - expInfo: contains information about all the experiments, such as name,
% animal, whether it is a volume
% - ind: which experiments do you want to run?

%% Step 0: Open the link to Fiji

Miji

%% Step 1: Set the correct folders and check if there are ROI files yet
for i = ind
    %set folder
    if analysisParams.server
        baseDir = ['Z:\Juliane\Data\2P_data\' char(expInfo.animal{i}) '\' char(expInfo.name{i}) '\Registered\'];
    else 
        baseDir = ['F:\Data\2P_data\' char(expInfo.animal{i}) '\' char(expInfo.name{i}) '\Registered\'];
    end

    % if it is a volume, combine the planes (if not already done)
    if expInfo.vol{i}
        RegTifsDir = [baseDir '\combined\'];
        if ~exist(RegTifsDir, 'dir')
            Suite2pSpineTifCombiner(analysisParams.server, char(expInfo.animal{i}), char(expInfo.name{i}))
        end
        ROIDir = [baseDir '\combined\Projection\'];
    else
        ROIDir = [baseDir '\slice1\Projection\'];
    end
    
    if ~exist(ROIDir, 'dir')% make new file directory 
        mkdir(ROIDir); 
    end 
    saveFile= [ROIDir 'Projection.tif'];
    
    %make ROIs if there is no ROI file yet
    ROIFiles = dir([ROIDir 'ROIs.mat']);
    if isempty(ROIFiles)
        %% Step 2: Load the tifFiles and make mean image
        %load tif files
        tifStack = [];
        tifFiles = dir([baseDir '\slice1\*.tif']);
        readFiles = min([3, length(tifFiles)]);
        for currentFile = 1:readFiles  
            % Specify stack name
            filePath = [baseDir '\slice1\' tifFiles(currentFile).name];

            % Read images into tifStack
            tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
        end
        
        %make mean image and save it
        cd(ROIDir)
        meanImg = uint16(mean(tifStack(:,:,1:1000),3));
        imwrite(meanImg, saveFile, 'tiff', 'writemode', 'overwrite', 'compression', 'none')
        
        %% Step 3: Run functions in Miji
        %read in the mean image
        avg = mijread([ROIDir 'Projection.tif']);
        MIJ.run('SpineROIs') %opens the ROImanager
        %add a figure and wait for user to add all ROIs
        f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
        h = uicontrol('Position',[10 10 190 30],'String','Save and Next Experiment?','Callback','uiresume(gcbf)');
        uiwait(gcf);
        %once the user has selected all the ROIs, continue
        MIJ.run('AxonROIs') %saves the ROIs in the folder
        MIJ.run('Close All'); %closes all figures in Miji
        close gcf %closes the waitbar
        
        %% Step 4: Convert the zip file to .mat file
        dim = size(meanImg);
        [ROIs, ~] = ROIconvert('RoiSet.zip', [dim(1) dim(2)]);
        save('ROIs.mat', 'ROIs', '-mat') 
    end
end