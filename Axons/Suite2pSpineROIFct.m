function  Suite2pSpineROIFct(server, animal, name, level)
%% loading experiment data
if server
    drive       = 'Z:\Juliane\';
else
    drive           = 'F:\';
end

baseDirectory   = [drive 'Data\2P_data\'];
filename        = 'Projection.tif';

if level
    dirName         = [baseDirectory animal '\' name '\Registered\combined\'];
    saveDir         = [baseDirectory animal '\' name '\Registered\combined\Projection\'];
else
    dirName         = [baseDirectory animal '\' name '\Registered\slice1\'];
    saveDir         = [baseDirectory animal '\' name '\Registered\slice1\Projection\'];
end

saveFile        = [saveDir filename];
if ~exist(saveDir, 'dir')
    % make new file directory
    mkdir(saveDir); 
end

%% read in the tiff files
tifStack = [];
tifFiles = dir([dirName '*.tif']);
readFiles = min([3, length(tifFiles)]);
for currentFile = 1:readFiles  
    % Specify stack name
    filePath = [dirName tifFiles(currentFile).name];

    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs(filePath,1, 50));
end

%% make projection for spine ROIing
cd(saveDir)
meanImg = uint16(mean(tifStack,3));
imwrite(meanImg, saveFile, 'tiff', 'writemode', 'overwrite', 'compression', 'none')
avg = mijread([saveDir filename]);
MIJ.run('SpineROIs')
f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
h = uicontrol('Position',[10 10 190 30],'String','Save and Next Experiment?','Callback','uiresume(gcbf)');
uiwait(gcf);
MIJ.run('AxonROIs')
MIJ.run('Close All');
close gcf

%% convert RoiSet to mat file
dim = size(meanImg);
[ROIs, ~] = ROIconvert('RoiSet.zip', [dim(1) dim(2)]);

%% extract ROI positions
data = [];
for m = 1:length(ROIs)    
    data.roi(m).xPos = ROIs(m).x;
    data.roi(m).yPos = ROIs(m).y;
    data.roi(m).mask = [ROIs(m).perimeter];
    data.roi(m).name = m;
    if level
        if data.roi(m).xPos > 512
            if data.roi(m).yPos > 512
                data.roi(m).plane = 4;
            else
                data.roi(m).plane = 3;
            end
        else
            if data.roi(m).yPos > 512
                data.roi(m).plane = 2; 
            else
                data.roi(m).plane = 1;
            end
        end
    end
end

%% save meanImg as template
data.template = double(meanImg);
save('ROIs.mat', 'data', '-mat') 
end