function Suite2pSpineTifCombiner(server, animal, name)
%combines tifs from multiple levels into one image

%% loading experiment folders
if server
    drive       = 'Z:\Juliane\';
else
    drive           = 'F:\';
end

baseDirectory   = [drive 'Data\2P_data\'];
dirName         = [baseDirectory animal '\' name '\Registered\'];
saveDir         = [baseDirectory animal '\' name '\Registered\combined\'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% look for folders and how many tifs are inside
folders = dir([dirName 'slice*']);
for planes = 1:length(folders)
    tifFiles{planes} = dir([dirName folders(planes).name  '\*.tif']);
end

% load files and write into new tif file
% assumes that we have 4 planes --> add if you have levels with 5 or
% more/less
if length(folders) == 4
    for tifNr = 1:size(tifFiles{1},1)
        tifStack = [];
        for planes = 1:length(folders)
            tifStack = cat(3, tifStack,read_Tiffs([dirName folders(planes).name filesep tifFiles{planes}(tifNr).name], 1, 50));
        end
        dimSize = size(tifStack);
        fileSize = dimSize(3)/4;
        catStack = [tifStack(:,:,1:fileSize) tifStack(:,:,fileSize+1:fileSize*2); tifStack(:,:,fileSize*2+1:fileSize*3) tifStack(:,:,fileSize*3+1:fileSize*4)];
        clear tifStack
        saveFile = [saveDir filesep 'file' num2str(tifNr) '.tif'];
        catStack = uint16(catStack);
        imwrite(catStack(:,:,1),saveFile,'tiff','writemode','overwrite','compression','none');
        for ind=2:size(catStack,3)
            imwrite(catStack(:,:,ind),saveFile,'tiff','writemode','append','compression', 'none');
        end
        disp(['Saving file ' saveFile])
        clear catStack
    end
elseif length(folders) == 5
    for tifNr = 1:size(tifFiles{1},1)
        tifStack = [];
        for planes = 1:length(folders)-1
            tifStack = cat(3, tifStack,read_Tiffs([dirName folders(planes).name filesep tifFiles{planes}(tifNr).name], 1, 50));
        end
        dimSize = size(tifStack);
        fileSize = dimSize(3)/4;
        catStack = [tifStack(:,:,1:fileSize) tifStack(:,:,fileSize+1:fileSize*2); tifStack(:,:,fileSize*2+1:fileSize*3) tifStack(:,:,fileSize*3+1:fileSize*4)];
        clear tifStack
        saveFile = [saveDir filesep 'file' num2str(tifNr) '.tif'];
        catStack = uint16(catStack);
        imwrite(catStack(:,:,1),saveFile,'tiff','writemode','overwrite','compression','none');
        for ind=2:size(catStack,3)
            imwrite(catStack(:,:,ind),saveFile,'tiff','writemode','append','compression', 'none');
        end
        disp(['Saving file ' saveFile])
        clear catStack
    end
else
    disp('There are more or less than four planes +/- flyback frame')
end

end