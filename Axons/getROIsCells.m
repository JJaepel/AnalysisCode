function getROIsCells(analysisParams, exp_info, ind)

computer = getenv('COMPUTERNAME');
switch computer
    case 'DF-LAB-WS38'
        RaidDir = 'F:\Data\2P_data\';
        ServerDir = 'Z:\Juliane\Data\2P_data\';
    case 'DF-LAB-WS22'
        RaidDir = 'C:\Data\2P_data\';
        ServerDir = 'Z:\Juliane\Data\2P_data\';
end

for i = ind
    %set folder
    if analysisParams.server
        baseDir = [ServerDir char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    else 
        baseDir = [RaidDir char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
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
    
    %make ROIs if there is no ROI file yet
    ROIFiles = dir([ROIDir 'ROIs.mat']);
    if isempty(ROIFiles)
        FijiSpineROIFct(analysisParams.server, char(exp_info.animal{i}), char(exp_info.name{i}), exp_info.vol{i})
    end
end