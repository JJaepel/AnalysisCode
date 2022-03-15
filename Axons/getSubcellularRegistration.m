function getSubcellularRegistration(analysisParams, exp_info, ind)

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
        baseDir = [ServerDirchar(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    else 
        baseDir = [RaidDir char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    end
    
        %if it is not registered, register the data
    if ~exist(baseDir, 'dir')
      subcellularRegistration(analysisParams.server, char(exp_info.animal{i}), char(exp_info.name{i}),exp_info.vol{i})
    end

end