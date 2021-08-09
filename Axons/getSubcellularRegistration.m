function getSubcellularRegistration(analysisParams, exp_info, ind)

for i = ind
    %set folder
    if analysisParams.server
        baseDir = ['Z:\Juliane\Data\2P_data\' char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    else 
        baseDir = ['F:\Data\2P_data\' char(exp_info.animal{i}) '\' char(exp_info.name{i}) '\Registered\'];
    end
    
        %if it is not registered, register the data
    if ~exist(baseDir, 'dir')
        subcellularRegistration(analysisParams.server, char(exp_info.animal{i}), char(exp_info.name{i}),exp_info.vol{i})
    end
end