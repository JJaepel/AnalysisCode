function getSubcellularRegistration(analysisParams, expInfo, ind)

%This function is the registration for all 2p experiments. It will only be
%performed if there are not registration files already existing or if
%analysisParams.reregister = 1.
%
%Input:
%- analysisParams: Info about the parameters of the analysis, like whether
%to work on the network or the drive, whether to reregister, etc.
%- expInfo: information about all the experiments
%- ind: which ones of the experiments do you want to write
%
%Output:
%- starts the registration, which will produce the registered tiffs, as
%well as the projection in the given folder

% let's check which computer this is being run on to get the RaidDir
computer = getenv('COMPUTERNAME');
switch computer
    case 'DF-LAB-WS38' %Julianes computer
        RaidDir = 'F:\Data\2P_data\';
        ServerDir = 'Z:\Juliane\Data\2P_data\';
    case 'DF-LAB-WS40' %Hannah's computer
        RaidDir = 'C:\Data\2P_data\';
        ServerDir = 'Z:\Juliane\Data\2P_data\';
end

%Now go through all the experiments that you want to work on
for i = ind
    %set folder
    if analysisParams.server
        baseDir = [ServerDir char(expInfo.animal{i}) '\' char(expInfo.name{i}) '\Registered\'];
    else 
        baseDir = [RaidDir char(expInfo.animal{i}) '\' char(expInfo.name{i}) '\Registered\'];
    end
    
    %if it is not registered, register the data
    if ~exist(baseDir, 'dir') || analysisParams.reregister
      subcellularRegistration(analysisParams.server, char(expInfo.animal{i}), char(expInfo.name{i}),expInfo.vol{i})
    end

end