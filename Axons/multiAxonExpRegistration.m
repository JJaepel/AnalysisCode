function multiAxonExpRegistration(animal, folder, server, level)

%% switch board for variables
plotROIs = 0;

%% define folders
if server == 1
    basedir = 'Z:\Juliane\Data\';
else
    basedir = 'F:\Data\';
end
Sp2Directory = [basedir 'Spike2Data\'];
TwoPDirectory  = [basedir '2P_data\' animal filesep folder filesep];
saveDirectory = [basedir '\ImageAnalysis\' animal filesep folder filesep];
ROIsaveDirectory = [saveDirectory 'ROIs' filesep];
ROIRespsaveDirectory = [saveDirectory 'ROIs_Responsive' filesep];
ROINonRespsaveDirectory = [saveDirectory 'ROIs_Nonresponsive' filesep];
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end
if ~exist(ROIsaveDirectory, 'dir')
    mkdir(ROIsaveDirectory);  
end
if ~exist(ROIRespsaveDirectory, 'dir')
    mkdir(ROIRespsaveDirectory);  
end
if ~exist(ROINonRespsaveDirectory, 'dir')
    mkdir(ROINonRespsaveDirectory);  
end

%% extract time series from suite2p registration
analysisParams.animal = animal;
analysisParams.name = folder;
analysisParams.server = server;
analysisParams.level = level;
Suite2pAxonExtractorFct(analysisParams)
axon_divider(animal, folder, server, level);

%% get all experiment ids, assumption that all spikeids = expids
files = dir(TwoPDirectory);
files =files(4:end);

for k = 1:length(files)
    expnumbs(k) = regexp({files(k).name},'t[0-9]*', 'match', 'once')';
end

[idxs, ~] = unique(expnumbs);

%% define stimtypes for all experiments in that folder and load data
for i = 1:length(idxs)
    exp(i) = char(idxs{i});
    Sp2Directory = ([basedir 'Spike2Data\' animal filesep exp(i) filesep]);
    files=dir(strcat(Sp2Directory, '\*.py'));
    files={files.name};
    for f = 1:length(files)
        if isempty(strfind(char(files{f}), 'serialTriggerDaqOut'))
            StimType(i)= strrep(char(files{f}), '.py', '');
        end
    end

    try
        %if reload the data, perform the stimulus analysis, but also
        %throw error if it is not one of the predefined stimulus types
        switch StimType(i)
            case 'driftingGrating'
                if reload_data
                    Ori_Grating_S2p_Axon(animal, exp(i), exp(i), folder, plotROIs, server) 
                end
            case 'driftingGrating_ori_sf'
                if reload_data
                    Ori_Sf_S2p_Axon(animal, expt_id, sp2id, name, reloadData, plotROIs, server)
                end
%                 case 'driftingGrating_ori_tf'
%                     if reload_data
%                     end
            otherwise
                error('Stimtype not valid for the analysis')
        end
        %load the analysis file into the master file
        datapath = [savedir animal filesep exp_numbers{exp} filesep];
        filename = [datapath filesep '*ana.mat'];
        files = dir(filename);
        master{exp} = load(fullfile(datapath, files(1).name), 'data', 'metadata', 'analysis');
        disp(['loaded Data from exp ' num2str(exp_numbers{exp})])
    catch
        disp(['Loading next experiment'])
    end       
end
%
