function axon_divider(name, folder, server, level)

if server
    drive       = 'Z:\Juliane\';
else
    drive           = 'F:\';
end

baseDirectory   = [drive 'Data\2P_data\'];
if level
    dirFileName         = [baseDirectory name '\' folder];
    DirROIName          = [baseDirectory name '\' folder '\suite2p\combined\'];
else
    dirFileName         = [baseDirectory name '\' folder];
    DirROIName          = [baseDirectory name '\' folder '\suite2p\plane0\'];
end

files = dir(dirFileName);
files =files(4:end);

for k = 1:length(files)
    expnumbs(k) = regexp({files(k).name},'t[0-9]*', 'match', 'once')';
    filesize(k) = size(imfinfo([dirFileName filesep '\' files(k).name]),1);
end

[idxs, ind] = unique(expnumbs);

for id = 1:length(idxs)-1
    exp_length{id} = sum(filesize(ind(id):ind(id+1)-1));
    exp_number{id} = expnumbs(ind(id));
end

exp_length{id+1} = sum(filesize(ind(id+1):end));
exp_number{id+1} = expnumbs(ind(id+1));

ROIfile = load([DirROIName filesep 'data.mat']);
ROIdata = ROIfile.data.roi;
data = struct;
data.template = ROIfile.data.template;

totalSize = size(ROIdata(1).rawF,2);
sizeFiles = 0;
for experiment = 1:length(exp_length)
    sizeFiles = sizeFiles+exp_length{experiment};
end

if totalSize == sizeFiles
     ind = 0;
     for exp = 1:length(exp_length)
        if level == 1
            saveDirectory = [baseDirectory name '\' char(exp_number{exp}) filesep 'suite2p' filesep 'combined'];
        else
            saveDirectory = [baseDirectory name '\' char(exp_number{exp}) filesep 'suite2p' filesep 'plane0'];
        end
        mkdir(saveDirectory)
        for m = 1:size(ROIdata,2)
            data.roi(m).xPos = ROIdata(m).xPos;
            data.roi(m).yPos = ROIdata(m).yPos;
            data.roi(m).mask = ROIdata(m).mask;
            data.roi(m).name = ROIdata(m).name;
            data.roi(m).rawF = ROIdata(m).rawF(:,ind+1:ind+exp_length{exp});
        end
        
%         if level == 1
%             Suite2p.F = F_all(:,ind+1:ind+ceil(exp_length{exp}/5));
%             Suite2p.Fneu = Fneu_all(:,ind+1:ind+ceil(exp_length{exp}/5));
%             Suite2p.spks = spks_all(:,ind+1:ind+ceil(exp_length{exp}/5));
%             save([saveDirectory filesep 'Fall.mat'], 'Suite2p')
%             disp(['saving experiment ' char(exp_number{exp})])
%             ind = ind + ceil(exp_length{exp}/5);
%         else
%             Suite2p.F = F_all(:,ind+1:ind+exp_length{exp});
%             Suite2p.Fneu = Fneu_all(:,ind+1:ind+exp_length{exp});
%             Suite2p.spks = spks_all(:,ind+1:ind+exp_length{exp});
        save([saveDirectory filesep 'data.mat'], 'data') %
        disp(['saving experiment ' char(exp_number{exp})]) %             
        ind = ind + exp_length{exp};
     end
else
     disp('Error: Not correct number of frames')
end