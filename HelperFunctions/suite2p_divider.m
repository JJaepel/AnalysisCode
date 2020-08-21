clear all
animal = 'F2425_2020-03-05\';
TwoPhontondir = 'Z:\Juliane\Data\2P_Data\';
base2pDirectory= [TwoPhontondir animal];
expt_id = 'Yellow';
level = 0;

files = dir([base2pDirectory filesep expt_id]);
files =files(4:end);

for k = 1:length(files)
    expnumbs(k) = regexp({files(k).name},'t[0-9]*', 'match', 'once')';
    filesize(k) = size(imfinfo([base2pDirectory filesep expt_id filesep '\' files(k).name]),1);
end

[idxs, ind] = unique(expnumbs);

for id = 1:length(idxs)-1
    exp_length{id} = sum(filesize(ind(id):ind(id+1)-1));
    exp_number{id} = expnumbs(ind(id));
end

exp_length{id+1} = sum(filesize(ind(id+1):end));
exp_number{id+1} = expnumbs(ind(id+1));


%cd Z:\\
if level == 1
    Suite2pDir = [base2pDirectory filesep expt_id filesep 'suite2p' filesep 'combined' filesep];
else
    Suite2pDir = [base2pDirectory filesep expt_id filesep 'suite2p' filesep 'plane0' filesep];
end
Suite2pFile = [Suite2pDir 'Fall.mat'];
s2p = load(Suite2pFile);
disp('loaded Suite2pFile');
F_all = s2p.F;
Fneu_all = s2p.Fneu;
spks_all = s2p.spks;
Suite2p.stat = s2p.stat;
Suite2p.ops = s2p.ops;
Suite2p.iscell = s2p.iscell;
clear s2p
totalSize = size(F_all,2);

sizeFiles = 0;
for experiment = 1:length(exp_length)
    sizeFiles = sizeFiles+exp_length{experiment};
end

if totalSize == sizeFiles
    ind = 0;
    for exp = 1:length(exp_length)
        if level == 1
            saveDirectory = [base2pDirectory char(exp_number{exp}) filesep 'suite2p' filesep 'combined'];
        else
            saveDirectory = [base2pDirectory char(exp_number{exp}) filesep 'suite2p' filesep 'plane0'];
        end
        mkdir(saveDirectory)
        if level == 1
            Suite2p.F = F_all(:,ind+1:ind+ceil(exp_length{exp}/5));
            Suite2p.Fneu = Fneu_all(:,ind+1:ind+ceil(exp_length{exp}/5));
            Suite2p.spks = spks_all(:,ind+1:ind+ceil(exp_length{exp}/5));
            save([saveDirectory filesep 'Fall.mat'], 'Suite2p')
            disp(['saving experiment ' char(exp_number{exp})])
            ind = ind + ceil(exp_length{exp}/5);
        else
            Suite2p.F = F_all(:,ind+1:ind+exp_length{exp});
            Suite2p.Fneu = Fneu_all(:,ind+1:ind+exp_length{exp});
            Suite2p.spks = spks_all(:,ind+1:ind+exp_length{exp});
            save([saveDirectory filesep 'Fall.mat'], 'Suite2p')
            disp(['saving experiment ' char(exp_number{exp})])
            ind = ind + exp_length{exp};
        end
    end
else
    disp('Error: Not correct number of frames')
end