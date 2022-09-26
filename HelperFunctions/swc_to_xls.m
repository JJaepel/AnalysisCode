
function [savefile] = swc_to_xls(path, filename)

file = [path filesep filename];
corefilename = char(filename);

b = loadfilelist(file);
a = b(4:end);

N = length(corefilename);
corefilename = corefilename(1:N-4);

ending = '.xls';
savefile = [corefilename ending];

c = cellfun(@str2num, a, 'UniformOutput', false);
d = cell2mat(c');

xlswrite(savefile, d);