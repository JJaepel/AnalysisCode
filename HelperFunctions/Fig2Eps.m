function Fig2Eps(strng)
if nargin < 1
    DI = dir('*.fig');
else 
    DI = dir(['*',strng,'*.fig']);
end
for i = 1:length(DI)
    nm = DI(i).name;
    if exist(fullfile(cd,strcat(nm(1:end-4),'.eps')),'file')~=2
        openfig(nm);
        fname = strcat(nm(1:end-4),'.eps');
        disp(fname)
        set(gcf,'Renderer','painters')
        saveas(gcf,fname,'epsc')
        close
    end
end
end
