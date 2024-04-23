function D = angdiffVNew(a)
D=min(abs(cat(3,a,a+(2*pi),a-(2*pi))),[],3);