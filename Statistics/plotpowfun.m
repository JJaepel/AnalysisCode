function f = plotpowfun(d,alpha,powers,indices)
close all;
f = plot(powers,indices,'o');
xl = xlabel('Power');
set(xl,'FontWeight','bold');
yl = ylabel('N');
set(yl,'FontWeight','bold');
titletext = sprintf('Power for detecting effect size d = %6.4f at alpha = %1.4f with N = %u:%u',d,alpha,min(indices),max(indices));
t = title(titletext);
set(t,'FontWeight','bold');