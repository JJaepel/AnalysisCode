function [dpower,nreq] = tpower(N, d, alpha, plot, despwr)
%TPOWER    Calculate the power of a t-test.
%   X=TPOWER(N,D,ALPHA) returns the power of a t-test at ALPHA with (N-1)
%   degrees of freedom. For a two-sample test N is the sample size per condition, NOT the total
%   sample size. Total sample size will equal 2*N.
%   
%   If N is an array then X will contain power values for min(N) through
%   max(N). 
%   
%   X=TPOWER(N,D,ALPHA,PLOT) returns power as above. If PLOT equals 1 a plot
%   of sample size against power is also generated. 
%   
%   [X,NREQ]=TPOWER(N,D,ALPHA,PLOT,DESPWR) returns power and plots as above. In addition, the function 
%   returns the minimum sample size needed to obtain a desired level of power. 
dpower = [];
if min(N) < 2
    error('Minimum N must be an integer greater than 1!');
end
% Calculate t_critical for desired alpha
tcrit = tinv( (1-alpha/2) , (N-1) );
% Calculate d_critical
dcrit = tcrit./sqrt(N);
% Standardize d
dstand = (d-dcrit)./(1./sqrt(N));
% Calculate power as probability of dstand
dpower = tcdf(dstand,(N-1));
% Determine whether a plot should be generated
if nargin >= 4
    if plot == 1
        f1 = plotpowfun(d,alpha,dpower,[min(N):max(N)]);    
        shg;
    end   
end
% Determine whether user wants as output the required sample size to reach desired power
if nargin >= 5
    nreq = max(find(((despwr-.01)<dpower)&(dpower<(despwr+.01))));
    if plot == 1
        hold on;
        l1 = line([despwr despwr],[0 nreq]);        
        l2 = line([0 despwr],[nreq nreq]);        
        xlabeltext = sprintf('Power\nRequired N (per condition) = %u',nreq);
        xl = xlabel(xlabeltext);
        set(xl,'FontWeight','bold');
        shg;
    end
end