function [ out ] = vonMisesLinSum(x, tdata)
%VONMISESFIT Linear Sum of two vonMises functions that's constrained to be
%pi radians apart
%   x is a vector of the parameters of the linear sum and vonMises
%   parameters
%   tdata is the x range over which to compute the vonMises

if sum(isnan(x)) > 0 || length(x) <6
    out = NaN * ones(size(tdata));
    return
end

A= x(1); % first direction
B=x(2);  % 2nd direction
D=x(3);  % a dc component
mu1=x(4);
kappa1=x(5);
kappa2=x(6);



mu2=mod(mu1 + pi, 2*pi);

out = A*vonMisesFunction(kappa1, mu1, tdata ) + B*vonMisesFunction(kappa2, mu2, tdata )+D;

end