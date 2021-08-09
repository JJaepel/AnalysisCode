function [ out ] = vonMisesFunction(kappa, mu, x )
%VONMISESFUNCTION Summary of this function goes here
%   Detailed explanation goes here

    out = exp(kappa * cos (x-mu))/ (2 * pi * besseli(0,kappa));
end
