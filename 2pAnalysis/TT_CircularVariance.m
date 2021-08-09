function [CircularVariance, ResultantLength, ResultantAngle] = TT_CircularVariance( TuningCurve )
% function [CircularVariance, ResultantLength, ResultantAngle] = ...
%                         tpd_circularVariance( TuningCurve )
%
% Calculates circular variance, resultant length and angle of the 
% resultant based on the tuning curve. Assumes a full 360 degree measured
% tuningcurve.
%
% Input:
% - TuningCurve: Array with mean response per direction stimulus
%
% Ouput:
% - CircularVariance: Circular variance of the tuning curve
% - ResultantLength:  Length of the tuning curve resultant
% - ResultantAngle:   Angle of the tuning curve resultant
%
% Written by Pieter Goltstein
% Version 1.0: July 25th, 2011
%

    % get number of datapoints on full circle
    [a, nDataPoints] = size(TuningCurve);
    
    % remove negatives from tuning curve
    TuningCurve( TuningCurve < 0 ) = 0;

    % convert theta to complex plane
    Directions = (1:nDataPoints) * (360/nDataPoints);
    Theta = (Directions/360) * 2 * pi;
    Theta = exp(1i*Theta);  

    % calculate complex resultant
    Resultant = sum(TuningCurve.*Theta)/sum(TuningCurve);

    % calculate resultant length
    ResultantLength = abs(Resultant);

    % calculate resultant angle
    ResultantAngle = (angle(Resultant)/(2*pi))*360;

    %imaginary numbers give negative angles, so convert to positive ones
    while ResultantAngle < 0
        ResultantAngle = ResultantAngle + 360;
    end

    % calculates circular variance
    CircularVariance = 1 - abs(Resultant);

end
