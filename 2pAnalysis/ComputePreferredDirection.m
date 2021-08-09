function [prefDir] = ComputePreferredDirection(medResponse, theta)
    bestCoeff=[NaN, NaN, NaN, NaN];
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    bestCoeff=fitVonMisesLinSumFunction(medResponse, theta, bestCoeff, options);
    xo=[0:pi/40:2*pi];
    fit=vonMisesLinSum(bestCoeff, xo);
    if bestCoeff(1) > bestCoeff(2)
        prefDir = rad2deg(bestCoeff(4));
    else
        prefDir = mod(rad2deg(bestCoeff(4))+180,360);
    end
end