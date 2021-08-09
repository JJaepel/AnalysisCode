function [prefOri, coeffOr, rsqPOr, rsqOr] = ComputePreferredOrientations(medResponse, theta)
% Calculates preferred orientation after fitting the data to a gaussian
% curve
%
% Input:
% - medResponse: a 1xorientation array containing the average response
% - theta: 
%
% Ouput:
% - prefOri: preferred orientation
% - coeffOr:  
% - rsqPOr: Rsquare for medResponse and optimized curve -> how well is fit
% using all parameters
% - rsqOr: Rsquare for medResponse and vonMises fit curve -> how well is fit
% fitting the data

    theta = theta *2;
    
    %do the fit
    options=optimoptions('lsqcurvefit', 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'Off');
    coeffOr=[NaN, NaN, NaN, NaN];
    [coeffOr, rsqPOr]=fitVonMisesFunction(medResponse, theta, coeffOr, options);
    
    %ca
    prefOri = rad2deg(coeffOr(3)/2);
    
    %plot the data
    hold off;
    plot(rad2deg(theta/2), medResponse, 'bo')
    hold on;
    plot(rad2deg([0:pi/32:pi]), vonMisesFit(coeffOr, [0:pi/16:2*pi]), 'r');
    title([num2str(prefOri),':',num2str(rsqPOr)]);
    drawnow()
    prefOri =mod(prefOri, 180);
    
    %calculate rsqOr
    rsqOr= Rsquared(medResponse, vonMisesFit(coeffOr, theta*2), true);
end

function [bestCoeff, bestrsq]=fitVonMisesFunction(ydata, tdata, bestCoeff, options)
    stims=length(ydata);
    fun=@(x, tdata) vonMisesFit(x,tdata);

    prefdir=mod(tdata(ydata == max(ydata(:))), pi);
    if length(prefdir) >1
        prefdir=prefdir(1);
    end
    %coeff0=[abs(max(ydata(:))), abs(max(ydata(:))), pi/2, pi/stims];
    
    lb=[0,0,0,0];
    ub=[abs(1.5*max(ydata(:))),abs(max(ydata(:))*2),2*pi,2.5*pi];
    coeff0=(ub-lb)./2 + lb;
    rsq=0;
    bestrsq=0;
    
    iterations=0;
    try
    if sum(isnan(lb))==0 && sum(isnan(ub)) ==0
        while rsq < 0.7 && iterations <5
            

            [currentCoeff, ~, ~,exitflag]=lsqcurvefit(fun, coeff0, tdata, ydata, lb, ub, options);

            [rsq, ~]=Rsquared(ydata, vonMisesFit(currentCoeff, tdata), true);

             
            if rsq > bestrsq
                bestCoeff= currentCoeff;
                bestrsq=rsq;
            end
            iterations=iterations+1;
            coeff0= ub-lb .* rand(1) + lb;
        end
        %disp(['Ran ', num2str(iterations),  ' iterations']);
    else
        bestCoeff=[NaN, NaN, NaN, NaN];
    end 
    catch
        bestCoeff=[NaN,NaN,NaN,NaN];
        rsq=0;
    end
        
end
function [ out ] = vonMisesFit(x, tdata)
%Linear Sum of two vonMises functions that's constrained to be
%pi radians apart

%Input:
% - x is a vector of the parameters of the linear sum and vonMises
%   parameters
% - tdata is the x range over which to compute the vonMises

%Output:
% - out: array containing the fit parameters  = [Alpha Beta mu kappa]

    if sum(isnan(x)) > 0 || length(x) <4
        out = NaN * ones(size(tdata));
        return
    end
    A= x(1);
    B= x(2);  % a dc component
    mu1=x(3); % center
    kappa1=x(4); % width
    out = A*vonMisesFunction(kappa1, mu1, tdata ) + B;
end
function [r2 rmse] = Rsquared(y,f,varargin)
    if isempty(varargin); c = true; 
    elseif length(varargin)>1; error 'Too many input arguments';
    elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
    else c = varargin{1}; 
    end

    % Compare inputs
    if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

    % Check for NaN
    tmp = ~or(isnan(y),isnan(f));
    y = y(tmp);
    f = f(tmp);

    if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
    else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
        if r2<0
        % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
            warning('Consider adding a constant term to your model') %#ok<WNTAG>
            r2 = 0;
        end
    end

    rmse = sqrt(mean((y(:) - f(:)).^2));
end
