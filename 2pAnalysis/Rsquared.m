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