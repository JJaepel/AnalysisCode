function y = percentileFilt1(x,percentile,n,blksz,DIM)
    narginchk(1,5);
    if nargin < 4, blksz = []; end
    if nargin < 5, DIM = []; end

    % Check the input data type. Single precision is not supported.
    % try
    %     chkinputdatatype(x,n,blksz,DIM);
    % catch ME
    %     throwAsCaller(ME);
    % end

    % Check if the input arguments are valid
    if isempty(n)
      n = 3;
    end

    if ~isempty(DIM) && DIM > ndims(x)
        error(message('signal:medfilt1:InvalidDimensions'))
    end

    % Reshape x into the right dimension.
    if isempty(DIM)
        % Work along the first non-singleton dimension
        [x, nshifts] = shiftdim(x);
    else
        % Put DIM in the first (row) dimension (this matches the order 
        % that the built-in filter function uses)
        perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
        x = permute(x,perm);
    end

    % Verify that the block size is valid.
    siz = size(x);
    if isempty(blksz)
        blksz = siz(1); % siz(1) is the number of rows of x (default)
    else
        blksz = blksz(:);
    end

    % Initialize y with the correct dimension
    y = zeros(siz); 

    % Call medfilt1D (vector)
    for i = 1:prod(siz(2:end))
        y(:,i) = prctilefilt1d(x(:,i),n,blksz,percentile);
    end

    % Convert y to the original shape of x
    if isempty(DIM)
        y = shiftdim(y, -nshifts);
    else
        y = ipermute(y,perm);
    end
end