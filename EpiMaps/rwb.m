function map = rwb(m)
    % Red-White-Blue LUT

    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    h = (0:m-1)'/max(m,1);
    if isempty(h)
      map = [];
    else
      map = zeros(m,3);
      for currentColorSide = 1:2
          switch currentColorSide
              case 1 % Blue
                  index = 1:floor(m/2);
                  map(index,3)     = 1;
                  map(index,[1 2]) = repmat(linspace(0,1,length(index))',[1 2]);
              case 2 % Red
                  index = (floor(m/2)+1):m;
                  map(index,1)     = 1;
                  map(index,[2 3]) = repmat(linspace(1,0,length(index))',[1 2]);
          end
      end
    end
end