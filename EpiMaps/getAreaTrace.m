function trace = getAreaTrace(AreaData,mask)

xVals = find(sum(mask,1)>0);
yVals = find(sum(mask,2)>0);
xRange = min(xVals):max(xVals);
yRange = min(yVals):max(yVals);
M = repmat(int16(mask(yRange,xRange)), [1 1 size(AreaData,3)] ); %make a mask over the whole range
trace = squeeze(squeeze(sum( sum( int16(AreaData(yRange,xRange,:)) .* M, 1 ), 2 ))) ./ nnz(mask);

