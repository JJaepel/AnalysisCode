function dff = NormalizeStack(metadata,rawF, binningFactor, percentileRank)
% The binning factor is used for subsampling the Fo. The approximation is 
% usually good and dramatically speeds up the computation
% The percentile cutoff for the rank order filter

if nargin < 3
    binningFactor = 10;
end
if nargin < 4
    percentileRank = 25;
end

stackSize = size(rawF);
imgStack = single(reshape(rawF,[stackSize(1)*stackSize(2) stackSize(3)]));
fps= metadata.Imaging.rate; 

parfor index = 1:(stackSize(1)*stackSize(2))
    rawSignalTrace = imgStack(index,:)';
    selectedFrames = 1:binningFactor:length(rawSignalTrace);
    baselineTrace  = baselinePercentileFilterEpi(rawSignalTrace(selectedFrames),fps/binningFactor,60,percentileRank);
    baselineTrace  = interp1(selectedFrames,baselineTrace,1:length(rawSignalTrace));
    baselineTrace  = reshape(baselineTrace,size(rawSignalTrace));
    imgStack(index,:) = (rawSignalTrace - baselineTrace)./baselineTrace;
end
dff = reshape(imgStack, stackSize);
    
end