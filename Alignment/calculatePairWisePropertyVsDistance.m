function [deltaDistBinned, binEdges] = plotPairWisePropertyVsDistance(distMatrix,deltaMatrix,varargin)

% Bins pairwise measurements based on the distances between the pairs and
% plots it as delta vs. pair-wise distance
%
% Input:
% - DistMatrix: spine x spine matrix of the distance between them
% - deltaMatrix: spine x spine matrix of their delta of a specific property
% - (type): Do we look at all pairs(0), inputs-only (1), non-inputs (2) or
% mixed pairs (3)
% - (edges): what are the bin edges?
% - (figureNr):what is the figure nuber
% 
% Output: 
% - deltaDistBinned = delta Data binned by underlying distances
% - binEdges = what are the edges of the bins
% (-figure: delta bins vs. edges)


if nargin > 2
    edges = varargin{1};
else
    edges = linspace(0,30,13);
end
if nargin > 3
    figureNr = varargin{2};
else
    figureNr = 1000;
end
if nargin > 4
    logInputs = varargin{3};
else
    logInputs = [];
end
if nargin > 3
    type = varargin{4};
else
    type = 0;
end


%% 1.) Reshape the pairs into a numPairsx1 vector

%do we look at all pairs or specific ones?
if type > 0 
    inputSpines = find(logInputs == 1);
    nonInputSpines = find(logInputs == 0);
else
    %logical matrix with only the upper triangle without the diagonal = true
    m = triu(true(size(confData.ROIs,2),size(confData.ROIs,2)),1); 
    deltaReshaped = deltaMatrix(m);
    distReshape = distMatrix(m);
    %remove the NaNs based on delta
    nonNanDelta = find(~isnan(deltaReshaped));
    deltaVector= deltaReshaped(nonNanDelta);
    distVector = distReshape(nonNanDelta);
end

%% 2.) Bin the delta vector based on distVector 

%sort the distVector into bins according to the edges
bins = discretize(distVector, edges, 'IncludedEdge', 'left');

%used these bins on the deltaVector and calculate mean and sem
binMeans = groupsummary(deltaVector, bins, @mean); %mean
binStd = groupsummary(deltaVector, bins, @std); %std

%get the bin center
diffEdges = diff(edges);
binCenters = edges(1:end-1) + diffEdges/2;

%% 3.) Plot

figure(figureNr)


