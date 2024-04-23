function [binMeans, binSEM, binCenters] = plotBinnedPairWisePropertyVsDistance(pwDistance,pwDelta,cocCode, varargin)

% Bins pairwise measurements based on the distances between the pairs and
% plots it as delta vs. pair-wise distance
%
% Input:
% - DistMatrix: spine x spine matrix of the distance between them
% - deltaMatrix: spine x spine matrix of their delta of a specific property
% - (figureNr):what is the figure nuber
% - (edges): what are the bin edges?

% 
% Output: 
% - binMeans = mean delta binned by underlying distances
% - binStd = std delta binned by underlying distances
% - binCenter = what are the centers of the bins
% (-figure: delta bins vs. edges)

% 
if nargin > 3
    edges = varargin{2};
else
    edges = linspace(0,30,13);
end


%% 1.) Bin the delta vector based on distVector 

%sort the distVector into bins according to the edges
bins = discretize(pwDistance,edges, 'IncludedEdge', 'left');

%if the bins contain a NaN element (outside of the edges), remove those
%elements from both the bins and the deltaVector
nanElements = find(isnan(bins));
if ~isempty(nanElements)
    bins(nanElements) = [];
    pwDelta(nanElements) = [];
    %pwDistance(nanElements) = [];
end

%used these bins on the deltaVector and calculate mean and sem
binMeans = groupsummary(pwDelta', bins', @mean); %mean
binSEM = groupsummary(pwDelta', bins', @std)./ groupsummary(pwDelta', bins', @numel);

%if there are no values for some bins, make sure we have nans for those
if length(binMeans) < 13
    %make a temp vector
    temp = NaN(13,1);
    temp(unique(bins)) = binMeans; %add the binMeans to the temp at the corresponding bins
    binMeans = temp;
    %repeat for SEM
    temp = NaN(13,1);
    temp(unique(bins)) = binSEM; %add the binMeans to the temp at the corresponding bins
    binSEM = temp;
end

%get the bin center
diffEdges = diff(edges);
binCenters = edges(1:end-1) + diffEdges/2;

%% 2.) Plot

%make sure the color is set right, otherwise set it to grey
if ~isequal(size(cocCode), [1, 3])
    cocCode = [.7, .7, .7];
end

h = errorbar(binCenters',binMeans(1:end-1),binSEM(1:end-1), 'o-');
%color it according to type
h.MarkerFaceColor = cocCode;
h.MarkerEdgeColor = cocCode;
h.Color = cocCode;
xlabel('Distance between spines in \mum')
box off
hold on
set(gcf, 'color', 'w');

