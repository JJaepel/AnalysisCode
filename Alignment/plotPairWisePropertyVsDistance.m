function [binMeans, binSEM, binCenters] = plotPairWisePropertyVsDistance(distMatrix,deltaMatrix,varargin)

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
% - binMeans = mean delta binned by underlying distances
% - binStd = std delta binned by underlying distances
% - binCenter = what are the centers of the bins
% (-figure: delta bins vs. edges)


if nargin > 2
    figureNr = varargin{1};
else
    figureNr = 1000;
end
if nargin > 3
    edges = varargin{2};
else
    edges = linspace(0,30,13);
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
    %make a new matrix that is specific for the type and multiple it with
    %the orginal matrices
    inputSpines = find(logInputs == 1);
    nonInputSpines = find(logInputs == 0);
    InputMatrix = NaN(size(deltaMatrix,1), size(deltaMatrix,2));
    
    %get the colors
    coc_prop = cbrewer('div', 'RdGy', 15);
    switch type
        case 1 %inputs only
            %set everything in input rows and columns to 1
            InputMatrix(inputSpines,inputSpines) = 1;
            
            %multiply with original matrices
            deltaMatrix = deltaMatrix .* InputMatrix;
            distMatrix = distMatrix .* InputMatrix;
            
            %logical matrix with only the upper triangle without the diagonal = true
            m = triu(true(size(distMatrix,2),size(distMatrix,2)),1); 
            
            %get the color
            coc = coc_prop(3,:);
        case 2 %non-inputs only
            %set everything in input rows and columns to 1
            InputMatrix(nonInputSpines,nonInputSpines) = 1;
            
            %multiply with original matrices
            deltaMatrix = deltaMatrix .* InputMatrix;
            distMatrix = distMatrix .* InputMatrix;
            
            %logical matrix with only the upper triangle without the diagonal = true
            m = triu(true(size(distMatrix,2),size(distMatrix,2)),1); 
            
            %get the color
            coc = coc_prop(13,:);
        case 3 %mixed of inputs
            %set everyting in input rows and non-input spines to 1
            InputMatrix(inputSpines,nonInputSpines) = 1;
            
            %multiply with original matrices
            deltaMatrix = deltaMatrix .* InputMatrix;
            distMatrix = distMatrix .* InputMatrix;
            
            %all pairs are valid as there are there only once, so take all of those
            m = true(size(distMatrix,2),size(distMatrix,2)); 
            
            %get the color
            coc = coc_prop(7,:);
        otherwise
            error('Correct the type');            
    end
else
    %logical matrix with only the upper triangle without the diagonal = true
    m = triu(true(size(distMatrix,2),size(distMatrix,2)),1); 
    
    %set the color
    coc = 'blue';
end

%reshape according to m
deltaReshaped = deltaMatrix(m);
distReshape = distMatrix(m);
%remove the NaNs based on delta
nonNanDelta = find(~isnan(deltaReshaped));
deltaVector= deltaReshaped(nonNanDelta);
distVector = distReshape(nonNanDelta);

%% 2.) Bin the delta vector based on distVector 

%sort the distVector into bins according to the edges
bins = discretize(distVector, edges, 'IncludedEdge', 'left');

%if the bins contain a NaN element (outside of the edges), remove those
%elements from both the bins and the deltaVector
nanElements = find(isnan(bins));
if ~isempty(nanElements)
    bins(nanElements) = [];
    deltaVector(nanElements) = [];
end

%used these bins on the deltaVector and calculate mean and sem
binMeans = groupsummary(deltaVector, bins, @mean); %mean
binSEM = groupsummary(deltaVector, bins, @std)./ groupsummary(deltaVector, bins, @numel);

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

%% 3.) Plot

figure(figureNr)
h = errorbar(binCenters,binMeans(1:end-1),binSEM(1:end-1), 'o-');
%color it according to type
h.MarkerFaceColor = coc;
h.MarkerEdgeColor = coc;
h.Color = coc;
xlabel('Distance between spines in \mum')
box off
hold on
set(gcf, 'color', 'w');

