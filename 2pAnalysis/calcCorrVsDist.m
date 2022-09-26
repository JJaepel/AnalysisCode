function [PopAnalysis] = calcCorrVsDist(PopAnalysis, analysisParams)

%% 1.) do it for all the data
%load the data and remove double pairs
L = triu(PopAnalysis.(analysisParams.field).distROIs);
L = reshape(L,1,length(L)^2);
L(L == 0) = NaN;
distROIs = L((~isnan(L)));

corrROIs = PopAnalysis.(analysisParams.field).corrROIs;
corrROIs = reshape(corrROIs,1,length(corrROIs)^2);
corrROIs = corrROIs((~isnan(L)));

%bin the data
maxDis = floor(max(distROIs));
edges = linspace(0, maxDis, 15);
[n, edges, binDist] = histcounts([distROIs],edges);
edge = edges(1:end-1)+25;
meanCorr = zeros(length(edge),1);
semCorr = zeros(length(edge),1);
for bin = 1:length(edge)
    meanCorr(bin) = nanmean(corrROIs(binDist == bin));
    semCorr(bin) = nanstd(corrROIs(binDist == bin))/sqrt(n(bin));
end

%shuffle and bin the data
meanCorrShuffle = zeros(analysisParams.shufflenum,bin);
semCorrShuffle = zeros(analysisParams.shufflenum,bin);
for rep = 1:analysisParams.shufflenum
    %shuffle the orientation preference differences onto the existing distances
    corrROIsShuffle = corrROIs(randperm(length(corrROIs)));
    for bin = 1:length(edge)
        meanCorrShuffle(rep,bin) = nanmean(corrROIsShuffle(binDist == bin));
        semCorrShuffle(rep,bin) = nanstd(corrROIsShuffle(binDist == bin))/sqrt(n(bin));
    end
end

%save in the structure
PopAnalysis.(analysisParams.field).edges = edges;
PopAnalysis.(analysisParams.field).meanCorr = meanCorr;
PopAnalysis.(analysisParams.field).semCorr = semCorr;
PopAnalysis.(analysisParams.field).meanCorrShuffle = meanCorrShuffle;
PopAnalysis.(analysisParams.field).semCorrShuffle = semCorrShuffle;

%% 2. Search for distances within certain limits
%look for long and short distances
shortDistPairs = find(distROIs < 100);
morethanDistPairs = find(distROIs >300);
removeDistPairs = find(distROIs < 400);
longDistPairs = intersect(morethanDistPairs, removeDistPairs);


%select into short and long distance correlations
PopAnalysis.(analysisParams.field).shortDistanceCorr = corrROIs(shortDistPairs);
PopAnalysis.(analysisParams.field).longDistanceCorr = corrROIs(longDistPairs);
