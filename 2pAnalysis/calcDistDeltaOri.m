function [PopAnalysis] = calcDistDeltaOri(analysisParams, data, analysis, umperpixel, level)    

%only do the distance calculation between ROIs that are active and ORi
%selective 
ori_sel = find([analysis.(analysisParams.field).roi.OSIFit] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1);

%calculate distances between all cell pairs
if isempty(ori_sel)
    edges = [];
    meanDeltaOri = [];
    semDeltaOri = [];
    meanDeltaOriShuffle = [];
    semDeltaOriShuffle = [];
else
    distRois = zeros(length(ori_sel), length(ori_sel));
    deltaOri = zeros(length(ori_sel), length(ori_sel));
    for A = 1:length(ori_sel)
        for B = 1:length(ori_sel)
            if level == 1
                if data.roi(A).plane == data.roi(B).plane
                    %calculate the distance
                    dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
                    dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
                    distRois(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel; %change it to um

                    %calculate the deltaOri preference
                    deltaOri(A,B) = abs(analysis.(analysisParams.field).roi(A).preferredOrientation - analysis.(analysisParams.field).roi(B).preferredOrientation);
                    if deltaOri(A,B) > 90
                        deltaOri(A,B) = 180 - deltaOri(A,B);
                    end
                else
                    distRois(A,B) = NaN;
                    deltaOri(A,B) = NaN; 
                end
            else
                dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
                dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
                distRois(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
                deltaOri(A,B) = abs(analysis.(analysisParams.field).roi(A).preferredOrientation - analysis.(analysisParams.field).roi(B).preferredOrientation);
                if deltaOri(A,B) > 90
                    deltaOri(A,B) = 180 - deltaOri(A,B);
                end
            end
        end
    end

    %remove double pairs
    L = triu(distRois);
    L = reshape(L,1,length(ori_sel)^2);
    L(L == 0) = NaN;
    distOri = L((~isnan(L)));
    deltaOri = reshape(deltaOri,1,length(ori_sel)^2);
    deltaOri = deltaOri((~isnan(L)));

    %save in one structure
    PopAnalysis.(analysisParams.field).distOriROIs = distRois;
    for i = 1:length(deltaOri)
        PopAnalysis.(analysisParams.field).ori_pair(i).distance = distOri(i);
        PopAnalysis.(analysisParams.field).ori_pair(i).deltaOri = deltaOri(i);
    end

    %bin the data
    maxDis = floor(max([PopAnalysis.(analysisParams.field).ori_pair.distance]));
    edges = linspace(0, maxDis, 15);
    [n, edges, binDist] = histcounts([PopAnalysis.(analysisParams.field).ori_pair.distance],edges);
    edge = edges(1:end-1)+25;
    meanDeltaOri = zeros(length(edge),1);
    semDeltaOri = zeros(length(edge),1);
    for bin = 1:length(edge)
        meanDeltaOri(bin) = nanmean(deltaOri(binDist == bin));
        semDeltaOri(bin) = nanstd(deltaOri(binDist == bin))/sqrt(n(bin));
    end

    %shuffle and bin the data
    meanDeltaOriShuffle = zeros(analysisParams.shufflenum,bin);
    semDeltaOriShuffle = zeros(analysisParams.shufflenum,bin);
    for rep = 1:analysisParams.shufflenum
        %shuffle the orientation preference differences onto the existing distances
        deltaOriShuffle = deltaOri(randperm(length(deltaOri)));
        for bin = 1:length(edge)
            meanDeltaOriShuffle(rep,bin) = nanmean(deltaOriShuffle(binDist == bin));
            semDeltaOriShuffle(rep,bin) = nanstd(deltaOriShuffle(binDist == bin))/sqrt(n(bin));
        end
    end
end

%save in the structure
PopAnalysis.(analysisParams.field).edges = edges;
PopAnalysis.(analysisParams.field).meanDeltaOri = meanDeltaOri;
PopAnalysis.(analysisParams.field).semDeltaOri = semDeltaOri;
PopAnalysis.(analysisParams.field).meanDeltaOriShuffle = meanDeltaOriShuffle;
PopAnalysis.(analysisParams.field).semDeltaOriShuffle = semDeltaOriShuffle;
