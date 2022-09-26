function [shortDistanceCorr, longDistanceCorr] = cellCorrelationDistance(analysis, data, metadata, field, umperpixel)
% Computs correlation of cells independent of the stimulus within certain
% distances
% Input:
% - metadata: structure containing StimParams information
% - data: containing ROI information, such as ROI position
% - analysis: structure containing for each ROI (analysis.roi) the
% stim response trace and whether they are significant responders
% - field: which kind of trace to use
%
% Ouput:
% - shortDistanceCorr: Correlation of cells within 100 um
% - longDistanceCorr: Correlation of cells more than 300 um apart

stimWindow=(analysis.(field).windowStart: analysis.(field).windowStop);

%% 1. Make structure containing distance between cells and correlation of responses
corrCoeff = zeros(length(data.roi),length(data.roi));
distROIs = zeros(length(data.roi),length(data.roi));

for A = 1:length(data.roi)
    %get normalized Response cell A
    Response = analysis.(field).roi(A).stimResponseTrace(1:end-1, :, :);
    Response = median(Response(:,:,stimWindow),3); %average over the stimWindow
    Response = reshape(Response',[(metadata.StimParams.uniqStims-1) * metadata.StimParams.numTrials,1]);
    
    %normalize Response
    xmin=min(Response(:));
    xmax=max(Response(:));
    normResponseA = (Response-xmin)/(xmax-xmin);
    
    for B = 1:length(data.roi)
        %get normalized Response cell B
        Response = analysis.(field).roi(B).stimResponseTrace(1:end-1, :, :);
        Response = median(Response(:,:,stimWindow),3); %average over the stimWindow
        Response = reshape(Response',[(metadata.StimParams.uniqStims-1) * metadata.StimParams.numTrials,1]);

        %normalize Response
        xmin=min(Response(:));
        xmax=max(Response(:));
        normResponseB = (Response-xmin)/(xmax-xmin);
        
        %caclulate correlation
        corrCoeff(A,B) = corr(normResponseA,normResponseB);
        
        %calculate distance
        dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
        dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
        distROIs(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;     
    end
end

%% 2. Search for distances within certain limits
%reshape both vectors
corrCoeff = reshape(corrCoeff, length(data.roi)*length(data.roi),1);
distROIs = reshape(distROIs, length(data.roi)*length(data.roi),1);

%look for long and short distances
shortDistPairs = find(distROIs < 100);
morethanDistPairs = find(distROIs >300);
removeDistPairs = find(distROIs < 400);
longDistPairs = intersect(morethanDistPairs, removeDistPairs);


%select into short and long distance correlations
shortDistanceCorr = corrCoeff(shortDistPairs);
longDistanceCorr = corrCoeff(longDistPairs);
end
