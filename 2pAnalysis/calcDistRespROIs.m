function [distRespROIs] = calcDistRespROIs(analysisParams, data, analysis, umperpixel, level)    

%only do the distance calculation between ROIs that are active and ORi
%selective 
respROIs = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);

%calculate distances between all cell pairs
distRois = zeros(length(respROIs), length(respROIs));
for A = 1:length(respROIs)
    for B = 1:length(respROIs)
        if level == 1
            if data.roi(A).plane == data.roi(B).plane
                %calculate the distance
                dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
                dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
                distRois(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel; %change it to um
            else
                distRois(A,B) = NaN;
            end
        else
            dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
            dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
            distRois(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
        end
    end
end

distRespROIs = distRois;


