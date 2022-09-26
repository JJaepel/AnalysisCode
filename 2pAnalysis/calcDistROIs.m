function [distROIs] = calcDistROIs(analysisParams, data, analysis, umperpixel, level)    

%calculate distances between all cell pairs
distRois = zeros(length(analysis.(analysisParams.field).roi), length(analysis.(analysisParams.field).roi));
for A = 1:length(analysis.(analysisParams.field).roi)
    for B = 1:length(analysis.(analysisParams.field).roi)
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

distROIs = distRois;


