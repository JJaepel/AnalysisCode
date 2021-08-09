function [distOri, deltaOri, distRois] = calcDistDeltaOri(data, analysis, ori_sel,umperpixel, level, field)    

distRois = zeros(length(ori_sel), length(ori_sel));
deltaOri = zeros(length(ori_sel), length(ori_sel));
for A = 1:length(ori_sel)
    for B = 1:length(ori_sel)
        if level == 1
            if data.roi(A).plane == data.roi(B).plane
                dist_x = abs(data.roi(A).xPos - data.roi(B).xPos);
                dist_y = abs(data.roi(A).yPos - data.roi(B).yPos);
                distRois(A,B) = sqrt(dist_x^2 + dist_y^2)*umperpixel;
                deltaOri(A,B) = abs(analysis.(field).roi(A).preferredOrientation - analysis.(field).roi(B).preferredOrientation);
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
            deltaOri(A,B) = abs(analysis.(field).roi(A).preferredOrientation - analysis.(field).roi(B).preferredOrientation);
            if deltaOri(A,B) > 90
                deltaOri(A,B) = 180 - deltaOri(A,B);
            end
        end
    end
end
L = triu(distRois);
L = reshape(L,1,length(ori_sel)^2);
L(L == 0) = NaN;
distOri = L((~isnan(L)));
deltaOri = reshape(deltaOri,1,length(ori_sel)^2);
deltaOri = deltaOri((~isnan(L)));
end