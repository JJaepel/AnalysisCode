function HI_SF = calculateHI_SF(PopAnalysis,analysis, metadata, level, data, field)
% Calculates homeogeneity index for a field of view for a variety of radiuses 
%
% Input:
% - umperpixel: 
% - level: volume imaging?
% - radius: a 1xnumber of radius array
%
% Ouput:
% - HI: a 1xnumRadius array stating homeogeneity index for different
% radiuses
%

respROIs = [analysis.(field).roi.isResponseSignificant] == 1;
radius = linspace(0,250,6);
if level == 1 
    plane_ori = [data.roi.plane];
end

HI_SF = zeros(length(respROIs),length(radius));
preferences = [analysis.(field).roi.prefSf];
sfs = metadata.StimParams.spatialFreq;

for dis = 1:length(radius)-1
    for cellID = 1:length(respROIs)-1
        distances = PopAnalysis.(field).distRespROIs(cellID, :);
        withinRadius1 = distances >= radius(dis);
        withinRadius2 = distances < radius (dis+1);
        withinRadius = withinRadius1&withinRadius2;
        withinRadius = withinRadius &respROIs;
        if level == 1
            sameplane = [plane_ori == data.roi(cellID).plane];
            withinRadius = withinRadius&sameplane;
        end
        withinRadius(cellID)=0;
        preferenceCell=preferences(cellID);
        samePreference = preferences == preferenceCell;
        samePrefWithinRadius = withinRadius & samePreference;
        fractPref = sum(samePrefWithinRadius)/sum(withinRadius);
        
        for drawing = 1:100
            randPreferences = randi([1,length(sfs)],1,length(respROIs));
            randSamePref = randPreferences == preferenceCell;
            randSamePrefWithinRadius = withinRadius & randSamePref;
            randFractAll(drawing) = sum(randSamePrefWithinRadius)/sum(withinRadius);
        end
        
        randFractPref = mean(randFractAll);
        
        HI_SF(cellID,dis) = fractPref - randFractPref;
       
        if HI_SF(cellID,dis) == 0
            HI_SF(cellID,dis) = NaN;
        end
    end
end
