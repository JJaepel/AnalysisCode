function HI = calculateHI(PopAnalysis,analysis, level, data, field)
% Calculates homeogeneity index for a field of view for a variety of radiuses 
%
% Input:
% - oriSel: ROIs that are orientation tuned and respond significantly
% - umperpixel: 
% - level: volume imaging?
% - radius: a 1xnumber of radius array
%
% Ouput:
% - HI: a 1xnumRadius array stating homeogeneity index for different
% radiuses
%

ori_sel = find([analysis.(field).roi.OSIFit] > 0.2 & [analysis.(field).roi.isResponseSignificant] == 1);

radius = linspace(0,250,6);
if level == 1 
    plane_ori = [data.roi(ori_sel).plane];
end
HI = zeros(length(ori_sel),length(radius));
preferences = [analysis.(field).roi(ori_sel).preferredOrientation];
for dis = 1:length(radius)-1
    for cellID = 1:length(ori_sel)-1
        distances = PopAnalysis.(field).distOriROIs(cellID, :);
        withinRadius1 = distances >= radius(dis);
        withinRadius2 = distances < radius (dis+1);
        withinRadius = withinRadius1&withinRadius2;
        if level == 1
            sameplane = [plane_ori == data.roi(cellID).plane];
            withinRadius = withinRadius&sameplane;
        end
        withinRadius(cellID)=0;
        thresholdOrientations=exp(2*1i*preferences(withinRadius));
        HI(cellID, dis) = 1-(abs(1/sum(withinRadius) * sum(thresholdOrientations)));
        if HI(cellID,dis) == 0
            HI(cellID,dis) = NaN;
        end
    end
end