function [Lat,LatON,LatOFF]=FindLatencyRotPatch(StimDur_corr,Significance,DisplayWindow,p)

LatON=StimDur_corr+find(Significance(1, find(DisplayWindow==1):end)<p,1,'first'); %we add 3 frames at the beginning corresponding to the previous patch presentation (frames -2to0)
LatOFF=StimDur_corr+find(Significance(2, find(DisplayWindow==1):length(DisplayWindow))<p,1,'first');


if isempty(LatON)&& isempty(LatOFF)
    LatON=nan;
    LatOFF=nan;
end

if isempty(LatON)
    LatON=nan;
end

if isempty(LatOFF)
    LatOFF=nan;
end
Lat=min(LatON,LatOFF);
end