function [Lat,LatON,LatOFF]=FindLatency(StimDur_corr,Significance,DisplayWindow,ind,p)

LatON=StimDur_corr+find(Significance(1, find(DisplayWindow==1):end,ind)<p,1,'first'); %we add 3 frames at the beginning corresponding to the previous patch presentation (frames -2to0)
LatOFF=StimDur_corr+find(Significance(2, find(DisplayWindow==1):length(DisplayWindow),ind)<p,1,'first');


if isempty(LatON)&& isempty(LatOFF)
    LatON=nan;
    LatOFF=nan;
end

if isempty(LatON)
    LatON=nan;
elseif LatON>(length(DisplayWindow)- find(DisplayWindow==1)-(2*StimDur_corr))
    LatON=nan;
end

if isempty(LatOFF)
    LatOFF=nan;
elseif LatOFF>(length(DisplayWindow)- find(DisplayWindow==1)-(2*StimDur_corr))   
LatOFF=nan;
end
Lat=min(LatON,LatOFF);
end