function HighData=HighNormalize(Frame, Mask, sigma_max)
%sigma_min: Lowpassfilter, removes local noise
%sigma_max: Highpass, removes global activation
%ferret data: sigma_min=2 (~26 micrometers), sigma_max=15 (~200
%Mikrometers)
if nargin<3    
    sigma_max=8;
end
    
Frame(~Mask(:))=0;
MaskF=1.*Mask;

HighData = imgaussfilt(Frame,sigma_max)./imgaussfilt(MaskF,sigma_max);
HighData(~Mask(:))=NaN;
end