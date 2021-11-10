function LowHighData=LowHighNormalize(Frame, Mask, varargin)
%sigma_min: Lowpassfilter, removes local noise
%sigma_max: Highpass, removes global activation
%ferret data: sigma_min=2 (~26 micrometers), sigma_max=15 (~200
%Mikrometers)
if nargin==2    
    sigma_min=1;
    sigma_max=8;
elseif nargin==3
    sigma_max=15;
    sigma_min=varargin{1};
elseif nargin==4
    sigma_min=varargin{1};
    sigma_max=varargin{2};
end
    
Frame(~Mask(:))=0;
MaskF=1.*Mask;

if sigma_min>0
    LowData = imgaussfilt(Frame,sigma_min)./imgaussfilt(MaskF,sigma_min);
    LowData(~Mask(:))=0;
else
    LowData = Frame;
end
LowHighData = LowData - imgaussfilt(Frame,sigma_max)./imgaussfilt(MaskF,sigma_max);
LowHighData(~Mask(:))=NaN;
end


