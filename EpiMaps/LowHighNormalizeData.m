function LowHighData=LowHighNormalizeData(data, Mask, sigma_min, sigma_max)
%sigma_min: Lowpassfilter, removes local noise

if nargin<4 
    sigma_max = 15;
end
if nargin < 3
    sigma_min = 1;
    sigma_max = 10;
end

LowHighData = zeros(size(data,1), size(data,2), size(data,3));
for frame = 1:size(data,3)
    LowHighData(:,:,frame) = LowHighNormalize(data(:,:,frame), Mask, sigma_min,sigma_max);
end

end