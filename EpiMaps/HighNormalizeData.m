function HighData=HighNormalizeData(data, Mask, sigma_max)
%sigma_min: Lowpassfilter, removes local noise

if nargin < 3
    sigma_max = 10;
end

HighData = zeros(size(data,1), size(data,2), size(data,3));
for frame = 1:size(data,3)
    HighData(:,:,frame) = HighNormalize(squeeze(data(:,:,frame)), Mask, sigma_max);
end

end