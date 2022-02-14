function data = downSampleData(data, downsampleFactor)

downsample = 1/downsampleFactor;
scaledX = round(size(data,1)*downsample);
scaledY = round(size(data,2)*downsample);
        
downsampledStack = zeros(scaledX,scaledY,size(data,3),'uint16');
for i = 1:size(data,3)
    downsampledStack(:,:,i) = imresize(data(:,:,i),[scaledX scaledY]);
end

data = downsampledStack;
clear downsampledStack