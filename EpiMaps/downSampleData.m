function data = downSampleData(data, downsampleFactor)

stackSize        = size(data.rawF);
stackSize     = floor(stackSize/downsampleFactor);
downsampledStack = zeros(stackSize,class(data.rawF));
for i=1:downsampleFactor
    downsampledStack = downsampledStack+data.rawF(i:downsampleFactor:downsampleFactor*stackSize(1),i:downsampleFactor:downsampleFactor*stackSize(2),i:downsampleFactor:downsampleFactor*stackSize(3))/downsampleFactor;
end
data.rawF = downsampledStack;
clear downsampledStack