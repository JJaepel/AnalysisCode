function RFfinal=PlotCleanRF(RF,metadata,threshold, ind, figSaveDirectory)

temp = abs(RF);
temp(RF < threshold) = 0;

if temp == 0
    threshold = mean(max(RF));
    temp = abs(RF);
    temp(RF < threshold) = 0;
end

%create binary mask
tempBW = bwlabel(temp,4);

%find max intensity for each supra-threshold area
props=regionprops(tempBW,temp,'MaxIntensity');
[~, I]=max([props.MaxIntensity]);
    
% keep only the area containing the max intensity as your RF
temp(tempBW~=I)=0;
RFfinal=temp;
    
x=[metadata.StimParams.minAzi metadata.StimParams.maxAzi];
y=[metadata.StimParams.maxElev metadata.StimParams.minElev];
    
f=figure;
subplot(1,2,1),imagesc(x,y,temp),axis image
set(gca, 'YDir','normal')
xlabel('azimuth');
ylabel('elevation');
title(['Max responsive RF areas. Threshold:' sprintf('%1d',threshold)])
subplot(1,2,2),imagesc(x,y,RFfinal),axis image
set(gca, 'YDir','normal')
xlabel('azimuth');
ylabel('elevation');
saveas(f, fullfile(figSaveDirectory, ['CleanRFs, thresholded, ROI ',num2str(ind) '.png']))
pause(0.1);
close(f);
