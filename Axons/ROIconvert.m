function [ROIs cvsROIs] = ROIconvert(zipname,imgsize)

[cvsROIs] = ReadImageJROI(zipname);

map = zeros(imgsize(:,1), imgsize(:,2));
finmap = zeros(imgsize(:,1), imgsize(:,2));
k = 1;

h = waitbar(0,'Converting ROIs...');

for i = 1:length(cvsROIs);
    try
        waitbar(i / length(cvsROIs));
        
        if sum(cvsROIs{i}.vnRectBounds)==0;
            continue
        end
        
        temp_r = poly2mask(cvsROIs{i}.mnCoordinates(:,1), cvsROIs{i}.mnCoordinates(:,2),imgsize(1),imgsize(2));
        map = map + temp_r;
        
        %overlap delete
        temp_r(map==2) = 0;
        
        %standard indices
        tempidc = find(temp_r==1);
        
        if length( tempidc )==0
            continue
        else
            ROIs(k).indices = tempidc;
        end
        
%         finmap(ROIs(k).indices) = 1;
        
        %type field
        ROIs(k).typen = 'd';
        
        %CHECK FOR DOUBLE REGIONPROPS!!! USE LARGEST ONE
        STATS = regionprops(temp_r,'Centroid', 'Area', 'PixelList', 'ConvexHull' );
        [~, idx] = max([STATS.Area]);

        ROIs(k).perimeter = cvsROIs{i}.mnCoordinates+1; %TR2014: '+1' corrects for xy offset introduced by the IJ converter...
%         ROIs(k).perimeter = STATS(idx).ConvexHull; 
        
        
        ROIs(k).body = STATS.PixelList;
        ROIs(k).type = 1;
        ROIs(k).typename = 'Neuron';
        ROIs(k).size = length(ROIs(k).body(:,2));
        ROIs(k).x = STATS(idx).Centroid(:,1)+1;
        ROIs(k).y = STATS(idx).Centroid(:,2)+1;
        
        k = k+1;
    end
end

close(h)

if ~isempty(strfind(zipname, 'RoiSet_bv'))
    ROIs = [];
    ROIs.indices = find(map==1);
end
end

