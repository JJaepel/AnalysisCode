function events = findLocations(events, areaData, AreaName, verbose, saveDir)
% finds the location of the hotspot of activation at the time of event
% onset, offset and peak
if nargin < 3
    AreaName = '';
end
if nargin < 4
    verbose = 0; %change this to 1 if you want to see all the peak locations
end
if nargin < 5
    saveDir = cd;
end

areaSize = nnz(mean(areaData,3));
minArea = areaSize*0.03; %shouldn't be smaller than 3 %
maxArea = areaSize*0.15; %but also not bigger than 15 %
peakThreshold = mean(areaData(:)) + std(double(areaData(:)));

for i = 1:length(events)
   %let's start with the onset location and then use a loop to find the
   %location of the offset by using type
   for type = 1:3
      if type == 1 %depending on the type, work on on- or offset location
          field = 'onset';
      elseif type == 2 
          field = 'offset';
      elseif type == 3
          field = 'peakTime';
      end
       
       %build a loop to make sure that the location is at least 100 pixel
       %big
       objectSize = 0;
       start = 0; %for on- and offset, slightly change the frames
       factor = 1; %for on- and offset, change the threshold
       numIt = 0;
       while objectSize < minArea && numIt < 20
           framesForLocation = areaData(:,:,events(i).(field)+start-1:events(i).(field)+start+1);
           stdValue = std(double(framesForLocation(:)),'omitnan');
           maxValue = max(double(framesForLocation),[],'all','omitnan');
           
           if type == 1 || type == 2
               threshold = maxValue - factor*stdValue;
           elseif type == 3
               threshold = peakThreshold;
           end
         

           %threshold by maximum of all location - one std and build a mask
           maxImage = max(framesForLocation,[],3);
           mask = maxImage;
           mask(mask<threshold) = 0;
           mask(mask>threshold) = 1;

           %find largest object
           objects = bwlabel(mask);
           objects(objects == 0) = NaN;

           [biggestObject, objectSize]= mode(objects(:)); %which object is the biggest one
           
           if isempty(biggestObject)
               objectSize = 0;
           end
           
           if objectSize > areaSize %this is to avoid super large objects in on- and offset
               objectSize = 0;
           end
           
           %in case the object is too small, go forward (for onset) or
           %backwards (for offset)
           if type == 1 
               start = start+1;
               if objectSize > maxArea %if it is quite big, adjust the threshold, but stay in the same frames
                   factor = factor/sqrt(2);
                   objectSize = 0;
                   start = start-1;
               end
           elseif type ==2 
               start = start-1;
               if objectSize > maxArea %if it is quite big, adjust the threshold, but stay in the same frames
                   factor = factor/sqrt(2);
                   objectSize = 0;
                   start = start+1;
               end
           elseif type == 3
               numIt = 20;
           end
           numIt = numIt+1;
       end
       if numIt >= 20
          [biggestObject, objectSize]= mode(objects(:)); 
       end
       %restrict the mask to the biggest Object
       objects(objects~=biggestObject) = 0;
       peakMask = mask;
       mask = mask .* int16(objects/biggestObject); %multiply mask with the mask from the biggest object and divide by the biggestObject to make sure the mask is only 1 and 0
       
       if verbose
           figure
           if type == 1 || type == 2
               imagesc(mask.*maxImage)
           elseif type == 3
               imagesc(peakMask.*maxImage)
           end
           saveas(gcf, fullfile(saveDir, [AreaName ' ' field ' Eventlocation Nr ' num2str(i)]))
       end
       
       %find the center of mass as the center location of it
       point = regionprops(mask, mat2gray(mask.*maxImage), 'WeightedCentroid');
       %save the mask and the center location
       name = [field 'Center'];
       events(i).(name) = round(point.WeightedCentroid);
       name = [field 'Area'];
       mask(mask>threshold) = 1;
       events(i).(name) = mask;
       if type == 3 %for peaks also save the percentage of involved area
           events(i).(name) = peakMask;
           events(i).percArea = nnz(peakMask)/areaSize;
       end
   end
end