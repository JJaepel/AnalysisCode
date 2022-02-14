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

for i = 1:length(events)
   %let's start with the onset location and then use a loop to find the
   %location of the offset by using type
   for type = 3%1:3
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
       factor = 1; %for peak, change the threshold

       while objectSize < minArea
          framesForLocation = areaData(:,:,events(i).(field)+start-1:events(i).(field)+start+1);
           stdValue = std(double(framesForLocation(:)),'omitnan');
           maxValue = max(double(framesForLocation),[],'all','omitnan');
           threshold = maxValue - factor*stdValue; 

           %threshold by maximum of all location - one std and build a mask
           maxImage = max(framesForLocation,[],3);
           mask = maxImage;
           mask(mask<threshold) = 0;
           mask(mask>threshold) = 1;

           %find largest object
           objects = bwlabel(mask);
           objects(objects == 0) = NaN;
           [biggestObject, objectSize]= mode(objects(:)); %which object is the biggest one
           
           if objectSize > areaSize %this is to avoid super large objects in on- and offset
               objectSize = 0;
           end
           
           %in case the object is too small, go forward (for onset) or
           %backwards (for offset)
           if type == 1 
               start = start+1;
               if ojectSize > maxArea %if it is quite big, adjust the threshold, but stay in the same frames
                   factor = factor/sqrt(2);
                   objectSize = 0;
                   start = start-1;
               end
           elseif type ==2 
               start = start-1;
               if ojectSize > maxArea %if it is quite big, adjust the threshold, but stay in the same frames
                   factor = factor/sqrt(2);
                   objectSize = 0;
                   start = start+1;
               end
           elseif type == 3 %for peaks, adjust the threshold
               if objectSize < minArea %when it is too small, make the threshold smaller by substracting more std
                   factor = factor*1.9; %it's not exactly 2 to make sure we don't get into a loop
               elseif objectSize > maxArea %when it is too big, make the threshold higher and set objectSize to 0 to restart the loop
                   factor = factor/sqrt(2);
                   objectSize = 0;
               end
           end
       end
       %restrict the mask to the biggest Object
       objects(objects~=biggestObject) = 0;
       mask = mask .* int16(objects/biggestObject); %multiply mask with the mask from the biggest object and divide by the biggestObject to make sure the mask is only 1 and 0
       
       if verbose
           figure
           imagesc(mask.*maxImage)
           saveas(gcf, fullfile(saveDir, [AreaName ' Eventlocation Nr ' num2str(i)]))
       end
       
       %find the center of mass as the center location of it
       point = regionprops(mask, mat2gray(mask.*maxImage), 'WeightedCentroid');
       %save the mask and the center location
       name = [field 'Center'];
       events(i).(name) = round(point.WeightedCentroid);
       name = [field 'Area'];
       mask(mask>threshold) = 1;
       events(i).(name) = mask;
   end
end