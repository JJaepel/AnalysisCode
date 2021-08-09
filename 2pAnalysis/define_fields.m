function [ONperimeter,ONx,ONy,ONsize,ONx_ax, ONy_ax, ONmajor_ax, ONminor_ax]= define_fields(ONfield,metadata)

UsedScreenWidth = abs(metadata.StimParams.minAzi - metadata.StimParams.maxAzi);
UsedScreenHeight = abs(metadata.StimParams.minElev - metadata.StimParams.maxElev);

[labeledObject,nPolygons] = bwlabel(ONfield,8);
 PolMax = [];
        for p = 1:nPolygons
            PolMax(p) = max(max(ONfield(labeledObject==p)));
        end
        PolNr = find(PolMax==max(PolMax));
Basic  = regionprops(labeledObject, 'Basic'); %calculates center
Pixels = regionprops(labeledObject, 'PixelList') ;%size

[y_min pos_min] = min(Pixels.PixelList(:,2));
[y_max pos_max] = max(Pixels.PixelList(:,2));
x_max = Pixels.PixelList(pos_max,1);
x_min = Pixels.PixelList(pos_min,1);
if x_min == x_max | x_min == x_max+1 | x_min+1 == x_max
    if x_min == 1 | x_max == 1
        mirror_labeledObject = flipdim(labeledObject,2);
        middle = zeros(size(labeledObject,1),1);
        middle(y_min,1) = 1;
        middle(y_max,1) = 1;
        all_Object = [mirror_labeledObject(:,1:1:size(labeledObject,2)-1) middle labeledObject(:,2:1:end)];
        xShift = (metadata.StimParams.minAzi)- size(labeledObject,2);
    elseif x_min == size(labeledObject,2) | x_max == size(labeledObject,2)
        mirror_labeledObject = flipdim(labeledObject,2);
        middle = zeros(size(labeledObject,1),1);
        middle(y_min,1) = 1;
        middle(y_max,1) = 1;
        all_Object = [labeledObject(:,1:1:size(labeledObject,2)-1) middle mirror_labeledObject(:,2:1:end)];
        xShift = (metadata.StimParams.minAzi);
    else
        all_Object = labeledObject;
        xShift = (metadata.StimParams.minAzi);   
    end
else
    all_Object = labeledObject;
    xShift = (metadata.StimParams.minAzi);
end
Pixels = regionprops(all_Object, 'PixelList'); 
[x_min2 pos_min_x] = min(Pixels.PixelList(:,1));
[x_max2 pos_max_x] = max(Pixels.PixelList(:,1));
y_max2 = Pixels.PixelList(pos_max_x,2);
y_min2 = Pixels.PixelList(pos_min_x,2);
if y_min2 == y_max2 | y_min2 == y_max2+1 | y_min2+1 == y_max2
    if y_min2 == 1 | y_max2 == 1
        mirror_all_Object = flipdim(all_Object,1);
        middle = zeros(1, size(all_Object,2));
        middle(1,x_min2) = 1;
        middle(1,x_max2) = 1;
        all_Object_final = [mirror_all_Object(1:1:size(labeledObject,1)-1,:); middle; all_Object(2:1:end,:)];
        yShift = (metadata.StimParams.minElev) - size(labeledObject,1);
    elseif y_min2 == size(labeledObject,1) | y_max2 == size(labeledObject,1)
        mirror_all_Object = flipdim(all_Object,1);
        middle = zeros(1, size(all_Object,2));
        middle(1,x_min2) = 1;
        middle(1,x_max2) = 1;
        all_Object_final = [all_Object(1:1:size(labeledObject,1)-1,:); middle; mirror_all_Object(2:1:end,:)];
        yShift = (metadata.StimParams.maxElev);
    else 
        all_Object_final = all_Object;
        yShift = (metadata.StimParams.maxElev);
    end
else
    all_Object_final = all_Object;
    yShift =(metadata.StimParams.maxElev);
end
xShift2 = (metadata.StimParams.minAzi);
yShift2 = (metadata.StimParams.maxElev);
imagesc(all_Object_final)
Area = regionprops(all_Object_final, 'FilledArea');
Perimeter = regionprops(all_Object_final, 'ConvexHull');
ONsize = Area.FilledArea;

ONperimeter =  Perimeter(PolNr).ConvexHull;
ONperimeter(:,1) = ONperimeter(:,1)+xShift;
ONperimeter(:,2) = yShift - ONperimeter(:,2);
ONx= Basic(PolNr).Centroid(1)+xShift2;
ONy = yShift2-Basic(PolNr).Centroid(2);

if metadata.StimParams.minElev > metadata.StimParams.maxElev
end

ONx_ax = abs(min(ONperimeter(:,1)) - max(ONperimeter(:,1)));
ONy_ax = abs(min(ONperimeter(:,2)) - max(ONperimeter(:,2)));

if ONx_ax > ONy_ax
    ONmajor_ax = ONx_ax;
    ONminor_ax = ONy_ax;
else
    ONmajor_ax = ONy_ax;
    ONminor_ax = ONx_ax;
end
% ONpixels = Pixels(PolNr).PixelList;
% ONpixels(:,1) = ONpixels(:,1)+xShift;
% ONpixels(:,2) = ONpixels(:,2)+yShift;
% ONsize=size(ONpixels,1);
