function [umperpixel] = getzoom(tifDirectory)
    filename = [tifDirectory filesep '*.tif'];
    files = dir(filename);
    filepath = [tifDirectory filesep files(1).name];
    InfoImage = imfinfo(filepath);
    try 
        a = InfoImage(1).Software;
        setup = 1;
    catch 
        a = InfoImage(1).ImageDescription;
        setup = 2;
    end
    zoom = regexp(a,'(?<=scanZoomFactor = )\d+\.?\d*', 'match');
    zoom = str2num(zoom{1});
    
    if setup == 1
        fieldofview = 1000/zoom;
        umperpixel = fieldofview/512;
        %disp('Setup Ben')
    elseif setup == 2
        umperpixel = 2.73/zoom;
    end
end