function zoom = getZoomSpines(tifDir)
    filename = [tifDir filesep '*.tif'];
    files = dir(filename);
    filepath = [tifDir filesep files(1).name];
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
end