%modify!
imageDirectory = 'F:\Data\2P_Data\'; 
metaDirectory = 'F:\Data\Spike2Data\';
imgAnalysisDirectory = 'F:\Data\ImageAnalysis\';
name = 'F2444_2020-09-07';
Sp2name = 'F2444_2020-09-07'; 
exptID = 10;
spk2ID = 10;

%switchpanel
spatialDownsamplingFactor = 1;
ori = 1;
tune_max = 0.2;
colorShift = 330/360; %adjust!
filter = 1;
resize = 0;
close all;
savedir = '\Pixelmaps\';

%define location of imaging data and meta data
imagingDataDirectory = sprintf('%s%s\\t%05d\\',imageDirectory,name,exptID);
imagingDataDirectory = [imagingDataDirectory '\Registered\slice1']; %make sure there are no other tifs in the dir
imagingDataDirectory = [imageDirectory name '\t00010\']
metaDataDirectory    = sprintf('%s%s\\t%05d\\',metaDirectory,Sp2name,spk2ID);
imgAnalysDataDirectory = sprintf('%s%s\\t%05d\\',imgAnalysisDirectory,name,exptID);

%Load imaging data
data = LoadImagingData(imagingDataDirectory, spatialDownsamplingFactor);

%load metadata
stim_matrix=[];
disp(strcat('Loading....',metaDataDirectory, '\stimontimes.txt'))
stimtimes=load(strcat(metaDataDirectory, '\stimontimes.txt'));
i=1;
while i <length(stimtimes)
    idx= stimtimes(i);
    stim_matrix(1, floor(i/2)+1)= stimtimes(i);
    stim_matrix(2, floor(i/2)+1)= stimtimes(i+1);

    i=i+2;
end

StimParams = struct;
StimParams.StimOnTimes = stim_matrix;
StimParams.uniqStimIds= unique(StimParams.StimOnTimes(1,1:end));
StimParams.uniqStims = length(StimParams.uniqStimIds);

files=dir(strcat(metaDataDirectory, '\*.py'));
files={files.name};
file = [metaDataDirectory filesep char(files{1})];
StimParams.stimDuration = FindStimulusParam(file, 'stimDuration'); 

ImagingTime = load(strcat(metaDataDirectory,'\frametrigger.txt')); %alternatively: '\frametrigger.txt'

% filter calculation
lowpass = 2; % um
pixel_in_x = size(data,2);
fov_in_um = 200;
lowpass_pix = lowpass * pixel_in_x / fov_in_um;
kernel_size = ceil(lowpass_pix * 5);
sp_filter = fspecial('gaussian', kernel_size, lowpass_pix);

% create stimIDs
n = length(StimParams.StimOnTimes);
stimStartIndex = zeros(n,1,'double'); 
stimStopIndex  = zeros(n,1,'double');

for i=1:n
    stimStartIndex(i) = find(ImagingTime>=StimParams.StimOnTimes(2,i),1,'first');
    stimStopIndex(i) = find(ImagingTime>=StimParams.StimOnTimes(2,i)+StimParams.stimDuration,1,'first');
end

stim_on = logical(zeros(1, size(data,3)));
stim_off = logical(ones(1, size(data,3)));
stim_id = uint8(zeros(1, size(data,3)));

for stimulus = 1:length(stimStartIndex)
    stim_on(1, stimStartIndex(stimulus):stimStopIndex(stimulus)) = 1;
    stim_off(1, stimStartIndex(stimulus):stimStopIndex(stimulus)) = 0;
    stim_id(1, stimStartIndex(stimulus):stimStopIndex(stimulus)) = StimParams.StimOnTimes(1,stimulus);
end

outputpath = [imgAnalysDataDirectory, savedir];
if ~exist(outputpath, 'dir')
    mkdir(outputpath);  
end
if ori
    orioutputpath = [outputpath filesep 'ori'];
    if ~exist(orioutputpath, 'dir')
        mkdir(orioutputpath);
    end
end

%calc base image
base = mean(data(:,:,stim_off'),3);
if filter
   base = Filter2Modified(sp_filter, base);
end

%calc F, dF and dFoverF images
no_conditions = StimParams.uniqStims - 1;
dim = size(data);
F_img = zeros(dim(1), dim(2), no_conditions);
dF_img = zeros(dim(1), dim(2), no_conditions);
dF_F_img = zeros(dim(1), dim(2), no_conditions);

sts = sort(unique(stim_id));
sts = sts(2:end);
for ind = 1:no_conditions
    st_ind = sts(ind);
    F_img(:,:,ind) = mean(data(:,:,stim_id == st_ind),3);
    dF_img(:,:,ind) = F_img(:,:,ind) - base;
    if filter
       dF_img(:,:,ind) = Filter2Modified(sp_filter, dF_img(:,:,ind));
    end
    dF_F_img(:,:,ind) = dF_img(:,:,ind) ./base;
end
max_dF = max(max(max(dF_img)));
max_dF_F = max(max(max(dF_F_img)));


for ind = 1:no_conditions
    imwrite(dF_img(:,:,ind) ./ max_dF, fullfile(outputpath,[ 'dF_', num2str(floor(ind)), '.tif']));
    imwrite(dF_F_img(:,:,ind) ./ max_dF_F, fullfile(outputpath,[ 'dF_F_', num2str(floor(ind)), '.tif']));
end

%calc map parameters for each pixel
params = CalcMapParams(dF_F_img);
[HLS, HLS_hc] = GenerateHLSMap(params, max_dF_F, 'maxTuning', tune_max, 'colorshift', colorShift);

%show HLS maps
figure
imshow(HLS)
imwrite(HLS, fullfile(outputpath,'HLS.tif'));
figure
imshow(HLS_hc)
imwrite(HLS_hc, fullfile(outputpath,'HLS_hc.tif'));

if ori %if orientation map wanted, do the same for orientation maps
    ori_F_img = zeros(dim(1), dim(2), no_conditions/2);
    ori_dF_img = zeros(dim(1), dim(2), no_conditions/2);
    ori_dF_F_img = zeros(dim(1), dim(2), no_conditions/2);
    for nori = 1:no_conditions/2
        st_ind = sts(nori);
        st2_ind = sts(nori+no_conditions/2);
        ori1 = stim_id == st_ind; 
        ori2 = stim_id == st2_ind;
        oris = logical(ori1 + ori2);
        ori_F_img(:,:,nori) = mean(data(:,:,oris),3);
        ori_dF_img(:,:,nori) = ori_F_img(:,:,nori) - base;
        if filter
           ori_dF_img(:,:,nori) = Filter2Modified(sp_filter, ori_dF_img(:,:,nori));
        end
        ori_dF_F_img(:,:,nori) = ori_dF_img(:,:,nori) ./base;
    end
    ori_max_dF = max(max(max(dF_img)));
    ori_max_dF_F = max(max(max(dF_F_img)));
    for nori = 1:no_conditions/2
        imwrite(dF_img(:,:,nori) ./ori_max_dF, fullfile(orioutputpath,[ 'ori_dF_', num2str(floor(nori)), '.tif']));
        imwrite(dF_F_img(:,:,nori) ./ori_max_dF_F, fullfile(orioutputpath,[ 'dF_F_', num2str(floor(nori)), '.tif']));
    end
    ori_params = CalcMapParams(ori_dF_F_img);
    [ori_HLS, ori_HLS_hc] = GenerateHLSMap(ori_params, ori_max_dF_F, 'maxTuning', tune_max, 'colorshift', colorShift);
    
    figure
    imshow(ori_HLS)
    imwrite(ori_HLS, fullfile(orioutputpath,'ori_HLS.tif'));
    figure
    imshow(ori_HLS_hc)
    imwrite(ori_HLS_hc, fullfile(orioutputpath,'ori_HLS_hc.tif'));
    angles = mod(ori_params.th + colorShift, 1 );
    filename = [orioutputpath '\anglemap.mat'];
    save(filename, 'angles')
    if resize
        centerPoint = [ceil(size(ori_HLS,1)/2),ceil(size(ori_HLS,2)/2)];
        dimx = ceil(size(ori_HLS,1)/4);
        dimy = ceil(size(ori_HLS,2)/4);
        ori_HLS_rs = ori_HLS(centerPoint(1)-dimx+1:centerPoint(1)+dimx,centerPoint(2)-dimy+1:centerPoint(2)+dimy,:);
        ori_HLS_hc_rs = ori_HLS_hc(centerPoint(1)-dimx+1:centerPoint(1)+dimx,centerPoint(2)-dimy+1:centerPoint(2)+dimy,:);
        angles_rs = angles(centerPoint(1)-dimx+1:centerPoint(1)+dimx,centerPoint(2)-dimy+1:centerPoint(2)+dimy);
        angles_rs = imresize(angles_rs, [size(angles,1), size(angles,2)]);
        imwrite(ori_HLS_rs, fullfile(orioutputpath,'ori_HLS_rs.tif'));
        imwrite(ori_HLS_hc_rs, fullfile(orioutputpath,'ori_HLS_hc_rs.tif'));
        filename = [orioutputpath '\anglemap_resized.mat'];
        save(filename, 'angles_rs')
    end
end

function data = LoadImagingData(imgPath, spatialDownsamplingFactor)
    downsample = 1/spatialDownsamplingFactor;
    % Determine number of multi-page tif images are present.    
    files= dir(strcat(imgPath,'\stack*.tif'));
    files = dir(strcat(imgPath,'t00010*.tif'));
    files={files.name}';
    numberOfFiles = size(files, 1);
    
    % Determine total number of frames present, and initialize rawF array
    framesPerFile = zeros(numberOfFiles,1);
    filesToCheck = 1:numberOfFiles;
    for n = filesToCheck 
        fileName  = char(strcat(imgPath,'\',files(n)));
        imageInfo = imfinfo(fileName);
        framesPerFile(n)=numel(imageInfo);
    end
    numberOfFrames = sum(framesPerFile);
    data = zeros([imageInfo(1).Height*downsample,imageInfo(1).Width*downsample,numberOfFrames],'uint16');
    
    % Read imaging frames into MATLAB 
    frameCounter = 0;
    for n = 1:numberOfFiles
        fileName = char(strcat(imgPath,'\',files(n)));
        data(:,:,frameCounter+[1:framesPerFile(n)])=read_Tiffs(fileName,downsample,100);
        frameCounter=frameCounter+framesPerFile(n);
    end

end
function tifStack = read_Tiffs(filePath,imgScaling,updateFrequency,useWaitBar)
% tifStack = read_Tiffs(filePath,updateFrequency)
%
% A fast method to read a tif stack: calls the tiff library directly.
% Used the same method described in this web-link: 
% http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/
% 
% filePath - fileDirectory and name of image stack
% imgScaling - Scales size of images by specified scaling factor (default is no scaling).
%              0.5 drops image size by half, while 2x doubles image size
% updateFrequency - Based on the number of images, how often in percentage
% of files read to inform user of current progress (0-100%)
% useWaitBar - displays a wait bar to denote current reading progress (default it
% false due to slowing down function)
%
% Written by David Whitney (10/3/2013)
% Updated by David Whitney (5/16/2016)
% David.Whitney@mpfi.org
% Max Planck Florida Institude
warning off;
if(nargin<2), imgScaling = 0.5;    end
if(nargin<3), updateFrequency = 5; end
if(nargin<4), useWaitBar = false;  end
tic;

disp(['Reading Image Stack - ' filePath]);

% Read TIF Header and Setup Information
InfoImage=imfinfo(filePath);
    xImage=InfoImage(1).Width;
    yImage=InfoImage(1).Height;
    NumberOfImages=length(InfoImage);
    disp(['Finished Reading Image Header - ' num2str(toc) ' seconds Elapsed']);

% use wait bar
if(useWaitBar)
    h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
%     currentPosition = get(h,'Position');
%     offset = 100;
%     set(h,'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
%     currentPosition = get(get(h,'Children'),'Position');
%     set(get(get(h,'Children')),'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
else
    updateFrequency = round((updateFrequency/100)*NumberOfImages);
end

% Initialize MATLAB array to contain tif stack
scaledX = round(xImage*imgScaling);
scaledY = round(yImage*imgScaling);
tifStack     = zeros(scaledY,scaledX,NumberOfImages,'uint16');
           
codeVersion = 'alternativeMethod'; % both methods seem pretty similar in performance, but original method is a bit faster
switch codeVersion
    case 'originalMethod'
        % uses the tifflib function to read images fast
        FileID = tifflib('open',filePath,'r');
        rps    = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
        rps    = min(rps,yImage);
        for i=1:(NumberOfImages)
            % display read progress
            if(useWaitBar)
                waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
            else
                if(mod(i+1,updateFrequency)==0)
                    disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                end
            end

            % turn off warnings
            warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            
            % Go through each strip of data.
            currentImage = zeros(yImage,xImage);
            tifflib('setDirectory',FileID,i-1);
            for r = 1:rps:yImage
              row_inds = r:min(yImage,r+rps-1);
              stripNum = tifflib('computeStrip',FileID,r)-1;
              currentImage(row_inds,:) = tifflib('readEncodedStrip',FileID,stripNum);
            end
            
            % Rescale data
            if(imgScaling ~= 1 && imgScaling>0)
                tifStack(:,:,i) = imresize(currentImage,[scaledY scaledX]); % Scales image size
            else
                tifStack(:,:,i) = currentImage;
            end
        end
        tifflib('close',FileID);
        disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
        warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
    case 'alternativeMethod'
        % Setup TIF object and Read-In Basic Information
        hTif = Tiff(filePath);
          
        warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
        for i=1:NumberOfImages
            if(useWaitBar)
                waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
            else
                if(mod(i+1,updateFrequency)==0)
                    disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                end
            end

            if(imgScaling ~= 1 && imgScaling>0)
                tifStack(:,:,i) = imresize(hTif.read(),[scaledY scaledX]); % Scales image size
            else
                tifStack(:,:,i) = hTif.read();
            end
            if(i == NumberOfImages)
                hTif.close();
                warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
                disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
            else
                hTif.nextDirectory();                
            end
        end
end

if(useWaitBar),close(h); warning on; end
end
function [val]= FindStimulusParam(fname,paramKey)
    fid=fopen(fname);
    try
        % get the rows for the paramKey
        fileText= textscan(fid, '%s', 'delimiter', '\n');
        fileText= fileText{1};
        fileText= strrep(fileText, ' ', ''); % delete all whitespace
        keys = strfind(fileText, strcat(paramKey, '='));
        keys= find(~cellfun(@isempty, keys));
        line = fileText(keys(1));
        line = strsplit(char(line), '#');
        line = strsplit(char(line(1)), '=');
        val= str2num(char(line(2)));
        if isempty(val)
            val = char(line(2));
        end
    catch ME
        val= '';
    end
    fclose(fid);
end
function out = Filter2Modified (ker, in)

% ker: filter kernel
% in: input 2-D image
% our: output 2-D image
%
% This 2-dimensional filter routine is a modification of filter2.
% Around edges, filter2 calculates values assuming all values outside image are zero.
% So, if the pixel values does not appoach zero around the edge, 
% filter2 will underestimate the value.
%
% Kenichi Ohki 08/11/04

[xDim, yDim] = size(in);
temp = ones(xDim, yDim);
out = filter2(ker, in);
temp = filter2(ker, temp);
out = out./temp;
end
function params = CalcMapParams(dir)
    % a set of parameters of direction (orientation) selectivity are obtained by vector averaging.
    % dir: input stack of direction or orientation dF images (one image per condition)
    % th: preferred angle obtained by vector averaging, normalized to 0-1.
    % mag: vector magnitude
    % ave_change: average signal change to all directions
    % max_change: max signal change to any directions
    % tune: vector sum / scalar sum of signal changes to all directions, according to Bonhoeffer et al. (1995)
    % this function can be used to estimate parameters of both orientation & direction selectivity.
    % if input is dF, mag, ave_change and max_change will be absolute change.
    % if input is ratio(dF/F), mag, ave_change and max_change will be ratio change.
    % th & tune will be the same, regardless the input is dF or ratio (dF/F).
    % when the signal change to one direction is negative, it is replaced by zero, because tune becomes strange.
    % thus, tune parameter is between 0 and 1.
    dim = size(dir);
    xsize = dim(1);
    ysize = dim(2);
    ndir = dim(3);

    params.th = zeros(xsize, ysize);
    params.mag = zeros(xsize, ysize);
    params.ave_change = zeros(xsize, ysize);
    params.max_change = zeros(xsize, ysize);
    params.tune = zeros(xsize, ysize);

    a = zeros(ndir,1);
    for x = 1:xsize
        for y = 1:ysize
            Vx = 0;
            Vy = 0;
            sum = 0;
            % vector averageing
            for i = 1:ndir
                if dir(x,y,i) < 0 
                    a(i) = 0;
                else
                    a(i) = dir(x,y,i); 
                end
                Vx = Vx + a(i) * cos(2 * (i - 1) * pi / ndir);
                Vy = Vy + a(i) * sin(2 * (i - 1) * pi / ndir);
                sum = sum + a(i);
            end
            params.th(x,y) = atan2(Vy, Vx);
            params.mag(x,y) = (Vx ^ 2 + Vy ^2 ) ^ 0.5;
            params.ave_change(x,y) = sum ./ ndir;
            params.max_change(x,y) = max(a);
            if sum~=0 
                params.tune(x,y) = params.mag(x,y) / sum; 
            else
                params.tune(x,y) = 0;
            end
        end
    end
    params.th = mod(params.th, 2 * pi)/2/pi; 
end
function [HLS, HLS_hc] = GenerateHLSMap(params, max, varargin)
    % Generate a HLS map
    % params.th: preferred angle obtained by vector averaging, normalized to 0-1 
    %       (result of e.g. CalcMapParams).
    % params.mag: vector magnitude
    % max: maximum value of dF images
    %
    % optional parameters
    % 'ColorShift', value
    % 'OutputPath', 'pathname' (in this case angle map will be saved)
    % 'maxtuning', value - scale to max tuning value [0, 1]
    %
    % output is a HLS map whereas
    % hue: preferred angle
    % intensity: max_change
    % saturation: tune (0-tune_max)
    %
    % HLS_hc: high contrast HLS map

    colorShift = 60/360; % this gives blue as horizontal orientation, assuming that the hue_table starts with red by default
    tune_max = 1;

    % process optional parameters
    if ~isempty(varargin)
        numIndex = find(cellfun('isclass', varargin(1:end-1), 'char'));
        for ind = 1:length(numIndex)
            switch lower(varargin{numIndex(ind)})
                case 'colorshift'
                    colorShift = varargin{numIndex(ind) + 1};
                case 'maxtuning' % should be [0, 1]
                    tune_max = varargin{numIndex(ind) + 1};
            end
        end
    end

    dim = size(params.th);
    HLS = zeros(dim(1),dim(2),3);
    HLS_hc = zeros(dim(1),dim(2),3);

    params.th = mod(params.th + colorShift, 1 );
    if (max > 0)
        for x = 1:dim(1)
            for y = 1:dim(2)
                m1 = params.max_change(x,y) ./ max;
                m2 = params.max_change(x,y) * 2 ./ max;
                s = params.tune(x,y) ./ tune_max;
                if m1 > 1; m1 = 1; end
                if m2 > 1; m2 = 1; end
                if s > 1; s = 1; end
                HLS(x,y,:) = HsvToRgb(params.th(x,y),s,m1)';
                HLS_hc(x,y,:) = HsvToRgb(params.th(x,y),s,m2)';
            end
        end
    end
end
function rgb = HsvToRgb (hue, saturation, lightness)

    if hue < 0
        hue = 0;
    else
        if hue > 1;
            hue = 1;
        end
    end

    if saturation < 0
        saturation = 0;
    else
        if saturation > 1;
            saturation = 1;
        end
    end

    if lightness < 0
        lightness = 0;
    else
        if lightness > 1;
            lightness = 1;
        end
    end


    white = [1;1;1];

    red = [1;0;0];
    yellow = [1;1;0];
    green= [0;1;0];
    cyan = [0;1;1];
    blue = [0;0;1];
    violet = [1;0;1];

    nColors = 6;

    if hue < 1/nColors
        rgb = red * (1-hue*nColors) + yellow * hue*nColors;  %red yellow
    else if hue < 2/nColors
            hue = hue - 1/nColors;
            rgb = yellow * (1-hue*nColors) + green * hue*nColors; % yellow green
        else if hue < 3/nColors
                hue = hue - 2/nColors;
                rgb = green * (1-hue*nColors) + cyan * hue*nColors; % green cyan
            else if hue < 4/nColors;
                    hue = hue - 3/nColors;
                    rgb = cyan * (1-hue*nColors) + blue * hue*nColors; % cyan blue
                else if hue < 5/nColors
                        hue = hue - 4/nColors;
                        rgb = blue * (1-hue*nColors) + violet * hue*nColors; % blue violet
                    else
                        hue = hue - 5/nColors;
                        rgb = violet * (1-hue*nColors) + red * hue*nColors; % blue red
                    end
                end
            end
        end
    end
    rgb= rgb.*saturation + white.*(1-saturation);
    rgb= rgb.*lightness;
end
