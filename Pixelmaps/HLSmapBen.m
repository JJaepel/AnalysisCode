colorShift = 45/360; % this gives blue as horizontal orientation, assuming that the hue_table starts with red by default
tune_max = 1;
filter = 1;
load('C:\Users\jaepelj\Dropbox\Work\oriMap.mat');
load('C:\Users\jaepelj\Dropbox\Work\angMap.mat');
dim = size(oriMap);
HLS = zeros(dim(1),dim(2),3);
HLS_hc = zeros(dim(1),dim(2),3);

%calculate magnitudeMap
magnitudeMap = abs(oriMap);
max_ori = max(max(magnitudeMap));
excludingNaNs = ~isnan(magnitudeMap(:));
highClipVal = mean(magnitudeMap(excludingNaNs(:)))+3*std(magnitudeMap(excludingNaNs(:)));
lowClipVal  = mean(magnitudeMap(excludingNaNs(:)))-3*std(magnitudeMap(excludingNaNs(:)));
magnitudeMap(magnitudeMap > highClipVal) = highClipVal;
magnitudeMap(magnitudeMap < lowClipVal ) = lowClipVal;
offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));

%calculate Phase map
phaseMap = angle(oriMap)*180/pi;
normalizedPhaseMap = wrapTo360(phaseMap)./360;
normalizedPhaseMap = mod(normalizedPhaseMap + colorShift,1);
if filter
    lowpass = 5; % um
    pixel_in_x = size(HLS,2);
    fov_in_um = 100;
    lowpass_pix = lowpass * pixel_in_x / fov_in_um;
    kernel_size = ceil(lowpass_pix * 5);
    sp_filter = fspecial('gaussian', kernel_size, lowpass_pix);
    normalizedMagnitudeMap = Filter2Modified(sp_filter, normalizedMagnitudeMap);
    normalizedPhaseMap = Filter2Modified(sp_filter, normalizedPhaseMap);
    angMap = Filter2Modified(sp_filter, angMap);
end
if (max_ori > 0)
    for x = 1:dim(1)
        for y = 1:dim(2)
            m1 = normalizedMagnitudeMap(x,y) ./ max_ori;
            m2 = normalizedMagnitudeMap(x,y) * 2 ./ max_ori;
            s = 1; %correct for later!
            if m1 > 1; m1 = 1; end
            if m2 > 1; m2 = 1; end
            if s > 1; s = 1; end
            HLS(x,y,:) = HsvToRgb(normalizedPhaseMap(x,y),s,m1)';
            HLS_hc(x,y,:) = HsvToRgb(normalizedPhaseMap(x,y),s,m2)';
        end
    end
end


    

figure
imshow(HLS)
figure
imshow(HLS_hc)

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