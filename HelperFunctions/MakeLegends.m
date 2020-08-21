function MakeLegends(type, dirName)

% make the legends for the map and cell maps
% modified pie chart function
format = '-djpeg';
resolution = [5 5];
fillfactor = 0;

x = 400; % resolution of colours = no. of segments
h = figure;
set(h, 'PaperUnits','centimeters');
set(h, 'Position', [0 0 600 600]);
set(h, 'PaperPosition', [0 0 resolution]);

maxpts = 1000;    % resolution of single segments

cax = axescheck(x);
cax = newplot(cax);
next = lower(get(cax, 'NextPlot'));
hold_state = ishold(cax);

switch type
    case 'dir'
        theta0 = -pi / 4;  % starting point
        hue_table = makeHueTable(x,210/360);
    case 'ori'
        theta0 = pi / 4;  % starting point        
        hue_table = repmat(makeHueTable(x/2,210/360), 2, 1);
end

legend = [];
n = max(1, ceil(maxpts * 1/x));
for i = 1:x
    r = [0; ones(n + 1, 1); 0];
    r2 = [0; repmat(fillfactor, n + 1, 1); 0];
    theta = theta0 - [0; 1/x * (0:n)' / n; 0] * 2 * pi; % the - sign makes the pie chart get filled clockwise
    [xx,yy] = pol2cart(theta, r);
    [xx2,yy2] = pol2cart(theta, r2);
    xx = [xx(2:n+2);xx2(n+2:-1:2)]; yy = [yy(2:n+2);yy2(n+2:-1:2)];
    theta0 = min(theta); % the "outer" bound is the minimum, since we are going clockwise
    legend = [legend, patch('XData', xx, 'YData', yy, 'CData', i * ones(size(xx)), ...
        'FaceColor', hue_table(i,:), 'LineStyle', 'none', 'parent', cax)];
end

if ~hold_state,
    view(cax, 2); set(cax, 'NextPlot', next);
    axis(cax, 'equal', 'off')
end

if ~isempty(dirName)
    switch type
        case 'dir'
            print(format, fullfile(dirName, 'legend_dir'));
        case 'ori'
            print(format, fullfile(dirName, 'legend_ori'));
    end
end

function hue_table = makeHueTable (N, colorShift)

if nargin < 2
    colorShift = 0;
end

shift = N * colorShift;

hue_table = [];

for i=1:N
    hue_table = [hue_table; HsvToRgb(mod(i-1+shift,N)/N,1,1)'];
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