function polarPlot( Degrees, Data, varargin )
% function PolarPlot( Degrees, Data, maxXY )
% 
% This function draws a polar plot.
%
% Inputs
% - Degrees:    Array with the angles of the responses in degrees
% - Data:       Array with the responses
% - (maxXY):    Optional 3rd argument to indicate the maximum value of 
%               the polar plot. In case of [], will remain unset
% - (DataLineColor): Color of the data line
% - (Style): LineStyle of the data line
% - (AngTicks): Do you want to plot the angle degrees

    
    if nargin > 2
        maxXY = varargin{1};
    else
        maxXY = [];
    end
    if nargin > 3
        DataLineColor = varargin{2};
    else
        DataLineColor = [1 0 0];
    end
    if nargin > 4
        Style = varargin{3};
    else
        Style = '-';
    end
    if nargin > 5
        AngTicks = varargin{4};
    else
        AngTicks = 1;
    end

    
    AxesColor = [0.5 0.5 0.5];
    TickMarkColor = [0.5 0.5 0.5];
    nTickMarks = 3;
    
    % convert degrees to radians
    Degrees = abs( Degrees - 360 );
    Degrees = mod( Degrees - 270,360 );
    Radians = ang2rad(Degrees);
    
    % get angular coordinates
    for i = 1:length(Degrees)
        x(i) = cos(Radians(i))*Data(i);
        y(i) = sin(Radians(i))*Data(i);
    end
    
    % get axes limits if necessary
    if isempty(maxXY)
        maxXY = max(abs([x y]))*1.3;
    end
    
    % get tickmark locations
    TickStep = (maxXY*0.8) / nTickMarks;
    UpRatio = ceil(-1*log10(TickStep));
    if UpRatio ~= 0
        TickStep = round(TickStep * (10^UpRatio)) / (10^UpRatio);
    end
    TickValues = TickStep:TickStep:(TickStep*nTickMarks);
    
    % draw axes
    hold on;
    axis equal;
    line([-maxXY maxXY],[0 0], 'color', AxesColor);
    line([0 0],[-maxXY maxXY], 'color', AxesColor);
    line([cosd(45)*maxXY cosd(225)*maxXY],[sind(45)*maxXY sind(225)*maxXY], 'color', AxesColor);
    line([cosd(135)*maxXY cosd(315)*maxXY],[sind(135)*maxXY sind(315)*maxXY], 'color', AxesColor);
    
    % draw tickmarks
    tCoords = linspace(0,2*pi,360);
    for v = 1:length(TickValues)
        line( cos(tCoords)*TickValues(v), sin(tCoords)*TickValues(v), ...
            'color', TickMarkColor, 'linewidth', 1);
    end
    
    %add degs
    if AngTicks
        angleCoord = linspace(0.5*pi,2.5*pi,17);
        angleNum = linspace(360,0,17);
        angleNum(1) = 0;
        for a = 1:length(angleCoord)-1
            text(cos(angleCoord(a))*maxXY*1.2, sin(angleCoord(a))*maxXY*1.2, [num2str(angleNum(a)) char(176)])
        end
    end
    
    % draw polar response patcu
    patch( x, y, [0 0 0], 'facecolor', DataLineColor, 'edgecolor', ...
        DataLineColor, 'linestyle', '-', 'linewidth', 2); %,'FaceAlpha',0.8);

    % draw polar response curve
    patch( x, y, [0 0 0], 'facecolor', 'none', 'edgecolor', ...
        DataLineColor, 'linestyle', Style, 'linewidth', 2 );
    
    % set axes
    set(gca,'xlim',[-maxXY maxXY]*1.2);
    set(gca,'ylim',[-maxXY maxXY]*1.2);
    axis off;
 
    set(gcf, 'color', 'w');
end