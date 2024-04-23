function polarPlotOri(Angles, Data, varargin )
% function PolarPlot( Degrees, Data, maxXY )
% 
% This function draws a polar plot in the shape of a halfplot from 0 to 180
% or -45 to 45
%
% Inputs
% - Degrees:    Array with the angles of the responses in degrees, if the
% minimum is < 0, it assumes that it is drawing from -90 to 90, otherwise
% from 0 to 180
% - Data:       Array with the responses, make sure that 0 response is
% repeated at the end as 180 deg
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
    Degrees = abs( Angles - 360 );
    Degrees = mod( Degrees -270,360 );
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
    line([0 maxXY],[0 0], 'color', AxesColor);
    line([0 0],[-maxXY maxXY], 'color', AxesColor);
    line([cosd(45)*maxXY 0],[sind(45)*maxXY 0], 'color', AxesColor);
    line([sind(45)*maxXY 0],[cosd(135)*maxXY 0], 'color', AxesColor);
    
    % draw tickmarks
    tCoords = linspace(-pi/2,pi/2, 180);
    for v = 1:length(TickValues)
        line( cos(tCoords)*TickValues(v), sin(tCoords)*TickValues(v), ...
            'color', TickMarkColor, 'linewidth', 1);
    end
    
    %add degs
    if AngTicks
        angleNum = linspace(180,0,9);
        angleCoord = linspace(-pi/2,pi/2,9);
        %angleCoordY = linspace(0,pi,9);
        for a = 1:length(angleCoord)
            text(cos(angleCoord(a))*maxXY*1.2, sin(angleCoord(a))*maxXY*1.2, [num2str(angleNum(a)) char(176)])
        end
    end
    % draw polar response patcu
    patch( x, y, [0 0 0], 'facecolor', DataLineColor, 'edgecolor', ...
        DataLineColor, 'linestyle', '-', 'linewidth', 2,  'FaceAlpha',0.8);
    
    % draw polar response curve
    patch( x, y, [0 0 0], 'facecolor', 'none', 'edgecolor', ...
        DataLineColor, 'linestyle', Style, 'linewidth', 2 );
    
    %set axes
    set(gca,'xlim',[0 maxXY]*1.2);
    set(gca,'ylim',[-maxXY maxXY]*1.2);
    axis off;
 
    set(gcf, 'color', 'w');
    
    
    