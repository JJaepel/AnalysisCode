function PlotConfInts( x, y, ci, color, w, Type, Direction )
    
    w = w./2;
    if isempty(ci)
        disp(['Cannot display confidence intervals. Content=' num2str(ci)]);
        return;
    end
    
    if strcmpi( Type, 'stderr' )
        if strcmpi( Direction, 'Both' )
            for i = 1:length(x)
                line([x(i)   x(i)  ],[y(i)-ci(i) y(i)+ci(i)], 'color', color );
                line([x(i)-w x(i)+w],[y(i)-ci(i) y(i)-ci(i)],'color', color );
                line([x(i)-w x(i)+w],[y(i)+ci(i) y(i)+ci(i)],'color', color );
            end
        elseif strcmpi( Direction, 'Patch' )
            patch( [x x(end:-1:1)], [y-ci y(end:-1:1)+ci(end:-1:1)], ...
                color, 'EdgeColor', 'None' );
        elseif strcmpi( Direction, 'Up' )
            for i = 1:length(x)
                line([x(i)   x(i)  ],[y(i)       y(i)+ci(i)], 'color', color );
                line([x(i)-w x(i)+w],[y(i)+ci(i) y(i)+ci(i)],'color', color );
            end
        elseif strcmpi( Direction, 'Down' )
            for i = 1:length(x)
                line([x(i)   x(i)  ],[y(i)-ci(i) y(i)      ], 'color', color );
                line([x(i)-w x(i)+w],[y(i)-ci(i) y(i)-ci(i)],'color', color );
            end
        end
    elseif strcmpi( Type, 'ci' )
        for i = 1:length(x)
            line([x(i)   x(i)  ],[ci(i,1) ci(i,2)], 'color', color );
            line([x(i)-w x(i)+w],[ci(i,1) ci(i,1)],'color', color );
            line([x(i)-w x(i)+w],[ci(i,2) ci(i,2)],'color', color );
        end
    end
    
end

