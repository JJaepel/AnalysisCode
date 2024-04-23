function plotFieldsOnCellTracing(dendrites, ROIs, IDs, figNR, color, plotOther)

if nargin < 6
    plotOther = 1;
    if nargin < 5
        color = 'red';
        if nargin < 4
            figNR = 1;
        end
    end
end

figure(figNR)
set(gca, 'YDir','reverse')
hold on

%plot cell tracing
for allD = 1:length(dendrites)
    plot(dendrites(allD).pixelCoord(:,1) ,dendrites(allD).pixelCoord(:,2), 'LineWidth',1.5, 'color', 'black');
    hold on
end

%plot ROIs of spines
for r = 1:length(ROIs)
    xpos= ROIs(r).xPos;
    ypos= ROIs(r).yPos;
    if IDs(r) == 1
        plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', color);
    else
        if plotOther
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end

axis off
set(gcf, 'color', 'w');