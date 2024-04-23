function plotFieldsonCell(template, ROIs, IDs, field, circ, LUT, figNR, plotOther)

if nargin < 7
    plotOther = 1;
    if nargin < 6
        figNR = 1;
        if nargin < 5
            LUT = jet(100);
            if nargin < 4
                circ = 1;
            end
        end
    end
end

figure(figNR)
imagesc(template)
axis image
colormap('gray')
hold on

for r = 1:length(ROIs)
    xpos= ROIs(r).xPos;
    ypos= ROIs(r).yPos;
    if IDs(r) == 1
        fieldContent = ROIs(r).funcData.(field);
        if circ
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT(1+floor(fieldContent),:));
        else
            try 
                plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:));
            catch
                disp(['Fig. ' num2str(figNR) ' , ROI ' num2str(r) ': Value out of range'])
            end
        end
    else
        if plotOther
            plot(xpos,ypos,'ok','MarkerSize',5,'MarkerFaceColor', 'white');
        end
    end
    hold on
end
axis off
set(gcf, 'color', 'w');