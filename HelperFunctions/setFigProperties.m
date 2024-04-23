function setFigProperties(n, fig)

%automatically sets figure properties
%Input:
% - n: Figure number
% - Fig: structure with figure properties
%   - fsz: FontSize
%   - alw: LineWidth
%   - width: Figure width
%   - height: Figure height
% Output:
% - pre-formated figure

fig.fig_hdl(n) = figure(n);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', fig.fsz, 'LineWidth', fig.alw, 'visible', 'on'); %<- Set properties
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) fig.width*100 fig.width*100 fig.height*100]); %<- Set size
box off;