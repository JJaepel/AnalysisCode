coc_prop = cbrewer('qual', 'Paired', 12);
shufflenum = 100;

ferret_V1{1}  = 'F2363_2019-09-13'; expID_V1{1}  = 't00021';
ferret_V1{2}  = 'F2350_2019-07-10'; expID_V1{2}  = 't00014';
ferret_V1{3}  = 'F2349_2019-07-08'; expID_V1{3}  = 't00018';
ferret_V1{4}  = 'F2336_2019-05-18'; expID_V1{4}  = 't00009';
ferret_V1{5}  = 'F2290_2019-01-17'; expID_V1{5}  = 't00029';
ferret_V1{6}  = 'F2289_2019-01-08'; expID_V1{6}  = 't00020';
ferret_V1{7}  = 'F2289_2019-01-08'; expID_V1{7}  = 't00021';
ferret_V1{8}  = 'F2276_2018-12-14'; expID_V1{8}  = 't00040';

ferret_V3{1}  = 'F2363_2019-09-13'; expID_V3{1}  = 't00037';
ferret_V3{2}  = 'F2363_2019-09-13'; expID_V3{2}  = 't00038';
ferret_V3{3}  = 'F2360_2019-09-10'; expID_V3{3}  = 't00065';
ferret_V3{4}  = 'F2360_2019-09-10'; expID_V3{4}  = 't00072';
ferret_V3{5}  = 'F2359_2019-09-06'; expID_V3{5}  = 't00001';
ferret_V3{6}  = 'F2350_2019-07-10'; expID_V3{6}  = 't00025';
ferret_V3{7}  = 'F2349_2019-07-08'; expID_V3{7}  = 't00013';
ferret_V3{8}  = 'F2308_2019-02-08'; expID_V3{8}  = 't00028';
ferret_V3{9}  = 'F2291_2019-01-15'; expID_V3{9}  = 't00021';
ferret_V3{10} = 'F2290_2019-01-17'; expID_V3{10} = 't00015';
ferret_V3{11} = 'F2289_2019-01-08'; expID_V3{11} = 't00019';
ferret_V3{12} = 'F2278_2018-12-11'; expID_V3{12} = 't00019';
ferret_V3{13} = 'F2276_2018-12-14'; expID_V3{13} = 't00025';
ferret_V3{14} = 'F2275_2018-12-05'; expID_V3{14} = 't00024';
ferret_V3{15} = 'F2244_2018-09-05'; expID_V3{15} = 't00014';
ferret_V3{16} = 'F2243_2018-08-27'; expID_V3{16} = 't00012';

adata_dir = 'F:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep 'V1_V3_Ori_Grating' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

%% load all data into two master files
for il = 1:length(ferret_V1)
    datapath = [adata_dir ferret_V1{il} filesep expID_V1{il} filesep];
    master_V1{il} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
    master_V1{il}.metadata.ferret = ferret_V1{il};
    master_V1{il}.metadata.expID = expID_V1{il};
    resp_V1{il} = find([master_V1{il}.analysis.dff.roi.isResponseSignificant] == 1);
    
end
for il = 1:length(ferret_V3)
    datapath = [adata_dir ferret_V3{il} filesep expID_V3{il} filesep];
    master_V3{il} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
    master_V3{il}.metadata.ferret = ferret_V3{il};
    master_V3{il}.metadata.expID = expID_V3{il};
    resp_V3{il} = find([master_V3{il}.analysis.dff.roi.isResponseSignificant] == 1);
end


%% data consolidation
OSI_V1 = [];
DSI_V1 = [];
ori_pair_dist_V1 = [];
ori_pair_deltaOri_V1 = [];
isResp_V1 = [];
HI_V1 = [];
for ferret =1:length(ferret_V1)
    OSI_V1 = [OSI_V1 master_V1{ferret}.analysis.dff.roi(resp_V1{ferret}).OSIFit];
    DSI_V1 = [DSI_V1 master_V1{ferret}.analysis.dff.roi(resp_V1{ferret}).DSI];
    ori_pair_dist_V1 = [ori_pair_dist_V1 master_V1{ferret}.analysis.dff.ori_pair.distance];
    ori_pair_deltaOri_V1 = [ori_pair_deltaOri_V1 master_V1{ferret}.analysis.dff.ori_pair.deltaOri];
    isResp_V1 = [isResp_V1 master_V1{ferret}.analysis.dff.roi.isResponseSignificant];
    HI_V1 = [HI_V1 master_V1{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,1)' master_V1{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,2)'];
end

OSI_V3 = [];
DSI_V3 = [];
ori_pair_dist_V3 = [];
ori_pair_deltaOri_V3 = [];
isResp_V3 = [];
HI_V3 = [];
for ferret =1:length(ferret_V3)
    OSI_V3 = [OSI_V3 master_V3{ferret}.analysis.dff.roi(resp_V3{ferret}).OSIFit];
    DSI_V3 = [DSI_V3 master_V3{ferret}.analysis.dff.roi(resp_V3{ferret}).DSI];
    ori_pair_dist_V3 = [ori_pair_dist_V3 master_V3{ferret}.analysis.dff.ori_pair.distance];
    ori_pair_deltaOri_V3 = [ori_pair_deltaOri_V3 master_V3{ferret}.analysis.dff.ori_pair.deltaOri];
    isResp_V3 = [isResp_V3 master_V3{ferret}.analysis.dff.roi.isResponseSignificant];
    HI_V3 = [HI_V3 master_V3{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,1)' master_V3{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,2)'];
end

%% plot responsiveness and selectivity portions
figure
subplot(2, 3, 1)
all_V1 = length(isResp_V1);
non_resp_V1 = length(find([isResp_V1] == 0)) ./all_V1;
perc_resp_V1 = length(find([isResp_V1] == 1)) ./all_V1;
h = pie([non_resp_V1 perc_resp_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(7,:));
title('Responsive V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 4)
all_V3 = length(isResp_V3);
non_resp_V3 = length(find([isResp_V3] == 0)) ./all_V3;
perc_resp_V3 = length(find([isResp_V3] == 1)) ./all_V3;
h = pie([non_resp_V3 perc_resp_V3]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(8,:));
title('Responsive V3')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 2)
ori_V1 = length(find([OSI_V1] > 0.2)) ./ length(OSI_V1);
non_ori_V1 = 1- ori_V1;
h = pie([non_ori_V1 ori_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(1,:));
title('V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 5)
ori_V3 = length(find([OSI_V3] > 0.2)) ./length(OSI_V3);
non_ori_V3 = 1- ori_V3;
h = pie([non_ori_V3 ori_V3]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(2,:));
title('V3')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 3)
dir_V1 = length(find([DSI_V1] > 0.2)) ./length(DSI_V1);
non_dir_V1 = 1- dir_V1;
h = pie([non_dir_V1 dir_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(3,:));
title('Direction-selective V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 6)
dir_V3 = length(find([DSI_V3] > 0.2)) ./length(DSI_V3);
non_dir_V3 = 1- dir_V3;
h = pie([non_dir_V3 dir_V3]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(4,:));
title('V3')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'resp_cells.png'))

%% plot OSI and DSI distribution
figure
subplot(2,2,1)
plot(nanmedian(OSI_V1),1.05 * max(histcounts(OSI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(1,:),'MarkerFaceColor',coc_prop(1,:)); hold on
text(nanmedian(OSI_V1),1.2 * max(histcounts(OSI_V1)),num2str(round(100*nanmedian(OSI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(1,:), 'FontSize', 12')
histogram(OSI_V1, 10, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(1,:));
ylabel('Cells');
xlabel('OSI');
title('V1')
set(gca,'Box','off');
subplot(2,2,3)
plot(nanmedian(OSI_V3),1.05 * max(histcounts(OSI_V3)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(2,:),'MarkerFaceColor',coc_prop(2,:)); hold on
text(nanmedian(OSI_V3),1.2 * max(histcounts(OSI_V3)),num2str(round(100*nanmedian(OSI_V3))/100),'HorizontalAlignment','center', 'Color', coc_prop(2,:), 'FontSize', 12')
histogram(OSI_V3, 10, 'FaceColor', coc_prop(2,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlabel('OSI');
title('V3')
set(gca,'Box','off');
subplot(2,2,2)
plot(nanmedian(DSI_V1),1.05 * max(histcounts(DSI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(3,:),'MarkerFaceColor',coc_prop(3,:)); hold on
text(nanmedian(DSI_V1),1.2 * max(histcounts(DSI_V1)),num2str(round(100*nanmedian(DSI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(3,:), 'FontSize', 12')
histogram(DSI_V1, 10, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(3,:));
ylabel('Cells');
xlabel('DSI');
title('V1')
set(gca,'Box','off');
subplot(2,2,4)
plot(nanmedian(DSI_V3),1.05 * max(histcounts(DSI_V3)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(4,:),'MarkerFaceColor',coc_prop(4,:)); hold on
text(nanmedian(DSI_V3),1.2 * max(histcounts(DSI_V3)),num2str(round(100*nanmedian(DSI_V3))/100),'HorizontalAlignment','center', 'Color', coc_prop(4,:), 'FontSize', 12')
histogram(DSI_V3, 10, 'FaceColor', coc_prop(4,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlabel('DSI');
title('V3')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_DSI_distribution.png'))

%% deltaOri vs distance
edges = linspace(0, 500, 11);
[n_V1, edges, bin_dist_V1] = histcounts(ori_pair_dist_V1,edges);
[n_V3, edges, bin_dist_V3] = histcounts(ori_pair_dist_V3,edges);
edge = edges(1:end-1)+25;
mean_deltaOri_V1 = zeros(length(edge),1);
SEM_deltaOri_V1 = zeros(length(edge),1);
mean_deltaOri_V3 = zeros(length(edge),1);
SEM_deltaOri_V3= zeros(length(edge),1);
for bin = 1:length(edge)
    mean_deltaOri_V1(bin) = nanmean(ori_pair_deltaOri_V1(bin_dist_V1 == bin));
    SEM_deltaOri_V1(bin) = nanstd(ori_pair_deltaOri_V1(bin_dist_V1 == bin))/sqrt(n_V1(bin));
    mean_deltaOri_V3(bin) = nanmean(ori_pair_deltaOri_V3(bin_dist_V3 == bin));
    SEM_deltaOri_V3(bin) = nanstd(ori_pair_deltaOri_V3(bin_dist_V3 == bin))/sqrt(n_V3(bin));
end

bin_deltaOri_V1_shuffle = zeros(shufflenum,bin);
bin_deltaOri_V3_shuffle = zeros(shufflenum,bin);
sem_deltaOri_V1_shuffle = zeros(shufflenum,bin);
sem_deltaOri_V3_shuffle = zeros(shufflenum,bin);
for rep = 1:shufflenum
    delta_ori_shuffle_V1 = ori_pair_deltaOri_V1(randperm(length(ori_pair_deltaOri_V1)));
    delta_ori_shuffle_V3 = ori_pair_deltaOri_V3(randperm(length(ori_pair_deltaOri_V3)));
    for bin = 1:length(edge)
        bin_deltaOri_V1_shuffle(rep,bin) = nanmean(delta_ori_shuffle_V1(bin_dist_V1 == bin));
        sem_deltaOri_V1_shuffle(rep,bin) = nanstd(delta_ori_shuffle_V1(bin_dist_V1 == bin))/sqrt(n_V1(bin));
        bin_deltaOri_V3_shuffle(rep,bin) = nanmean(delta_ori_shuffle_V3(bin_dist_V3 == bin));
        sem_deltaOri_V3_shuffle(rep,bin) = nanstd(delta_ori_shuffle_V3(bin_dist_V3 == bin))/sqrt(n_V3(bin));
    end
end
mean_deltaOri_V1_shuffle = nanmean(bin_deltaOri_V1_shuffle);
SEM_deltaOri_V1_shuffle = nanmean(sem_deltaOri_V1_shuffle);
mean_deltaOri_V3_shuffle = nanmean(bin_deltaOri_V3_shuffle);
SEM_deltaOri_V3_shuffle = nanmean(sem_deltaOri_V3_shuffle);

figure
errorbar(edge,mean_deltaOri_V1,SEM_deltaOri_V1, 'o-', 'Color', coc_prop(1,:), 'MarkerFaceColor', coc_prop(1,:))
hold all
errorbar(edge,mean_deltaOri_V3,SEM_deltaOri_V3, 'o-', 'Color', coc_prop(2,:), 'MarkerFaceColor', coc_prop(2,:))
hold all
errorbar(edge,mean_deltaOri_V1_shuffle,SEM_deltaOri_V1_shuffle, 'o-', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7])
hold all
errorbar(edge,mean_deltaOri_V3_shuffle,SEM_deltaOri_V3_shuffle, 'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel('Distance in \mum')
ylabel('\DeltaOrientation preference (\circ)')
ylim([0 90])
xlim([0 500])
legend('V1 Data', 'V3 Data', 'V1 Shuffle', 'V3 Shuffle')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_distance.png'))

%% plot HI 
figure
subplot(2,2,1)
plot(nanmedian(HI_V1), 1.05 * max(histcounts(HI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(5,:),'MarkerFaceColor',coc_prop(5,:)); hold on
text(nanmedian(HI_V1), 1.15 * max(histcounts(HI_V1)),num2str(round(100*nanmedian(HI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(5,:), 'FontSize', 12')
histogram(HI_V1,10,'FaceColor', coc_prop(5,:), 'EdgeColor', coc_prop(5,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (100 \mum)');
title('V1')
set(gca,'Box','off');

subplot(2,2,3)
plot(nanmedian(HI_V3), 1.05 * max(histcounts(HI_V3)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(6,:),'MarkerFaceColor',coc_prop(6,:)); hold on
text(nanmedian(HI_V3), 1.15 * max(histcounts(HI_V3)),num2str(round(100*nanmedian(HI_V3))/100),'HorizontalAlignment','center', 'Color', coc_prop(6,:), 'FontSize', 12')
histogram(HI_V3,10,'FaceColor', coc_prop(6,:), 'EdgeColor', coc_prop(5,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (100 \mum)');
title('V3')
set(gca,'Box','off');

subplot(2,2,2)
hdl(1) = cdfplot(HI_V1); hold all
hdl(2) = cdfplot(HI_V3); hold all
set(hdl(1), 'Color', coc_prop(5,:), 'LineWidth', 2)
set(hdl(2), 'Color', coc_prop(6,:), 'LineWidth', 2)
grid off
xlabel('Homogeneity Index')
xlim([0 1])
ylabel('Cumulative fraction of cells')
legend('V1', 'V3', 'Location', 'SouthEast'); legend('boxoff')
set(gca, 'Box', 'off')
set(gcf, 'color', 'w');
title('')
saveas(gcf, fullfile(save_dir, 'HomeogeneityIndex.png'));


function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% 
% [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging), 'qual' (qualitative)
%   - cname: name of colortable. It changes depending on ctype.
%   - ncol:  number of color in the table. It changes according to ctype and
%            cname
%   - interp_method: interpolation method (see interp1.m). Default is "cubic" )
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of 
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types and
%           names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
% ------------------------------
% 18.09.2015  Minor fixes, fixed a bug where the 'spectral' color table did not appear in the preview


    % load colorbrewer data
    load('C:\Users\jaepelj\Dropbox\Work\colorbrewer.mat')
    % initialise the colormap is there are any problems
    colormap=[];
    if (~exist('interp_method', 'var'))
        interp_method='cubic';
    end

    % If no arguments
    if (~exist('ctype', 'var') | ~exist('cname', 'var') | ~exist('ncol', 'var'))
        disp(' ')
        disp('[colormap] = cbrewer(ctype, cname, ncol [, interp_method])')
        disp(' ')
        disp('INPUT:')
        disp('  - ctype: type of color table *seq* (sequential), *div* (divergent), *qual* (qualitative)')
        disp('  - cname: name of colortable. It changes depending on ctype.')
        disp('  - ncol:  number of color in the table. It changes according to ctype and cname')
        disp('  - interp_method:  interpolation method  (see interp1.m). Default is "cubic" )')

        disp(' ')
        disp('Sequential tables:')
        z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
                 'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
        disp(z')     

        disp('Divergent tables:')
        z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
        disp(z')

        disp(' ')
        disp('Qualitative tables:')
        %getfield(colorbrewer, 'qual')
        z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
        disp(z')

        plot_brewer_cmap
        return
    end

    % Verify that the input is appropriate
    ctype_names={'div', 'seq', 'qual'};
    if (~ismember(ctype,ctype_names))
        disp('ctype must be either: *div*, *seq* or *qual*')
        colormap=[];
        return
    end

    if (~isfield(colorbrewer.(ctype),cname))
        disp(['The name of the colortable of type *' ctype '* must be one of the following:'])
        getfield(colorbrewer, ctype)
        colormap=[];
        return
    end

    if (ncol>length(colorbrewer.(ctype).(cname)))
    %     disp(' ')
    %     disp('----------------------------------------------------------------------')
    %     disp(['The maximum number of colors for table *' cname '* is ' num2str(length(colorbrewer.(ctype).(cname)))])
    %     disp(['The new colormap will be extrapolated from these ' num2str(length(colorbrewer.(ctype).(cname))) ' values'])
    %     disp('----------------------------------------------------------------------')
    %     disp(' ')
        cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
        colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
        colormap=colormap./255;
        return
    end

    if (isempty(colorbrewer.(ctype).(cname){ncol}))

        while(isempty(colorbrewer.(ctype).(cname){ncol}))
            ncol=ncol+1;
        end        
        disp(' ')
        disp('----------------------------------------------------------------------')
        disp(['The minimum number of colors for table *' cname '* is ' num2str(ncol)])
        disp('This minimum value shall be defined as ncol instead')
        disp('----------------------------------------------------------------------')
        disp(' ')
    end

    colormap=(colorbrewer.(ctype).(cname){ncol})./255;
 end