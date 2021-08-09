%% switch board for analysis variable
analysisParams = struct;
%which type of stimulus should it run
analysisParams.dataType = 1; %data type: 1 = cells, 2 = axons, 3 = spines
analysisParams.stimType = 1;

%what should it do?
analysisParams.reloadData = 1; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse =1; %should you reanalyse the data or just plot?
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.plotROIs = 0;   %should you plot traces for all resp ROIs?
analysisParams.plotRespROIsOnly = 0; %should you also plot traces for all non-resp ROIs?
analysisParams.server = 0; %load from the server (1) or the raid (0)
analysisParams.makeROIs = 1;

%analysisParameters
analysisParams.zThresh = 4;
analysisParams.fraction = 0.5;
analysisParams.predictor = 0;
analysisParams.shufflenum = 100;
analysisParams.field = 'dff';
analysisParams.windowStart = 0;
analysisParams.windowStop = 2;
analysisParams.pre = 1;

%% 0.) Set folders, list all experiments and reanalyze if necessary

adata_dir = 'F:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep 'V1_V3_Ori_Dir' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

filePath = 'F:\Organization\Animals\';
file = '2pExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
exp_info = findExpInfo(xls_txt, xls_all);

allExpInd = find(exp_info.run); %all experiments that need to be reanalyzed
V1Ind =find(exp_info.region == 1); %all V1 exp
V3Ind = find(exp_info.region == 3); %all V3 exp

%reanalyze if necessary
for i = allExpInd
    disp(['Currently analyzing: Ferret ' char(exp_info.animal{i}) ', Experiment ' char(exp_info.exp_id{i})])
    if exp_info.vol{i} == 1
        analysisParams.level = 1;
    else
        analysisParams.level = 0;
    end
    analysisParams.animal = char(exp_info.animal{i});
    analysisParams.expID = char(exp_info.exp_id{i});
    analysisParams.sp2ID = char(exp_info.sp2_id{i});
    analysisParams.name = char(exp_info.name{i});
    GratingAnalysis(analysisParams);
end

%% 1.) Load all data into two master files and do data consolidation
for ilV1 = 1:length(V1Ind)
    datapath = [adata_dir char(exp_info.animal{V1Ind(ilV1)}) filesep char(exp_info.exp_id{V1Ind(ilV1)}) filesep];
    try
        masterV1{ilV1} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterV1{ilV1} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'metadata','analysis');
    end
    masterV1{ilV1}.metadata.ferret = char(exp_info.animal{V1Ind(ilV1)});
    masterV1{ilV1}.metadata.expID = exp_info.exp_id{V1Ind(ilV1)};
    resp_V1{ilV1} = find([masterV1{ilV1}.analysis.dff.roi.isResponseSignificant] == 1);
    
end
for ilV3 = 1:length(V3Ind)
    datapath = [adata_dir char(exp_info.animal{V3Ind(ilV3)}) filesep char(exp_info.exp_id{V3Ind(ilV3)}) filesep];
    try
        masterV3{ilV3} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        masterV3{ilV3} = load(fullfile(datapath, 's1_ori_Grating.mat'), 'metadata','analysis');
    end
    masterV3{ilV3}.metadata.ferret = char(exp_info.animal{V3Ind(ilV3)});
    masterV3{ilV3}.metadata.expID = exp_info.exp_id{V3Ind(ilV3)};
    resp_V3{ilV3} = find([masterV3{ilV3}.analysis.dff.roi.isResponseSignificant] == 1);
end
%%
%data consolidation
OSIFit_V1 = [];
OriCircVar_V1 = [];
cohensD_V1 = [];
DSI_V1 = [];
DirCircVar_V1 = [];
isResp_V1 = [];
for ferret =1:length(V1Ind)
    OSIFit_V1 = [OSIFit_V1 masterV1{ferret}.analysis.dff.roi(resp_V1{ferret}).OSIFit];
    OriCircVar_V1 = [OriCircVar_V1 masterV1{ferret}.analysis.dff.roi(resp_V1{ferret}).OriCircVar];
    cohensD_V1 = [cohensD_V1 masterV1{ferret}.analysis.dff.roi(resp_V1{ferret}).cohensD];
    DSI_V1 = [DSI_V1 masterV1{ferret}.analysis.dff.roi(resp_V1{ferret}).DSI];
    DirCircVar_V1 = [DirCircVar_V1 masterV1{ferret}.analysis.dff.roi(resp_V1{ferret}).DirCircVar];
    isResp_V1 = [isResp_V1 masterV1{ferret}.analysis.dff.roi.isResponseSignificant];
end

OSIFit_V3 = [];
OriCircVar_V3 = [];
cohensD_V3 = [];
DSI_V3= [];
DirCircVar_V3 = [];
isResp_V3 = [];
for ferret =1:length(V3Ind)
    OSIFit_V3 = [OSIFit_V3 masterV3{ferret}.analysis.dff.roi(resp_V3{ferret}).OSIFit];
    OriCircVar_V3 = [OriCircVar_V3 masterV3{ferret}.analysis.dff.roi(resp_V3{ferret}).OriCircVar];
    cohensD_V3 = [cohensD_V3 masterV3{ferret}.analysis.dff.roi(resp_V3{ferret}).cohensD];
    DSI_V3 = [DSI_V3 masterV3{ferret}.analysis.dff.roi(resp_V3{ferret}).DSI];
    DirCircVar_V3 = [DirCircVar_V3 masterV3{ferret}.analysis.dff.roi(resp_V3{ferret}).DirCircVar];
    isResp_V3 = [isResp_V3 masterV3{ferret}.analysis.dff.roi.isResponseSignificant];
end

%% 2.) Calculate population indices
%% 3.) Plot results
coc_prop = cbrewer('qual', 'Paired', 12);

%% a)responsiveness and selectivity portions
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
ori_V1 = length(find([OSIFit_V1] > 0.2)) ./ length(OSIFit_V1);
% ori_V1 = length(find([OriCircVar_V1] > 0.2)) ./ length(OriCircVar_V1);
non_ori_V1 = 1- ori_V1;
h = pie([non_ori_V1 ori_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(1,:));
title('V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(2, 3, 5)
ori_V3 = length(find([OSIFit_V3] > 0.2)) ./length(OSIFit_V3);
%ori_V3 = length(find([OriCircVar_V3] > 0.2)) ./length(OriCircVar_V3);
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

%% b) OSIFit and DSI distribution
figure
subplot(2,2,1)
plot(nanmedian(OSIFit_V1),1.05 * max(histcounts(OSIFit_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(1,:),'MarkerFaceColor',coc_prop(1,:)); hold on
text(nanmedian(OSIFit_V1),1.2 * max(histcounts(OSIFit_V1)),num2str(round(100*nanmedian(OSIFit_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(1,:), 'FontSize', 12')
histogram(OSIFit_V1, 20, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(1,:));
ylabel('Cells');
xlabel('OSI');
title('V1')
set(gca,'Box','off');
subplot(2,2,3)
plot(nanmedian(OSIFit_V3),1.05 * max(histcounts(OSIFit_V3)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(2,:),'MarkerFaceColor',coc_prop(2,:)); hold on
text(nanmedian(OSIFit_V3),1.2 * max(histcounts(OSIFit_V3)),num2str(round(100*nanmedian(OSIFit_V3))/100),'HorizontalAlignment','center', 'Color', coc_prop(2,:), 'FontSize', 12')
histogram(OSIFit_V3, 20, 'FaceColor', coc_prop(2,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlabel('OSI');
title('V3')
set(gca,'Box','off');
subplot(2,2,2)
plot(nanmedian(DSI_V1),1.05 * max(histcounts(DSI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(3,:),'MarkerFaceColor',coc_prop(3,:)); hold on
text(nanmedian(DSI_V1),1.2 * max(histcounts(DSI_V1)),num2str(round(100*nanmedian(DSI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(3,:), 'FontSize', 12')
histogram(DSI_V1, 20, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(3,:));
ylabel('Cells');
xlabel('DSI');
title('V1')
set(gca,'Box','off');
subplot(2,2,4)
plot(nanmedian(DSI_V3),1.05 * max(histcounts(DSI_V3)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(4,:),'MarkerFaceColor',coc_prop(4,:)); hold on
text(nanmedian(DSI_V3),1.2 * max(histcounts(DSI_V3)),num2str(round(100*nanmedian(DSI_V3))/100),'HorizontalAlignment','center', 'Color', coc_prop(4,:), 'FontSize', 12')
histogram(DSI_V3, 20, 'FaceColor', coc_prop(4,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlabel('DSI');
xlim([0 1])
title('V3')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_DSI_distribution.png'))

%% c) circVar & dirCircVar
figure
subplot(1,4,1)
distributionPlot(OSIFit_V1','color', coc_prop(1,:)); hold on
boxplot(OSIFit_V1,'Label', {'V1'})
% distributionPlot(OriCircVar_V1','color', coc_prop(1,:)); hold on
% boxplot(OriCircVar_V1,'Label', {'V1'})
ylim([0 1])
xlabel('OSI')
set(gca,'Box','off');

subplot(1,4,2)
% distributionPlot(OriCircVar_V3','color', coc_prop(2,:)); hold on
% boxplot(OriCircVar_V3, 'Label', {'V3'})
distributionPlot(OSIFit_V3','color', coc_prop(2,:)); hold on
boxplot(OSIFit_V3, 'Label', {'19'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

subplot(1,4,3)
% distributionPlot(DirCircVar_V1','color', coc_prop(3,:)); hold on
% boxplot(DirCircVar_V1, 'Label', {'V1'})
distributionPlot(DSI_V1','color', coc_prop(3,:)); hold on
boxplot(DSI_V1, 'Label', {'V1'})
ylim([0 1])
xlabel('DSI')
set(gca,'Box','off');

subplot(1,4, 4)
% distributionPlot(DirCircVar_V3','color', coc_prop(4,:)); hold on
% boxplot(DirCircVar_V3, 'Label', {'V3'})
distributionPlot(DSI_V3','color', coc_prop(4,:)); hold on
boxplot(DSI_V3, 'Label', {'19'})
ylim([0 1])
set(gca,'box','off','ycolor','w')

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSIDSIViolin.png'))

[pOriCircVar, ~, stats_pOriCircVar] = ranksum(OriCircVar_V1, OriCircVar_V3);
[pDirCircVar, ~, stats_pDirCircVar] = ranksum(DirCircVar_V1, DirCircVar_V3);
[pOSIFit, ~, stats_pOSIFit] = ranksum(OSIFit_V1, OSIFit_V3);
[pDSI, ~, stats_pDSI] = ranksum(DSI_V1, DSI_V3);

%% d) cohensD
figure
hdl(1) = cdfplot(cohensD_V1); hold all
hdl(2) = cdfplot(cohensD_V3); hold all
set(hdl(1), 'Color', coc_prop(1,:), 'LineWidth', 4);
set(hdl(2), 'Color', coc_prop(2,:), 'LineWidth', 4);
hold all
grid off;
xlim([0 5]);
legend('V1', '19', 'Location', 'SouthEast'); legend('boxoff')
ylabel('Cumulative probability', 'FontSize', 12)
xlabel('Cohens D', 'FontSize', 12);

h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', 'w');
set(gcf, 'color', 'w');
set(gca, 'box', 'off')

saveas(gcf, fullfile(save_dir, 'cohensD.png'))
%% additional functions
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
    load('C:\Users\jaepelj\Documents\GIT\AnalysisCode\HelperFunctions\colorbrewer.mat')
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