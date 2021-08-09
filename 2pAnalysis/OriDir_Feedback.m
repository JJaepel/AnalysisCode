%% switch board for analysis variable
analysisParams = struct;
%which type of stimulus should it run
analysisParams.dataType = 2; %data type: 1 = cells, 2 = axons, 3 = spines
analysisParams.stimType = 1;

%what should it do?
analysisParams.reloadData = 1; %should you reload from suite2p/Miji and do baselining?
analysisParams.reanalyse =1; %should you reanalyse the data or just plot?
analysisParams.select = 1; %load only selected data (1, marked in column run) or all data (0)?
analysisParams.plotROIs = 0;   %should you plot traces for all resp ROIs?
analysisParams.plotRespROIsOnly = 0; %should you also plot traces for all non-resp ROIs?
analysisParams.server = 1; %load from the server (1) or the raid (0)
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
save_dir = [adata_dir filesep 'Feedbac_Axons_Ori_Dir' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

filePath = 'F:\Organization\Animals\';
file = '2pExpByStimulusAxon.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
exp_info = findExpInfo(xls_txt, xls_all);

allExpInd = find(exp_info.run); %all experiments that need to be reanalyzed
temp =strcmp(exp_info.region, 'feedback'); %all feedback exp
Ind = find(temp);

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
for il = 1:length(Ind)
    datapath = [adata_dir char(exp_info.animal{Ind(il)}) filesep char(exp_info.exp_id{Ind(il)}) filesep];
    try
        master{il} = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    catch
        master{il} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'metadata','analysis');
    end
    master{il}.metadata.ferret = char(exp_info.animal{Ind(il)});
    master{il}.metadata.expID = exp_info.exp_id{Ind(il)};
    resp{il} = find([master{il}.analysis.dff.roi.isResponseSignificant] == 1);
    
end
%%
%data consolidation
OSIFit = [];
OriCircVar = [];
cohensD = [];
DSI = [];
DirCircVar = [];
isResp = [];
for ferret =1:length(Ind)
    OSIFit = [OSIFit master{ferret}.analysis.dff.roi(resp{ferret}).OSIFit];
    OriCircVar = [OriCircVar master{ferret}.analysis.dff.roi(resp{ferret}).OriCircVar];
    cohensD = [cohensD master{ferret}.analysis.dff.roi(resp{ferret}).cohensD];
    DSI = [DSI master{ferret}.analysis.dff.roi(resp{ferret}).DSI];
    DirCircVar = [DirCircVar master{ferret}.analysis.dff.roi(resp{ferret}).DirCircVar];
    isResp = [isResp master{ferret}.analysis.dff.roi.isResponseSignificant];
end


%% 2.) Calculate population indices
%% 3.) Plot results
coc_prop = cbrewer('qual', 'Paired', 12);

%% a)responsiveness and selectivity portions
figure
subplot(1, 3, 1)
all = length(isResp);
non_resp = length(find([isResp] == 0)) ./all;
perc_resp = length(find([isResp] == 1)) ./all;
h = pie([non_resp perc_resp]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(7,:));
title('Responsiveness')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1, 3, 2)
ori = length(find([OSIFit] > 0.2)) ./ length(OSIFit);
% ori_V1 = length(find([OriCircVar_V1] > 0.2)) ./ length(OriCircVar_V1);
non_ori = 1- ori;
h = pie([non_ori ori]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(1,:));
title('Ori-selective')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1, 3, 3)
dir = length(find([DSI] > 0.2)) ./length(DSI);
non_dir = 1- dir;
h = pie([non_dir dir]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(3,:));
title('Dir-selective')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'pie_charts.png'))

%% b) OSIFit and DSI distribution
figure
subplot(1,2,1)
plot(nanmedian(OSIFit),1.05 * max(histcounts(OSIFit)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(1,:),'MarkerFaceColor',coc_prop(1,:)); hold on
text(nanmedian(OSIFit),1.1 * max(histcounts(OSIFit)),num2str(round(100*nanmedian(OSIFit))/100),'HorizontalAlignment','center', 'Color', coc_prop(1,:), 'FontSize', 12')
histogram(OSIFit, 20, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(1,:));
ylabel('Cells');
xlabel('OSI');
set(gca,'Box','off');


subplot(1,2,2)
plot(nanmedian(DSI),1.05 * max(histcounts(DSI)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(3,:),'MarkerFaceColor',coc_prop(3,:)); hold on
text(nanmedian(DSI),1.1 * max(histcounts(DSI)),num2str(round(100*nanmedian(DSI))/100),'HorizontalAlignment','center', 'Color', coc_prop(3,:), 'FontSize', 12')
histogram(DSI, 20, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(3,:));
ylabel('Cells');
xlabel('DSI');
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_DSI_distribution.png'))

%% c) circVar & dirCircVar
figure
subplot(1,2,1)
distributionPlot(OSIFit','color', coc_prop(1,:)); hold on
boxplot(OSIFit)
% distributionPlot(OriCircVar_V1','color', coc_prop(1,:)); hold on
% boxplot(OriCircVar_V1,'Label', {'V1'})
ylim([0 1])
xlabel('OSI')
set(gca,'Box','off');

subplot(1,2,2)
% distributionPlot(DirCircVar_V1','color', coc_prop(3,:)); hold on
% boxplot(DirCircVar_V1, 'Label', {'V1'})
distributionPlot(DSI','color', coc_prop(3,:)); hold on
boxplot(DSI)
ylim([0 1])
xlabel('DSI')
set(gca,'Box','off');

set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSIDSIViolin.png'))

%% d) cohensD
figure
hdl(1) = cdfplot(cohensD); hold all
set(hdl(1), 'Color', coc_prop(1,:), 'LineWidth', 4);
hold all
grid off;
xlim([0 5]);
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