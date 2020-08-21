close all
clear all

ferret{1} = 'F2425_2020-03-05';     expID_Elev{1} = 't00007';  expID_Azi{1} = 't00004';
ferret{2} = 'F2429_2020-02-27';     expID_Elev{2} = 't00005';  expID_Azi{2} = 't00007';
ferret{3} = 'F2363_2019-09-13';     expID_Elev{3} = 't00003';  expID_Azi{3} = 't00004';
ferret{4} = 'F2428_2020-02-20';     expID_Elev{4} = 't00008';  expID_Azi{4} = 't00020';

adata_dir = 'F:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep 'Retinotopy_adult' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

coc_V1 = cbrewer('seq', 'Blues', length(ferret)+1);
coc_V3 = cbrewer('seq', 'Greens', length(ferret)+1);

%% load all data into a master files
for il = 1:length(ferret)
    datapathAzi = [adata_dir ferret{il} filesep expID_Azi{il} filesep];
    datapathElev = [adata_dir ferret{il} filesep expID_Elev{il} filesep];
    master{il}.azimuth = load(fullfile(datapathAzi, 'continousEdge_ana.mat'), 'expParam', 'metadata', 'analysis');
    master{il}.elevation = load(fullfile(datapathElev, 'continousEdge_ana.mat'), 'expParam', 'metadata', 'analysis');
    master{il}.metadata.ferret = ferret{il};
    master{il}.metadata.expID_Azi = expID_Azi{il};
    master{il}.metadata.expID_Elev = expID_Elev{il};
end

%% calculate change in retinotopy position
for il = 1:length(ferret)
    azimuthMinV3 = prctile(master{il}.azimuth.analysis.quant.V3ProfileDeg,1);
    azimuthRangeV3 = prctile(master{il}.azimuth.analysis.quant.V3ProfileDeg,99) - prctile(master{il}.azimuth.analysis.quant.V3ProfileDeg,1);
    master{il}.azimuth.analysis.quant.V3ProfileDegChange = master{il}.azimuth.analysis.quant.V3ProfileDeg - azimuthMinV3;
    master{il}.azimuth.analysis.quant.V3ProfileDegNormalized = master{il}.azimuth.analysis.quant.V3ProfileDegChange ./ azimuthRangeV3;
    
    azimuthMinV1 = prctile(master{il}.azimuth.analysis.quant.V1ProfileDeg,1);
    azimuthRangeV1 = prctile(master{il}.azimuth.analysis.quant.V1ProfileDeg,99) - prctile(master{il}.azimuth.analysis.quant.V1ProfileDeg,1);
    master{il}.azimuth.analysis.quant.V1ProfileDegChange = master{il}.azimuth.analysis.quant.V1ProfileDeg - azimuthMinV1;
    master{il}.azimuth.analysis.quant.V1ProfileDegNormalized = master{il}.azimuth.analysis.quant.V1ProfileDegChange ./ azimuthRangeV1;
    
    elevationMinV3 = prctile(master{il}.elevation.analysis.quant.V3ProfileDeg,1);
    elevationRangeV3 = prctile(master{il}.elevation.analysis.quant.V3ProfileDeg,99) - prctile(master{il}.elevation.analysis.quant.V3ProfileDeg,1);
    master{il}.elevation.analysis.quant.V3ProfileDegChange = master{il}.elevation.analysis.quant.V3ProfileDeg - elevationMinV3;
    master{il}.elevation.analysis.quant.V3ProfileDegNormalized = master{il}.elevation.analysis.quant.V3ProfileDegChange ./ elevationRangeV3;
    
    elevationMinV1 = prctile(master{il}.elevation.analysis.quant.V1ProfileDeg,1);
    elevationRangeV1 = prctile(master{il}.elevation.analysis.quant.V1ProfileDeg,99) - prctile(master{il}.elevation.analysis.quant.V1ProfileDeg,1);
    master{il}.elevation.analysis.quant.V1ProfileDegChange = master{il}.elevation.analysis.quant.V1ProfileDeg - elevationMinV1;
    master{il}.elevation.analysis.quant.V1ProfileDegNormalized = master{il}.elevation.analysis.quant.V1ProfileDegChange ./ elevationRangeV1;
end

%% data plotting
% first azimuth data
figure
for a=1:length(ferret)
    plot(master{a}.azimuth.analysis.quant.lenghtProfile, smoothdata(master{a}.azimuth.analysis.quant.V3ProfileDeg','movmedian',10,'includenan'),'Color',coc_V3(a+1,:))
    hold on
    plot(master{a}.azimuth.analysis.quant.lenghtProfile, smoothdata(master{a}.azimuth.analysis.quant.V1ProfileDeg','movmedian',10,'includenan'),'Color', coc_V1(a+1,:))
    hold on
end

title('Retinotopic position - Azimuth')
ylabel('Position in Deg')
xlabel('Position in \mum')
set(gcf,'color','w');

% elevation data
figure
for a=1:length(ferret)
    plot(master{a}.elevation.analysis.quant.lenghtProfile, smoothdata(master{a}.elevation.analysis.quant.V3ProfileDeg','movmedian',10,'includenan'),'Color',coc_V3(a+1,:))
    hold on
    plot(master{a}.elevation.analysis.quant.lenghtProfile, smoothdata(master{a}.elevation.analysis.quant.V1ProfileDeg','movmedian',10,'includenan'),'Color', coc_V1(a+1,:))
    hold on
end

title('Retinotopic position - Elevation')
ylabel('Position in Deg')
xlabel('Position in \mum')
set(gcf,'color','w');

% now change in Retinotopy
% first azimuth data
figure
for a=1:length(ferret)
    subplot(2,1,2)
    plot(master{a}.azimuth.analysis.quant.lenghtProfile, smoothdata(master{a}.azimuth.analysis.quant.V3ProfileDegChange', 'movmedian',10,'includenan'),'Color',coc_V3(a+1,:))
    hold on
    subplot(2,1,1)
    plot(master{a}.azimuth.analysis.quant.lenghtProfile, smoothdata(master{a}.azimuth.analysis.quant.V1ProfileDegChange', 'movmedian',10,'includenan'),'Color', coc_V1(a+1,:))
    hold on
end

title('Azimuth')
ylabel('Change in Deg - 17/18')
xlabel('Position in \mum')

subplot(2,1,2)
ylabel('Change in Deg - 19')
xlabel('Position in \mum')
set(gcf,'color','w');

% elevation data
figure
for a=1:length(ferret)
    subplot(2,1,2)
    plot(master{a}.elevation.analysis.quant.lenghtProfile, smoothdata(master{a}.elevation.analysis.quant.V3ProfileDegChange', 'movmedian',10,'includenan'),'Color',coc_V3(a+1,:))
    hold on
    subplot(2,1,1)
    plot(master{a}.elevation.analysis.quant.lenghtProfile, smoothdata(master{a}.elevation.analysis.quant.V1ProfileDegChange', 'movmedian',10,'includenan'),'Color', coc_V1(a+1,:))
    hold on
end

title('Elevation')
ylabel('Change in Deg - 17/18')

subplot(2,1,2)
ylabel('Change in Deg - 19')
xlabel('Position in \mum')
set(gcf,'color','w');

% now normalized
% first azimuth data
figure
for a=1:length(ferret)
    subplot(2,1,2)
    plot(master{a}.azimuth.analysis.quant.lenghtProfile, smoothdata(master{a}.azimuth.analysis.quant.V3ProfileDegNormalized', 'movmedian',10,'includenan'),'Color',coc_V3(a+1,:))
    hold on
    subplot(2,1,1)
    plot(master{a}.azimuth.analysis.quant.lenghtProfile, smoothdata(master{a}.azimuth.analysis.quant.V1ProfileDegNormalized', 'movmedian',10,'includenan'),'Color', coc_V1(a+1,:))
    hold on
end

title('Azimuth')
ylabel('% in Deg - 17/18')
xlabel('Position in \mum')

subplot(2,1,2)
ylabel('% in Deg - 19')
xlabel('Position in \mum')
set(gcf,'color','w');

% elevation data
figure
for a=1:length(ferret)
    subplot(2,1,2)
    plot(master{a}.elevation.analysis.quant.lenghtProfile, smoothdata(master{a}.elevation.analysis.quant.V3ProfileDegNormalized', 'movmedian',10,'includenan'),'Color',coc_V3(a+1,:))
    hold on
    subplot(2,1,1)
    plot(master{a}.elevation.analysis.quant.lenghtProfile, smoothdata(master{a}.elevation.analysis.quant.V1ProfileDegNormalized', 'movmedian',10,'includenan'),'Color', coc_V1(a+1,:))
    hold on
end

title('Elevation')
ylabel('% in Deg - 17/18')

subplot(2,1,2)
ylabel('% in Deg - 19')
xlabel('Position in \mum')
set(gcf,'color','w');


%% add ons
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
    load('C:\Users\jaepelj\Dropbox\Work\Code\colorbrewer.mat')
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