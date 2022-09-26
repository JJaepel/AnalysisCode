function plotGratingMetrices(analysisParams, analysis, metadata,types, saveDirectory)
%plots summary of all the metrices:
%Figure 1: orietation selectivity - OSI, OSI fit, CircVar, Bandwidth
%Figure 2: direction selectivity - DSI, DirCircVar
%Figure 3: Variability - Cohen's D, Fano Factor, Variability Index
%Figure 4: other - SFSI, SFvar, TSI,...

if types == 1
    selectedCells = linspace(1,length(analysis.(analysisParams.field).roi),length(analysis.(analysisParams.field).roi));
    number = 0;
    ending = 'all';
elseif types == 2
    selectedCells = find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1);
    number =10;
    ending = 'resp';
else
    error('No ROIs specified')
end

if ~isempty(selectedCells)
    figure(201+number)
    subplot(1,4,1)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).OSI]','color', analysisParams.coc_prop(2,:)); hold all
    boxplot([analysis.(analysisParams.field).roi(selectedCells).OSI])
    title('OSI')
    ylim([0 1])
    subplot(1,4,2)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).OSIFit]','color', analysisParams.coc_prop(2,:)); hold all
    boxplot([analysis.(analysisParams.field).roi(selectedCells).OSI])
    title('OSI after fitting')
    ylim([0 1])
    subplot(1,4,3)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).OriCircVar]','color', analysisParams.coc_prop(2,:)); hold all
    boxplot([analysis.(analysisParams.field).roi(selectedCells).OriCircVar])
    title('CircVar')
    ylim([0 1])
    subplot(1,4,4)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).Bandwidth]','color', analysisParams.coc_prop(2,:)); hold all
    boxplot([analysis.(analysisParams.field).roi(selectedCells).Bandwidth])
    title('Bandwidth')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, [num2str(201+number) 'OSIdistribution_' ending '.png']))

    figure(202+number)
    subplot(1,2,1)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).DSI]','color', analysisParams.coc_prop(4,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(selectedCells).cohensD])
    title('DSI')
    subplot(1,2,2)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).DSI]','color', analysisParams.coc_prop(4,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(selectedCells).DirCircVar])
    title('DirCircVar')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, [num2str(201+number) 'DSI_distribution_' ending '.png']))

    figure(203+number)
    subplot(1,3,1)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).cohensD]','color', analysisParams.coc_prop(3,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(selectedCells).cohensD])
    title('cohensD')
    subplot(1,3,2)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).fanoFactor]','color', analysisParams.coc_prop(3,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(selectedCells).fanoFactor])
    title('fanoFactor')
    subplot(1,3,3)
    distributionPlot([analysis.(analysisParams.field).roi(selectedCells).VI]','color', analysisParams.coc_prop(3,:)); hold on
    boxplot([analysis.(analysisParams.field).roi(selectedCells).VI])
    title('Variability Index')
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, [num2str(201+number) 'Variability_' ending '.png']))
    
    if metadata.StimParams.numCon >1
        figure(204+number)
        if metadata.StimParams.numSf >1
            subplot(1,2,1)
            distributionPlot([analysis.(analysisParams.field).roi(selectedCells).SFSI]','color', analysisParams.coc_prop(5,:)); hold on
            boxplot([analysis.(analysisParams.field).roi(selectedCells).SFSI])
            title('SFSI')
            subplot(1,2,2)
            distributionPlot([analysis.(analysisParams.field).roi(selectedCells).SFVar]','color', analysisParams.coc_prop(5,:)); hold on
            boxplot([analysis.(analysisParams.field).roi(selectedCells).SFVar])
            title('SFVar')
            set(gcf, 'color', 'w');
            saveas(gcf, fullfile(saveDirectory, ['SF_' ending '.png']))
        elseif metadata.StimParams.numTf>1
            distributionPlot([analysis.(analysisParams.field).roi(selectedCells).TFSI]','color', analysisParams.coc_prop(5,:)); hold on
            boxplot([analysis.(analysisParams.field).roi(selectedCells).TFSI])
            title('TFSI')
            set(gcf, 'color', 'w');
            saveas(gcf, fullfile(saveDirectory, [num2str(201+number) 'TF_' ending '.png']))
        end
    end
end