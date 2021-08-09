coc_prop = cbrewer('qual', 'Paired', 12);

ferret_V1{1}  = 'F2290_2019-01-17'; expID_V1{1}  = 't00038';
ferret_V1{2}  = 'F2360_2019-09-10'; expID_V1{2}  = 't00051';

ferret_V3{1}  = 'F2275_2018-12-05'; expID_V3{1}  = 't00048';
ferret_V3{2}  = 'F2278_2018-12-11'; expID_V3{2}  = 't00039';
ferret_V3{3}  = 'F2290_2019-01-17'; expID_V3{3}  = 't00026';
ferret_V3{4}  = 'F2359_2019-09-06'; expID_V3{4}  = 't00002';

adata_dir = 'F:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep 'V1_V3_Retinotopty' filesep];
if ~exist(save_dir)
    mkdir(save_dir)
end

%% load all data into two master files
for il = 1:length(ferret_V1)
    datapath = [adata_dir ferret_V1{il} filesep expID_V1{il} filesep];
    master_V1{il} = load(fullfile(datapath, 'Patches.mat'));
    master_V1{il}.metadata.ferret = ferret_V1{il};
    master_V1{il}.metadata.expID = expID_V1{il};
end

for il = 1:length(ferret_V3)
    datapath = [adata_dir ferret_V3{il} filesep expID_V3{il} filesep];
    master_V3{il} = load(fullfile(datapath, 'Patches.mat'));
    master_V3{il}.metadata.ferret = ferret_V3{il};
    master_V3{il}.metadata.expID = expID_V3{il};
end

%% data consolidation
All_OFFsize_V1 = [];
All_ONsize_V1 = [];
All_OFFsize_V3 = [];
All_ONsize_V3 = [];

for ferret =1:length(ferret_V1)
    All_OFFsize_V1 = [All_OFFsize_V1 master_V1{ferret}.analysis.roi.OFFsize];
    All_ONsize_V1 = [All_ONsize_V1 master_V1{ferret}.analysis.roi.ONsize];
end

for ferret =1:length(ferret_V3)
    All_OFFsize_V3 = [All_OFFsize_V3 master_V3{ferret}.analysis.roi.OFFsize];
    All_ONsize_V3 = [All_ONsize_V3 master_V3{ferret}.analysis.roi.ONsize];
end

%% plot area size 
All_OFFsize_V1(All_OFFsize_V1 ==0) = NaN;
All_ONsize_V1(All_ONsize_V1 ==0) = NaN;
All_OFFsize_V3(All_OFFsize_V3 ==0) = NaN;
All_ONsize_V3(All_ONsize_V3 ==0) = NaN;
All_OFFsize_V1(All_OFFsize_V1 < 100) = 100;
All_ONsize_V1(All_ONsize_V1 < 100) = 100;
All_OFFsize_V3(All_OFFsize_V3 < 100) = 100;
All_ONsize_V3(All_ONsize_V3 < 100) = 100;

figure 
subplot(1,4,1) %first plot V1 sizes
boxplot([All_ONsize_V1], 'Labels', {'ON, V1'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),coc_prop(5,:));
hold all
ml = findobj(gca, 'Tag', 'Median');
line(get(ml(1),'XData'),get(ml(1),'YData'), 'Color', [0 0 0])
ylabel('Field size (deg^2)');
box off;
ylim([50 350])

subplot(1,4,3) %first plot V1 sizes
boxplot([All_OFFsize_V1], 'Labels', {'OFF, V1'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),coc_prop(1,:));
hold all
ml = findobj(gca, 'Tag', 'Median');
line(get(ml(1),'XData'),get(ml(1),'YData'), 'Color', [0 0 0])
set(gca,'box','off','ycolor','w')
ylim([50 350])

subplot(1,4,2) %first plot V1 sizes
boxplot([All_ONsize_V3], 'Labels', {'ON, V3'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),coc_prop(6,:));
hold all
ml = findobj(gca, 'Tag', 'Median');
line(get(ml(1),'XData'),get(ml(1),'YData'), 'Color', [0 0 0])
set(gca,'box','off','ycolor','w')
ylim([50 350])

subplot(1,4,4) %first plot V1 sizes
boxplot([All_OFFsize_V3], 'Labels', {'OFF, V3'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),coc_prop(2,:));
hold all
ml = findobj(gca, 'Tag', 'Median');
line(get(ml(1),'XData'),get(ml(1),'YData'), 'Color', [0 0 0])
set(gca,'box','off','ycolor','w')
ylim([50 350])
set(gcf, 'color', 'w');

saveas(gcf, fullfile(save_dir, 'Area size.png'))