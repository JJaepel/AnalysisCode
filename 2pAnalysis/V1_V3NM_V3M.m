dcoc_prop = cbrewer('qual', 'Paired', 12);
shufflenum = 100;

ferret_V1{1}  = 'F2363_2019-09-13'; expID_V1{1}  = 't00021';
ferret_V1{2}  = 'F2350_2019-07-10'; expID_V1{2}  = 't00014';
ferret_V1{3}  = 'F2349_2019-07-08'; expID_V1{3}  = 't00018';
ferret_V1{4}  = 'F2336_2019-05-18'; expID_V1{4}  = 't00009';
ferret_V1{5}  = 'F2290_2019-01-17'; expID_V1{5}  = 't00029';
ferret_V1{6}  = 'F2289_2019-01-08'; expID_V1{6}  = 't00020';
ferret_V1{7}  = 'F2289_2019-01-08'; expID_V1{7}  = 't00021';
ferret_V1{8}  = 'F2276_2018-12-14'; expID_V1{8}  = 't00040';

ferret_V3M{1}  = 'F2363_2019-09-13'; expID_V3M{1}  = 't00037';
ferret_V3M{2}  = 'F2363_2019-09-13'; expID_V3M{2}  = 't00038';
ferret_V3M{3}  = 'F2359_2019-09-06'; expID_V3M{3}  = 't00001';
ferret_V3M{4}  = 'F2350_2019-07-10'; expID_V3M{4}  = 't00025';
ferret_V3M{5}  = 'F2349_2019-07-08'; expID_V3M{5}  = 't00013';
ferret_V3M{6}  = 'F2278_2018-12-11'; expID_V3M{6}  = 't00019';
ferret_V3M{7}  = 'F2276_2018-12-14'; expID_V3M{7}  = 't00025';
ferret_V3M{8}  = 'F2275_2018-12-05'; expID_V3M{8}  = 't00024';

ferret_V3NM{1}  = 'F2360_2019-09-10'; expID_V3NM{1}  = 't00065';
ferret_V3NM{2}  = 'F2360_2019-09-10'; expID_V3NM{2}  = 't00072';
ferret_V3NM{3}  = 'F2308_2019-02-08'; expID_V3NM{3}  = 't00028';
ferret_V3NM{4}  = 'F2291_2019-01-15'; expID_V3NM{4}  = 't00021';
ferret_V3NM{5}  = 'F2290_2019-01-17'; expID_V3NM{5}  = 't00015';
ferret_V3NM{6}  = 'F2289_2019-01-08'; expID_V3NM{6}  = 't00019';
%ferret_V3NM{7}  = 'F2244_2018-09-05'; expID_V3NM{7}  = 't00014';
%ferret_V3NM{8}  = 'F2243_2018-08-27'; expID_V3NM{8}  = 't00012';

adata_dir = 'E:\Data\ImageAnalysis\';
save_dir = [adata_dir filesep 'V1_V3NM_V3M_Ori_Grating' filesep];
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

for il = 1:length(ferret_V3NM)
    datapath = [adata_dir ferret_V3NM{il} filesep expID_V3NM{il} filesep];
    master_V3NM{il} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
    master_V3NM{il}.metadata.ferret = ferret_V3NM{il};
    master_V3NM{il}.metadata.expID = expID_V3NM{il};
    resp_V3NM{il} = find([master_V3NM{il}.analysis.dff.roi.isResponseSignificant] == 1);
end
for il = 1:length(ferret_V3M)
    datapath = [adata_dir ferret_V3M{il} filesep expID_V3M{il} filesep];
    master_V3M{il} = load(fullfile(datapath, 's1_ori_Grating_ana.mat'), 'data', 'metadata', 'sliceparams', 'analysis');
    master_V3M{il}.metadata.ferret = ferret_V3M{il};
    master_V3M{il}.metadata.expID = expID_V3M{il};
    resp_V3M{il} = find([master_V3M{il}.analysis.dff.roi.isResponseSignificant] == 1);
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

OSI_V3NM = [];
DSI_V3NM = [];
ori_pair_dist_V3NM = [];
ori_pair_deltaOri_V3NM = [];
isResp_V3NM = [];
HI_V3NM = [];
for ferret =1:length(ferret_V3NM)
    OSI_V3NM = [OSI_V3NM master_V3NM{ferret}.analysis.dff.roi(resp_V3NM{ferret}).OSIFit];
    DSI_V3NM = [DSI_V3NM master_V3NM{ferret}.analysis.dff.roi(resp_V3NM{ferret}).DSI];
    ori_pair_dist_V3NM = [ori_pair_dist_V3NM master_V3NM{ferret}.analysis.dff.ori_pair.distance];
    ori_pair_deltaOri_V3NM = [ori_pair_deltaOri_V3NM master_V3NM{ferret}.analysis.dff.ori_pair.deltaOri];
    isResp_V3NM = [isResp_V3NM master_V3NM{ferret}.analysis.dff.roi.isResponseSignificant];
    HI_V3NM = [HI_V3NM master_V3NM{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,1)' master_V3NM{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,2)'];
end

OSI_V3M = [];
DSI_V3M = [];
ori_pair_dist_V3M = [];
ori_pair_deltaOri_V3M = [];
isResp_V3M = [];
HI_V3M = [];
for ferret =1:length(ferret_V3M)
    OSI_V3M = [OSI_V3M master_V3M{ferret}.analysis.dff.roi(resp_V3M{ferret}).OSIFit];
    DSI_V3M = [DSI_V3M master_V3M{ferret}.analysis.dff.roi(resp_V3M{ferret}).DSI];
    ori_pair_dist_V3M = [ori_pair_dist_V3M master_V3M{ferret}.analysis.dff.ori_pair.distance];
    ori_pair_deltaOri_V3M = [ori_pair_deltaOri_V3M master_V3M{ferret}.analysis.dff.ori_pair.deltaOri];
    isResp_V3M = [isResp_V3M master_V3M{ferret}.analysis.dff.roi.isResponseSignificant];
    HI_V3M = [HI_V3M master_V3M{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,1)' master_V3M{ferret}.analysis.dff.ori_cells.HomeogeneityIndex(:,2)'];
end

%% plot responsiveness and selectivity portions
figure
subplot(3, 3, 1)
all_V1 = length(isResp_V1);
non_resp_V1 = length(find([isResp_V1] == 0)) ./all_V1;
perc_resp_V1 = length(find([isResp_V1] == 1)) ./all_V1;
h = pie([non_resp_V1 perc_resp_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(2,:));
title('V1')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 4)
all_V3NM = length(isResp_V3NM);
non_resp_V3NM = length(find([isResp_V3NM] == 0)) ./all_V3NM;
perc_resp_V3NM = length(find([isResp_V3NM] == 1)) ./all_V3NM;
h = pie([non_resp_V3NM perc_resp_V3NM]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(4,:));
title('V3NM')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 7)
all_V3M = length(isResp_V3M);
non_resp_V3M = length(find([isResp_V3M] == 0)) ./all_V3M;
perc_resp_V3M = length(find([isResp_V3M] == 1)) ./all_V3M;
h = pie([non_resp_V3M perc_resp_V3M]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(10,:));
title('V3M')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 2)
ori_V1 = length(find([OSI_V1] > 0.2)) ./ length(OSI_V1);
non_ori_V1 = 1- ori_V1;
h = pie([non_ori_V1 ori_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(2,:));
title('Orientation-selective V1')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 5)
ori_V3NM = length(find([OSI_V3NM] > 0.2)) ./length(OSI_V3NM);
non_ori_V3NM = 1- ori_V3NM;
h = pie([non_ori_V3NM ori_V3NM]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(4,:));
title('V3NM')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 8)
ori_V3M = length(find([OSI_V3M] > 0.2)) ./length(OSI_V3M);
non_ori_V3M = 1- ori_V3M;
h = pie([non_ori_V3M ori_V3M]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(10,:));
title('V3M')
legend({'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 3)
dir_V1 = length(find([DSI_V1] > 0.2)) ./length(DSI_V1);
non_dir_V1 = 1- dir_V1;
h = pie([non_dir_V1 dir_V1]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(2,:));
title('Direction-selective V1')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(3, 3, 6)
dir_V3NM = length(find([DSI_V3NM] > 0.2)) ./length(DSI_V3NM);
non_dir_V3NM = 1- dir_V3NM;
h = pie([non_dir_V3NM dir_V3NM]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(4,:));
title('V3NM')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');

subplot(3, 3, 9)
dir_V3M = length(find([DSI_V3M] > 0.2)) ./length(DSI_V3M);
non_dir_V3M = 1- dir_V3M;
h = pie([non_dir_V3M dir_V3M]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', [0.5, 0.5, 0.5]);
set(hp(2), 'FaceColor', coc_prop(10,:));
title('V3M')
legend({'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'resp_cells.png'))

%% plot OSI and DSI distribution

edges = linspace(0,1,11);
figure
subplot(2,3,1)
plot(nanmedian(OSI_V1),1.05 * max(histcounts(OSI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(2,:),'MarkerFaceColor',coc_prop(1,:)); hold on
text(nanmedian(OSI_V1),1.2 * max(histcounts(OSI_V1)),num2str(round(100*nanmedian(OSI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(2,:), 'FontSize', 12')
histogram(OSI_V1, edges, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlabel('OSI');
title('V1')
set(gca,'Box','off');

subplot(2,3,2)
plot(nanmedian(OSI_V3NM),1.05 * max(histcounts(OSI_V3NM)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(4,:),'MarkerFaceColor',coc_prop(3,:)); hold on
text(nanmedian(OSI_V3NM),1.2 * max(histcounts(OSI_V3NM)),num2str(round(100*nanmedian(OSI_V3NM))/100),'HorizontalAlignment','center', 'Color', coc_prop(4,:), 'FontSize', 12')
histogram(OSI_V3NM, edges, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlabel('OSI');
title('V3NM')
set(gca,'Box','off');

subplot(2,3,3)
plot(nanmedian(OSI_V3M),1.05 * max(histcounts(OSI_V3M)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(10,:),'MarkerFaceColor',coc_prop(9,:)); hold on
text(nanmedian(OSI_V3M),1.2 * max(histcounts(OSI_V3M)),num2str(round(100*nanmedian(OSI_V3M))/100),'HorizontalAlignment','center', 'Color', coc_prop(10,:), 'FontSize', 12')
histogram(OSI_V3M, edges, 'FaceColor', coc_prop(9,:), 'EdgeColor', coc_prop(10,:));
ylabel('Cells');
xlabel('OSI');
title('V3M')
set(gca,'Box','off');

subplot(2,3,4)
plot(nanmedian(DSI_V1),1.05 * max(histcounts(DSI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(2,:),'MarkerFaceColor',coc_prop(1,:)); hold on
text(nanmedian(DSI_V1),1.2 * max(histcounts(DSI_V1)),num2str(round(100*nanmedian(DSI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(2,:), 'FontSize', 12')
histogram(DSI_V1, edges, 'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlabel('DSI');
title('V1')
set(gca,'Box','off');

subplot(2,3,5)
plot(nanmedian(DSI_V3NM),1.05 * max(histcounts(DSI_V3NM)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(4,:),'MarkerFaceColor',coc_prop(3,:)); hold on
text(nanmedian(DSI_V3NM),1.2 * max(histcounts(DSI_V3NM)),num2str(round(100*nanmedian(DSI_V3NM))/100),'HorizontalAlignment','center', 'Color', coc_prop(4,:), 'FontSize', 12')
histogram(DSI_V3NM, edges, 'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlabel('DSI');
title('V3NM')
set(gca,'Box','off');
set(gcf, 'color', 'w');

subplot(2,3,6)
plot(nanmedian(DSI_V3M),1.05 * max(histcounts(DSI_V3M)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(10,:),'MarkerFaceColor',coc_prop(9,:)); hold on
text(nanmedian(DSI_V3M),1.2 * max(histcounts(DSI_V3M)),num2str(round(100*nanmedian(DSI_V3M))/100),'HorizontalAlignment','center', 'Color', coc_prop(10,:), 'FontSize', 12')
histogram(DSI_V3M, edges, 'FaceColor', coc_prop(9,:), 'EdgeColor', coc_prop(10,:));
ylabel('Cells');
xlabel('DSI');
title('V3M')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_DSI_distribution.png'))

%% deltaOri vs distance
edges = linspace(0, 500, 11);
[n_V1, edges, bin_dist_V1] = histcounts(ori_pair_dist_V1,edges);
[n_V3NM, edges, bin_dist_V3NM] = histcounts(ori_pair_dist_V3NM,edges);
[n_V3M, edges, bin_dist_V3M] = histcounts(ori_pair_dist_V3M,edges);
edge = edges(1:end-1)+25;
mean_deltaOri_V1 = zeros(length(edge),1);
SEM_deltaOri_V1 = zeros(length(edge),1);
mean_deltaOri_V3NM = zeros(length(edge),1);
SEM_deltaOri_V3NM = zeros(length(edge),1);
mean_deltaOri_V3M = zeros(length(edge),1);
SEM_deltaOri_V3M = zeros(length(edge),1);
for bin = 1:length(edge)
    mean_deltaOri_V1(bin) = nanmean(ori_pair_deltaOri_V1(bin_dist_V1 == bin));
    SEM_deltaOri_V1(bin) = nanstd(ori_pair_deltaOri_V1(bin_dist_V1 == bin))/sqrt(n_V1(bin));
    mean_deltaOri_V3NM(bin) = nanmean(ori_pair_deltaOri_V3NM(bin_dist_V3NM == bin));
    SEM_deltaOri_V3NM(bin) = nanstd(ori_pair_deltaOri_V3NM(bin_dist_V3NM == bin))/sqrt(n_V3NM(bin));
    mean_deltaOri_V3M(bin) = nanmean(ori_pair_deltaOri_V3M(bin_dist_V3M == bin));
    SEM_deltaOri_V3M(bin) = nanstd(ori_pair_deltaOri_V3M(bin_dist_V3M == bin))/sqrt(n_V3M(bin));
end

bin_deltaOri_V1_shuffle = zeros(shufflenum,bin);
bin_deltaOri_V3NM_shuffle = zeros(shufflenum,bin);
bin_deltaOri_V3M_shuffle = zeros(shufflenum,bin);
sem_deltaOri_V1_shuffle = zeros(shufflenum,bin);
sem_deltaOri_V3NM_shuffle = zeros(shufflenum,bin);
sem_deltaOri_V3M_shuffle = zeros(shufflenum,bin);
for rep = 1:shufflenum
    delta_ori_shuffle_V1 = ori_pair_deltaOri_V1(randperm(length(ori_pair_deltaOri_V1)));
    delta_ori_shuffle_V3NM = ori_pair_deltaOri_V3NM(randperm(length(ori_pair_deltaOri_V3NM)));
    delta_ori_shuffle_V3M = ori_pair_deltaOri_V3M(randperm(length(ori_pair_deltaOri_V3M)));
    for bin = 1:length(edge)
        bin_deltaOri_V1_shuffle(rep,bin) = nanmean(delta_ori_shuffle_V1(bin_dist_V1 == bin));
        sem_deltaOri_V1_shuffle(rep,bin) = nanstd(delta_ori_shuffle_V1(bin_dist_V1 == bin))/sqrt(n_V1(bin));
        bin_deltaOri_V3NM_shuffle(rep,bin) = nanmean(delta_ori_shuffle_V3NM(bin_dist_V3NM == bin));
        sem_deltaOri_V3NM_shuffle(rep,bin) = nanstd(delta_ori_shuffle_V3NM(bin_dist_V3NM == bin))/sqrt(n_V3NM(bin));
        bin_deltaOri_V3M_shuffle(rep,bin) = nanmean(delta_ori_shuffle_V3M(bin_dist_V3M == bin));
        sem_deltaOri_V3M_shuffle(rep,bin) = nanstd(delta_ori_shuffle_V3M(bin_dist_V3M == bin))/sqrt(n_V3M(bin));
    end
end
mean_deltaOri_V1_shuffle = nanmean(bin_deltaOri_V1_shuffle);
SEM_deltaOri_V1_shuffle = nanmean(sem_deltaOri_V1_shuffle);
mean_deltaOri_V3NM_shuffle = nanmean(bin_deltaOri_V3NM_shuffle);
SEM_deltaOri_V3NM_shuffle = nanmean(sem_deltaOri_V3NM_shuffle);
mean_deltaOri_V3M_shuffle = nanmean(bin_deltaOri_V3M_shuffle);
SEM_deltaOri_V3M_shuffle = nanmean(sem_deltaOri_V3M_shuffle);

figure
errorbar(edge,mean_deltaOri_V1,SEM_deltaOri_V1, 'o-', 'Color', coc_prop(2,:), 'MarkerFaceColor', coc_prop(1,:))
hold all
errorbar(edge,mean_deltaOri_V3NM,SEM_deltaOri_V3NM, 'o-', 'Color', coc_prop(4,:), 'MarkerFaceColor', coc_prop(3,:))
hold all
errorbar(edge,mean_deltaOri_V3M,SEM_deltaOri_V3M, 'o-', 'Color', coc_prop(10,:), 'MarkerFaceColor', coc_prop(9,:))
hold all
errorbar(edge,mean_deltaOri_V1_shuffle,SEM_deltaOri_V1_shuffle, 'o-', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7])
hold all
errorbar(edge,mean_deltaOri_V3NM_shuffle,SEM_deltaOri_V3NM_shuffle, 'o-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
hold all
errorbar(edge,mean_deltaOri_V3M_shuffle,SEM_deltaOri_V3M_shuffle, 'o-', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3])
xlabel('Distance in \mum')
ylabel('\DeltaOrientation preference (\circ)')
ylim([0 90])
xlim([0 500])
legend('V1 Data', 'V3NM Data', 'V3M Data', 'V1 Shuffle', 'V3NM Shuffle', 'V3M Shuffle')
legend('boxoff')
set(gca,'Box','off');
set(gcf, 'color', 'w');
saveas(gcf, fullfile(save_dir, 'OSI_distance.png'))

%% plot HI 
figure
subplot(2,2,1)
plot(nanmedian(HI_V1), 1.05 * max(histcounts(HI_V1)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(2,:),'MarkerFaceColor',coc_prop(1,:)); hold on
text(nanmedian(HI_V1), 1.15 * max(histcounts(HI_V1)),num2str(round(100*nanmedian(HI_V1))/100),'HorizontalAlignment','center', 'Color', coc_prop(2,:), 'FontSize', 12')
histogram(HI_V1,10,'FaceColor', coc_prop(1,:), 'EdgeColor', coc_prop(2,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (100 \mum)');
title('V1')
set(gca,'Box','off');

subplot(2,2,2)
plot(nanmedian(HI_V3NM), 1.05 * max(histcounts(HI_V3NM)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(4,:),'MarkerFaceColor',coc_prop(3,:)); hold on
text(nanmedian(HI_V3NM), 1.15 * max(histcounts(HI_V3NM)),num2str(round(100*nanmedian(HI_V3NM))/100),'HorizontalAlignment','center', 'Color', coc_prop(4,:), 'FontSize', 12')
histogram(HI_V3NM,10,'FaceColor', coc_prop(3,:), 'EdgeColor', coc_prop(4,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (100 \mum)');
title('V3NM')
set(gca,'Box','off');

subplot(2,2,3)
plot(nanmedian(HI_V3M), 1.05 * max(histcounts(HI_V3M)),'v','MarkerSize', 8','MarkerEdgeColor',coc_prop(10,:),'MarkerFaceColor',coc_prop(9,:)); hold on
text(nanmedian(HI_V3M), 1.15 * max(histcounts(HI_V3M)),num2str(round(100*nanmedian(HI_V3M))/100),'HorizontalAlignment','center', 'Color', coc_prop(10,:), 'FontSize', 12')
histogram(HI_V3M,10,'FaceColor', coc_prop(9,:), 'EdgeColor', coc_prop(10,:));
ylabel('Cells');
xlim([0 1])
xlabel('Homeogeneity Index (100 \mum)');
title('V3M')
set(gca,'Box','off');

subplot(2,2,4)
hdl(1) = cdfplot(HI_V1); hold all
hdl(2) = cdfplot(HI_V3NM); hold all
hdl(3) = cdfplot(HI_V3M); hold all
set(hdl(1), 'Color', coc_prop(2,:), 'LineWidth', 2)
set(hdl(2), 'Color', coc_prop(4,:), 'LineWidth', 2)
set(hdl(3), 'Color', coc_prop(10,:), 'LineWidth', 2)
grid off
xlabel('Homogeneity Index')
xlim([0 1])
ylabel('Cumulative fraction of cells')
legend('V1', 'V3NM', 'V3M','Location', 'SouthEast'); legend('boxoff')
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