function plotPercentagesPiesOri(analysisParams, analysis, saveDirectory)

figure(200)
subplot(1,3,1)
all = length(analysis.(analysisParams.field).roi);
non_resp = length(find([analysis.(analysisParams.field).roi.isResponseSignificant] == 0)) ./all;
resp = length(find([analysis.(analysisParams.field).roi.isResponseSignificant] == 1)) ./all;
h = pie([non_resp resp]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', analysisParams.coc_prop(7,:));
try set(hp(2), 'FaceColor', analysisParams.coc_prop(8,:)); end
title('Responsive')
legend({'Non-resp', 'Resp'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,2)
ori = length(find([analysis.(analysisParams.field).roi.OSIFit] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1)) ./all;
non_ori = 1- ori - non_resp;
h = pie([non_resp non_ori ori]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', analysisParams.coc_prop(7,:));
try set(hp(3), 'FaceColor', analysisParams.coc_prop(2,:)); end
try set(hp(2), 'FaceColor', analysisParams.coc_prop(1,:)); end
title('Orientation-selective')
legend({'Non-resp', 'Non-selective', 'Ori-selective'}, 'Location', 'southoutside')
legend('boxoff')

subplot(1,3,3)
dirResp = length(find([analysis.(analysisParams.field).roi.DSI] > 0.2 & [analysis.(analysisParams.field).roi.isResponseSignificant] == 1)) ./all;
non_dir = 1-dirResp-non_resp;
h = pie([non_resp non_dir dirResp]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', analysisParams.coc_prop(7,:));
try set(hp(3), 'FaceColor', analysisParams.coc_prop(4,:)); end
try set(hp(2), 'FaceColor', analysisParams.coc_prop(3,:)); end
title('Direction-selective')
set(gca, 'box', 'off')
legend({'Non-resp', 'Non-selective', 'Dir-selective'}, 'Location', 'southoutside')
legend('boxoff')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, '200_Piecharts.png'))