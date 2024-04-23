baseDir = 'Z:\Juliane\InputAnalysis\';
saveDir = [baseDir filesep 'Pooled Data\'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%% Quick multi-cell analysis
cocA19 = cbrewer('seq', 'RdPu', 12);
cocV1 = cbrewer('seq', 'YlGn', 12);
cocAll = cbrewer('seq', 'Greys', 12);

%1) Load all pooled spine data
load('Z:/Juliane/InputAnalysis/Pooled Data/SpineData.mat');
%add in a new field input
allSpines = arrayfun(@(x) setfield(x, 'input', zeros(1,1)), allSpines);
cells = {'F2688', 'F2691/Cell right', 'F2712', 'F2715', 'F2726', 'F2815'};
cellNrs = [1, 3, 4, 5, 6, 7];
inputType = [1, 1, 2, 2, 2, 1];
%current cell numbers
% 1 - F2688, 2 - F2691 left, 3 - F2691 right, 4 - F2712, 5 - F2715, 6 - F2726, 7 - F2815

for c = 1:length(cellNrs)
    %load the multimodaldata
    temp = load(['Z:/Juliane/InputAnalysis/' char(cells{c}) '/E - Analysis/04_FunctionalInputReconstruction.mat']);
    InputSpines = temp.Spines;
    %find the inputspines
    inputs = find([InputSpines.Input] == 1);
    %which of allSpines are coming from this cell?
    cellSpines = find([allSpines.CellNr] == cellNrs(c));
    %go through all of them and mark those as inputs
    cellInputSpines = cellSpines(inputs);
    for ci = 1:length(cellInputSpines)
        allSpines(cellInputSpines(ci)).input = inputType(c);
    end
end

%2) Selectors
%which ones are 2p spines
TwoPROIsNR = find([allSpines.TwoPMatch] ==1);

%which of those are responsive
goodROIs = find([allSpines.good] == 1);

%which of those are selective
oriTwoPROIs = find([allSpines.OSI] > 0.1); %all oriselectROIs, size = allfuncSpines
oriSelect = TwoPROIsNR(oriTwoPROIs);
oriGood = intersect(oriSelect, goodROIs);

dirTwoPROIs = find([allSpines.DSIvect] > 0.1); %all oriselectROIs, size = allfuncSpines
dirSelect = TwoPROIsNR(dirTwoPROIs);
dirGood = intersect(dirSelect, goodROIs);

%which ones are A19 spines
A19ROIs = find([allSpines.input] == 1);
V1ROIs = find([allSpines.input] == 2);

%intersect for all 2p inputs
A19Funct = intersect(TwoPROIsNR, A19ROIs);
V1Funct = intersect(TwoPROIsNR, V1ROIs);

%intersect for responsive inputs
goodA19 = intersect(goodROIs, A19ROIs);
goodV1 = intersect(goodROIs, V1ROIs);

%intersect for ori select input
origoodA19 = intersect(oriGood, A19ROIs);
origoodV1 = intersect(oriGood, V1ROIs);

%intersect for dir select input
dirgoodA19 = intersect(dirGood, A19ROIs);
dirgoodV1 = intersect(dirGood, V1ROIs);

%% Plots
%4) Plot deltaOri of all the sum spines
getDeltaOri = @(x) x.funcData.deltaOri;
goodDeltaOri = cell2mat(arrayfun(getDeltaOri, allSpines(oriGood), 'UniformOutput', false))';
A19DeltaOri = cell2mat(arrayfun(getDeltaOri, allSpines(origoodA19), 'UniformOutput', false))';
V1DeltaOri = cell2mat(arrayfun(getDeltaOri, allSpines(origoodV1), 'UniformOutput', false))';

edges = linspace(0, 90, 10);

figure
plot(histcounts(V1DeltaOri,edges,'Normalization', 'probability'), 'color', cocV1(9,:))
hold on
plot(histcounts(A19DeltaOri,edges,'Normalization', 'probability'), 'color', cocA19(9,:))
plot(histcounts(goodDeltaOri,edges,'Normalization', 'probability'), 'color',cocAll(9,:))
xticks(linspace(0, 9, 10))
xticklabels(edges)
xlabel('deltaOri soma - spine')
ylabel('Probability')
legend(['all Spines', 'V1 inputs', 'A19 inputs'])
set(gcf, 'color', 'w');

% Chi2 test for independendce of the first datapoint 
% Step 0: Get counts first and make vectors
AllCounts = histcounts(goodDeltaOri,edges,'Normalization', 'count');
A19Counts = histcounts(A19DeltaOri,edges,'Normalization', 'count');
V1Counts = histcounts(V1DeltaOri,edges,'Normalization', 'count');
AllFirstVsAll = [repmat('1',AllCounts(1),1); repmat('0',sum(AllCounts)-AllCounts(1),1)];
A19FirstVsAll = [repmat('1',A19Counts(1),1); repmat('0',sum(A19Counts)-A19Counts(1),1)];
V1FirstVsAll = [repmat('1',V1Counts(1),1); repmat('0',sum(V1Counts)-V1Counts(1),1)];

% Step 1: do all vs. A19 
typeOfInputAllA19 = [repmat('N', sum(AllCounts), 1); repmat('A', sum(A19Counts),1)];
typeOfOriDeltaAllA19 = [AllFirstVsAll; A19FirstVsAll];
[~,~,pvalDeltaOriAllA19] = crosstab(typeOfInputAllA19,typeOfOriDeltaAllA19);

% Step 2:  all vs. v1
typeOfInputAllV1 = [repmat('N', sum(AllCounts), 1); repmat('V', sum(V1Counts),1)];
typeOfOriDeltaAllV1 = [AllFirstVsAll; V1FirstVsAll];
[~,~,pvalDeltaOriAllV1] = crosstab(typeOfInputAllV1,typeOfOriDeltaAllV1);

% Step 3: a19 vs. v1
typeOfInputA19V1 = [repmat('A', sum(A19Counts),1); repmat('V', sum(V1Counts),1)];
typeOfOriDeltaA19V1 = [A19FirstVsAll;V1FirstVsAll];
[~,~,pvalDeltaOriA19V1] = crosstab(typeOfInputA19V1,typeOfOriDeltaA19V1);

% Step 4: Do Bonferoni correction for 2 comparisons
numPopulations = 3;
pvalDeltaOriAllA19 = pvalDeltaOriAllA19/numPopulations;
pvalDeltaOriAllV1 = pvalDeltaOriAllV1/numPopulations;
pvalDeltaOriA19V1 = pvalDeltaOriA19V1/numPopulations;

% Step 5: Display result
disp(['P value deltaOri All vs. A19: ' num2str(pvalDeltaOriAllA19)]);
disp(['P value deltaOri All vs. V1: ' num2str(pvalDeltaOriAllV1)]);
disp(['P value deltaOri A19 vs. V1: ' num2str(pvalDeltaOriA19V1)]);

%5) Distribution of responsiveness/ori select
figure
subplot(1,3,1)
h = pie([length(oriGood)/length(TwoPROIsNR), length(goodROIs)/length(TwoPROIsNR)-length(oriGood)/length(TwoPROIsNR),1-length(goodROIs)/length(TwoPROIsNR)]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocAll(9,:));
set(hp(2), 'FaceColor', cocAll(7,:));
set(hp(3), 'FaceColor', cocAll(4,:));
title('All spines')
legend({'Ori-select', 'Responsive', 'Non-responsive'}, 'Location', 'southoutside')
legend('boxoff')
subplot(1,3,2)
h = pie([length(origoodA19)/length(A19Funct), length(goodA19)/length(A19Funct)-length(origoodA19)/length(A19Funct),1-length(goodA19)/length(A19Funct)]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocA19(9,:));
set(hp(2), 'FaceColor', cocA19(7,:));
set(hp(3), 'FaceColor', cocA19(4,:));
title('A19 inputs')
legend({'Ori-select', 'Responsive', 'Non-responsive'}, 'Location', 'southoutside')
legend('boxoff')
subplot(1,3,3)
h = pie([length(origoodV1)/length(V1Funct), length(goodV1)/length(V1Funct)-length(origoodV1)/length(V1Funct),1-length(goodV1)/length(V1Funct)]);
hp = findobj(h, 'Type', 'patch');
set(hp(1), 'FaceColor', cocV1(9,:));
set(hp(2), 'FaceColor', cocV1(7,:));
set(hp(3), 'FaceColor', cocV1(4,:));
title('V1 inputs')
legend({'Ori-select', 'Responsive', 'Non-responsive'}, 'Location', 'southoutside')
legend('boxoff')
set(gcf, 'color', 'w');

% Chi2 test for independendce
% Step 0: Get numbers for the individual fractions
noInput = [repmat('N', length(TwoPROIsNR)-length(goodROIs), 1); repmat('G', length(goodROIs)-length(oriGood),1); repmat('O', length(oriGood),1)];
A19Input = [repmat('N', length(A19Funct)-length(goodA19), 1); repmat('G', length(goodA19)-length(origoodA19),1); repmat('O', length(origoodA19),1)];
V1Input = [repmat('N', length(V1Funct)-length(goodV1), 1); repmat('G', length(goodV1)-length(origoodV1),1); repmat('O', length(origoodV1),1)];

% Step 1: do all vs. A19 
typeOfInputAllA19 = [repmat('N', length(TwoPROIsNR), 1); repmat('A', length(A19Funct),1)];
typeOfRespAllA19 = [noInput; A19Input];
[~,~,pvalAllA19] = crosstab(typeOfInputAllA19,typeOfRespAllA19);

% Step 2:  all vs. v1
typeOfInputAllV1 = [repmat('N', length(TwoPROIsNR), 1); repmat('V', length(V1Funct),1)];
typeOfRespAllV1 = [noInput; V1Input];
[~,~,pvalAllV1] = crosstab(typeOfInputAllV1,typeOfRespAllV1);

% Step 3: a19 vs. v1
typeOfInputA19V1 = [repmat('A', length(A19Funct),1); repmat('V', length(V1Funct),1)];
typeOfRespA19V1 = [A19Input;V1Input];
[~,~,pvalA19V1] = crosstab(typeOfInputA19V1,typeOfRespA19V1);

% Step 4: Do Bonferoni correction for 2 comparisons
numPopulations = 3;
pvalAllA19 = pvalAllA19/numPopulations;
pvalAllV1 = pvalAllV1/numPopulations;
pvalA19V1 = pvalA19V1/numPopulations;

% Step 5: Display result
disp(['P value for responder type All vs. A19: ' num2str(pvalAllA19)]);
disp(['P value for responder type All vs. V1: ' num2str(pvalAllV1)]);
disp(['P value for responder type A19 vs. V1: ' num2str(pvalA19V1)]);

%6) Distribution OSI 
getOSI = @(x) x.OSI;
goodOSI = cell2mat(arrayfun(getOSI, allSpines(goodROIs), 'UniformOutput', false))';
A19OSI = cell2mat(arrayfun(getOSI, allSpines(goodA19), 'UniformOutput', false))';
V1OSI = cell2mat(arrayfun(getOSI, allSpines(goodV1), 'UniformOutput', false))';

figure
subplot(1,3,1)
distributionPlot(goodOSI,'color', cocAll(9,:)); hold on
boxplot(goodOSI','Label', {'All spines'})
ylim([0 1])
ylabel('OSI')
set(gca,'Box','off');
subplot(1,3,2)
distributionPlot(A19OSI,'color', cocA19(9,:)); hold on
boxplot(A19OSI','Label', {'A19 Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
subplot(1,3,3)
distributionPlot(V1OSI,'color', cocV1(9,:)); hold on
boxplot(V1OSI','Label', {'V1 Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');

%let's do statistics
maxNGood = max([length(goodV1) length(goodA19) length(goodROIs)]);
OSIall = NaN(maxNGood,3);
OSIall(1:1:length(goodROIs),1) = goodOSI;
OSIall(1:1:length(goodA19),2) = A19OSI;
OSIall(1:1:length(goodV1),3) = V1OSI;
group = {'all', 'A19 inputs', 'V1 inputs'};
[pOSI,~,statsOSI] = kruskalwallis(OSIall, group, 'on');

if pOSI < 0.05
    OSIStatsMulti = multcompare(statsOSI, 'CType','bonferroni', 'display', 'off');
    OSI_AllvsA19 = OSIStatsMulti(1,6);
    OSI_AllvsV1 = OSIStatsMulti(2,6);
    OSI_A19vsV1 = OSIStatsMulti(3,6);
    disp(['OSI: All vs. A19: p = ' num2str(OSI_AllvsA19) ', All vs. V1: p = ' num2str(OSI_AllvsV1) ', A19 vs. V1: p = ' num2str(OSI_A19vsV1)]);
else
    disp('No significant difference in OSI distribution')
end


%7) Distribution DSI
getDSI = @(x) x.DSIvect;
goodDSI = cell2mat(arrayfun(getDSI, allSpines(goodROIs), 'UniformOutput', false))';
A19DSI = cell2mat(arrayfun(getDSI, allSpines(goodA19), 'UniformOutput', false))';
V1DSI = cell2mat(arrayfun(getDSI, allSpines(goodV1), 'UniformOutput', false))';

figure
subplot(1,3,1)
distributionPlot(goodDSI,'color', cocAll(9,:)); hold on
boxplot(goodDSI','Label', {'All spines'})
ylim([0 1])
ylabel('DSI')
set(gca,'Box','off');
subplot(1,3,2)
distributionPlot(A19DSI,'color', cocA19(9,:)); hold on
boxplot(A19DSI','Label', {'A19 Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
subplot(1,3,3)
distributionPlot(V1DSI,'color', cocV1(9,:)); hold on
boxplot(V1DSI','Label', {'Recurrent Inputs'})
ylim([0 1])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');

%let's do statistics
DSIall = NaN(maxNGood,3);
DSIall(1:1:length(goodROIs),1) = goodDSI;
DSIall(1:1:length(goodA19),2) = A19DSI;
DSIall(1:1:length(goodV1),3) = V1DSI;
[pDSI,~,statsDSI] = kruskalwallis(DSIall, group, 'on');

if pDSI < 0.05
    DSIStatsMulti = multcompare(statsDSI, 'CType','bonferroni', 'display', 'off');
    DSI_AllvsA19 = DSIStatsMulti(1,6);
    DSI_AllvsV1 = DSIStatsMulti(2,6);
    DSI_A19vsV1 = DSIStatsMulti(3,6);
    disp(['DSI: All vs. A19: p = ' num2str(DSI_AllvsA19) ', All vs. V1: p = ' num2str(DSI_AllvsV1) ', A19 vs. V1: p = ' num2str(DSI_A19vsV1)]);
else
    disp('No significant difference in DSI distribution')
end

%7) Distribution Bandwidth
getBW = @(x) x.Bandwidth;
goodBW = cell2mat(arrayfun(getBW, allSpines(goodROIs), 'UniformOutput', false))';
A19BW = cell2mat(arrayfun(getBW, allSpines(goodA19), 'UniformOutput', false))';
V1BW = cell2mat(arrayfun(getBW, allSpines(goodV1), 'UniformOutput', false))';

figure
subplot(1,3,1)
distributionPlot(goodBW,'color', cocAll(9,:)); hold on
boxplot(goodBW','Label', {'All spines'})
ylim([0 180])
ylabel('Bandwidth')
set(gca,'Box','off');
subplot(1,3,2)
distributionPlot(A19BW,'color', cocA19(9,:)); hold on
boxplot(A19BW','Label', {'A19 Inputs'})
ylim([0 180])
set(gca,'box','off','ycolor','w')
subplot(1,3,3)
distributionPlot(V1BW,'color', cocV1(9,:)); hold on
boxplot(V1BW','Label', {'Recurrent Inputs'})
ylim([0 180])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');

%let's do statistics
BWall = NaN(maxNGood,3);
BWall(1:1:length(goodROIs),1) = goodBW;
BWall(1:1:length(goodA19),2) = A19BW;
BWall(1:1:length(goodV1),3) = V1BW;
[pBW,~,statsBW] = kruskalwallis(BWall, group, 'on');

if pBW < 0.05
    BWStatsMulti = multcompare(statsBW, 'CType','bonferroni', 'display', 'off');
    BW_AllvsA19 = BWStatsMulti(1,6);
    BW_AllvsV1 = BWStatsMulti(2,6);
    BW_A19vsV1 = BWStatsMulti(3,6);
    disp(['BW: All vs. A19: p = ' num2str(BW_AllvsA19) ', All vs. V1: p = ' num2str(BW_AllvsV1) ', A19 vs. V1: p = ' num2str(BW_A19vsV1)]);
else
    disp('No significant difference in BW distribution')
end

%7) nearest Neighbor with similar pref
getDirNeighbor = @(x) x.nearestSimilarPrefDir;
allDirNeigh = cell2mat(arrayfun(getDirNeighbor, allSpines(dirGood), 'UniformOutput', false))';
A19DirNeigh = cell2mat(arrayfun(getDirNeighbor, allSpines(dirgoodA19), 'UniformOutput', false))';
V1DirNeigh = cell2mat(arrayfun(getDirNeighbor, allSpines(dirgoodV1), 'UniformOutput', false))';

getOriNeighbor = @(x) x.nearestSimilarPrefOri;
allOriNeigh = cell2mat(arrayfun(getOriNeighbor, allSpines(oriGood), 'UniformOutput', false))';
A19OriNeigh = cell2mat(arrayfun(getOriNeighbor, allSpines(origoodA19), 'UniformOutput', false))';
V1OriNeigh = cell2mat(arrayfun(getOriNeighbor, allSpines(origoodV1), 'UniformOutput', false))';


%% Summary
disp(['Number of inputs A19: ' num2str(length(A19ROIs), '%d')])
disp(['Number of 2p matched inputs A19: ' num2str(length(A19Funct), '%d')])
disp(['Number of responsive A19 inputs: ' num2str(length(goodA19), '%d')])
disp(['Number of ori-selective A19 spines: ' num2str(length(origoodA19), '%d')])

disp(['Number of inputs V1: ' num2str(length(V1ROIs), '%d')])
disp(['Number of 2p matched inputs V1: ' num2str(length(V1Funct), '%d')])
disp(['Number of responsive V1 inputs: ' num2str(length(goodV1), '%d')])
disp(['Number of ori-selective V1 spines: ' num2str(length(origoodV1), '%d')])

%% Save
%write to text file
fid = fopen([saveDir filesep 'V1InputsSummary.txt'], 'w');
formatSpec = ['Number of V1 inputs: %d \n'...
    'Number of 2p matched V1 inputs: %d \n' ...
    'Number of responsive V1 inputs: %d \n'...
    'Number of ori-selective V1 spines:: %d'];
fprintf(fid, formatSpec, length(V1ROIs), length(V1Funct),length(goodV1), length(origoodV1));
fclose(fid);

fid = fopen([saveDir filesep 'A19InputsSummary.txt'], 'w');
formatSpec = ['Number of A19 inputs: %d \n'...
    'Number of 2p matched A19 inputs: %d \n' ...
    'Number of responsive A19 inputs: %d \n'...
    'Number of ori-selective A19 spines:: %d'];
fprintf(fid, formatSpec, length(A19ROIs), length(A19Funct),length(goodA19), length(origoodA19));
fclose(fid);

%save variables
%save([saveDir filesep 'SpineData.mat'], 'allSpines', 'allBranches', 'pwMeasures', '-mat') 
