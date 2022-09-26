%% switch board for analysis variable
analysisParams = struct;

analysisParams.stimType = 1;  % 1 = driftingGrating, 2 = OriSf, 3 = OriTf
analysisParams.server = 1;
analysisParams.select = 1;
analysisParams.downsample = 2;
analysisParams.bandpass = 1;
analysisParams.intrinsic = 0;
analysisParams.pre = 1;
analysisParams.changeThreshold = 0;
analysisParams.clean = 0;
analysisParams.reload = 1;
analysisParams.register = 0;

analysisParams.field = 'rawF';

cocV3 = cbrewer('seq', 'RdPu',30);
cocNaive = cocV3(13:17,:);
cocEarly = cocV3(19:24,:);
cocAdult = cocV3(26:30,:);
cocRdn = cocV3(6:10,:);

cocV1 = cbrewer('seq', 'PuBuGn',30);
cocNaiveV1 = cocV1(13:17,:);
cocEarlyV1 = cocV1(19:24,:);
cocAdultV1 = cocV1(26:30,:);
cocRdnV1 = cocV1(6:10,:);

%% 0.) Set folders, list all experiments and reanalyze if necessary
adata_dir = 'Z:\Juliane\Data\ImageAnalysis\';
saveDirectory = [adata_dir filesep 'MapDevelopment' filesep];
if ~exist(saveDirectory)
    mkdir(saveDirectory)
end

filePath = 'F:\Organization\Animals\';
file = 'EpiExpByStimulus.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
exp_info = findExpInfo(xls_txt, xls_all);

allExpIncluded = find(cell2mat(exp_info.included));
allExpInd = find(exp_info.run);

NaiveAge =find(cell2mat(exp_info.EO) == 0); %EO+0
EarlyAge = find(cell2mat(exp_info.EO) == 2 | cell2mat(exp_info.EO)== 3); %EO2-3
AdultAge = find(cell2mat(exp_info.EO) > 4); %EO>8

EarlyInd = intersect(allExpIncluded, EarlyAge);
NaiveInd = intersect(allExpIncluded, NaiveAge);
AdultInd = intersect(allExpIncluded, AdultAge);

%reanalyze everything if necessary
if ~analysisParams.select  
    allExpInd = [AdultIndV1];
end

%reanalyze if necessary
for i = allExpInd
    disp(['Currently analyzing: Ferret ' char(exp_info.animal{i}) ', Experiment ' char(exp_info.exp_id{i})])
    analysisParams.animal = char(exp_info.animal{i});
    analysisParams.expID = char(exp_info.exp_series{i});
    analysisParams.sp2ID = char(exp_info.sp2_id{i});
    driftingGrating_Fct(analysisParams);
end

%% 1.) Load all data
disp('loading naive data');
overlapA19Naive = []; overlapV1Naive = []; overlapA19RndNaive = []; overlapV1RndNaive = [];
corrA19Naive = []; corrV1Naive = []; corrA19RndNaive = []; corrV1RndNaive = [];
convFactA19Naive = []; convFactV1Naive = [];
for ilNaive = 1:length(NaiveInd)
    disp(['Loading animal ' char(exp_info.animal{NaiveInd(ilNaive)})]);
    datapath = [adata_dir char(exp_info.animal{NaiveInd(ilNaive)}) filesep char(exp_info.exp_series{NaiveInd(ilNaive)}) filesep];
    animal = load(fullfile(datapath, 'AnaData.mat'), 'analysis');
    overlapA19Naive = [overlapA19Naive; animal.analysis.(analysisParams.field).perOverlapA19(:)];
    overlapA19Naive_animal(ilNaive) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19(:));
    overlapV1Naive = [overlapV1Naive; animal.analysis.(analysisParams.field).perOverlapV1(:)];
    overlapV1Naive_animal(ilNaive) = nanmean(animal.analysis.(analysisParams.field).perOverlapV1(:));    
    overlapA19RndNaive = [overlapA19RndNaive; animal.analysis.(analysisParams.field).perOverlapA19Rnd(:)];
    overlapA19RndNaive_animal(ilNaive) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19Rnd(:));    
    overlapV1RndNaive = [overlapV1RndNaive; animal.analysis.(analysisParams.field).perOverlapV1Rnd(:)];    overlapA19RndNaive_animal(ilNaive) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19Rnd(:));
    overlapV1RndNaive_animal(ilNaive) = nanmean(animal.analysis.(analysisParams.field).perOverlapV1Rnd(:));

    corrA19Naive = [corrA19Naive; animal.analysis.(analysisParams.field).corrMapsA19(:)];
    corrA19Naive_animal(ilNaive) =  nanmean(animal.analysis.(analysisParams.field).corrMapsA19(:));
    corrV1Naive = [corrV1Naive; animal.analysis.(analysisParams.field).corrMapsV1(:)];
    corrV1Naive_animal(ilNaive) =  nanmean(animal.analysis.(analysisParams.field).corrMapsV1(:));
    corrA19RndNaive = [corrA19RndNaive; animal.analysis.(analysisParams.field).corrMapsA19Rnd(:)];
    corrA19RndNaive_animal(ilNaive) =  nanmean(animal.analysis.(analysisParams.field).corrMapsA19Rnd(:));
    corrV1RndNaive = [corrV1RndNaive; animal.analysis.(analysisParams.field).corrMapsV1Rnd(:)];
    corrV1RndNaive_animal(ilNaive) =  nanmean(animal.analysis.(analysisParams.field).corrMapsV1Rnd(:));
    
    convFactA19Naive = [convFactA19Naive; animal.analysis.(analysisParams.field).convFactA19(:)];
    convFactV1Naive = [convFactV1Naive; animal.analysis.(analysisParams.field).convFactV1(:)];
    clear animal
end

disp('loading early data');
overlapA19Early = []; overlapV1Early = []; overlapA19RndEarly = []; overlapV1RndEarly = [];
corrA19Early = []; corrV1Early = []; corrA19RndEarly = []; corrV1RndEarly = [];
convFactA19Early = []; convFactV1Early = [];
for ilEarly = 1:length(EarlyInd)
    disp(['Loading animal ' char(exp_info.animal{EarlyInd(ilEarly)})]);
    datapath = [adata_dir char(exp_info.animal{EarlyInd(ilEarly)}) filesep char(exp_info.exp_series{EarlyInd(ilEarly)}) filesep];
    animal= load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    overlapA19Early = [overlapA19Early; animal.analysis.(analysisParams.field).perOverlapA19(:)];
    overlapA19Early_animal(ilEarly) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19(:));
    overlapV1Early = [overlapV1Early; animal.analysis.(analysisParams.field).perOverlapV1(:)];
    overlapV1Early_animal(ilEarly) = nanmean(animal.analysis.(analysisParams.field).perOverlapV1(:));   
    overlapA19RndEarly = [overlapA19RndEarly; animal.analysis.(analysisParams.field).perOverlapA19Rnd(:)];
    overlapA19RndEarly_animal(ilEarly) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19Rnd(:));    
    overlapV1RndEarly = [overlapV1RndEarly; animal.analysis.(analysisParams.field).perOverlapV1Rnd(:)];
    overlapV1RndEarly_animal(ilEarly) = nanmean(animal.analysis.(analysisParams.field).perOverlapV1Rnd(:));
        
    corrA19Early = [corrA19Early; animal.analysis.(analysisParams.field).corrMapsA19(:)];
    corrA19Early_animal(ilEarly) =  nanmean(animal.analysis.(analysisParams.field).corrMapsA19(:));
    corrV1Early = [corrV1Early; animal.analysis.(analysisParams.field).corrMapsA19(:)];
    corrV1Early_animal(ilEarly) =  nanmean(animal.analysis.(analysisParams.field).corrMapsV1(:));
    corrA19RndEarly = [corrA19RndEarly;animal.analysis.(analysisParams.field).corrMapsA19Rnd(:)];
    corrA19RndEarly_animal(ilEarly) =  nanmean(animal.analysis.(analysisParams.field).corrMapsA19Rnd(:));
    corrV1RndEarly = [corrV1RndEarly; animal.analysis.(analysisParams.field).corrMapsA19Rnd(:)];
    corrV1RndEarly_animal(ilEarly) =  nanmean(animal.analysis.(analysisParams.field).corrMapsV1Rnd(:));

    convFactA19Early = [convFactA19Early; animal.analysis.(analysisParams.field).convFactA19(:)];
    convFactV1Early = [convFactV1Early; animal.analysis.(analysisParams.field).convFactV1(:)];
    clear animal
end

disp('loading adult data');
overlapA19Adult = []; overlapV1Adult = []; overlapA19RndAdult = []; overlapV1RndAdult = [];
corrA19Adult = []; corrV1Adult = []; corrA19RndAdult = []; corrV1RndAdult = [];
convFactA19Adult = []; convFactV1Adult = [];
for ilAdult= 1:length(AdultInd)
    disp(['Loading animal ' char(exp_info.animal{AdultInd(ilAdult)})]);
    datapath = [adata_dir char(exp_info.animal{AdultInd(ilAdult)}) filesep char(exp_info.exp_series{AdultInd(ilAdult)}) filesep];
    animal = load(fullfile(datapath, 'AnaData.mat'), 'metadata','analysis');
    overlapA19Adult = [overlapA19Adult; animal.analysis.(analysisParams.field).perOverlapA19(:)];
    overlapA19Adult_animal(ilAdult) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19(:));
    overlapV1Adult= [overlapV1Adult; animal.analysis.(analysisParams.field).perOverlapV1(:)];
    overlapV1Adult_animal(ilAdult) = nanmean(animal.analysis.(analysisParams.field).perOverlapV1(:));      
    overlapA19RndAdult = [overlapA19RndAdult; animal.analysis.(analysisParams.field).perOverlapA19Rnd(:)];
    overlapA19RndAdult_animal(ilAdult) = nanmean(animal.analysis.(analysisParams.field).perOverlapA19Rnd(:));    
    overlapV1RndAdult= [overlapV1RndAdult; animal.analysis.(analysisParams.field).perOverlapV1Rnd(:)];
    overlapV1RndAdult_animal(ilAdult) = nanmean(animal.analysis.(analysisParams.field).perOverlapV1Rnd(:));

    corrA19Adult = [corrA19Adult; animal.analysis.(analysisParams.field).corrMapsA19(:)];
    corrA19Adult_animal(ilAdult) =  nanmean(animal.analysis.(analysisParams.field).corrMapsA19(:));
    corrV1Adult = [corrV1Adult; animal.analysis.(analysisParams.field).corrMapsV1(:)];
    corrV1Adult_animal(ilAdult) =  nanmean(animal.analysis.(analysisParams.field).corrMapsV1(:));
    corrA19RndAdult = [corrA19RndAdult; animal.analysis.(analysisParams.field).corrMapsA19Rnd(:)];
    corrA19RndAdult_animal(ilAdult) =  nanmean(animal.analysis.(analysisParams.field).corrMapsA19Rnd(:));
    corrV1RndAdult = [corrV1RndAdult; animal.analysis.(analysisParams.field).corrMapsV1Rnd(:)];
    corrV1RndAdult_animal(ilAdult) =  nanmean(animal.analysis.(analysisParams.field).corrMapsV1Rnd(:));

    convFactA19Adult = [convFactA19Adult; animal.analysis.(analysisParams.field).convFactA19(:)];
    convFactV1Adult = [convFactV1Adult; animal.analysis.(analysisParams.field).convFactV1(:)];
end

%% 2.) Plot results
%overlap
figure
subplot(1,2,1)
allOverlap = [overlapA19Naive(:); overlapA19RndNaive(:); overlapA19Early(:); overlapA19RndEarly(:);overlapA19Adult(:); overlapA19RndAdult(:)];
boxHelp = [zeros(length(overlapA19Naive(:)), 1); ones(length(overlapA19RndNaive(:)), 1); 2*ones(length(overlapA19Early(:)), 1); 3*ones(length(overlapA19RndEarly(:)), 1); 4*ones(length(overlapA19Adult(:)), 1);  5*ones(length(overlapA19RndAdult(:)), 1)];
boxplot(allOverlap, boxHelp, 'Labels',{'Naive','Rnd','Early', 'Rnd','Experienced', 'Rnd'})
h = findobj(gca,'Tag','Box');
patch(get(h(6),'XData'),get(h(6),'YData'),cocNaive(4,:),'FaceAlpha',.5);
patch(get(h(5),'XData'),get(h(5),'YData'),cocRdn(4,:),'FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),cocEarly(4,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),cocRdn(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocAdult(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocRdn(4,:),'FaceAlpha',.5);
box off
ylabel('percentage overlap')
title('A19')

subplot(1,2,2)
allOverlap = [overlapV1Naive(:); overlapV1RndNaive(:); overlapV1Early(:); overlapV1RndEarly(:); overlapV1Adult(:); overlapV1RndAdult(:)];
boxHelp = [zeros(length(overlapV1Naive(:)), 1); ones(length(overlapV1RndNaive(:)), 1); 2*ones(length(overlapV1Early(:)), 1); 3*ones(length(overlapV1RndEarly(:)), 1); 4*ones(length(overlapV1Adult(:)), 1);  5*ones(length(overlapV1RndAdult(:)), 1)];
boxplot(allOverlap, boxHelp, 'Labels',{'Naive','Rnd','Early', 'Rnd','Experienced', 'Rnd'})
h = findobj(gca,'Tag','Box');
patch(get(h(6),'XData'),get(h(6),'YData'),cocNaiveV1(4,:),'FaceAlpha',.5);
patch(get(h(5),'XData'),get(h(5),'YData'),cocRdnV1(4,:),'FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),cocEarlyV1(4,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),cocRdnV1(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocAdultV1(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocRdnV1(4,:),'FaceAlpha',.5);
box off
ylabel('percentage overlap')
title('V1')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'ContourOverlap_Development.png'))

%per animal, substracting the random
overlapAnimalA19NaiveNorm = overlapA19Naive_animal-overlapA19RndNaive_animal;
overlapAnimalA19EarlyNorm = overlapA19Early_animal-overlapA19RndEarly_animal;
overlapAnimalA19AdultNorm = overlapA19Adult_animal-overlapA19RndAdult_animal;
overlapAnimalV1NaiveNorm = overlapV1Naive_animal-overlapV1RndNaive_animal;
overlapAnimalV1EarlyNorm = overlapV1Early_animal-overlapV1RndEarly_animal;
overlapAnimalV1AdultNorm = overlapV1Adult_animal-overlapV1RndAdult_animal;

figure
subplot(1,2,1)
allOverlapAnimal = [overlapAnimalA19NaiveNorm(:); overlapAnimalA19EarlyNorm(:); overlapAnimalA19AdultNorm(:)]; 
boxHelp = [zeros(length(overlapAnimalA19NaiveNorm(:)), 1); ones(length(overlapAnimalA19EarlyNorm(:)), 1); 2*ones(length(overlapAnimalA19AdultNorm(:)),1)];
boxplot(allOverlapAnimal, boxHelp, 'Labels',{'Naive','Early', 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(3),'XData'),get(h(3),'YData'),cocNaive(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarly(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocAdult(4,:),'FaceAlpha',.5);
box off
ylabel('Random corrected percentage overlap')
ylim([0 0.25])
title('A19')

subplot(1,2,2)
allOverlapAnimal = [overlapAnimalV1NaiveNorm(:); overlapAnimalV1EarlyNorm(:); overlapAnimalV1AdultNorm(:)]; 
boxHelp = [zeros(length(overlapAnimalV1NaiveNorm(:)), 1); ones(length(overlapAnimalV1EarlyNorm(:)), 1); 2*ones(length(overlapAnimalV1AdultNorm(:)),1)];
boxplot(allOverlapAnimal, boxHelp, 'Labels',{'Naive','Early', 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(3),'XData'),get(h(3),'YData'),cocNaiveV1(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarlyV1(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocAdultV1(4,:),'FaceAlpha',.5);
box off
ylabel('Random corrected percentage overlap')
ylim([0 0.25])
title('V1')
saveas(gcf, fullfile(saveDirectory, 'OverlapCorrected_Development_Animal.png'))

%correlation
figure
subplot(1,2,1)
allCorr = [corrA19Naive(:); corrA19RndNaive(:); corrA19Early(:); corrA19RndEarly(:); corrA19Adult(:); corrA19RndAdult(:)];
boxHelp = [zeros(length(corrA19Naive(:)), 1); ones(length(corrA19RndNaive(:)), 1); 2*ones(length(corrA19Early(:)), 1); 3*ones(length(corrA19RndEarly(:)), 1); 4*ones(length(corrA19Adult(:)), 1);  5*ones(length(corrA19RndAdult(:)), 1)];
boxplot(allCorr, boxHelp, 'Labels',{'Naive','Rnd','Early', 'Rnd','Experienced', 'Rnd'})
h = findobj(gca,'Tag','Box');
patch(get(h(6),'XData'),get(h(6),'YData'),cocNaive(4,:),'FaceAlpha',.5);
patch(get(h(5),'XData'),get(h(5),'YData'),cocRdn(4,:),'FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),cocEarly(4,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),cocRdn(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocAdult(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocRdn(4,:),'FaceAlpha',.5);
box off
ylabel('correlation')
ylim([0 1])
title('A19')

subplot(1,2,2)
allCorr = [corrV1Naive(:); corrV1RndNaive(:); corrV1Early(:);corrV1RndEarly(:); corrV1Adult(:); corrV1RndAdult(:)];
boxHelp = [zeros(length(corrV1Naive(:)), 1); ones(length(corrV1RndNaive(:)), 1); 2*ones(length(corrV1Early(:)), 1); 3*ones(length(corrV1RndEarly(:)), 1); 4*ones(length(corrV1Adult(:)), 1);  5*ones(length(corrV1RndAdult(:)), 1)];
boxplot(allCorr, boxHelp, 'Labels',{'Naive','Rnd','Early', 'Rnd','Experienced', 'Rnd'})
h = findobj(gca,'Tag','Box');
patch(get(h(6),'XData'),get(h(6),'YData'),cocNaiveV1(4,:),'FaceAlpha',.5);
patch(get(h(5),'XData'),get(h(5),'YData'),cocRdnV1(4,:),'FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),cocEarlyV1(4,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),cocRdnV1(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocAdultV1(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocRdnV1(4,:),'FaceAlpha',.5);
box off
ylabel('correlation')
ylim([0 1])
title('V1')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'Correlation_Development.png'))

%correlation as distribution plots
figure
subplot(1,6,1)
distributionPlot(corrA19Naive,'color', cocNaive(4,:)); hold on
boxplot(corrA19Naive,'Label', {'Naive'})
ylim([-0.2 1.2])
set(gca,'Box','off');
ylabel('correlation')
subplot(1,6,2)
distributionPlot(corrA19RndNaive,'color', cocRdn(4,:)); hold on
boxplot(corrA19RndNaive,'Label', {'Rnd'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,3)
distributionPlot(corrA19Early,'color', cocEarly(4,:)); hold on
boxplot(corrA19Early,'Label', {'Early'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,4)
distributionPlot(corrA19RndEarly,'color', cocRdn(4,:)); hold on
boxplot(corrA19RndEarly,'Label', {'Rnd'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,5)
distributionPlot(corrA19Adult,'color', cocAdult(4,:)); hold on
boxplot(corrA19Adult,'Label', {'Adult'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,6)
distributionPlot(corrA19RndAdult,'color', cocRdn(4,:)); hold on
boxplot(corrA19RndAdult,'Label', {'Rnd'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'CorrelationA19_Violins.png'))

figure
subplot(1,6,1)
distributionPlot(corrV1Naive,'color', cocNaiveV1(4,:)); hold on
boxplot(corrV1Naive,'Label', {'Naive'})
ylim([-0.2 1.2])
set(gca,'Box','off');
ylabel('correlation')
subplot(1,6,2)
distributionPlot(corrV1RndNaive,'color', cocRdnV1(4,:)); hold on
boxplot(corrV1RndNaive,'Label', {'Rnd'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,3)
distributionPlot(corrV1Early,'color', cocEarlyV1(4,:)); hold on
boxplot(corrV1Early,'Label', {'Early'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,4)
distributionPlot(corrV1RndEarly,'color', cocRdnV1(4,:)); hold on
boxplot(corrV1RndEarly,'Label', {'Rnd'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,5)
distributionPlot(corrV1Adult,'color', cocAdultV1(4,:)); hold on
boxplot(corrV1Adult,'Label', {'Adult'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
subplot(1,6,6)
distributionPlot(corrV1RndAdult,'color', cocRdnV1(4,:)); hold on
boxplot(corrV1RndAdult,'Label', {'Rnd'})
ylim([-0.2 1.2])
set(gca,'box','off','ycolor','w')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'CorrelationV1_Violins.png'))

% correlation as a difference between random and stim-driven
corrAnimalA19NaiveNorm = corrA19Naive_animal-corrA19RndNaive_animal;
corrAnimalA19EarlyNorm = corrA19Early_animal-corrA19RndEarly_animal;
corrAnimalA19AdultNorm = corrA19Adult_animal-corrA19RndAdult_animal;
corrAnimalV1NaiveNorm = corrV1Naive_animal-corrV1RndNaive_animal;
corrAnimalV1EarlyNorm = corrV1Early_animal-corrV1RndEarly_animal;
corrAnimalV1AdultNorm = corrV1Adult_animal-corrV1RndAdult_animal;

figure
subplot(1,2,1)
allcorrAnimal = [corrAnimalA19NaiveNorm(:); corrAnimalA19EarlyNorm(:); corrAnimalA19AdultNorm(:)]; 
boxHelp = [zeros(length(corrAnimalA19NaiveNorm(:)), 1); ones(length(corrAnimalA19EarlyNorm(:)), 1); 2*ones(length(corrAnimalA19AdultNorm(:)),1)];
boxplot(allcorrAnimal, boxHelp, 'Labels',{'Naive','Early', 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(3),'XData'),get(h(3),'YData'),cocNaive(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarly(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocAdult(4,:),'FaceAlpha',.5);
box off
ylabel('Random corrected correlation')
ylim([0 0.3])
title('A19')

subplot(1,2,2)
allcorrAnimal = [corrAnimalV1NaiveNorm(:); corrAnimalV1EarlyNorm(:); corrAnimalV1AdultNorm(:)]; 
boxHelp = [zeros(length(corrAnimalV1NaiveNorm(:)), 1); ones(length(corrAnimalV1EarlyNorm(:)), 1); 2*ones(length(corrAnimalV1AdultNorm(:)),1)];
boxplot(allcorrAnimal, boxHelp, 'Labels',{'Naive','Early', 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(3),'XData'),get(h(3),'YData'),cocNaiveV1(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarlyV1(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocAdultV1(4,:),'FaceAlpha',.5);
box off
ylabel('Random corrected correlation')
ylim([0 0.3])
title('V1')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'CorrCorrected_Development_Animal.png'))

% conversion factor
figure
subplot(1,2,1)
allCorr = [convFactA19Naive(:); convFactA19Early(:); convFactA19Adult(:)];
boxHelp = [zeros(length(convFactA19Naive(:)), 1); ones(length(convFactA19Early(:)), 1); 2*ones(length(convFactA19Adult(:)), 1)];
boxplot(allCorr, boxHelp, 'Labels',{'Naive','Early', 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(3),'XData'),get(h(3),'YData'),cocNaive(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarly(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocAdult(4,:),'FaceAlpha',.5);
box off
ylabel('Convergence Factor')
ylim([0 0.4])
title('A19')

subplot(1,2,2)
allCorr = [convFactV1Naive(:); convFactV1Early(:); convFactV1Adult(:)];
boxHelp = [zeros(length(convFactV1Naive(:)), 1); ones(length(convFactV1Early(:)), 1); 2*ones(length(convFactV1Adult(:)), 1)];
boxplot(allCorr, boxHelp, 'Labels',{'Naive','Early', 'Experienced'})
h = findobj(gca,'Tag','Box');
patch(get(h(3),'XData'),get(h(3),'YData'),cocNaiveV1(4,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),cocEarlyV1(4,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),cocAdultV1(4,:),'FaceAlpha',.5);
box off
ylabel('Convergence Factor')
ylim([0 0.4])
title('V1')
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'ConvergenceFactor_Development.png'))