clear all
%function CombineInputsFunctionalData(animalID)
animalID = 'F2688';
addInputsID{1} = '202'; %add inputs that are directly matched, but not with confocal data
CellNr = 1;

%% 1.)Find all experiments with that animal 

filePath =  'Z:\Juliane\Organization\Animals\';
file = 'SpinePerAnimal.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'bimodal');

exp_info = findExpInfo(xls_txt, xls_all);
ind = 1:1:length(exp_info.animal);

%select experiments containing the animal
onlyAnimal= cellfun(@(x) find(contains(x, char(animalID))),exp_info.animal,'UniformOutput',false);
onlyOthers = cellfun(@isempty,onlyAnimal); 
AnimalExp = setdiff(ind,find(onlyOthers));

%get the cellID
onlyCells = cellfun(@(x) find(contains(x, 'cells')),exp_info.region,'UniformOutput',false);
onlySpines = cellfun(@isempty,onlyCells); 
CellInd = setdiff(ind,find(onlySpines)); 
CellExp = intersect(AnimalExp, CellInd);
SpineExp = setdiff(AnimalExp, CellInd);
saveDir = ['Z:\Juliane\Data\ImageAnalysis\' char(exp_info.animal{CellExp}) filesep];

%% 2.) Load all grouping files and make a structure containing all functional characterized inputs
%load the STED groups
[STEDFile, filePathSTED] = uigetfile('*.mat', 'Please select the STED match file', saveDir);
STEDData = load([filePathSTED STEDFile]);
%find the input confocal IDs
confocalInputs = {STEDData.MatchData(:).ConfocalID};

%load the invivo groups
[InvivoFile, filePath2p] = uigetfile('*.mat', 'Please select the 2p match file', saveDir);
InVivoData = load([filePath2p InvivoFile]);

%look for all confocal Inputs if you find a matching group there as well
confocalGroupsInvivo = {InVivoData.MatchData.ConfocalID};
inputCounter =1;
for ci = 1:size(confocalInputs,2)
    Matches= cellfun(@(x) find(strcmp(x, confocalInputs{ci})),confocalGroupsInvivo,'UniformOutput',false);
    noMatches= cellfun(@isempty,Matches); 
    matchID = find(noMatches==0);
    if ~isempty(matchID)
        inputGroups(inputCounter) = InVivoData.MatchData(matchID);
        inputCounter = inputCounter+1;
    end
end

%add inputs that are matched across STED and 2p if applicable
if ~isempty(addInputsID)
    invivoIDs = {InVivoData.ROIsData2.id};
    for i = 1:length(addInputsID)
        idMatch= cellfun(@(x) find(strcmp(x, addInputsID{i})),invivoIDs,'UniformOutput',false);
        noidMatch= cellfun(@isempty,idMatch); 
        match = find(noidMatch==0);
        inputGroups(inputCounter).TwoP_ID = addInputsID{i};
        inputGroups(inputCounter).expID = InVivoData.ROIsData2(match).exp;
        inputGroups(inputCounter).isResp =  InVivoData.ROIsData2(match).isReps;
        inputGroups(inputCounter).OSI =  InVivoData.ROIsData2(match).OSI;
        inputGroups(inputCounter).prefOri =  InVivoData.ROIsData2(match).prefOri;
        inputGroups(inputCounter).DSI =  InVivoData.ROIsData2(match).DSI;
        inputGroups(inputCounter).prefDir =  InVivoData.ROIsData2(match).prefDir;
        inputCounter = inputCounter+1;
    end
end

%% 3.) Load the cell data and get soma pref
%cellData
loadDirCell = ['Z:\Juliane\Data\ImageAnalysis\' char(exp_info.animal{CellExp}) filesep char(exp_info.name{CellExp}) filesep];
load([loadDirCell 'ROIsAna.mat'], 'ce');
allCellData = ce; clear ce;
CellOriPref = allCellData(CellNr).prefOri;

%% 4.) Now load also the positional information of input ROIs and calculate deltaOri and add tuning curves to file

%load the cell summary containing all measurements
SummaryDir = ['Z:\Juliane\Data\ImageAnalysis\' char(exp_info.animal{CellExp}) filesep];
load([SummaryDir char(animalID) '_CellSummary.mat']);
a = {completeCell{1}.SpineData.expNr};

%input ROIs
for i = 1:size(inputGroups,2)
    idMatch= cellfun(@(x) find(strcmp(x, inputGroups(i).TwoP_ID)),invivoIDs,'UniformOutput',false);
    noidMatch= cellfun(@isempty,idMatch);
    match = find(noidMatch==0);
    inputGroups(i).AbsLeft = InVivoData.ROIsData2(match).AbsLeft;
    inputGroups(i).AbsUp = InVivoData.ROIsData2(match).AbsUp;
    deltaOriAbs = abs(InVivoData.ROIsData2(match).prefOri - CellOriPref);
    if deltaOriAbs >90
         inputGroups(i).deltaOriAbs = 180-deltaOriAbs;
    else
         inputGroups(i).deltaOriAbs = deltaOriAbs;
    end
    deltaOri = InVivoData.ROIsData2(match).prefOri - CellOriPref;
    if deltaOri > 90
        inputGroups(i).deltaOri = 180-deltaOri;
    elseif deltaOri < -90
        inputGroups(i).deltaOri = -180-deltaOri;
    else
        inputGroups(i).deltaOri = deltaOri;
    end    
    %get the tuning curves
    if ~isnan(inputGroups(i).isResp) %only for spines, no dendrites
        expMatch= cellfun(@(x) find(strcmp(x, inputGroups(i).expID)),a,'UniformOutput',false); %find all spines from that experiments
        noexpMatch= cellfun(@isempty,expMatch);
        expMatchID = find(noexpMatch==0);
        ROIID = str2double(inputGroups(i).TwoP_ID(end-1:end)); %find all experiments that have that ROI ID
        ROIIDMatch = find([completeCell{1}.SpineData(:).spineROINr] == ROIID);
        spineID = intersect(expMatchID,ROIIDMatch); %combine those to find the spineID 
        inputGroups(i).meanResp = completeCell{1}.SpineData(spineID(1)).meanResp;
        inputGroups(i).events = completeCell{1}.SpineData(spineID(1)).events;
        inputGroups(i).eventAmps = completeCell{1}.SpineData(spineID(1)).eventAmps;
    end
end

%allROIs
for sp = 1:length(InVivoData.ROIsData2)
    deltaOriAbs = abs(InVivoData.ROIsData2(sp).prefOri - CellOriPref); 
    if deltaOriAbs >90
        InVivoData.ROIsData2(sp).deltaOriAbs = 180-deltaOriAbs;
    else
        InVivoData.ROIsData2(sp).deltaOriAbs = deltaOriAbs;
    end
    
    deltaOri = InVivoData.ROIsData2(match).prefOri - CellOriPref;
    if deltaOri > 90
        InVivoData.ROIsData2(sp).deltaOri = 180-deltaOri;
    elseif deltaOri < -90
        InVivoData.ROIsData2(sp).deltaOri = -180-deltaOri;
    else
        InVivoData.ROIsData2(sp).deltaOri = deltaOri;
    end
    
    InVivoData.ROIsData2(sp).depth = completeCell{1}.SpineData(spineID(1)).depth;
    InVivoData.ROIsData2(sp).denType = completeCell{1}.SpineData(spineID(1)).denType;

    
    %get the tuning curves & other information that you need to sneek in
    if ~isnan(InVivoData.ROIsData2(sp).isReps) %only for spines, no dendrites
        expMatch= cellfun(@(x) find(strcmp(x, InVivoData.ROIsData2(sp).exp)),a,'UniformOutput',false); %find all spines from that experiments
        noexpMatch= cellfun(@isempty,expMatch);
        expMatchID = find(noexpMatch==0);
        ROIIDMatch = find([completeCell{1}.SpineData(:).spineROINr] == InVivoData.ROIsData2(sp).RoiId);%find all experiments that have that ROI ID
        spineID = intersect(expMatchID,ROIIDMatch); %combine those to find the spineID 
        InVivoData.ROIsData2(sp).meanResp = completeCell{1}.SpineData(spineID(1)).meanResp;
        InVivoData.ROIsData2(sp).events = completeCell{1}.SpineData(spineID(1)).events;
        InVivoData.ROIsData2(sp).eventAmps = completeCell{1}.SpineData(spineID(1)).eventAmps;
        
    end
end

%% 5.) Plot the data
%load the cell image
CellImageDir = ['Z:\Juliane\Data\2P_Data\' char(exp_info.animal{CellExp}) filesep];
[CellFile, filePathImage] = uigetfile('*.mat', 'Please select the file containing the cell data', CellImageDir);
load([filePathImage CellFile], 'InVivoImage');
InVivoImage = mat2gray(InVivoImage);

%Ori pref of all functional spines
figure
imagesc(InVivoImage)
colormap('gray')
axis off
hold on
LUT = hsv(180);
spineCounter = 0;
apiSpineCounter = 0;
basalSpineCounter = 0;
for sp = 1:length(InVivoData.ROIsData2)
    switch InVivoData.ROIsData2(sp).type
        case 'Spine'
            if ~isnan(InVivoData.ROIsData2(sp).isReps)
                if  InVivoData.ROIsData2(sp).isReps
                   xPos = InVivoData.ROIsData2(sp).AbsLeft;
                   yPos = InVivoData.ROIsData2(sp).AbsUp;
                   plot(xPos, yPos, 'ok', 'MarkerSize', 5', 'MarkerFaceColor', LUT(1+ floor(InVivoData.ROIsData2(sp).prefOri),:));
                   spineCounter = spineCounter +1;
                   denType = InVivoData.ROIsData2(sp).denType;
                   switch denType
                       case 'basal'
                           basalSpineCounter = basalSpineCounter +1;
                       case 'apical'
                           apiSpineCounter = apiSpineCounter+1;
                   end
                end
            end
            if isnan(InVivoData.ROIsData2(sp).isReps)
                InVivoData.ROIsData2(sp).isRep = 0;
            end
    end
end

todayDate = today('datetime');
F = getframe(gcf);
[X, ~] = frame2im(F);
imwrite(X, [saveDir filesep 'All ' num2str(spineCounter) ' Functional Spines' char(todayDate) '.tif'], 'tiff');

%Ori prefs of all input spines
figure
imagesc(InVivoImage)
colormap('gray')
axis off
hold on
LUT = hsv(180);
inputCounter = 0;
for i = 1:length(inputGroups)
    if ~isnan(inputGroups(i).isResp)
        if inputGroups(i).isResp
           xPos = inputGroups(i).AbsLeft;
           yPos = inputGroups(i).AbsUp;
           plot(xPos, yPos, 'ok', 'MarkerSize', 5', 'MarkerFaceColor', LUT(1+ floor(inputGroups(i).prefOri),:));
           inputCounter = inputCounter +1;
        end
    else
        inputGroups(i).isResp = 0;
    end
end

F = getframe(gcf);
[X, ~] = frame2im(F);
imwrite(X, [saveDir filesep 'All ' num2str(inputCounter) ' FunctionalInputs' char(todayDate) '.tif'], 'tiff');

%all non-inputs white, other colored
figure
imagesc(InVivoImage)
colormap('gray')
axis off
hold on
for sp = 1:length(InVivoData.ROIsData2)
    switch InVivoData.ROIsData2(sp).type
        case 'Spine'
            if ~isnan(InVivoData.ROIsData2(sp).isReps)
                if  InVivoData.ROIsData2(sp).isReps
                   xPos = InVivoData.ROIsData2(sp).AbsLeft;
                   yPos = InVivoData.ROIsData2(sp).AbsUp;
                   plot(xPos, yPos, 'ok', 'MarkerSize', 5', 'MarkerFaceColor', 'white');
                end
            end
            if isnan(InVivoData.ROIsData2(sp).isReps)
                InVivoData.ROIsData2(sp).isRep = 0;
            end
    end
end
hold on 
LUT = hsv(180);
for i = 1:length(inputGroups)
    if ~isnan(inputGroups(i).isResp)
        if inputGroups(i).isResp
           xPos = inputGroups(i).AbsLeft;
           yPos = inputGroups(i).AbsUp;
           plot(xPos, yPos, 'ok', 'MarkerSize', 5', 'MarkerFaceColor', LUT(1+ floor(inputGroups(i).prefOri),:));
           inputCounter = inputCounter +1;
        end
    else
        inputGroups(i).isResp = 0;
    end
end
F = getframe(gcf);
[X, ~] = frame2im(F);
imwrite(X, [saveDir filesep 'All ' num2str(inputCounter) ' FunctionalInputsofAllInputs' char(todayDate) '.tif'], 'tiff');

%% 6.a) Calculate cell aggregate for all spines and all input spines as circular mean of all preferences
%all spines
allPrefDeg = [InVivoData.ROIsData2(find([InVivoData.ROIsData2(:).isReps] == 1)).prefOri];
allPrefAng = deg2rad(allPrefDeg);
allAggregate = rad2deg(circ_mean(allPrefAng'*2)/2);

%input spines
inputPrefDeg = [inputGroups(logical([inputGroups(:).isResp])).prefOri];
inputPrefAng = deg2rad(inputPrefDeg);
inputAggregate = rad2deg(circ_mean(inputPrefAng'*2)/2);

%% 6.b) Calculate summed spine responses 
coc_prop = cbrewer('qual', 'Set1', 9);
%all spines
allSpinesResp = zeros(size(inputGroups(i).meanResp,2)/2, spineCounter);
apiSpinesResp = zeros(size(inputGroups(i).meanResp,2)/2, apiSpineCounter);
basalSpinesResp = zeros(size(inputGroups(i).meanResp,2)/2, basalSpineCounter);
allSpinesEvents = allSpinesResp;
counter=1;
apiCounter = 1;
basCounter = 1;
for sp = 1:length(InVivoData.ROIsData2)
    switch InVivoData.ROIsData2(sp).type
        case 'Spine'
            if ~isnan(InVivoData.ROIsData2(sp).isReps)
                if  InVivoData.ROIsData2(sp).isReps
                    %summed based on mean Resp
                    meanResp = InVivoData.ROIsData2(sp).meanResp; %get the mean resp
                    meanRespFolded = (meanResp(1:8)+meanResp(9:16))/2; %fold it into orientation space
                    meanRespNorm = (meanRespFolded-min(meanRespFolded))/(max(meanRespFolded)-min(meanRespFolded)); %normalize the resp
                    allSpinesResp(:,counter) = meanRespNorm; %add it to the array
                    
                    %summed responses by type
                    denType = InVivoData.ROIsData2(sp).denType;
                    switch denType
                       case 'basal'
                           basalSpinesResp(:,basCounter) = meanRespNorm;
                           basCounter = basCounter+1;
                       case 'apical'
                           apiSpinesResp(:,apiCounter) = meanRespNorm;
                           apiCounter = apiCounter+1;
                    end
                    
                    %summed based on events
                    events = InVivoData.ROIsData2(sp).events; %get all events
                    eventsAllTrials = sum(events,2);            %sum over all trials
                    eventsFolded = (eventsAllTrials(1:8)'+eventsAllTrials(9:16)')/2; %fold it into orientation space
                    eventsNorm = (eventsFolded-min(eventsFolded))/(max(eventsFolded)-min(eventsFolded)); %normalize the events
                    allSpinesEvents(:,counter) = eventsNorm; %add it to the array
                    counter = counter+1;
                end
            end
    end
end
meanSpineResp = mean(allSpinesResp,2); %get the mean resp of all inputs
meanSpineRespNorm2 = (meanSpineResp)/(max(meanSpineResp)); %normalize the resp

meanApiSpineResp = mean(apiSpinesResp,2); %get the mean resp of all inputs
meanApiSpineRespNorm2 = (meanApiSpineResp)/(max(meanApiSpineResp)); %normalize the resp

meanBasalSpineResp = mean(basalSpinesResp,2); %get the mean resp of all inputs
meanBasalSpineRespNorm2 = (meanBasalSpineResp)/(max(meanBasalSpineResp)); %normalize the resp

meanSpineEvents = mean(allSpinesEvents,2); %get the mean events of all inputs
meanSpineEventsNorm2 = (meanSpineEvents)/(max(meanSpineEvents)); %normalize the resp

%inputs spines
allInputResp = zeros(size(inputGroups(i).meanResp,2)/2, inputCounter);
allInputEvents = allInputResp;
counter=1;
for i = 1:length(inputGroups)
    if ~isnan(inputGroups(i).isResp)
        if inputGroups(i).isResp
            %summed based on mean Resp
            meanResp = inputGroups(i).meanResp; %get the mean resp
            meanRespFolded = (meanResp(1:8)+meanResp(9:16))/2; %fold it into orientation space
            meanRespNorm = (meanRespFolded-min(meanRespFolded))/(max(meanRespFolded)-min(meanRespFolded)); %normalize the resp
            allInputResp(:,counter) = meanRespNorm; %add it to the array
            
            %summed based on events
            eventsInputs = inputGroups(i).events; %get all events
            eventsAllTrials = sum(eventsInputs,2);            %sum over all trials
            eventsFolded = (eventsAllTrials(1:8)'+eventsAllTrials(9:16)')/2; %fold it into orientation space
            eventsNorm = (eventsFolded-min(eventsFolded))/(max(eventsFolded)-min(eventsFolded)); %normalize the events
            allInputEvents(:,counter) = eventsNorm; %add it to the array
            counter = counter+1;
        end
    end
end
meanInputResp = mean(allInputResp,2); %get the mean resp of all inputs
meanInputRespNorm = (meanInputResp-min(meanInputResp))/(max(meanInputResp)-min(meanInputResp)); %normalize the resp
meanInputRespNorm2 = (meanInputResp)/(max(meanInputResp)); %normalize the resp

meanInputEvents = mean(allInputEvents,2); %get the mean events of all inputs
meanInputEventsNorm = (meanInputEvents-min(meanInputEvents))/(max(meanInputEvents)-min(meanInputEvents)); %normalize the events
meanInputEventsNorm2 = (meanInputEvents)/(max(meanInputEvents)); %normalize the events

%cell normalized response
cellResp = completeCell{1}.meanResp;
cellRespFolded = (cellResp(1:8)+cellResp(9:16))/2; %fold it into orientation space
cellRespNorm = (cellRespFolded-min(cellRespFolded))/(max(cellRespFolded)-min(cellRespFolded)); %normalize the resp
cellRespNorm2 = (cellRespFolded)/(max(cellRespFolded)); %normalize the resp

%get y axes
angles = linspace(0,180,9);
angles = angles(1:end-1);

% mean resp vs. input resp in comparison to cell tuning
figure
plot(angles, cellRespNorm, 'Color', 'black', 'LineWidth', 5);
hold on
plot(angles, meanSpineRespNorm2, 'Color', coc_prop(9,:), 'LineWidth', 2);
hold on
plot(angles, meanInputRespNorm2, 'Color', coc_prop(1,:), 'LineWidth', 2);
xlabel('Angles in deg')
xticks([0 45 90 135])
ylabel('Normalized responses')
yticks([0 0.25 0.5 0.75 1])
ylim([0 1]);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'SummedResponses.png'))

%shift it to the prefered angle
[~, prefInd] = max(cellRespNorm);
anglesAligned = linspace(-90,90,9);
if prefInd < 5
    addBefore = 5-prefInd-1;
    addUntil =  8-prefInd;
    cellNew = [cellRespNorm(end-addBefore:end), cellRespNorm(1:addUntil)];
    ApiNew = [meanApiSpineRespNorm2(end-addBefore:end)', meanApiSpineRespNorm2(1:addUntil)'];
    BasalNew = [meanBasalSpineRespNorm2(end-addBefore:end)', meanBasalSpineRespNorm2(1:addUntil)'];
    inputNew = [meanInputRespNorm2(end-addBefore:end)', meanInputRespNorm2(1:addUntil)'];
    inputEventsNew = [meanInputEventsNorm2(end-addBefore:end)', meanInputEventsNorm2(1:addUntil)'];
    meanNew = [meanSpineRespNorm2(end-addBefore:end)', meanSpineRespNorm2(1:addUntil)'];
    meanEventsNew = [meanSpineEventsNorm2(end-addBefore:end)', meanSpineEventsNorm2(1:addUntil)'];
else
end

% mean resp vs. input resp in comparison to cell tuning aligned to soma
% respponse
figure
plot(anglesAligned, cellNew, 'Color', 'black', 'LineWidth', 5);
hold on
plot(anglesAligned, meanNew, 'Color', coc_prop(9,:), 'LineWidth', 2);
hold on
plot(anglesAligned, inputNew, 'Color', coc_prop(1,:), 'LineWidth', 2);
xlabel('Angles in deg')
xticks([-90 -45 0 45 90])
ylabel('Normalized responses')
yticks([0 0.25 0.5 0.75 1])
ylim([0 1]);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'SummedResponsesAligned.png'))

% events vs. input event in comparison to cell tuning aligned to soma
% respponse
figure
plot(anglesAligned, cellNew, 'Color', 'black', 'LineWidth', 5);
hold on
plot(anglesAligned, meanEventsNew, 'Color', coc_prop(9,:), 'LineWidth', 2);
hold on
plot(anglesAligned, inputEventsNew, 'Color', coc_prop(1,:), 'LineWidth', 2);
xlabel('Angles in deg')
xticks([-90 -45 0 45 90])
ylabel('Normalized events')
yticks([0 0.25 0.5 0.75 1])
ylim([0 1]);
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'SummedEventsAligned.png'))

% mean resp split by type in comparison to cell tuning
figure
plot(anglesAligned, cellNew, 'Color', 'black', 'LineWidth', 5);
hold on
plot(anglesAligned, meanNew, 'Color', coc_prop(9,:), 'LineWidth', 2);
hold on
plot(anglesAligned, ApiNew, 'Color', coc_prop(2,:), 'LineWidth', 2);
hold on
plot(anglesAligned, BasalNew, 'Color', coc_prop(3,:), 'LineWidth', 2);
xlabel('Angles in deg')
xticks([0 45 90 135])
ylabel('Normalized responses')
yticks([0 0.25 0.5 0.75 1])
ylim([0 1]);
set(gcf, 'color', 'w');
legend('soma', 'allSpines', 'apical', 'basal', 'Location', 'SouthEast')
saveas(gcf, fullfile(saveDir,'SummedResponsesAlignedByType.png'))

%% 7.) Show histogram of delta ori/ori pref for all spines and all input spines
figure
histogram([InVivoData.ROIsData2(find([InVivoData.ROIsData2(:).isReps] == 1)).prefOri],linspace(0,180,8), 'Normalization','probability', 'DisplayStyle', 'stairs', 'EdgeColor', coc_prop(9,:))
hold on 
histogram([inputGroups(logical([inputGroups(:).isResp])).prefOri], linspace(0,180,8), 'Normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', coc_prop(1,:))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'OriPrefHistogram.png'))

figure
subplot(1,2,1)
distributionPlot([InVivoData.ROIsData2(find([InVivoData.ROIsData2(:).isReps] == 1)).prefOri]', 'color', coc_prop(9,:))
ylim([0 180])
subplot(1,2,2)
distributionPlot([inputGroups(logical([inputGroups(:).isResp])).prefOri]', 'color', coc_prop(1,:))
ylim([0 180])
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'OriPrefDistribution.png'))

figure
h(1) = cdfplot([InVivoData.ROIsData2(find([InVivoData.ROIsData2(:).isReps] == 1)).prefOri]);
hold on
h(2) = cdfplot([inputGroups(logical([inputGroups(:).isResp])).prefOri]);
set(h(1),'Color', coc_prop(9,:));
set(h(2),'Color', coc_prop(1,:));
set(gcf, 'color', 'w');
grid off
saveas(gcf, fullfile(saveDir,'OriPrefCumulative.png'))

figure
histogram([InVivoData.ROIsData2(find([InVivoData.ROIsData2(:).isReps] == 1)).deltaOriAbs],linspace(0,90,6), 'Normalization','probability', 'DisplayStyle', 'stairs', 'EdgeColor', coc_prop(9,:))
hold on 
histogram([inputGroups(logical([inputGroups(:).isResp])).deltaOriAbs], linspace(0,90,6), 'Normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', coc_prop(1,:))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'DeltaOriAbsHistogram.png'))

figure
histogram([InVivoData.ROIsData2(find([InVivoData.ROIsData2(:).isReps] == 1)).deltaOri],linspace(-90,90,8), 'Normalization','probability', 'DisplayStyle', 'stairs', 'EdgeColor', coc_prop(9,:))
hold on 
histogram([inputGroups(logical([inputGroups(:).isResp])).deltaOri], linspace(-90,90,8), 'Normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', coc_prop(1,:))
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir,'DeltaOriHistogram.png'))



