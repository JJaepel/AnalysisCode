function [analysis, metadata] = calculateGratingMetrices(analysisParams, metadata, analysis, i)

theta= 0:2*pi/(metadata.StimParams.numDirections): 2*pi-2*pi/(metadata.StimParams.numDirections);%theta is needed for the fitting
thetaOri = theta(1:metadata.StimParams.numOrientations);
analysis.(analysisParams.field).stimWindow=(analysis.(analysisParams.field).stimStart: analysis.(analysisParams.field).stimStop);

%for multiple sf or tfs or ipsi/contra/binocular responses, find the prefered stimulus
if metadata.StimParams.numTf * metadata.StimParams.numSf > 1 || analysisParams.stimType == 5
    [~, ind] = max(analysis.dff.roi(i).avgStimResponse(1:end-1));
    prefCon = ceil(ind/metadata.StimParams.numDirections);
    dirResponse = zeros(metadata.StimParams.numDirections,metadata.StimParams.numTrials,size(analysis.(analysisParams.field).roi(i).stimResponseTrace,3));
    j = 1;
    for direction = (prefCon-1)*metadata.StimParams.numDirections+1:prefCon*metadata.StimParams.numDirections
        dirResponse(j,:,:) = analysis.(analysisParams.field).roi(i).stimResponseTrace(direction,:, :);
        j = j+1; 
    end
    switch analysisParams.stimType
        case 2
            analysis.(analysisParams.field).roi(i).prefSfStim = prefCon;
            analysis.(analysisParams.field).roi(i).prefSf = metadata.StimParams.spatialFreq(prefCon);
        case 4
            analysis.(analysisParams.field).roi(i).prefTfStim = prefCon;
            analysis.(analysisParams.field).roi(i).prefTf = metadata.StimParams.temporalFreq(prefCon);
        case 5
            analysis.(analysisParams.field).roi(i).prefCon = prefCon;
    end
else
    dirResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace(1:end-1, :, :);
end

%calculate the median response curve for each direction and orientation
dirResponse = mean(dirResponse(:,:,analysis.(analysisParams.field).stimWindow),3); %average over the stimWindow
Responses = dirResponse; %this contains the variability over the trials and is needed for Cohen's D

medResponse(1:size(dirResponse,1)/2, 1:size(dirResponse,2)) = dirResponse(1:size(dirResponse,1)/2,:); %rearrange so that the two corresponding directions can be summed
medResponse(1:size(dirResponse,1)/2, 1+size(dirResponse,2):2*size(dirResponse,2))= dirResponse(1+size(dirResponse,1)/2:size(dirResponse,1),:);
allOriResponses = medResponse;
medResponse= median(medResponse,2);
oriResponse = medResponse(:)';  %this is the median response for each orientation
dirResponse = median(dirResponse,2)'; %this is the median response for each direction

%OSI/DSI: (Pref Ori - Orth Ori)/ (Pref Ori + Orth Ori)
[analysis.(analysisParams.field).roi(i).OSI] = computeOSI(metadata.StimParams.numOrientations ,oriResponse);
[analysis.(analysisParams.field).roi(i).DSI] = computeDSI(metadata.StimParams.numOrientations *2, dirResponse, 'preferred');

%CircVar, DirCircVar
[analysis.(analysisParams.field).roi(i).OriCircVar] = 1-OriCircularVariance(dirResponse);
[analysis.(analysisParams.field).roi(i).DirCircVar] = 1-DirCircularVariance(dirResponse);

%calculate bandwidth
[analysis.(analysisParams.field).roi(i).Bandwidth] = computeBandWidth(dirResponse);

%Cohen's D
[analysis.(analysisParams.field).roi(i).cohensD] = computeCohensD(metadata.StimParams.numOrientations ,Responses);

%Calculate Fano factor for preferred stimulus
[analysis.(analysisParams.field).roi(i).fanoFactor] = computeFanoFactor(Responses);

%calculate Variability index
[analysis.(analysisParams.field).roi(i).VI] = computeVI(Responses);

%Find preferred orientation based on gaussian fit
[analysis.(analysisParams.field).roi(i).preferredOrientation,analysis.(analysisParams.field).roi(i).coeffOr,analysis.(analysisParams.field).roi(i).rsqOr] = ComputePreferredOrientations(oriResponse, thetaOri);
analysis.(analysisParams.field).roi(i).OSIFit = computeOSIFit(analysis.(analysisParams.field).roi(i).coeffOr); %compute OSI from fit


%Find preferred direction based on double gaussian fit
try
    [analysis.(analysisParams.field).roi(i).preferredDirection] = ComputePreferredDirection(dirResponse, theta);
catch
    [analysis.(analysisParams.field).roi(i).preferredDirection] = NaN;
end

interval = metadata.StimParams.directions(2) - metadata.StimParams.directions(1);
ind = round(analysis.(analysisParams.field).roi(i).preferredOrientation / interval);
ind = ind + 1;
if ind > metadata.StimParams.numOrientations || isnan(ind)
    ind = 1;
end
analysis.(analysisParams.field).roi(i).prefOriStimInd = ind;

%check if significantly ori/dir-selective
vectorRespOri = (allOriResponses'*transpose(exp(sqrt(-1)*2*mod(metadata.StimParams.orientations*pi/180,pi))));
try
    [analysis.(analysisParams.field).roi(i).OriSignH, analysis.(analysisParams.field).roi(i).pOriSign] =  hotellingt2test(vectorRespOri, [0 0]);
    [analysis.(analysisParams.field).roi(i).DirSignH, analysis.(analysisParams.field).roi(i).pDirSign] =  computeDirSignificanceDotProduct(metadata.StimParams.directions, Responses');
catch
    analysis.(analysisParams.field).roi(i).OriSignH = NaN; analysis.(analysisParams.field).roi(i).pOriSign = NaN;
    analysis.(analysisParams.field).roi(i).DirSignH = NaN; analysis.(analysisParams.field).roi(i).pDirSign = NaN;
end

if metadata.StimParams.numTf * metadata.StimParams.numSf > 1 || analysisParams.stimType == 5
    if analysisParams.stimType == 5
        metadata.StimParams.numCon = 3;
    else
        metadata.StimParams.numCon = max([metadata.StimParams.numSf metadata.StimParams.numTf]);
    end

    % compute ori preference at each condition
    for con = 1:metadata.StimParams.numCon
        dirResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace((con-1)*metadata.StimParams.numDirections+1:con*metadata.StimParams.numDirections, :,:);
        dirResponse = mean(dirResponse(:,:,analysis.(analysisParams.field).stimWindow),3); %average over the stimWindow
        medResponse(1:size(dirResponse,1)/2, 1:size(dirResponse,2)) = dirResponse(1:size(dirResponse,1)/2,:); %rearrange so that the two corresponding directions can be summed
        medResponse(1:size(dirResponse,1)/2, 1+size(dirResponse,2):2*size(dirResponse,2))= dirResponse(1+size(dirResponse,1)/2:size(dirResponse,1),:);
        medResponse= median(medResponse,2);
        oriResponse = medResponse(:)';  %this is the median response for each orientation
        dirResponse = median(dirResponse,2)';

        %find max response per condition
        maxResponse(con)=max(dirResponse);

        %Find preferred orientation based on gaussian fit
        [prefOri_allCon(con),coeffOr,~] = ComputePreferredOrientations(oriResponse, thetaOri);
        OSIFit_allCon(con) = computeOSIFit(coeffOr); %compute OSI from fit

        %Find preferred direction based on double gaussian fit
        [prefDir_allCon(con)] = ComputePreferredDirection(dirResponse, theta);
        [DSI_allCon(con)] = computeDSI(metadata.StimParams.numDirections, dirResponse, 'preferred');
    end

    % calculate indices for spatial/temporal preference by looking
    % at the invidivual conditions - make array with responses at
    % preferred direction

    ConSI = (max(maxResponse) - min(maxResponse)) ./ min(maxResponse);
    switch analysisParams.stimType
        case 2
            analysis.(analysisParams.field).roi(i).preferredOrientation_allSf = prefOri_allCon;
            analysis.(analysisParams.field).roi(i).OSIFit_allSf = OSIFit_allCon;
            analysis.(analysisParams.field).roi(i).preferredDirection_allSf = prefDir_allCon;
            analysis.(analysisParams.field).roi(i).DSI_allSf = DSI_allCon;
            analysis.(analysisParams.field).roi(i).SFSI = ConSI;
            analysis.(analysisParams.field).roi(i).SFVar = computeSFVar(maxResponse,metadata.StimParams.spatialFreq);
        case 3
            analysis.(analysisParams.field).roi(i).preferredOrientation_allTf = prefOri_allCon;
            analysis.(analysisParams.field).roi(i).OSIFit_allTf = OSIFit_allCon;
            analysis.(analysisParams.field).roi(i).preferredDirection_allTf = prefDir_allCon;
            analysis.(analysisParams.field).roi(i).DSI_allTf = DSI_allCon;
            analysis.(analysisParams.field).roi(i).TFSI = ConSI;
        case 5
            analysis.(analysisParams.field).roi(i).preferredOrientation_allEyes = prefOri_allCon;
            analysis.(analysisParams.field).roi(i).OSIFit_allEyes = OSIFit_allCon;
            analysis.(analysisParams.field).roi(i).preferredDirection_allEyes = prefDir_allCon;
            analysis.(analysisParams.field).roi(i).DSI_allEyes = DSI_allCon;analysis.(analysisParams.field).roi(i).ODI = ConSI;
    end
else
    metadata.StimParams.numCon = 1;
end