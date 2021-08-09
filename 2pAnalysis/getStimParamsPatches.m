function [analysisParams, metadata] = getStimParamsPatches(analysisParams, metadata)
    Stimtype = metadata.StimParams.type;
    switch Stimtype
        case 'Retinotopy_2D'
            metadata.StimParams.isi = metadata.StimParams.ISI;
            metadata.StimParams.numTrials = metadata.StimParams.numberOfTrials;
            analysisParams.windowStop=metadata.StimParams.stimDuration+0.5;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            stimSize = strsplit(metadata.StimParams.stimSize(2:end-1), ',');
            metadata.StimParams.stimSize = [str2double(stimSize{1}), str2double(stimSize{2})];
            if analysisParams.special
                metadata.StimParams.stimSize(1) = metadata.StimParams.stimSize(1)/2;
                metadata.StimParams.stimSize(2) = metadata.StimParams.stimSize(2)/2;
            end

            if mod((metadata.StimParams.endPointy - metadata.StimParams.startPointy), (str2double(metadata.StimParams.stimSize(2)))) == 0
                metadata.StimParams.numElevation = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(metadata.StimParams.stimSize(2)));
                metadata.StimParams.numAzimuth = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(metadata.StimParams.stimSize(1)));
            else
                metadata.StimParams.numElevation = floor((metadata.StimParams.endPointy - metadata.StimParams.startPointy)/(metadata.StimParams.stimSize(2)))+1;
                metadata.StimParams.numAzimuth = floor((metadata.StimParams.endPointx - metadata.StimParams.startPointx)/(metadata.StimParams.stimSize(1)))+1;
            end

            metadata.StimParams.stimPosX = linspace(metadata.StimParams.startPointx,(metadata.StimParams.numAzimuth-1)*str2double(metadata.StimParams.stimSize(2))+metadata.StimParams.startPointx,metadata.StimParams.numAzimuth);
            metadata.StimParams.stimPosY = linspace(metadata.StimParams.startPointy,(metadata.StimParams.numElevation-1)*str2double(metadata.StimParams.stimSize(2))+metadata.StimParams.startPointy,metadata.StimParams.numElevation);

            metadata.StimParams.minAzi = metadata.StimParams.startPointx;
            metadata.StimParams.maxElev = metadata.StimParams.startPointy;
            if min(metadata.StimParams.stimSize) < 2
                metadata.StimParams.maxAzi = floor(metadata.StimParams.startPointx + (metadata.StimParams.numAzimuth) * metadata.StimParams.stimSize(1) - 0.5);
                metadata.StimParams.minElev = floor(metadata.StimParams.startPointy + (metadata.StimParams.numElevation) * metadata.StimParams.stimSize(2) - 0.5);
            else
                metadata.StimParams.maxAzi = floor(metadata.StimParams.startPointx + (metadata.StimParams.numAzimuth) * metadata.StimParams.stimSize(1) - 1);
                metadata.StimParams.minElev = floor(metadata.StimParams.startPointy + (metadata.StimParams.numElevation) * metadata.StimParams.stimSize(2) - 1);
            end

            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
            metadata.StimParams.PatchX = repmat(metadata.StimParams.stimPosX',1,metadata.StimParams.numAzimuth);
            metadata.StimParams.PatchY = repmat(metadata.StimParams.stimPosY,metadata.StimParams.numElevation,1);
            metadata.StimParams.TwoDStim = 1;
            metadata.StimParams.numCol = 1;
            metadata.StimParams.numX = repmat(linspace(1, metadata.StimParams.numAzimuth, metadata.StimParams.numAzimuth), 1, metadata.StimParams.numElevation);
            metadata.StimParams.numY = reshape(repmat(linspace(1, metadata.StimParams.numElevation, metadata.StimParams.numElevation), metadata.StimParams.numAzimuth, 1), metadata.StimParams.numPatches, [])';

        case 'PatchGrating'
            analysisParams.windowStop=0.5;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
            metadata.StimParams.numElevation = metadata.StimParams.numElevation;
            metadata.StimParams.numAzimuth = metadata.StimParams.numAzimuth;
            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
            minAzim = 10 - ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*metadata.StimParams.numAzimuth+minAzim,metadata.StimParams.numAzimuth);
            maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*metadata.StimParams.numElevation,maxElev,metadata.StimParams.numElevation);
            metadata.StimParams.TwoDStim = 1;
        case 'RF_localbar'
            analysisParams.windowStop=0.5;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
            metadata.StimParams.numElevation = metadata.StimParams.numElevation;
            metadata.StimParams.numAzimuth = metadata.StimParams.numAzimuth;
            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
            minAzim = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*metadata.StimParams.numAzimuth+minAzim,metadata.StimParams.numAzimuth);
            maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*metadata.StimParams.numElevation,maxElev,metadata.StimParams.numElevation);
            metadata.StimParams.TwoDStim = 0;
        case 'Patch'
            analysisParams.windowStop=1;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
            metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth*2;
            minAzim = -((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.stimPosX = linspace(minAzim,metadata.StimParams.stimSize(1)*(metadata.StimParams.numAzimuth-1)+minAzim,metadata.StimParams.numAzimuth);
            maxElev = ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.stimPosY = linspace(maxElev-metadata.StimParams.stimSize(2)*(metadata.StimParams.numElevation-1),maxElev,metadata.StimParams.numElevation);
            metadata.StimParams.TwoDStim = 1;
            metadata.StimParams.numCol = 2;
            numUniqStims = metadata.StimParams.numElevation * metadata.StimParams.numAzimuth;
            metadata.StimParams.minAzi = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.maxAzi = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) - metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.numX = repmat(linspace(1, metadata.StimParams.numAzimuth, metadata.StimParams.numAzimuth), 1, metadata.StimParams.numElevation);
            metadata.StimParams.minElev = metadata.StimParams.centerPoint(2) - ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) + metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.maxElev = metadata.StimParams.centerPoint(2) + ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.numY = reshape(repmat(linspace(1, metadata.StimParams.numElevation, metadata.StimParams.numElevation), metadata.StimParams.numAzimuth, 1), numUniqStims, [])';
        case 'azimuthBars'
            analysisParams.windowStop=1;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.numElevation = 1;
            metadata.StimParams.numAzimuth = metadata.StimParams.numAzimuth;
            if metadata.StimParams.minAzim == 0
                metadata.StimParams.startAzim = 0;
            else
                metadata.StimParams.startAzim = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) +metadata.StimParams. stimSize(1)/2;
            end
            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
            metadata.StimParams.stimPosX = linspace(metadata.StimParams.startAzim,metadata.StimParams.stimSize(2)*(metadata.StimParams.numAzimuth-1)+metadata.StimParams.startAzim,metadata.StimParams.numAzimuth);
            metadata.StimParams.stimPosY = 1;
            metadata.StimParams.TwoDStim = 0;
        case 'RF_local'
            analysisParams.windowStop=1;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.numElevation = 1;
            metadata.StimParams.numElevation = 1;
            metadata.StimParams.numAzimuth = metadata.StimParams.numAzimuth;
            if metadata.StimParams.minAzim == 0
                metadata.StimParams.startAzim = 0;
            else
                metadata.StimParams.startAzim = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(2))/2) +metadata.StimParams. stimSize(2)/2;
            end
            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth;
            metadata.StimParams.stimPosX = linspace(metadata.StimParams.startAzim,metadata.StimParams.stimSize(2)*(metadata.StimParams.numAzimuth-1)+metadata.StimParams.startAzim,metadata.StimParams.numAzimuth);
            metadata.StimParams.stimPosY = 1;
            metadata.StimParams.TwoDStim = 0;
        case 'RetWedge'
            metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
            analysisParams.windowStop=metadata.StimParams.stimDuration;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.numPatches = metadata.StimParams.numWedges;
            metadata.StimParams.TwoDStim = 0;
        case 'Ret_Annulus'
            metadata.StimParams.TwoDStim = 0;
            metadata.StimParams.stimDuration = metadata.StimParams.stimDuration * metadata.StimParams.numOrientations;
            analysisParams.windowStop=metadata.StimParams.stimDuration;
            analysisParams.windowStart=0;
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            metadata.StimParams.numPatches = metadata.StimParams.numSizes;
        case 'rotatingGratingPatch'
            analysisParams.pre= 0.5;
            analysisParams.post = 0.5;
            analysisParams.windowStart=0;
            analysisParams.windowStop=metadata.StimParams.stimDuration;
            metadata.StimParams.doBlank = 0;
            metadata.StimParams.numElevation = metadata.StimParams.numStimElev;
            metadata.StimParams.numAzimuth = metadata.StimParams.numStimAzi;
            metadata.StimParams.numPatches = metadata.StimParams.numAzimuth*metadata.StimParams.numElevation;
            metadata.StimParams.TwoDStim = 1;
            metadata.StimParams.stimPosX = linspace(metadata.StimParams.centerPoint(1),metadata.StimParams.stimSize(1)*(metadata.StimParams.numAzimuth-1)+metadata.StimParams.centerPoint(1),metadata.StimParams.numAzimuth);
            metadata.StimParams.stimPosY = linspace(metadata.StimParams.centerPoint(2),metadata.StimParams.stimSize(2)*(metadata.StimParams.numElevation-1)+metadata.StimParams.centerPoint(2),metadata.StimParams.numElevation);
            metadata.StimParams.numCol = 1;
            numUniqStims = metadata.StimParams.numElevation * metadata.StimParams.numAzimuth;
            metadata.StimParams.minAzi = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.maxAzi = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) - metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.numX = repmat(linspace(1, metadata.StimParams.numAzimuth, metadata.StimParams.numAzimuth), 1, metadata.StimParams.numElevation);
            metadata.StimParams.minElev = metadata.StimParams.centerPoint(2) - ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) + metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.maxElev = metadata.StimParams.centerPoint(2) + ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.numY = reshape(repmat(linspace(1, metadata.StimParams.numElevation, metadata.StimParams.numElevation), metadata.StimParams.numAzimuth, 1), numUniqStims, [])';
        case 'SparseNoise'
            %read information from csv File
            cvsFile = [analysisParams.Sp2dir analysisParams.animal filesep analysisParams.cvsFile];
            metadata.StimParams.stimOrder = csvread(cvsFile);
            if size(metadata.StimParams.StimOnTimes,2) < metadata.StimParams.stimOrder
                missingTriggerSites = find(diff(metadata.StimParams.StimOnTimes(2,:)) > 0.65);
                addIns = 0;
                medianTriggerTime = median(diff(metadata.StimParams.StimOnTimes(2,:)));
                for m = 1:length(missingTriggerSites)
                    lastTrigger = missingTriggerSites(m)+addIns;
                    timeDiff = metadata.StimParams.StimOnTimes(2,lastTrigger+1)-metadata.StimParams.StimOnTimes(2,lastTrigger);
                    triggersToAdd = round(timeDiff/medianTriggerTime)-1;
                    matrixToAdd = ones(2,triggersToAdd);
                    for t=1:triggersToAdd
                        matrixToAdd(2,t)=metadata.StimParams.StimOnTimes(2,lastTrigger)+t*medianTriggerTime;
                    end
                    metadata.StimParams.StimOnTimes = [metadata.StimParams.StimOnTimes(:,1:lastTrigger) matrixToAdd metadata.StimParams.StimOnTimes(:,lastTrigger+1:end)];
                    addIns = addIns + tiggersToAdd;
                end
            end

            %read information from .mat file
            timeStampSplit = split(cvsFile,'Noise_');
            timeStamp = timeStampSplit{2};
            StimInfoFile = [analysisParams.Sp2dir analysisParams.animal filesep 'StimulusInformation_' timeStamp(1:end-4) '.mat'];
            imageInformation = load(StimInfoFile);
            try
                metadata.StimParams.stimSize = imageInformation.stimSize; 
                metadata.StimParams.numElevation = imageInformation.PresArea(1) / imageInformation.stimSize(1); 
                metadata.StimParams.numAzimuth = imageInformation.PresArea(2) / imageInformation.stimSize(2);
            catch
                metadata.StimParams.stimSize = [5 5]; 
                metadata.StimParams.numElevation = 5; 
                metadata.StimParams.numAzimuth = 5;
            end
            metadata.StimParams.numCol = 2;
            metadata.StimParams.numPatches = metadata.StimParams.numElevation*metadata.StimParams.numAzimuth*metadata.StimParams.numCol;

            stimArray = reshape(imageInformation.allUsedStims, 4, size(imageInformation.selectedFrames,2));
            orderedStimArray = stimArray(:,metadata.StimParams.stimOrder);
            allStims = reshape(orderedStimArray,1,size(orderedStimArray,1)*size(orderedStimArray,2));
            for P = 1:metadata.StimParams.numPatches
                temp= ceil(find(allStims == P)./4);
                stimTimesAll{P} = temp(temp <= size(metadata.StimParams.StimOnTimes,2));
            end
            metadata.StimParams.stimTimesAll = stimTimesAll;

            metadata.StimParams.minAzi = metadata.StimParams.centerPoint(1) - ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) + metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.maxAzi = metadata.StimParams.centerPoint(1) + ((metadata.StimParams.numAzimuth * metadata.StimParams.stimSize(1))/2) - metadata.StimParams.stimSize(1)/2;
            metadata.StimParams.numX = repmat(linspace(1, metadata.StimParams.numAzimuth, metadata.StimParams.numAzimuth), 1, metadata.StimParams.numElevation);
            metadata.StimParams.minElev = metadata.StimParams.centerPoint(2) - ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) + metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.maxElev = metadata.StimParams.centerPoint(2) + ((metadata.StimParams.numElevation * metadata.StimParams.stimSize(2))/2) - metadata.StimParams.stimSize(2)/2;
            metadata.StimParams.numY = reshape(repmat(linspace(1, metadata.StimParams.numElevation, metadata.StimParams.numElevation), metadata.StimParams.numAzimuth, 1), metadata.StimParams.numPatches/2, [])';

            metadata.StimParams.numTrials = min(cellfun('size', stimTimesAll, 2));
    end
end