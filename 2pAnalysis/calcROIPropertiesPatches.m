function [analysisParams, metadata, data, analysis] = calcROIPropertiesPatches(analysisParams, metadata, data, analysis)
    
    Stimtype = metadata.StimParams.type;

    %% 1.) Chop traces
    disp('Chopping traces')
    [analysis, metadata, data] = ChopStimulusTraceRet(analysis,metadata,data,analysisParams.level, analysisParams.field, 'pre', analysisParams.pre, 'post',analysisParams.post,'windowStart',analysisParams.windowStart, 'windowStop',analysisParams.windowStop);

    
    %% 2.) find maxResponses and significant responses
    disp('Calculating significant responses')
    for i = 1:length(analysis.(analysisParams.field).roi)
        analysis.(analysisParams.field).roi(i).avgRespNormalized = analysis.(analysisParams.field).roi(i).avgStimResponse ./  max(analysis.(analysisParams.field).roi(i).avgStimResponse); 


        pretrialTime= analysis.(analysisParams.field).preTrialTime;
        preTrialIndex= (1:floor(pretrialTime * metadata.TwoPhoton.rate));
        stimWindow=(analysis.(analysisParams.field).windowStart: analysis.(analysisParams.field).windowStop);
        %collect our pretrial interval
        analysis.(analysisParams.field).roi(i).isRespSignificant = false;
        analysis.(analysisParams.field).roi(i).respThreshold = [];
        baselines=analysis.(analysisParams.field).roi(i).stimResponseTrace(:,:,preTrialIndex);
        analysis.(analysisParams.field).roi(i).baselineSD = std(baselines,[],3);
        analysis.(analysisParams.field).roi(i).baselineMean = mean(baselines,3);
        analysis.(analysisParams.field).roi(i).baselineMean(analysis.(analysisParams.field).roi(i).baselineMean < 0) = mean(mean(analysis.(analysisParams.field).roi(i).baselineMean,2),1);
        analysisPeriod=(analysis.(analysisParams.field).windowStart:analysis.(analysisParams.field).windowStop);
        stimResp = analysis.(analysisParams.field).roi(i).stimResponseTrace(:,:,analysisPeriod);
        analysis.(analysisParams.field).roi(i).peaks = max(stimResp,[],3);
        analysis.(analysisParams.field).roi(i).zscore = ([analysis.(analysisParams.field).roi(i).peaks]-[analysis.(analysisParams.field).roi(i).baselineMean])./[analysis.(analysisParams.field).roi(i).baselineSD];
        analysis.(analysisParams.field).roi(i).crosser = sum(analysis.(analysisParams.field).roi(i).zscore > analysisParams.zThresh,2);
        analysis.(analysisParams.field).roi(i).respStim = analysis.(analysisParams.field).roi(i).crosser >= ((metadata.StimParams.numTrials)*analysisParams.fraction);
        if sum(analysis.(analysisParams.field).roi(i).respStim) > 0
            analysis.(analysisParams.field).roi(i).isResponseSignificant = 1;
        else 
            analysis.(analysisParams.field).roi(i).isResponseSignificant = 0;
        end
    end
    
    %% 3.) Calculate ROI properties
    disp('Calculating ROI properties')
    for i = 1:length(data.roi)
        %find preferred Patch for all cells
        switch Stimtype
            case 'Patch'
                medResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace;
                medResponse = mean(medResponse(:,:,stimWindow),3);
                medResponseA = zeros(metadata.StimParams.numPatches/2, size(medResponse,2)*2);
                medResponseA(:,1:size(medResponse,2)) = medResponse(1:metadata.StimParams.numPatches/2,:);
                medResponseA(:,size(medResponse,2)+1:size(medResponse,2)*2) = medResponse(metadata.StimParams.numPatches/2+1:metadata.StimParams.numPatches,:);
                medResponse = medResponseA;
            case 'Retinotopy_2D'
                medResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace;
                medResponse = mean(medResponse(:,:,stimWindow),3);    
            case 'rotatingGratingPatch'
                medResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace;
                medResponse = mean(medResponse(:,:,stimWindow),3);
            otherwise
                medResponse = analysis.(analysisParams.field).roi(i).stimResponseTrace(1:end-1, :, :);
                medResponse = mean(medResponse(:,:,stimWindow),3);
        end

        medResponsePatch= median(medResponse,2);
        [~, prefPatch] = max(medResponsePatch);
        analysis.(analysisParams.field).roi(i).prefPatch = prefPatch;

        %find preferred Elevation for all cells for all 2D stimuli
        if metadata.StimParams.TwoDStim == 1
            switch Stimtype
                case 'Patch'
                    medResponseElev = zeros(metadata.StimParams.numElevation, metadata.StimParams.numAzimuth*metadata.StimParams.numTrials*2);
                    for elev = 1:metadata.StimParams.numElevation
                        medResponseElev(elev,:) = reshape(medResponse((elev-1)*metadata.StimParams.numAzimuth+1:elev*metadata.StimParams.numAzimuth,:),1,metadata.StimParams.numTrials*2*metadata.StimParams.numAzimuth);
                    end
                otherwise
                    medResponseElev = zeros(metadata.StimParams.numElevation, metadata.StimParams.numAzimuth*metadata.StimParams.numTrials);
                    for elev = 1:metadata.StimParams.numElevation
                        medResponseElev(elev,:) = reshape(medResponse((elev-1)*metadata.StimParams.numAzimuth+1:elev*metadata.StimParams.numAzimuth,:),1,metadata.StimParams.numTrials*metadata.StimParams.numAzimuth);
                    end
            end

            medResponseElev= median(medResponseElev,2);
            [~, prefElev] = max(medResponseElev);

            %find preferred Azimuth for all cells
            switch Stimtype
                case 'Patch'
                    medResponseAzi = zeros(metadata.StimParams.numAzimuth, metadata.StimParams.numElevation*metadata.StimParams.numTrials*2);
                    for azi = 1:metadata.StimParams.numAzimuth
                        medResponseAzi(azi,:) = reshape(medResponse(azi:metadata.StimParams.numAzimuth:metadata.StimParams.numPatches/2,:),1,metadata.StimParams.numTrials*2*metadata.StimParams.numElevation);
                    end
                otherwise
                    medResponseAzi = zeros(metadata.StimParams.numAzimuth, metadata.StimParams.numElevation*metadata.StimParams.numTrials);
                    for azi = 1:metadata.StimParams.numAzimuth
                        medResponseAzi(azi,:) = reshape(medResponse(azi:metadata.StimParams.numAzimuth:end,:),1,metadata.StimParams.numTrials*metadata.StimParams.numElevation);
                    end
            end
            medResponseAzi= median(medResponseAzi,2);
            [~, prefAzi] = max(medResponseAzi);

            %write to struct
            analysis.(analysisParams.field).roi(i).prefPatch = prefPatch;
            analysis.(analysisParams.field).roi(i).prefElev = prefElev;
            analysis.(analysisParams.field).roi(i).prefAzi = prefAzi;
            analysis.(analysisParams.field).roi(i).prefElevDeg = metadata.StimParams.stimPosY(prefElev);
            analysis.(analysisParams.field).roi(i).prefAziDeg = metadata.StimParams.stimPosX(prefAzi);
        end
    end
end